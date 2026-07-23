# R/wgs_dragen_functions.R
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(gridExtra)
  library(grid)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

StandardizeChrom <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x <- ifelse(x %in% c("23"), "X", x)
  x <- ifelse(x %in% c("24"), "Y", x)
  x
}

ReadDragenAnnotationRaw <- function(anno_file) {
  data.table::fread(anno_file, data.table = FALSE)
}

ReadDragenCoverageRaw <- function(coverage_file, gender) {
  cov <- data.table::fread(coverage_file, data.table = FALSE)

  if (ncol(cov) < 5) {
    stop("Coverage file must have at least five columns. The 5th column should be coverage.")
  }

  colnames(cov)[5] <- "cov"

  required_cov <- c("contig", "start", "stop", "cov")
  missing_cov <- setdiff(required_cov, names(cov))
  if (length(missing_cov) > 0) {
    stop("Coverage file is missing columns: ", paste(missing_cov, collapse = ", "))
  }

  cov <- cov %>%
    dplyr::select(contig, start, stop, cov) %>%
    dplyr::mutate(
      contig = StandardizeChrom(contig),
      start = as.numeric(start),
      stop = as.numeric(stop),
      cov = as.numeric(cov)
    )

  if (gender == "female") {
    cov <- cov %>% dplyr::filter(contig != "Y")
  }

  cov
}

PreprocessDragenAnnotation <- function(df, gender, user_purity, user_diploid_coverage) {
  required_cols <- c(
    "Chrom", "QUAL", "SubClone", "Start", "End", "GT", "SD", "SM",
    "CN", "CNF", "MAF", "Call", "FILTER", "EstimatedTumorPurity",
    "DiploidCoverage", "OverallPloidy"
  )

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Annotation file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- df %>%
    dplyr::select(dplyr::all_of(required_cols)) %>%
    dplyr::mutate(
      EstimatedTumorPurity = user_purity,
      DiploidCoverage = user_diploid_coverage,
      Chrom = StandardizeChrom(Chrom),
      SM = 2^as.numeric(SM),
      CNF = ifelse(is.na(CN), SM * 2, CNF),
      CN = ifelse(is.na(SubClone), round(CNF), CNF),
      FILTER = gsub("segmentMean|;", "", FILTER),
      FILTER = ifelse(FILTER == "", "PASS", FILTER),
      QUAL = ifelse(FILTER == "PASS", 1, 0),
      MAF = as.numeric(MAF),
      FILTER = ifelse(FILTER == "PASS", "PASS", "FAILED")
    )

  colnames(df)[1] <- "seqnames"
  df$seqnames <- StandardizeChrom(df$seqnames)

  if (gender == "female") {
    df <- df %>% dplyr::filter(seqnames != "Y")
  }

  df
}

BinRawCoverageFast <- function(cov, cov_binsize = 200000) {
  # Pre-bin raw coverage once. Slider changes can reuse this table.
  cov %>%
    dplyr::mutate(
      midpoint = floor((start + stop) / 2),
      bin_end = ceiling(midpoint / cov_binsize) * cov_binsize,
      bin_start = bin_end - cov_binsize + 1
    ) %>%
    dplyr::group_by(contig, bin_start, bin_end) %>%
    dplyr::summarize(
      median_cov = median(cov, na.rm = TRUE),
      .groups = "drop"
    )
}

NormalizeCoverageBinsByUserModel <- function(cov_bins_raw, gender, purity, diploid_coverage) {
  if (is.na(diploid_coverage) || diploid_coverage <= 0) {
    stop("Diploid coverage must be greater than 0.")
  }

  if (is.na(purity) || purity < 0 || purity > 1) {
    stop("Purity must be between 0 and 1.")
  }

  cov_bins <- cov_bins_raw %>%
    dplyr::mutate(
      contig = StandardizeChrom(contig),
      median_cov = as.numeric(median_cov)
    )

  if (gender == "female") {
    cov_bins <- cov_bins %>% dplyr::filter(contig != "Y")
  }

  if (purity <= 0) {
    cov_bins <- cov_bins %>%
      dplyr::mutate(cov_to_cnf = median_cov / diploid_coverage * 2)
  } else {
    cov_bins <- cov_bins %>%
      dplyr::mutate(
        normal_cn = dplyr::case_when(
          gender == "male" & contig %in% c("X", "Y") ~ 1,
          TRUE ~ 2
        ),
        cov_to_cnf = (
          median_cov / diploid_coverage * 2 -
            normal_cn * (1 - purity)
        ) / purity
      )
  }

  cov_bins %>%
    dplyr::mutate(
      cov_to_cnf = ifelse(is.finite(cov_to_cnf), cov_to_cnf, NA_real_),
      smoothed_cnf = cov_to_cnf,
      smoothed_bin_cnf = round(smoothed_cnf * 10) / 10
    )
}

ReadPrepareRoundDragenBAllele <- function(ballele_file, gender) {
  ai <- data.table::fread(ballele_file, data.table = FALSE)

  required_cols <- c(
    "contig", "start", "stop", "refAllele", "allele1", "allele2",
    "allele1Count", "allele2Count"
  )

  missing_cols <- setdiff(required_cols, names(ai))
  if (length(missing_cols) > 0) {
    stop("B-allele file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  ai <- ai %>%
    dplyr::filter(allele1Count + allele2Count >= 10) %>%
    dplyr::select(dplyr::all_of(required_cols)) %>%
    dplyr::filter(!grepl(",", allele1), !grepl(",", allele2)) %>%
    dplyr::mutate(
      contig = StandardizeChrom(contig),
      seqnames = contig,
      start = as.numeric(start),
      stop = as.numeric(stop)
    )

  if (gender == "female") {
    ai <- ai %>% dplyr::filter(contig != "Y")
  }

  allelePreference <- c(0:3)
  names(allelePreference) <- c("A", "T", "G", "C")

  ai <- ai %>%
    dplyr::mutate(
      af = ifelse(
        allelePreference[allele1] > allelePreference[allele2],
        allele1Count / (allele1Count + allele2Count),
        allele2Count / (allele1Count + allele2Count)
      ),
      norm_af = round(af * 20) / 20
    ) %>%
    dplyr::select(seqnames, contig, start, stop, allele1Count, allele2Count, norm_af)

  ai
}

CleanHomAlt <- function(clean_ai, call_seg) {
  call_seg <- call_seg %>% dplyr::filter(MAF > 0.7 | MAF < 0.3)

  if (nrow(call_seg) == 0) {
    return(clean_ai %>% dplyr::filter(allele1Count > 1, allele2Count > 1))
  }

  noloh_range <- GenomicRanges::makeGRangesFromDataFrame(
    call_seg,
    seqnames.field = "seqnames",
    start.field = "Start",
    end.field = "End",
    keep.extra.columns = FALSE
  )

  ai_range <- GenomicRanges::makeGRangesFromDataFrame(
    clean_ai,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "stop",
    keep.extra.columns = TRUE
  )

  overlaps <- GenomicRanges::findOverlaps(ai_range, noloh_range)

  clean_ai$matched_seqnames <- NA
  clean_ai$matched_start <- NA
  clean_ai$matched_end <- NA

  if (length(overlaps) > 0) {
    overlapping_noloh <- noloh_range[S4Vectors::subjectHits(overlaps)]
    overlapping_noloh_df <- as.data.frame(overlapping_noloh)
    clean_ai[
      S4Vectors::queryHits(overlaps),
      c("matched_seqnames", "matched_start", "matched_end")
    ] <- overlapping_noloh_df[, c(1:3)]
  }

  clean_ai %>%
    dplyr::filter(allele1Count > 1, allele2Count > 1) %>%
    dplyr::filter(
      !(is.na(matched_seqnames) & (allele1Count < 4 | allele2Count < 4)) |
        seqnames %in% c("X", "Y")
    )
}

AiBinSmooth <- function(x) {
  counts <- data.frame(table(x))
  mode_value <- counts[order(counts$Freq, decreasing = TRUE), ]
  colnames(mode_value) <- c("value", "Freq")
  total_dots <- length(x)
  mode_value$prop <- mode_value$Freq / total_dots

  if (total_dots > 30) {
    re <- rep(mode_value$value, round(mode_value$prop * 30))
  } else {
    re <- x
  }

  as.numeric(as.character(re))
}

SmoothAi <- function(df, ai_binsize = 100000, gender) {
  if (gender == "female") {
    df <- df %>% dplyr::filter(seqnames != "Y")
  }

  ai_bin_chr <- split(df, f = df$seqnames)
  ai_bin_chr <- Filter(function(x) nrow(x) > 0, ai_bin_chr)

  sep_bin <- lapply(ai_bin_chr, function(chr_ai) {
    chr_ai_smoothed <- chr_ai %>%
      dplyr::mutate(bin_norm_id = ceiling(start / ai_binsize) * ai_binsize) %>%
      dplyr::group_by(seqnames, bin_norm_id) %>%
      dplyr::summarize(smoothed_ai = list(AiBinSmooth(norm_af)), .groups = "drop") %>%
      tidyr::unnest(smoothed_ai) %>%
      dplyr::mutate(bin_start = bin_norm_id - ai_binsize + 1) %>%
      dplyr::select(seqnames, bin_start, bin_norm_id, dplyr::everything()) %>%
      dplyr::ungroup()

    colnames(chr_ai_smoothed)[3] <- "bin_end"
    chr_ai_smoothed
  })

  do.call(rbind, sep_bin)
}

PlotWgsCnvGrid <- function(df_cov, df_ai, call_seg, gender, prefix) {
  if (gender == "male") {
    chrom_levels <- c(as.character(1:22), "X", "Y")
  } else {
    chrom_levels <- c(as.character(1:22), "X")
  }
  
  color <- rep(
    c("darkblue", "darkgreen", "darkred", "darkorchid4", "darkgoldenrod"),
    length.out = length(chrom_levels)
  )
  names(color) <- chrom_levels
  
  df_cov <- df_cov %>%
    dplyr::mutate(seqnames = StandardizeChrom(contig))
  
  df_cov$seqnames <- factor(StandardizeChrom(df_cov$seqnames), levels = chrom_levels)
  df_ai$seqnames <- factor(StandardizeChrom(df_ai$seqnames), levels = chrom_levels)
  call_seg$seqnames <- factor(StandardizeChrom(call_seg$seqnames), levels = chrom_levels)
  
  y_lim <- 8
  
  p_margin <- ggplot2::margin(t = 1, r = 6, b = 0.5, l = 38, unit = "pt")
  q_margin <- ggplot2::margin(t = 0, r = 6, b = 0.5, l = 38, unit = "pt")
  qc_margin <- ggplot2::margin(t = 0, r = 6, b = 0.5, l = 38, unit = "pt")
  
  line_pos <- c(1, 3, 4, 5, 6, 7, 8)
  
  df_cov <- df_cov %>%
    dplyr::mutate(
      smoothed_bin_cnf = ifelse(smoothed_bin_cnf >= y_lim, y_lim, smoothed_bin_cnf),
      smoothed_bin_cnf = ifelse(smoothed_bin_cnf < 0, 0, smoothed_bin_cnf)
    )
  
  call_seg_maf <- call_seg %>%
    dplyr::mutate(
      MAF = as.numeric(MAF),
      MAF = ifelse(MAF > 1, MAF / 100, MAF),
      baf_lower = pmin(MAF, 1 - MAF),
      baf_upper = 1 - baf_lower
    ) %>%
    dplyr::filter(
      !is.na(baf_lower),
      !is.na(baf_upper),
      baf_lower >= 0,
      baf_upper <= 1
    )
  
  y_axis_theme <- ggplot2::theme(
    axis.title.y = ggplot2::element_text(size = 10, color = "black"),
    axis.text.y = ggplot2::element_text(size = 9, color = "black"),
    axis.text.y.left = ggplot2::element_text(size = 9, color = "black"),
    axis.ticks.y = ggplot2::element_line(color = "black"),
    axis.line.y.left = ggplot2::element_line(color = "black")
  )
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = 2,
      color = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = line_pos,
      color = "grey",
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      data = df_cov,
      ggplot2::aes(
        x = bin_start,
        xend = bin_end,
        y = smoothed_bin_cnf - 0.1,
        yend = smoothed_bin_cnf + 0.1,
        color = seqnames
      ),
      alpha = 0.2
    ) +
    ggplot2::scale_color_manual(values = color, drop = FALSE) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(seqnames),
      scales = "free_x",
      space = "free_x",
      drop = FALSE
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_continuous(
      limits = c(0, y_lim),
      breaks = seq(0, y_lim, by = 1),
      labels = sprintf("%.2f", seq(0, y_lim, by = 1)),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title = prefix,
      y = "Copy Number"
    ) +
    ggplot2::theme(
      legend.position = "none",
      panel.spacing = grid::unit(0, "lines"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 13),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 13),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = p_margin
    ) +
    y_axis_theme
  
  q <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = 0.5,
      color = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = c(0, 0.25, 0.75, 1),
      color = "grey",
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      data = df_ai,
      ggplot2::aes(
        x = bin_start,
        xend = bin_start,
        y = smoothed_ai - 0.05,
        yend = smoothed_ai + 0.05,
        color = seqnames
      ),
      alpha = 0.01
    ) +
    ggplot2::scale_color_manual(values = color, drop = FALSE) +
    ggplot2::geom_segment(
      data = call_seg_maf,
      ggplot2::aes(
        x = Start,
        xend = End,
        y = baf_lower,
        yend = baf_lower
      ),
      linewidth = 1,
      color = "cyan3"
    ) +
    ggplot2::geom_segment(
      data = call_seg_maf,
      ggplot2::aes(
        x = Start,
        xend = End,
        y = baf_upper,
        yend = baf_upper
      ),
      linewidth = 1,
      color = "cyan3"
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(seqnames),
      scales = "free_x",
      space = "free_x",
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-0.05, 1.05),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = sprintf("%.2f", c(0, 0.25, 0.5, 0.75, 1)),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.05, 1.05), clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = "BAF") +
    ggplot2::theme(
      legend.position = "none",
      panel.spacing = grid::unit(0, "lines"),
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = q_margin
    ) +
    y_axis_theme
  
  quality <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = call_seg,
      ggplot2::aes(
        x = Start,
        xend = End,
        y = 0.5,
        yend = 0.5,
        color = FILTER
      ),
      linewidth = 4
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "PASS" = "grey20",
        "FAILED" = "grey50"
      ),
      na.value = "grey20"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0.5),
      labels = c("QC"),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(seqnames),
      scales = "free_x",
      space = "free_x",
      drop = FALSE
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = "") +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_blank(),
      panel.spacing = grid::unit(0, "lines"),
      strip.background = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 12, color = "black"),
      axis.text.y.left = ggplot2::element_text(size = 12, color = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y.left = ggplot2::element_line(color = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = qc_margin
    )
  
  gridExtra::arrangeGrob(
    p,
    q,
    quality,
    ncol = 1,
    heights = c(5, 3, 0.7)
  )
}

BuildWgsDragenCnvPlotFromCachedData <- function(
  call_seg_static,
  cov_bins_raw,
  final_ai,
  gender = "female",
  diploid_coverage,
  purity,
  prefix = "Sample"
) {
  call_seg <- call_seg_static %>%
    dplyr::mutate(
      EstimatedTumorPurity = purity,
      DiploidCoverage = diploid_coverage
    )

  smooth_cov <- NormalizeCoverageBinsByUserModel(
    cov_bins_raw = cov_bins_raw,
    gender = gender,
    purity = purity,
    diploid_coverage = diploid_coverage
  )

  PlotWgsCnvGrid(
    df_cov = smooth_cov,
    df_ai = final_ai,
    call_seg = call_seg,
    gender = gender,
    prefix = prefix
  )
}

BuildWgsDragenCnvPlot <- function(
  anno_file,
  coverage_file,
  ballele_file,
  gender = "female",
  diploid_coverage,
  purity,
  cov_binsize = 200000,
  ai_binsize = 100000,
  prefix = "Sample"
) {
  anno_raw <- ReadDragenAnnotationRaw(anno_file)
  call_seg_static <- PreprocessDragenAnnotation(
    df = anno_raw,
    gender = gender,
    user_purity = purity,
    user_diploid_coverage = diploid_coverage
  )

  cov_raw <- ReadDragenCoverageRaw(coverage_file, gender = gender)
  cov_bins_raw <- BinRawCoverageFast(cov_raw, cov_binsize = cov_binsize)

  ai_rounded <- ReadPrepareRoundDragenBAllele(ballele_file, gender = gender)
  clean_ai <- CleanHomAlt(clean_ai = ai_rounded, call_seg = call_seg_static)
  final_ai <- SmoothAi(df = clean_ai, ai_binsize = ai_binsize, gender = gender)

  BuildWgsDragenCnvPlotFromCachedData(
    call_seg_static = call_seg_static,
    cov_bins_raw = cov_bins_raw,
    final_ai = final_ai,
    gender = gender,
    diploid_coverage = diploid_coverage,
    purity = purity,
    prefix = prefix
  )
}

ReadDragenModelGrid <- function(model_file) {
  models <- data.table::fread(model_file, data.table = FALSE)

  if (ncol(models) < 3) {
    stop("Model grid table must have at least 3 columns: Purity, Coverage, logL.")
  }

  names(models)[1:3] <- c("Purity", "Coverage", "logL")

  models <- models %>%
    dplyr::mutate(
      Purity = as.numeric(Purity),
      Coverage = as.numeric(Coverage),
      logL = as.numeric(logL)
    ) %>%
    dplyr::filter(
      !is.na(Purity),
      !is.na(Coverage),
      !is.na(logL),
      is.finite(Purity),
      is.finite(Coverage),
      is.finite(logL)
    )

  if (nrow(models) == 0) {
    stop("No valid rows were found in the model grid table.")
  }

  models
}

BuildWgsDragenModelPlotFromGrid <- function(
    models,
    user_purity,
    user_coverage,
    sample_name = "Sample"
) {
  max_model <- models[which.max(models$logL), , drop = FALSE]
  
  coverage_breaks <- pretty(
    range(models$Coverage, user_coverage, na.rm = TRUE),
    n = 6
  )
  
  ggplot2::ggplot(
    models,
    ggplot2::aes(
      x = Purity,
      y = Coverage,
      fill = logL
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_point(
      data = max_model,
      ggplot2::aes(
        x = Purity,
        y = Coverage
      ),
      color = "red",
      size = 6,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = data.frame(
        Purity = user_purity,
        Coverage = user_coverage
      ),
      ggplot2::aes(
        x = Purity,
        y = Coverage
      ),
      color = "blue",
      size = 4,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_gradient(
      low = "navy",
      high = "yellow"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, by = 0.1),
      labels = sprintf("%.2f", seq(0, 1, by = 0.1)),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      breaks = coverage_breaks,
      labels = sprintf("%.0f", coverage_breaks),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::labs(
      title = paste0(sample_name, " model plot"),
      x = "Purity",
      y = "Diploid Coverage",
      fill = "logL"
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 18,
        face = "bold",
        color = "black",
        margin = ggplot2::margin(b = 8)
      ),
      axis.title.x = ggplot2::element_text(
        size = 16,
        color = "black",
        margin = ggplot2::margin(t = 8)
      ),
      axis.title.y = ggplot2::element_text(
        size = 16,
        color = "black",
        margin = ggplot2::margin(r = 8)
      ),
      axis.text.x = ggplot2::element_text(
        size = 13,
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        size = 13,
        color = "black"
      ),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.line = ggplot2::element_line(color = "black"),
      legend.title = ggplot2::element_text(size = 13, color = "black"),
      legend.text = ggplot2::element_text(size = 12, color = "black"),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(
        t = 18,
        r = 20,
        b = 18,
        l = 20,
        unit = "pt"
      )
    )
}

BuildWgsDragenModelPlot <- function(model_file, user_purity, user_coverage, sample_name = "Sample") {
  models <- ReadDragenModelGrid(model_file)
  BuildWgsDragenModelPlotFromGrid(
    models = models,
    user_purity = user_purity,
    user_coverage = user_coverage,
    sample_name = sample_name
  )
}
