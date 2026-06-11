# app.R
library(shiny)
library(shinyFiles)
library(ggplot2)
library(dplyr)
library(data.table)
library(gridExtra)
library(grid)

source("R/wgs_dragen_functions.R")

clean_path <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return("")
  }
  x <- trimws(x)
  x <- gsub('^["\']|["\']$', "", x)
  path.expand(x)
}

validate_local_file <- function(path, label = "file", required = TRUE) {
  path <- clean_path(path)

  if (path == "") {
    if (required) {
      stop("Please provide a path for ", label, ".")
    }
    return(NULL)
  }

  if (!file.exists(path)) {
    stop("The selected ", label, " does not exist:\n", path)
  }

  normalizePath(path, mustWork = TRUE)
}

local_only_check <- function(session) {
  remote_addr <- session$request$REMOTE_ADDR
  forwarded_for <- session$request$HTTP_X_FORWARDED_FOR

  allowed <- is.null(remote_addr) || remote_addr %in% c("127.0.0.1", "::1", "localhost")
  no_proxy <- is.null(forwarded_for) || forwarded_for == ""

  if (!allowed || !no_proxy) {
    stop("This app is configured for local use only. Please run from localhost / 127.0.0.1.")
  }
}

volumes <- c(
  Home = normalizePath("~"),
  Working_Directory = normalizePath(getwd())
)

if (dir.exists("/Volumes")) {
  volumes <- c(volumes, Volumes = normalizePath("/Volumes", mustWork = TRUE))
}

if (dir.exists("/")) {
  volumes <- c(volumes, Root = normalizePath("/", mustWork = TRUE))
}

ui <- fluidPage(
  tags$head(
    tags$style(
      HTML(
        "
        .container-fluid {
          max-width: 1450px;
          margin-left: auto;
          margin-right: auto;
        }
        .main-plot-box {
          max-width: 1180px;
          margin-left: auto;
          margin-right: auto;
        }
        .model-plot-box {
          max-width: 650px;
          margin-left: auto;
          margin-right: auto;
        }
        .file-block {
          margin-bottom: 14px;
          padding-bottom: 10px;
          border-bottom: 1px solid #eeeeee;
        }
        .small-note {
          color: #666666;
          font-size: 12px;
          margin-top: -5px;
          margin-bottom: 8px;
        }
        "
      )
    )
  ),

  titlePanel("Interactive DRAGEN WGS CNV Plot"),

  fluidRow(
    column(
      width = 12,
      sliderInput(
        inputId = "diploid_coverage",
        label = "Diploid coverage",
        min = 0,
        max = 1000,
        value = 100,
        step = 1,
        round = 0,
        ticks = TRUE,
        width = "100%"
      ),
      sliderInput(
        inputId = "purity",
        label = "Purity",
        min = 0,
        max = 1,
        value = 0.90,
        step = 0.01,
        round = 2,
        ticks = TRUE,
        width = "100%"
      )
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      textInput(
        inputId = "sample_name",
        label = "Sample name",
        value = "Sample"
      ),

      selectInput(
        inputId = "gender",
        label = "Gender",
        choices = c("female", "male"),
        selected = "female"
      ),

      tags$hr(),

      tags$div(
        class = "file-block",
        tags$b("Annotation file"),
        tags$div(class = "small-note", "*_Annotation_full.txt"),
        textInput(
          inputId = "anno_path_manual",
          label = NULL,
          value = "",
          placeholder = "/path/to/sample_Annotation_full.txt",
          width = "100%"
        ),
        shinyFilesButton(
          id = "anno_file",
          label = "Browse",
          title = "Select annotation file",
          multiple = FALSE
        )
      ),

      tags$div(
        class = "file-block",
        tags$b("Coverage file"),
        tags$div(class = "small-note", "*.tumor.target.counts.gz or *.tumor.target.counts.gc-corrected.gz"),
        textInput(
          inputId = "coverage_path_manual",
          label = NULL,
          value = "",
          placeholder = "/path/to/sample.tumor.target.counts.gc-corrected.gz",
          width = "100%"
        ),
        shinyFilesButton(
          id = "coverage_file",
          label = "Browse",
          title = "Select coverage file",
          multiple = FALSE
        )
      ),

      tags$div(
        class = "file-block",
        tags$b("B-allele file"),
        tags$div(class = "small-note", "*.tumor.ballele.counts.gz"),
        textInput(
          inputId = "ballele_path_manual",
          label = NULL,
          value = "",
          placeholder = "/path/to/sample.tumor.ballele.counts.gz",
          width = "100%"
        ),
        shinyFilesButton(
          id = "ballele_file",
          label = "Browse",
          title = "Select B-allele file",
          multiple = FALSE
        )
      ),

      tags$div(
        class = "file-block",
        tags$b("Optional model grid table"),
        tags$div(class = "small-note", "*.cnv.purity.coverage.models.tsv"),
        textInput(
          inputId = "model_path_manual",
          label = NULL,
          value = "",
          placeholder = "/path/to/sample.cnv.purity.coverage.models.tsv",
          width = "100%"
        ),
        shinyFilesButton(
          id = "model_file",
          label = "Browse",
          title = "Select model grid table",
          multiple = FALSE
        )
      ),

      actionButton(
        inputId = "load_files",
        label = "Load files",
        class = "btn-primary"
      ),

      tags$div(
        class = "small-note",
        "No whitelist filtering is performed. After loading, slider changes reuse cached data."
      ),

      tags$br(),

      downloadButton("download_cnv_png", "Save CNV plot PNG"),
      downloadButton("download_model_png", "Save model plot PNG")
    ),

    mainPanel(
      width = 9,
      div(
        class = "main-plot-box",
        plotOutput("cnv_plot", height = "640px", width = "100%")
      ),
      tags$br(),
      div(
        class = "model-plot-box",
        plotOutput("model_plot", height = "450px", width = "100%")
      )
    )
  )
)

server <- function(input, output, session) {
  local_only_check(session)

  shinyFileChoose(input, "anno_file", roots = volumes, session = session)
  shinyFileChoose(input, "coverage_file", roots = volumes, session = session)
  shinyFileChoose(input, "ballele_file", roots = volumes, session = session)
  shinyFileChoose(input, "model_file", roots = volumes, session = session)

  selected_path_from_browser <- function(file_input) {
    if (is.null(file_input) || length(file_input) == 0) {
      return(NULL)
    }
    parsed <- shinyFiles::parseFilePaths(volumes, file_input)
    if (nrow(parsed) == 0) {
      return(NULL)
    }
    normalizePath(parsed$datapath[1], mustWork = TRUE)
  }

  observeEvent(input$anno_file, {
    path <- selected_path_from_browser(input$anno_file)
    if (!is.null(path)) updateTextInput(session, "anno_path_manual", value = path)
  })

  observeEvent(input$coverage_file, {
    path <- selected_path_from_browser(input$coverage_file)
    if (!is.null(path)) updateTextInput(session, "coverage_path_manual", value = path)
  })

  observeEvent(input$ballele_file, {
    path <- selected_path_from_browser(input$ballele_file)
    if (!is.null(path)) updateTextInput(session, "ballele_path_manual", value = path)
  })

  observeEvent(input$model_file, {
    path <- selected_path_from_browser(input$model_file)
    if (!is.null(path)) updateTextInput(session, "model_path_manual", value = path)
  })

  loaded_data <- eventReactive(input$load_files, {
    anno_final <- validate_local_file(input$anno_path_manual, "annotation file")
    coverage_final <- validate_local_file(input$coverage_path_manual, "coverage file")
    ballele_final <- validate_local_file(input$ballele_path_manual, "B-allele file")
    model_final <- validate_local_file(
      input$model_path_manual,
      "model grid table",
      required = FALSE
    )

    withProgress(message = "Reading files once...", value = 0, {
      incProgress(0.15, detail = "Reading annotation")
      anno_raw <- ReadDragenAnnotationRaw(anno_final)

      call_seg_static <- PreprocessDragenAnnotation(
        df = anno_raw,
        gender = input$gender,
        user_purity = input$purity,
        user_diploid_coverage = input$diploid_coverage
      )

      incProgress(0.35, detail = "Reading and pre-binning coverage")
      cov_raw <- ReadDragenCoverageRaw(coverage_final, gender = input$gender)
      cov_bins_raw <- BinRawCoverageFast(cov_raw, cov_binsize = 200000)

      incProgress(0.55, detail = "Reading BAF")
      ai_rounded <- ReadPrepareRoundDragenBAllele(ballele_final, gender = input$gender)

      incProgress(0.72, detail = "Smoothing BAF")
      clean_ai <- CleanHomAlt(clean_ai = ai_rounded, call_seg = call_seg_static)
      final_ai <- SmoothAi(df = clean_ai, ai_binsize = 100000, gender = input$gender)

      incProgress(0.90, detail = "Reading model grid, if provided")
      models <- NULL
      if (!is.null(model_final)) {
        models <- ReadDragenModelGrid(model_final)
      }

      incProgress(1, detail = "Done")

      list(
        anno_raw = anno_raw,
        call_seg_static = call_seg_static,
        cov_bins_raw = cov_bins_raw,
        final_ai = final_ai,
        models = models,
        gender = input$gender,
        loaded_at = Sys.time()
      )
    })
  })

  cnv_plot_obj <- reactive({
    data <- loaded_data()
    req(data)

    sample_name <- input$sample_name
    if (is.null(sample_name) || sample_name == "") {
      sample_name <- "Sample"
    }

    BuildWgsDragenCnvPlotFromCachedData(
      call_seg_static = data$call_seg_static,
      cov_bins_raw = data$cov_bins_raw,
      final_ai = data$final_ai,
      gender = data$gender,
      diploid_coverage = input$diploid_coverage,
      purity = input$purity,
      prefix = paste0(
        sample_name,
        " | purity = ", formatC(input$purity, format = "f", digits = 2),
        " | diploid coverage = ", formatC(input$diploid_coverage, format = "f", digits = 0)
      )
    )
  })

  model_plot_obj <- reactive({
    data <- loaded_data()
    req(data)

    if (is.null(data$models)) {
      return(NULL)
    }

    sample_name <- input$sample_name
    if (is.null(sample_name) || sample_name == "") {
      sample_name <- "Sample"
    }

    BuildWgsDragenModelPlotFromGrid(
      models = data$models,
      user_purity = input$purity,
      user_coverage = input$diploid_coverage,
      sample_name = sample_name
    )
  })

  output$cnv_plot <- renderPlot({
    req(cnv_plot_obj())
    grid::grid.newpage()
    grid::grid.draw(cnv_plot_obj())
  })

  output$model_plot <- renderPlot({
    p <- model_plot_obj()
    if (is.null(p)) {
      plot.new()
      text(
        x = 0.5,
        y = 0.5,
        labels = "Optional model plot: provide a model grid table, then click Load files.",
        cex = 0.9
      )
    } else {
      print(p)
    }
  })

  output$download_cnv_png <- downloadHandler(
    filename = function() {
      sample_name <- input$sample_name
      if (is.null(sample_name) || sample_name == "") sample_name <- "Sample"
      sample_name <- gsub("[^A-Za-z0-9_\\-]", "_", sample_name)
      paste0(sample_name, "_WGS_CNV_plot.png")
    },
    content = function(file) {
      req(cnv_plot_obj())
      png(file, width = 11.8, height = 6.4, units = "in", res = 300)
      grid::grid.newpage()
      grid::grid.draw(cnv_plot_obj())
      dev.off()
    }
  )

  output$download_model_png <- downloadHandler(
    filename = function() {
      sample_name <- input$sample_name
      if (is.null(sample_name) || sample_name == "") sample_name <- "Sample"
      sample_name <- gsub("[^A-Za-z0-9_\\-]", "_", sample_name)
      paste0(sample_name, "_WGS_model_plot.png")
    },
    content = function(file) {
      p <- model_plot_obj()
      if (is.null(p)) stop("No model plot is available. Please provide a model grid table.")
      png(file, width = 6.5, height = 4.5, units = "in", res = 300)
      print(p)
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)
