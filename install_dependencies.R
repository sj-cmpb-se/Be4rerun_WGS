# install_dependencies.R
# Run once:
# source("install_dependencies.R")

cran_packages <- c(
  "shiny",
  "shinyFiles",
  "data.table",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "scales",
  "gridExtra",
  "R.utils"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "GenomicRanges",
  "IRanges",
  "S4Vectors"
)

missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("All required packages are installed.")
