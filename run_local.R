# run_local.R
# Run from the project root:
# source("run_local.R")

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("The shiny package is not installed. Please run source('install_dependencies.R') first.")
}

shiny::runApp(
  appDir = ".",
  host = "127.0.0.1",
  port = 7879,
  launch.browser = TRUE
)
