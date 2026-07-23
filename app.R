# app.R

library(shiny)
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

local_only_check <- function(session) {
  remote_addr <- session$request$REMOTE_ADDR
  forwarded_for <- session$request$HTTP_X_FORWARDED_FOR
  
  allowed <- is.null(remote_addr) ||
    remote_addr %in% c("127.0.0.1", "::1", "localhost")
  
  no_proxy <- is.null(forwarded_for) || forwarded_for == ""
  
  if (!allowed || !no_proxy) {
    stop(
      "This app is configured for local use only. ",
      "Please run from localhost / 127.0.0.1."
    )
  }
}

BuildInputPaths <- function(base_path, sample_id) {
  base_path <- clean_path(base_path)
  sample_id <- trimws(sample_id)
  
  if (base_path == "") {
    stop("Please enter the input folder path.")
  }
  
  if (sample_id == "") {
    stop("Please enter the sample ID.")
  }
  
  list(
    anno = file.path(base_path, paste0(sample_id, "_Annotation_full.txt")),
    coverage = file.path(base_path, paste0(sample_id, ".tumor.target.counts.gc-corrected.gz")),
    ballele = file.path(base_path, paste0(sample_id, ".tumor.ballele.counts.gz")),
    model = file.path(base_path, paste0(sample_id, ".cnv.purity.coverage.models.tsv"))
  )
}

validate_local_file <- function(path, label = "file", required = TRUE) {
  path <- clean_path(path)
  
  if (path == "") {
    if (required) {
      stop("Please provide a path for ", label, ".")
    } else {
      return(NULL)
    }
  }
  
  if (!file.exists(path)) {
    stop("The selected ", label, " does not exist:\n", path)
  }
  
  normalizePath(path, mustWork = TRUE)
}

GetParameterValue <- function(text_value, slider_value, label) {
  if (!is.null(text_value) &&
      length(text_value) > 0 &&
      !is.na(text_value)) {
    value <- as.numeric(text_value)
  } else {
    value <- as.numeric(slider_value)
  }
  
  if (!is.finite(value)) {
    stop("Invalid ", label, " value.")
  }
  
  if (label == "purity" && (value < 0 || value > 1)) {
    stop("Purity must be between 0 and 1.")
  }
  
  if (label == "diploid coverage" && value <= 0) {
    stop("Diploid coverage must be greater than 0.")
  }
  
  value
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

        pre {
          font-size: 11px;
          white-space: pre-wrap;
          word-break: break-all;
        }
        "
      )
    )
  ),
  
  titlePanel("Interactive DRAGEN WGS CNV Plot"),
  
  tags$div(
    style = "padding: 8px; background-color: #f7f7f7; border-left: 4px solid #555; margin-bottom: 10px;",
    tags$b("Local-only mode: "),
    "Data are loaded only when you click Load data. ",
    "Plots are regenerated only when you click Generate plot."
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      textInput(
        inputId = "base_path",
        label = "Input folder path",
        value = "",
        placeholder = "/path/to/folder"
      ),
      
      textInput(
        inputId = "sample_id",
        label = "Sample ID",
        value = "",
        placeholder = "SampleID"
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
        tags$b("Diploid coverage"),
        tags$div(
          class = "small-note",
          "Typed value overrides slider. Leave blank to use slider."
        ),
        numericInput(
          inputId = "diploid_coverage_value",
          label = "Diploid coverage value",
          value = NA_real_,
          min = 1,
          max = 1000,
          step = 1
        ),
        sliderInput(
          inputId = "diploid_coverage_slider",
          label = "Diploid coverage slider",
          min = 0,
          max = 1000,
          value = 100,
          step = 1,
          round = 0,
          ticks = TRUE,
          width = "100%"
        )
      ),
      
      tags$div(
        class = "file-block",
        tags$b("Purity"),
        tags$div(
          class = "small-note",
          "Typed value overrides slider. Leave blank to use slider."
        ),
        numericInput(
          inputId = "purity_value",
          label = "Purity value",
          value = NA_real_,
          min = 0,
          max = 1,
          step = 0.01
        ),
        sliderInput(
          inputId = "purity_slider",
          label = "Purity slider",
          min = 0,
          max = 1,
          value = 0.90,
          step = 0.01,
          round = 2,
          ticks = TRUE,
          width = "100%"
        )
      ),
      
      actionButton(
        inputId = "load_data",
        label = "Load data",
        class = "btn-success"
      ),
      
      tags$br(),
      tags$br(),
      
      actionButton(
        inputId = "generate_plot",
        label = "Generate plot",
        class = "btn-primary"
      ),
      
      tags$br(),
      tags$br(),
      
      downloadButton("download_cnv_png", "Save CNV plot PNG"),
      downloadButton("download_model_png", "Save model plot PNG"),
      
      tags$hr(),
      
      tags$b("Auto-generated input files"),
      verbatimTextOutput("auto_paths"),
      
      tags$hr(),
      
      tags$b("Data status"),
      verbatimTextOutput("data_status")
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
  session$allowReconnect(TRUE)
  
  loaded_data <- reactiveVal(NULL)
  loaded_paths <- reactiveVal(NULL)
  loaded_status <- reactiveVal("No data loaded yet.")
  
  current_paths <- reactive({
    tryCatch(
      BuildInputPaths(
        base_path = input$base_path,
        sample_id = input$sample_id
      ),
      error = function(e) {
        NULL
      }
    )
  })
  
  output$auto_paths <- renderPrint({
    paths <- current_paths()
    
    if (is.null(paths)) {
      cat("Enter input folder path and sample ID.")
      return()
    }
    
    cat("Annotation:\n", paths$anno, "\n\n", sep = "")
    cat("Coverage:\n", paths$coverage, "\n\n", sep = "")
    cat("B-allele:\n", paths$ballele, "\n\n", sep = "")
    cat("Model grid:\n", paths$model, "\n", sep = "")
  })
  
  output$data_status <- renderPrint({
    cat(loaded_status())
  })
  
  observeEvent(input$load_data, {
    loaded_status("Loading data... please wait.")
    loaded_data(NULL)
    
    tryCatch({
      paths <- BuildInputPaths(
        base_path = input$base_path,
        sample_id = input$sample_id
      )
      
      anno_final <- validate_local_file(paths$anno, "annotation file")
      coverage_final <- validate_local_file(paths$coverage, "coverage file")
      ballele_final <- validate_local_file(paths$ballele, "B-allele file")
      
      model_final <- NULL
      if (file.exists(paths$model)) {
        model_final <- normalizePath(paths$model, mustWork = TRUE)
      }
      
      purity_for_loading <- GetParameterValue(
        text_value = input$purity_value,
        slider_value = input$purity_slider,
        label = "purity"
      )
      
      diploid_coverage_for_loading <- GetParameterValue(
        text_value = input$diploid_coverage_value,
        slider_value = input$diploid_coverage_slider,
        label = "diploid coverage"
      )
      
      withProgress(message = "Reading files once...", value = 0, {
        incProgress(0.10, detail = "Reading annotation")
        
        anno_raw <- ReadDragenAnnotationRaw(anno_final)
        
        call_seg_static <- PreprocessDragenAnnotation(
          df = anno_raw,
          gender = input$gender,
          user_purity = purity_for_loading,
          user_diploid_coverage = diploid_coverage_for_loading
        )
        
        incProgress(0.30, detail = "Reading and pre-binning coverage")
        
        cov_raw <- ReadDragenCoverageRaw(
          coverage_file = coverage_final,
          gender = input$gender
        )
        
        cov_bins_raw <- BinRawCoverageFast(
          cov = cov_raw,
          cov_binsize = 200000
        )
        
        incProgress(0.50, detail = "Reading BAF")
        
        ai_rounded <- ReadPrepareRoundDragenBAllele(
          ballele_file = ballele_final,
          gender = input$gender
        )
        
        incProgress(0.70, detail = "Cleaning and smoothing BAF")
        
        clean_ai <- CleanHomAlt(
          clean_ai = ai_rounded,
          call_seg = call_seg_static
        )
        
        final_ai <- SmoothAi(
          df = clean_ai,
          ai_binsize = 100000,
          gender = input$gender
        )
        
        incProgress(0.90, detail = "Reading model grid, if available")
        
        models <- NULL
        if (!is.null(model_final)) {
          models <- ReadDragenModelGrid(model_final)
        }
        
        incProgress(1, detail = "Done")
        
        loaded_data(
          list(
            anno_raw = anno_raw,
            call_seg_static = call_seg_static,
            cov_bins_raw = cov_bins_raw,
            final_ai = final_ai,
            models = models,
            gender = input$gender,
            loaded_at = Sys.time()
          )
        )
        
        loaded_paths(paths)
      })
      
      loaded_status(
        paste0(
          "Data loaded successfully.\n\n",
          "Sample: ", input$sample_id, "\n",
          "Gender: ", input$gender, "\n",
          "Loaded at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
          "Annotation:\n", paths$anno, "\n\n",
          "Coverage:\n", paths$coverage, "\n\n",
          "B-allele:\n", paths$ballele, "\n\n",
          "Model grid:\n",
          ifelse(file.exists(paths$model), paths$model, "Not found / skipped")
        )
      )
      
      showNotification(
        "Data loaded successfully.",
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      loaded_data(NULL)
      
      loaded_status(
        paste0(
          "Data loading failed.\n\n",
          conditionMessage(e)
        )
      )
      
      showNotification(
        paste0("Data loading failed: ", conditionMessage(e)),
        type = "error",
        duration = 10
      )
    })
  })
  
  cnv_plot_obj <- eventReactive(input$generate_plot, {
    data <- loaded_data()
    
    if (is.null(data)) {
      stop("Please click Load data before Generate plot.")
    }
    
    sample_name <- trimws(input$sample_id)
    if (sample_name == "") {
      sample_name <- "Sample"
    }
    
    purity <- GetParameterValue(
      text_value = input$purity_value,
      slider_value = input$purity_slider,
      label = "purity"
    )
    
    diploid_coverage <- GetParameterValue(
      text_value = input$diploid_coverage_value,
      slider_value = input$diploid_coverage_slider,
      label = "diploid coverage"
    )
    
    BuildWgsDragenCnvPlotFromCachedData(
      call_seg_static = data$call_seg_static,
      cov_bins_raw = data$cov_bins_raw,
      final_ai = data$final_ai,
      gender = data$gender,
      diploid_coverage = diploid_coverage,
      purity = purity,
      prefix = paste0(
        sample_name,
        " | purity = ",
        formatC(purity, format = "f", digits = 2),
        " | diploid coverage = ",
        formatC(diploid_coverage, format = "f", digits = 0)
      )
    )
  })
  
  model_plot_obj <- eventReactive(input$generate_plot, {
    data <- loaded_data()
    
    if (is.null(data)) {
      return(NULL)
    }
    
    if (is.null(data$models)) {
      return(NULL)
    }
    
    sample_name <- trimws(input$sample_id)
    if (sample_name == "") {
      sample_name <- "Sample"
    }
    
    purity <- GetParameterValue(
      text_value = input$purity_value,
      slider_value = input$purity_slider,
      label = "purity"
    )
    
    diploid_coverage <- GetParameterValue(
      text_value = input$diploid_coverage_value,
      slider_value = input$diploid_coverage_slider,
      label = "diploid coverage"
    )
    
    BuildWgsDragenModelPlotFromGrid(
      models = data$models,
      user_purity = purity,
      user_coverage = diploid_coverage,
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
        labels = "Optional model plot: model grid was not loaded.",
        cex = 0.9
      )
    } else {
      print(p)
    }
  })
  
  output$download_cnv_png <- downloadHandler(
    filename = function() {
      sample_name <- trimws(input$sample_id)
      if (sample_name == "") {
        sample_name <- "Sample"
      }
      
      sample_name <- gsub("[^A-Za-z0-9_\\-]", "_", sample_name)
      paste0(sample_name, "_WGS_CNV_plot.png")
    },
    content = function(file) {
      req(cnv_plot_obj())
      
      png(
        file,
        width = 11.8,
        height = 6.4,
        units = "in",
        res = 300
      )
      
      grid::grid.newpage()
      grid::grid.draw(cnv_plot_obj())
      
      dev.off()
    }
  )
  
  output$download_model_png <- downloadHandler(
    filename = function() {
      sample_name <- trimws(input$sample_id)
      if (sample_name == "") {
        sample_name <- "Sample"
      }
      
      sample_name <- gsub("[^A-Za-z0-9_\\-]", "_", sample_name)
      paste0(sample_name, "_WGS_model_plot.png")
    },
    content = function(file) {
      p <- model_plot_obj()
      
      if (is.null(p)) {
        stop("No model plot is available. Please provide a model grid table.")
      }
      
      png(
        file,
        width = 6.5,
        height = 4.5,
        units = "in",
        res = 300
      )
      
      print(p)
      
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)