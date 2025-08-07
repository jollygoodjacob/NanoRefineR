library(shiny)
library(terra)
library(signal)
library(doParallel)
library(plotly)

# === UI ===
ui <- fluidPage(
  titlePanel("NanoRefineR: A toolbox for dealing with noisy Nano-Hyperspec imagery"),
  
  tags$head(
    tags$style(HTML(".boxed-section {
      border: 1px solid #ccc;
      border-radius: 5px;
      padding: 10px;
      margin-bottom: 10px;
      background-color: #f9f9f9;
    }"))
  ),
  
  fluidRow(
    column(3,
           div(class = "boxed-section",
               h4("Processing Settings"),
               textInput("folder", "Path to Folder with Hyperspectral Cubes", value = ""),
               textInput("file_pattern", "File Pattern (Regex)", value = "flightline_[0-9]+$"),
               textInput("output_suffix", "Suffix for Output Files", value = "_processed"),
               
               checkboxInput("do_spatial", "Enable Spatial Downscaling", value = TRUE),
               textInput("target_crs", "Target CRS (e.g., EPSG:32610)", value = "EPSG:32610"),
               numericInput("target_res", "Target Resolution (meters)", value = 0.25, min = 0.01),
               
               checkboxInput("do_binning", "Enable Spectral Binning", value = TRUE),
               numericInput("bin_width", "Spectral Bin Width (nm)", value = 5, min = 1),
               
               checkboxInput("do_sg", "Enable Savitzky-Golay Filtering", value = TRUE),
               numericInput("sgolay_p", "SG Polynomial Order (p)", value = 3, min = 1),
               numericInput("sgolay_n", "SG Filter Length (n, odd)", value = 5, min = 3),
               
               numericInput("num_cores", "Number of Cores to Use", value = max(1, parallel::detectCores() - 2), min = 1),
               actionButton("run", "Run Processing"),
               verbatimTextOutput("status")
           )
    ),
    
    column(6,
           div(class = "boxed-section",
               h4("RGB Plots"),
               fluidRow(
                 column(6, plotOutput("rgb_input", click = "click_input", height = "435px")),
                 column(6, plotOutput("rgb_output", click = "click_output", height = "435px"))
               )
           ),
           div(class = "boxed-section",
               h4("Spectral Profiles"),
               fluidRow(
                 column(6, plotlyOutput("spectral_input", height = "300px")),
                 column(6, plotlyOutput("spectral_output", height = "300px"))
               )
           )
    ),
    
    column(3,
           div(class = "boxed-section",
               h4("Visualization Controls"),
               selectInput("preview_file", "Preview Processed File", choices = NULL),
               numericInput("r_band", "Red Wavelength (nm)", value = 670),
               numericInput("g_band", "Green Wavelength (nm)", value = 550),
               numericInput("b_band", "Blue Wavelength (nm)", value = 470)
           )
    )
  )
)


# === Server ===
server <- function(input, output, session) {
  values <- reactiveValues(input_rast = NULL, output_rast = NULL, wavelengths_in = NULL, wavelengths_out = NULL)
  
  observeEvent(input$run, {
    req(input$folder)
    paths <- list.files(input$folder, pattern = input$file_pattern, full.names = TRUE)
    
    if (length(paths) == 0) {
      output$status <- renderText("No matching files found.")
      return(NULL)
    }
    
    # Extract Shiny inputs to pass into foreach
    suffix <- isolate(input$output_suffix)
    do_spatial <- isolate(input$do_spatial)
    target_crs <- isolate(input$target_crs)
    target_res <- isolate(input$target_res)
    do_binning <- isolate(input$do_binning)
    bin_width <- isolate(input$bin_width)
    do_sg <- isolate(input$do_sg)
    sgolay_p <- isolate(input$sgolay_p)
    sgolay_n <- isolate(input$sgolay_n)
    
    cl <- makeCluster(input$num_cores)
    registerDoParallel(cl)
    
    output$status <- renderText("Processing started...")
    
    foreach(i = seq_along(paths), .packages = c("terra", "signal")) %dopar% {
      img <- rast(paths[[i]])
      
      if (do_spatial) {
        img <- project(img, target_crs, res = target_res)
      }
      
      if (do_binning) {
        wl <- names(img)
        options(digits = 6)
        wl <- round(as.double(gsub(" .*", "", wl)), 0)
        
        spectral_binning <- function(image, wavelengths, bin_width) {
          bins <- seq(min(wavelengths), max(wavelengths), by = bin_width)
          binned <- rast()
          band_centers <- c()
          for (j in 1:(length(bins) - 1)) {
            idx <- which(wavelengths >= bins[j] & wavelengths < bins[j + 1])
            if (length(idx) > 0) {
              mean_band <- mean(image[[idx]])
              binned <- c(binned, mean_band)
              band_centers <- c(band_centers, round((bins[j] + bins[j + 1]) / 2, 1))
            }
          }
          names(binned) <- as.character(band_centers)
          return(binned)
        }
        img <- spectral_binning(img, wl, bin_width)
      }
      
      if (do_sg) {
        apply_sg_filter <- function(x) {
          if (all(is.na(x))) return(rep(NA, length(x)))
          sgolayfilt(x, p = sgolay_p, n = sgolay_n)
        }
        old_names <- names(img)
        img <- app(img, apply_sg_filter)
        names(img) <- old_names  # ✅ Restore band center names with 1 decimal
      }
      
      base_name <- tools::file_path_sans_ext(basename(paths[[i]]))
      out_path <- file.path(dirname(paths[[i]]), paste0(base_name, suffix, ".tif"))
      writeRaster(img, out_path, overwrite = TRUE)
    }
    
    stopCluster(cl)
    
    # Update preview file list
    processed <- list.files(input$folder, pattern = paste0(suffix, ".tif$"), full.names = TRUE)
    updateSelectInput(session, "preview_file", choices = processed, selected = processed[1])
    output$status <- renderText("Processing complete! Select a file to preview.")
  })
  
  observeEvent(input$preview_file, {
    req(input$preview_file)
    out_rast <- rast(input$preview_file)
    
    base_name <- gsub(input$output_suffix, "", tools::file_path_sans_ext(basename(input$preview_file)))
    orig_candidates <- list.files(input$folder, full.names = TRUE)
    orig_file <- orig_candidates[grepl(base_name, orig_candidates) & !grepl(input$output_suffix, orig_candidates)]
    if (length(orig_file) > 0) {
      in_rast <- rast(orig_file[1])
      values$input_rast <- in_rast
      values$output_rast <- out_rast
      values$wavelengths_in <- round(as.double(gsub(" .*", "", names(in_rast))), 0)
      values$wavelengths_out <- round(as.double(gsub(" .*", "", names(out_rast))), 1)
      output$status <- renderText("Input and output rasters loaded. Click a pixel to view spectra.")
    } else {
      output$status <- renderText("Could not find matching input file.")
    }
  })
  
  get_closest_band <- function(wavelengths, target) {
    which.min(abs(wavelengths - target))
  }
  
  render_rgb <- function(rast_obj, wavelengths, input) {
    r_idx <- get_closest_band(wavelengths, input$r_band)
    g_idx <- get_closest_band(wavelengths, input$g_band)
    b_idx <- get_closest_band(wavelengths, input$b_band)
    terra::plotRGB(rast_obj[[c(r_idx, g_idx, b_idx)]], stretch = "lin")
  }
  
  output$rgb_input <- renderPlot({
    req(values$input_rast, values$wavelengths_in)
    render_rgb(values$input_rast, values$wavelengths_in, input)
  })
  
  output$rgb_output <- renderPlot({
    req(values$output_rast, values$wavelengths_out)
    render_rgb(values$output_rast, values$wavelengths_out, input)
  })
  
  observeEvent(input$click_input, {
    req(values$input_rast)
    xy <- c(input$click_input$x, input$click_input$y)
    spectra <- terra::extract(values$input_rast, vect(matrix(xy, ncol = 2), crs = crs(values$input_rast)))
    if (!is.null(spectra)) {
      sp <- as.numeric(spectra[1, -1])
      output$spectral_input <- renderPlotly({
        plot_ly(x = values$wavelengths_in, y = sp, type = 'scatter', mode = 'lines', name = "Input") %>%
          layout(title = "Input Spectral Profile", xaxis = list(title = "Wavelength (nm)"), yaxis = list(title = "Reflectance"))
      })
    }
  })
  
  observeEvent(input$click_output, {
    req(values$output_rast)
    xy <- c(input$click_output$x, input$click_output$y)
    spectra <- terra::extract(values$output_rast, vect(matrix(xy, ncol = 2), crs = crs(values$output_rast)))
    if (!is.null(spectra)) {
      sp <- as.numeric(spectra[1, -1])
      output$spectral_output <- renderPlotly({
        plot_ly(x = values$wavelengths_out, y = sp, type = 'scatter', mode = 'lines', name = "Output") %>%
          layout(title = "Output Spectral Profile", xaxis = list(title = "Wavelength (nm)"), yaxis = list(title = "Reflectance"))
      })
    }
  })
}

# === Run App ===
shinyApp(ui = ui, server = server)
