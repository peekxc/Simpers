#' dgm_var_UI
#' @description Variables UI control selection for persistence diagrams.
#' @param id, character used to specify namespace, see \code{shiny::\link[shiny]{NS}}
#' @return a \code{shiny::\link[shiny]{tagList}} containing UI elements
#' @details \code{shiny::\link[shiny]{tagList}} returned includes a selectizeInput, a numericInput, and
#' a dynamically generated uiOutput that becomes a range slider.
#' @import shiny
#' @export
dgm_var_UI <- function(id) {
  ns <- NS(id)
  tagList(
    selectizeInput(ns("dim_select"), label = "Dimensions", choices = c("H0", "H1", "H2"), selected = list("H0", "H1"), multiple = TRUE),
    numericInput(ns("noise_threshold"), label = "Noise threshold", value = 0, min = 0),
    uiOutput(ns("time_window_ui")),
    selectInput(ns("plot_type"), label = "Plot type", choices = c("Barcodes", "Diagrams", "Landscapes"), selected = "Barcodes")
  )
}

#' dgm_var
#' @param dgm a data.frame with names 'appear', 'disappear', and 'dim'.
#' @return list with following components
#' \describe{
#'   \item{time_window}{reactive character string indicating x variable selection}
#'   \item{dims_selected}{reactive character string indicating y variable selection}
#'   \item{noise_threshold}{ }
#' }
#' @export
dgm_var <- function(input, output, session, dgm){
  if (any(c("appear", "disappear", "dim") %in% names(dgm))){
    dgm_names <- names(dgm)
    names(dgm) <- sapply(dgm_names, function(nm) switch(nm, "dim"="dimension", "appear"="birth", "disappear"="death", nm))
  }
  stopifnot(all(c("birth", "death", "dimension") %in% names(dgm)), is.data.frame(dgm), is.factor(dgm$dimension))

  ## Time window
  output$time_window_ui <- renderUI({
    ns <- session$ns
    min_time <- min(dgm$birth)
    max_time <- max(dgm$death[dgm$death != Inf])
    sliderInput(ns("time_window"), label = "Time window",
                min = round(min_time, 3), max = round(max_time, 3),
                value = c(min_time, max_time), dragRange = TRUE)
  })

  ## Update noise threshold step value
  updateNumericInput(session, inputId = "noise_threshold", step = min(dgm$death - dgm$birth))

  ## Return the selected diagram
  return(list(
    time_window = reactive({ input$time_window }),
    dims_selected = reactive({ (match(input$dim_select, c("H0", "H1", "H2"))-1L) }),
    noise_threshold = reactive({ input$noise_threshold })
  ))
}


#' linkedBarcodeOutput
#' @description Output UI module for creating a barcode plot.
#' @export
linkedBarcodeOutput <- function(id, ...){
  ns <- NS(id)
  plotOutput(ns("barcode_plot"), click = ns("bar_click"), ...)
}

#' linkedBarcode
#' @param dgm a (non-reactive) data.frame with names 'appear', 'disappear', and 'dim'.
#' @param dgm_vars list containing reactive expressions 'time_window', 'dims_selected', and 'noise_threshold'
#' @param plot_params reactive expression of ggplot2 elements to add
#' @return list with following components
#' \describe{
#'   \item{xvar}{reactive character string indicating x variable selection}
#'   \item{yvar}{reactive character string indicating y variable selection}
#' }
#' @import ggtda
#' @export
linkedBarcode <- function(input, output, session, dgm, dgm_vars, plot_params=reactive({ NULL })){
  if (any(c("appear", "disappear", "dim") %in% names(dgm))){
    dgm_names <- names(dgm)
    names(dgm) <- sapply(dgm_names, function(nm) switch(nm, "dim"="dimension", "appear"="birth", "disappear"="death", nm))
  }
  stopifnot(all(c("birth", "death", "dimension") %in% names(dgm)), is.data.frame(dgm), is.factor(dgm$dimension))

  ## Reactive that chooses the subset of the diagram based on the 'dgm_vars'
  lifetimes <- abs(dgm$death - dgm$birth)
  sub_dgm <- reactive({
    dim_idx <- as.integer(as.character(dgm$dimension)) %in% dgm_vars$dims_selected()
    thresh_idx <- lifetimes >= dgm_vars$noise_threshold()
    dgm[(dim_idx & thresh_idx),,drop=FALSE]
  })

  ## The output barcodes
  output$barcode_plot <- renderPlot({
    aes_map <- ggplot2::aes(start = birth, end = death, colour = dimension)
    ggplot2::ggplot(sub_dgm(), aes_map) +
      ggtda::geom_barcode(size=1.5) +
      ggtda::theme_tda() +
      ggplot2::coord_cartesian(xlim = dgm_vars$time_window()) +
      plot_params()
  })

  ## Return the current diagram subset in view, and the index of the selected feature
  return(list(
    dgm_subset = sub_dgm,
    selected_feature_idx = reactive({
      if (!is.null(input$bar_click)) {
        dim_idx <- as.integer(as.character(dgm$dimension)) %in% dgm_vars$dims_selected()
        thresh_idx <- lifetimes >= dgm_vars$noise_threshold()
        which(dim_idx & thresh_idx)[round(input$bar_click$y)]
      } else { NULL } })
  ))
}



