library(shiny)
library(emba)
library(usefun)
library(dplyr)
library(DT)
library(plotly)

source('fun.R')

# load initial example files data
input_dir = paste0(getwd(), '/examples')
original_file = paste0(input_dir, '/ensemble_synergies')
random_file = paste0(input_dir, '/ensemble_synergies_random')
observed_synergies_file = paste0(input_dir, '/observed_synergies')

observed_synergies = emba::get_observed_synergies(observed_synergies_file)
orig_res = emba::get_synergy_scores(original_file)
rand_res = emba::get_synergy_scores(random_file)

validate_input(orig_res, rand_res, observed_synergies)
predictions_example = create_predictions_df(orig_res, rand_res, observed_synergies)

ui = fluidPage(
  # App title
  titlePanel("Druglogics ROC"),

  sidebarLayout(
    # Selecting files menu
    sidebarPanel(width = 4,
      selectInput(inputId = 'method', label = 'Select Output',
                  choices = c('original', 'random', 'original-random',
                              'original+random')),
      fileInput(inputId = 'observed', label = 'Select Observed Synergies File',
                buttonLabel = 'Browse...', placeholder = 'No file selected'),
      fileInput(inputId = 'orig', label = 'Select Original Ensemble Synergies File',
                buttonLabel = 'Browse...', placeholder = 'No file selected'),
      fileInput(inputId = 'rand', label = 'Select Random Ensemble Synergies File',
                buttonLabel = 'Browse...', placeholder = 'No file selected'),
      htmlOutput(outputId = "message")
      ),
    # Plotting the  ROC curve
    mainPanel(
      plotlyOutput(outputId = 'rocPlot')
    )
  ),
  basicPage(DT::DTOutput(outputId = 'rocStatsTable'))
)

server = function(input, output) {

  get_message = reactiveVal(value = "Init")

  get_roc_stats = reactive({
    print("GET ROC STATISTICS")

    method = input$method
    observed_synergies_file = input$observed$datapath
    original_file = input$orig$datapath
    random_file = input$rand$datapath

    print(paste0("Method: ", method))
    print(paste0("observed_synergies_file: ", observed_synergies_file))
    print(paste0("original_file: ", original_file))
    print(paste0("random_file: ", random_file))

    # which `predictions` object to use?
    if (!is.null(observed_synergies_file) & (!is.null(original_file) | !is.null(random_file))) {
      observed_synergies = emba::get_observed_synergies(observed_synergies_file)
      orig_res = NULL
      rand_res = NULL
      if (!is.null(original_file) & !is.null(random_file)) {
        orig_res = emba::get_synergy_scores(original_file)
        rand_res = emba::get_synergy_scores(random_file)
      } else if (!is.null(original_file) & is.null(random_file)) {
        orig_res = emba::get_synergy_scores(original_file)
      } else if (is.null(original_file) & !is.null(random_file)) {
        rand_res = emba::get_synergy_scores(random_file)
      }

      validate_input(orig_res, rand_res, observed_synergies)
      get_message("Using uploaded files")
      predictions = create_predictions_df(orig_res, rand_res, observed_synergies)
    }
    else {
      # use example `predictions`
      get_message("Using example files")
      predictions = predictions_example
    }

    # Add columns (combine biomarkers)
    predictions = combine_biomarkers(predictions)

    # check if method is compatible with input data
    validate(
      need(method %in% colnames(predictions), "Please select an appropriate method!")
    )

    # generate ROC stats
    roc_stats = gen_roc_stats(predictions, method)
  })

  output$rocPlot = renderPlotly({
    roc_stats = get_roc_stats()

    x = roc_stats$FPR
    y = roc_stats$TPR
    auc = sum(diff(x) * (head(y,-1)+tail(y,-1)))/2

    plotly_roc(roc_stats, auc, input$method)
  })

  output$rocStatsTable = DT::renderDT({
    datatable(data = get_roc_stats(), options = list(
      searching = FALSE,
      caption = htmltools::tags$caption('ROC Statistics', style="color:#dd4814; font-size: 18px")
      )) %>% formatRound(c(1,6,7,8,9), digits = 3)
  })

  output$message = renderUI({
    htmlText = paste("<b><div style='background-color:lightblue; font-size:2em; width:fit-content'>",
      get_message(), "</div></b>")
    HTML(htmlText)
  })
}

shinyApp(ui, server)
