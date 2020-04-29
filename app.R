library(shiny)
library(emba)
library(usefun)
library(dplyr)
library(DT)
library(plotly)
library(PRROC)

########################
### Useful Functions ###
########################
plotly_roc = function(roc_stats, AUC, method) {
  # find optimal cutoffs
  youden_index_df = roc_stats %>%
    filter(dist_from_chance == max(dist_from_chance))
  min_dist_from_0_1 = roc_stats %>%
    filter(dist_from_0_1 == min(dist_from_0_1))

  chance_line = rbind.data.frame(c(0,0), c(1,1)) # y = x
  colnames(chance_line) = c('x', 'y')

  plot_ly() %>%
    config(displayModeBar = FALSE) %>%
    layout(title = 'ROC curve', showlegend = FALSE,
      xaxis = list(title = 'False Positive Rate (FPR)'),
      yaxis = list(title = 'True Positive Rate (TPR)'),
      margin = list(l = 50, r = 50, b = 30, t = 50, pad = 4)) %>%
    add_annotations(x = 0.8, y = 0.25, showarrow = FALSE,
      text = paste0('AUC: ', specify_decimal(AUC, digits.to.keep = 3)),
      font = list(color = 'blue', size = 18)) %>%
    add_annotations(x = 0.75, y = 0.15, showarrow = FALSE,
      text = paste0('Method: ', method),
      font = list(color = '#264E86', size = 18)) %>%
    add_annotations(x = youden_index_df$FPR, y = youden_index_df$TPR,
      ax = 55, ay = 60, showarrow = TRUE,
      text = 'Max Youden Index',
      font = list(color = 'green')) %>%
    add_annotations(x = min_dist_from_0_1$FPR, y = min_dist_from_0_1$TPR,
      ax = 0, ay = -60, showarrow = TRUE,
      text = 'Min distance from (0,1)',
      font = list(color = 'purple')) %>%
    add_trace(data = roc_stats, x = ~FPR, y = ~TPR, name = '',
      text = paste0('Threshold: ', specify_decimal(
        roc_stats$threshold, digits.to.keep = 5)),
      color = I('blue'), type = 'scatter', mode = 'lines+markers') %>%
    add_trace(data = chance_line, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
      line = list(color = 'lightgrey', width = 2, dash = 'dot'))
}

# adds 2 further columns to the `predictions` data.frame
# 2 linear combinations of biomarkers (+,-)
combine_biomarkers = function(predictions) {
  if(all(c('original', 'random') %in% colnames(predictions))) {
    predictions = mutate(predictions, 'original+random' = original + random)
    predictions = mutate(predictions, 'original-random' = original - random)
  }
  return(predictions)
}

validate_input = function(orig_res = NULL, rand_res = NULL, observed.synergies) {
  if (!is.null(orig_res) & !(is.null(rand_res)))
    check_perturbations(orig_res$perturbation, rand_res$perturbation)

  # check that the observed synergies are all part of the perturbations tested
  if (!is.null(orig_res))
    stopifnot(all(observed.synergies %in% orig_res$perturbation))
  if (!is.null(rand_res))
    stopifnot(all(observed.synergies %in% rand_res$perturbation))
}

# checks that the same perturbations were tested in the
# original and random files and in the same order
check_perturbations = function(perturbations.orig, perturbations.rand) {
  stopifnot(all(perturbations.orig == perturbations.rand))
}

create_predictions_tbl = function(orig_res = NULL, rand_res = NULL, observed_synergies) {
  predictions = NULL
  if (is.null(rand_res) & !is.null(orig_res)) {
    perturbations = orig_res$perturbation
    observed = sapply(perturbations %in% observed_synergies, as.numeric)
    predictions = as_tibble(cbind.data.frame(perturbations, observed, orig_res$score,
      stringsAsFactors = FALSE))
    colnames(predictions) = c('perturbation', 'observed', 'original')
  } else if (is.null(orig_res) & !is.null(rand_res)) {
    perturbations = rand_res$perturbation
    observed = sapply(perturbations %in% observed_synergies, as.numeric)
    predictions = as_tibble(cbind.data.frame(perturbations, observed, rand_res$score,
      stringsAsFactors = FALSE))
    colnames(predictions) = c('perturbation', 'observed', 'random')
  } else if (!is.null(orig_res) & !is.null(rand_res)) { # we have 2 drabme files
    perturbations = orig_res$perturbation
    observed = sapply(perturbations %in% observed_synergies, as.numeric)
    predictions = as_tibble(cbind.data.frame(perturbations, observed, orig_res$score,
      rand_res$score, stringsAsFactors = FALSE))
    colnames(predictions) = c('perturbation', 'observed', 'original', 'random')
  }
  return(predictions)
}

#########################
### Load example data ###
#########################
input_dir = paste0(getwd(), '/examples')
original_file = paste0(input_dir, '/ensemble_synergies')
random_file = paste0(input_dir, '/ensemble_synergies_random')
observed_synergies_file = paste0(input_dir, '/observed_synergies')

observed_synergies = emba::get_observed_synergies(observed_synergies_file)
orig_res = emba::get_synergy_scores(original_file)
rand_res = emba::get_synergy_scores(random_file)

validate_input(orig_res, rand_res, observed_synergies)
predictions_example = create_predictions_tbl(orig_res, rand_res, observed_synergies)

##########
### UI ###
##########
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
    # ROC & PR curves
    mainPanel(
      fluidRow(
        splitLayout(
          style = "border: 1px solid silver;",
          cellWidths = 500,
          cellArgs = list(style = "padding: 6px"),
          plotlyOutput(outputId = 'rocPlot'),
          plotOutput(outputId = 'prPlot')
          )
        )
    )
  ),
  basicPage(DT::DTOutput(outputId = 'rocStatsTable'))
)

##############
### SERVER ###
##############
server = function(input, output) {

  get_message = reactiveVal(value = "Init")

  get_data = reactive({
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
      predictions = create_predictions_tbl(orig_res, rand_res, observed_synergies)
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

    data_list = list()
    data_list[["pred"]] = predictions
    data_list[["method"]] = method
    return(data_list)
  })

  get_data_roc = reactive({
    data_list = get_data()
    predictions = data_list$pred
    method = data_list$method

    # generate ROC stats
    usefun::get_roc_stats(df = predictions, pred_col = method, label_col = "observed", direction = "<")
  })

  get_data_pr = reactive({
    data_list = get_data()
    predictions = data_list$pred
    method = data_list$method

    PRROC::pr.curve(scores.class0 = predictions %>% pull(method) %>% (function(x) {-x}),
      weights.class0 = predictions %>% pull(observed), curve = TRUE, rand.compute = TRUE)
  })

  output$rocPlot = renderPlotly({
    res = get_data_roc()
    plotly_roc(res$roc_stats, res$AUC, input$method)
  })

  output$prPlot = renderPlot({
    res = get_data_pr()
    plot(res, color = "#377EB8", rand.plot = TRUE)
  })

  output$rocStatsTable = DT::renderDT({
    datatable(data = get_data_roc()[[1]], options = list(
      searching = FALSE,
      caption = htmltools::tags$caption('ROC Statistics', style="color:#dd4814; font-size: 18px")
      )) %>% formatRound(c(1,6,7,8,9), digits = 3)
  }, server = FALSE)

  output$message = renderUI({
    htmlText = paste("<b><div style='background-color:lightblue; font-size:2em; width:fit-content'>",
      get_message(), "</div></b>")
    HTML(htmlText)
  })
}

shinyApp(ui, server)
