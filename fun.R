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

create_predictions_df = function(orig_res = NULL, rand_res = NULL, observed_synergies) {
    predictions = NULL
    if (is.null(rand_res) & !is.null(orig_res)) {
      perturbations = orig_res$perturbation
      observed = sapply(perturbations %in% observed_synergies, as.numeric)
      predictions = cbind.data.frame(perturbations, observed, orig_res$score,
                                     stringsAsFactors = FALSE)
      colnames(predictions) = c('perturbation', 'observed', 'original')
    } else if (is.null(orig_res) & !is.null(rand_res)) {
      perturbations = rand_res$perturbation
      observed = sapply(perturbations %in% observed_synergies, as.numeric)
      predictions = cbind.data.frame(perturbations, observed, rand_res$score,
                                     stringsAsFactors = FALSE)
      colnames(predictions) = c('perturbation', 'observed', 'random')
    } else if (!is.null(orig_res) & !is.null(rand_res)) { # we have 2 drabme files
      perturbations = orig_res$perturbation
      observed = sapply(perturbations %in% observed_synergies, as.numeric)
      predictions = cbind.data.frame(perturbations, observed, orig_res$score,
                                     rand_res$score, stringsAsFactors = FALSE)
      colnames(predictions) = c('perturbation', 'observed', 'original', 'random')
    }
    return(predictions)
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

gen_roc_stats = function(predictions, method) {
  thresholds = sort(unique(predictions[,method]))

  stats = list()
  index = 1
  for(thres in thresholds) {
    stats[[index]] = c(thres, get_conf_mat_for_thres(predictions, method, thres))
    index = index + 1
  }

  # get the (1,1) point in the ROC curve!
  stats[[index]] = c(thres + 0.001, get_conf_mat_for_thres(predictions, method, thres + 0.001))

  roc_stats = as.data.frame(do.call(rbind, stats))
  colnames(roc_stats)[1] = 'threshold'

  return(roc_stats)
}

# get the confusion matrix values (TP, FN, TN, FP) + TPR, FPR from a
# `predictions` data.frame by comparing the score values from the column
# `method` to the `thres.value`
get_conf_mat_for_thres = function(predictions, method, thres.value) {
  tp = 0
  fn = 0
  tn = 0
  fp = 0

  for(i in 1:nrow(predictions)) {
    value = predictions[i, method]
    obs   = predictions[i, 'observed']
    if (value < thres.value & obs == 1) {
      tp = tp + 1
    } else if (value < thres.value & obs == 0) {
      fp = fp + 1
    } else if (value >= thres.value & obs == 1) {
      fn = fn + 1
    } else if (value >= thres.value & obs == 0) {
      tn = tn + 1
    }
  }

  tpr = tp / (tp + fn)
  fpr = fp / (fp + tn)
  youden.index = tpr - fpr
  dist.from.0.1 = (fpr - 0)^2 + (tpr - 1)^2

  res = c(tp, fn, tn, fp, fpr, tpr, youden.index, dist.from.0.1)
  names(res) = c('TP', 'FN', 'TN', 'FP', 'FPR', 'TPR', 'YoudenIndex', 'distFrom0.1')

  return(res)
}

# plot_roc = function(x, y, auc, method) {
#   plot(x, y, type = 'b', col = 'blue',
#        main = 'ROC curve', xlab = 'False Positive Rate (FPR)',
#        ylab = 'True Positive Rate (TPR)')
#   legend('bottomright', legend = specify_decimal(auc, digits.to.keep = 3),
#          title = paste0('AUC (method: ', method, ')'), col = 'blue', pch = 19)
#   grid()
#   abline(a = 0, b = 1, col = 'lightgray', lty = 2) # y=bx+a
# }

plotly_roc = function(roc_stats, auc, method) {
  # find optimal cutoffs
  max.youden.index = roc_stats %>%
    filter(YoudenIndex == max(YoudenIndex))
  min.dist.from.0.1 = roc_stats %>%
    filter(distFrom0.1 == min(distFrom0.1))

  chance.line = rbind.data.frame(c(0,0), c(1,1)) # y = x
  colnames(chance.line) = c('x', 'y')

  plot_ly(height = 500, width = 550) %>%
    config(displayModeBar = FALSE) %>%
    layout(title = 'ROC curve', showlegend = FALSE,
           xaxis = list(title = 'False Positive Rate (FPR)'),
           yaxis = list(title = 'True Positive Rate (TPR)'),
           margin = list(l = 50, r = 50, b = 30, t = 50, pad = 4)) %>%
    add_annotations(x = 0.8, y = 0.25, showarrow = FALSE,
                    text = paste0('AUC: ', specify_decimal(auc, digits.to.keep = 3)),
                    font = list(color = 'blue', size = 18)) %>%
    add_annotations(x = 0.75, y = 0.15, showarrow = FALSE,
                    text = paste0('Method: ', method),
                    font = list(color = '#264E86', size = 18)) %>%
    add_annotations(x = max.youden.index$FPR, y = max.youden.index$TPR,
                    ax = 55, ay = 60, showarrow = TRUE,
                    text = 'Max Youden Index',
                    font = list(color = 'green')) %>%
    add_annotations(x = min.dist.from.0.1$FPR, y = min.dist.from.0.1$TPR,
                    ax = 0, ay = -60, showarrow = TRUE,
                    text = 'Min distance from (0,1)',
                    font = list(color = 'purple')) %>%
    add_trace(data = roc_stats, x = ~FPR, y = ~TPR, name = '',
              text = paste0('Threshold: ', specify_decimal(
                roc_stats$threshold, digits.to.keep = 5)),
              color = I('blue'), type = 'scatter', mode = 'lines+markers') %>%
    add_trace(data = chance.line, x = ~x, y = ~y, type = 'scatter', mode = 'lines',
              line = list(color = 'lightgrey', width = 2, dash = 'dot'))
}
