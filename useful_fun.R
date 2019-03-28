library(dplyr) # 0.8.0.1

# load observed synergies file
get.observed.synergies = function(observed.synergies.file) {
  lines = readLines(observed.synergies.file)
  observed.synergies = gsub("~", "-", lines)
  return(observed.synergies)
}

# load ensemble-wise synergy file(s)
read.ensemble.synergies.drame.file = function(ensemble.file) {
  res = NULL
  if (!is.null(ensemble.file)) {
    res = read.table(file = ensemble.file, sep = "\t",
                     col.names = c("perturbation", "score"),
                     stringsAsFactors = FALSE)
    res$perturbation = gsub("\\[|\\]", "", res$perturbation)
  }

  return(res)
}

data.validation = function(orig.res = NULL, rand.res = NULL, observed.synergies) {
  if (!is.null(orig.res) & !(is.null(rand.res)))
    check.perturbations(orig.res$perturbation, rand.res$perturbation)

  # check that the observed synergies are all part of the perturbations tested
  if (!is.null(orig.res))
    stopifnot(all(observed.synergies %in% orig.res$perturbation))
  if (!is.null(rand.res))
    stopifnot(all(observed.synergies %in% rand.res$perturbation))
}

# checks that the same perturbations were tested in the
# original and random files and in the same order
check.perturbations = function(perturbations.orig, perturbations.rand) {
  stopifnot(all(perturbations.orig == perturbations.rand))
}

# checks that you can perform the desired scaling/normalization `method`
# with the `predictions` given
method.check = function(method, predictions) {
  if ((!"original" %in% colnames(predictions) &
      (method %in% c("original", "exp-prod-norm", "exp-fold-change-norm"))) |
      (!"random" %in% colnames(predictions) &
       (method %in% c("random", "exp-prod-norm", "exp-fold-change-norm"))))
     stop("method selected not compatible with predictions data")
}

create.predictions.data.frame =
  function(orig.res = NULL, rand.res = NULL, observed.synergies) {
    predictions = NULL
    if (is.null(rand.res) & !is.null(orig.res)) {
      perturbations = orig.res$perturbation
      observed = sapply(perturbations %in% observed.synergies, as.numeric)
      predictions = cbind.data.frame(perturbations, observed, orig.res$score,
                                     stringsAsFactors = FALSE)
      colnames(predictions) = c("perturbation", "observed", "original")
    } else if (is.null(orig.res) & !is.null(rand.res)) {
      perturbations = rand.res$perturbation
      observed = sapply(perturbations %in% observed.synergies, as.numeric)
      predictions = cbind.data.frame(perturbations, observed, rand.res$score,
                                     stringsAsFactors = FALSE)
      colnames(predictions) = c("perturbation", "observed", "random")
    } else if (!is.null(orig.res) & !is.null(rand.res)) { # we have 2 drabme files
      perturbations = orig.res$perturbation
      observed = sapply(perturbations %in% observed.synergies, as.numeric)
      predictions = cbind.data.frame(perturbations, observed, orig.res$score,
                                     rand.res$score, stringsAsFactors = FALSE)
      colnames(predictions) = c("perturbation", "observed", "original", "random")
    }
    return(predictions)
  }

# add further columns to the `predictions` data.frame
normalize = function(predictions) {
  if(all(c("original", "random") %in% colnames(predictions))) {
    predictions = mutate(predictions, exp.prod.norm = exp(original + random))
    predictions = mutate(predictions, exp.fold.change.norm = exp(original - random))
  }
  return(predictions)
}

# `method` is the string defined in the UI whereas the
# returned value is the corresponding column name for that method
map.method.to.column = function(method) {
  if      (method == "exp-orig") return("exp.orig")
  else if (method == "exp-rand") return("exp.rand")
  else if (method == "exp-prod-norm") return("exp.prod.norm")
  else if (method == "exp-fold-change-norm") return("exp.fold.change.norm")
  else return(method)
}

# get the confusion matrix values (TP, FN, TN, FP) + TPR, FPR from a
# `predictions` data.frame by comparing the score values from the column
# `method.column` to the `thres.value`
get.conf.mat.for.thres = function(predictions, method.column, thres.value) {
  tp = 0
  fn = 0
  tn = 0
  fp = 0

  for(i in 1:nrow(predictions)) {
    value = predictions[i, method.column]
    obs   = predictions[i, "observed"]
    if (value < thres.value & obs == 1) {
      # print(paste0("Value: ", value))
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

  res = c(tp, fn, tn, fp, tpr, fpr, youden.index, dist.from.0.1)
  names(res) = c("TP", "FN", "TN", "FP", "TPR", "FPR", "YoudenIndex", "distFrom0.1")

  return(res)
}

gen.roc.stats = function(predictions, method.column) {
  thres.values = sort(unique(predictions[,method.column]))

  stats = list()
  index = 1
  for(thres.value in thres.values) {
    stats[[index]] =
      c(thres.value, get.conf.mat.for.thres(predictions, method.column, thres.value))
    index = index + 1
  }

  roc_stats = as.data.frame(do.call(rbind, stats))
  colnames(roc_stats)[1] = "threshold"

  return(roc_stats)
}

# `digits.to.keep` refers to digits after decimal point '.'
specify.decimal = function(number, digits.to.keep) {
  trimws(format(round(number, digits.to.keep), nsmall = digits.to.keep))
}

plot.roc = function(x, y, auc, method) {
  plot(x, y, type = "b", col = "blue",
       main = "ROC curve", xlab = "False Positive Rate (FPR)",
       ylab = "True Positive Rate (TPR)")
  legend("bottomright", legend = specify.decimal(AUC, digits.to.keep = 3),
         title = paste0("AUC (method: ", method, ")"), col = "blue", pch = 19)
  grid()
  abline(a = 0, b = 1, col = "lightgray", lty = 2) # y=bx+a
}

plotly.roc = function(roc.stats, auc, method) {
  # find optimal cutoffs
  max.youden.index = roc.stats %>%
    filter(YoudenIndex == max(YoudenIndex))
  min.dist.from.0.1 = roc.stats %>%
    filter(distFrom0.1 == min(distFrom0.1))

  chance.line = rbind.data.frame(c(0,0), c(1,1)) # y = x
  colnames(chance.line) = c("x", "y")

  plot_ly(height = 500, width = 500) %>%
    config(displayModeBar = FALSE) %>%
    layout(title = 'ROC curve', showlegend = FALSE,
           xaxis = list(title = 'False Positive Rate (FPR)'),
           yaxis = list(title = 'True Positive Rate (TPR)')) %>%
    add_annotations(x = 0.8, y = 0.25, showarrow = FALSE,
                    text = paste0("AUC: ", specify.decimal(auc, digits.to.keep = 3)),
                    font = list(color = "blue", size = 18)) %>%
    add_annotations(x = 0.75, y = 0.15, showarrow = FALSE,
                    text = paste0("Method: ", method),
                    font = list(color = "#264E86", size = 18)) %>%
    add_annotations(x = max.youden.index$FPR, y = max.youden.index$TPR,
                    ax = 55, ay = 60, showarrow = TRUE,
                    text = "Max Youden Index",
                    font = list(color = "green")) %>%
    add_annotations(x = min.dist.from.0.1$FPR, y = min.dist.from.0.1$TPR,
                    ax = 0, ay = -60, showarrow = TRUE,
                    text = "Min distance from (0,1)",
                    font = list(color = "purple")) %>%
    add_trace(data = roc.stats, x = ~FPR, y = ~TPR, name = "",
              text = paste0("Threshold: ", specify.decimal(
                roc.stats$threshold, digits.to.keep = 5)),
              color = I("blue"), type = "scatter", mode = "lines+markers") %>%
    add_trace(data = chance.line, x = ~x, y = ~y, type = "scatter", mode = "lines",
              line = list(color = 'lightgrey', width = 2, dash = 'dot'))
}
