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
      (method %in% c("original", "exp-orig", "exp-prod-norm", "exp-fold-change-norm"))) |
      (!"random" %in% colnames(predictions) &
       (method %in% c("random", "exp-rand", "exp-prod-norm", "exp-fold-change-norm"))))
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
normalize.predictions = function(predictions) {
  if ("original" %in% colnames(predictions)) {
    predictions = mutate(predictions, exp.orig = exp(original))
  }
  if ("random" %in% colnames(predictions)) {
    predictions = mutate(predictions, exp.rand = exp(random))
  }
  if(all(c("original", "random") %in% colnames(predictions))) {
    predictions = mutate(predictions, exp.prod.norm = exp(original + random))
    predictions = mutate(predictions, exp.fold.change.norm = exp(original - random))
  }
  return(predictions)
}