library(dplyr) # 0.8.0.1

# load ensemble-wise synergy file(s)
read.ensemble.synergies.drame.file = function(ensemble.file) {
  res = read.table(file = ensemble.file, sep = "\t",
                   col.names = c("perturbation", "score"),
                   stringsAsFactors = FALSE)
  res$perturbation = gsub("\\[|\\]", "", res$perturbation)
  return(res)
}

# checks that the same perturbations were tested in the
# original and random files and in the same order
check.perturbations = function(perturbations.orig, perturbations.rand) {
  stopifnot(all(perturbations.orig == perturbations.rand))
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