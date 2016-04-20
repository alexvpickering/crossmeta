


same_features <- function(esets) {

  #same number?
  n_features <- sapply(esets, function(eset) nrow(eset))
  if (diff(range(n_features) >= 1)) return(FALSE)

  #same names?
  feature_names <- lapply(esets, function(eset) featureNames(eset))
  com <- Reduce(intersect, feature_names)
  if (length(com) != length(feature_names[[1]])) return(FALSE)

  return(TRUE)
}
