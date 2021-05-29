#' Filter genes in RNA-seq ExpressionSet
#'
#' Uses \link[edgeR]{filterByExpr} to filter based on 'counts' assay or 'exprs'
#' assay if 'counts' isn't available (for ARCHS4 data).
#'
#'
#' @param eset ExpressionSet with 'counts' assayDataElement and group column in
#'  pData
#'
#' @return filtered \code{eset}
#' @export
#' @seealso \link[edgeR]{filterByExpr}
#' @examples
#'
#' # example ExpressionSet
#' eset <- makeExampleCountsEset()
#' eset <- filter_genes(eset)
#' 
filter_genes <- function(eset) {
  els <- Biobase::assayDataElementNames(eset)
  els <- ifelse("counts" %in% els, "counts", "exprs")
  counts <- Biobase::assayDataElement(eset, els)
  
  keep <- edgeR::filterByExpr(counts, group = eset$group)
  eset <- eset[keep, ]
  if (!nrow(eset)) stop("No genes with reads after filtering")
  return(eset)
}

# adapted from DESeq2::makeExampleDESeqDataSet
makeExampleCountsEset <- function (n = 1000, m = 12, betaSD = 0, interceptMean = 4, interceptSD = 2, dispMeanRel = function(x) 4/x + 0.1, sizeFactors = rep(1, m)) {
  beta <- cbind(rnorm(n, interceptMean, interceptSD), rnorm(n, 0, betaSD))
  dispersion <- dispMeanRel(2^(beta[, 1]))
  colData <- data.frame(condition = factor(rep(c("A", "B"), 
                                              times = c(ceiling(m/2), floor(m/2)))))
  x <- if (m > 1) {
    stats::model.matrix.default(~colData$condition)
  }
  else {
    cbind(rep(1, m), rep(0, m))
  }
  mu <- t(2^(x %*% t(beta)) * sizeFactors)
  countData <- matrix(rnbinom(m * n, mu = mu, size = 1/dispersion), 
                      ncol = m)
  
  colnames(countData) <- paste("sample", 1:m, sep = "")
  row.names(countData) <- paste0("gene", 1:n)
  
  eset <- Biobase::ExpressionSet(countData)
  eset$group <- colData$condition
  return(eset)

}

dispMeanRel <- function(x) {
  4/x + 0.1
}
