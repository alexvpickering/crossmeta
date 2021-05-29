#' Get variance stabilized data for exploratory data analysis
#'
#' @param eset ExpressionSet loaded with \link{load_raw}.
#'
#' @return \code{matrix} with variance stabilized expression data.
#' 
#' @export
#' 
get_vsd <- function(eset) {
  pdata <- Biobase::pData(eset)
  lib.size <- pdata$lib.size * pdata$norm.factors
  
  # GK Smyth advice for variance stabilization for plotting
  # https://support.bioconductor.org/p/114133/
  vsd <- edgeR::cpm(Biobase::exprs(eset),
                    lib.size = lib.size,
                    log = TRUE, prior.count = 5)
  
  return(vsd)
}