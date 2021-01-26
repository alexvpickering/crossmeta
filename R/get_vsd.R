#' Get variance stabilized data for exploratory data analysis
#'
#' @param eset ExpressionSet loaded with \link{load_raw}.
#'   Requires group column in \code{pData(eset)} specifying sample groupings.
#' @param rlog_cutoff Sample number above which will use
#'   \code{\link[DESeq2]{vst}} instead of \code{\link[DESeq2]{rlog}}.
#'   Default is 50.
#'
#' @return \code{DESeqTransform} with variance stabilized expression data.
#' 
#' @keywords internal
#' 
get_vsd <- function(eset, rlog_cutoff = 50) {
  trans_fun <- if (ncol(eset) > rlog_cutoff) DESeq2::vst else DESeq2::rlog
  els <- Biobase::assayDataElementNames(eset)
  pdata <- Biobase::pData(eset)
  
  if (all(c("abundance", "counts", "length") %in% els)) {
    txi.deseq <- list(
      countsFromAbundance = "no",
      abundance = Biobase::assayDataElement(eset, "abundance"),
      counts = Biobase::assayDataElement(eset, "counts"),
      length = Biobase::assayDataElement(eset, "length")
    )
    
    dds <- DESeq2::DESeqDataSetFromTximport(txi.deseq,
                                            pdata,
                                            design = ~group
    )
  } else {
    # this is e.g. for eset from load_archs4_seq
    dds <- DESeq2::DESeqDataSetFromMatrix(Biobase::exprs(eset),
                                          pdata,
                                          design = ~group
    )
  }
  dds <- DESeq2::estimateSizeFactors(dds)
  vsd <- trans_fun(dds, blind = FALSE)
  
  return(vsd)
}