#' Opens non-normalized Illumina txt and xls files.
#'
#' Helper utility opens non-normalized Illumina files. User must check and
#' possibly edit the files.
#'
#' To help, I recommend Sublime Text 2 (a text editor with regular expression
#' capabilities). See package vignette for example of process.
#'
#' @param gse_names Character vector of GSE names to open raw files for.
#' @param data_dir String specifying directory with GSE folders.
#'
#' @export
#' @seealso \code{\link{get_raw}}, \code{\link{load_raw}}.
#' @return Character vector specifying successfully formatted GSEs. Only these
#'   Illumina GSEs can be loaded by \code{load_raw}.
#' @examples \dontrun{
#' library(lydata)
#' illum_names <- c("GSE50841", "GSE34817", "GSE29689")
#'
#' #this will open raw data using default program
#' illum_names <- open_raw_illum(illum_names, data_dir)
#'}

open_raw_illum <- function (gse_names, data_dir) {

    out_names <- gse_names
    for (i in seq_along(gse_names)) {
        #get data paths
        gse_dir <- paste(data_dir, gse_names[i], sep="/")
        data_paths <- list.files(gse_dir, pattern="non.norm.*txt", full.names=T)
        data_paths <- c(data_paths, list.files(gse_dir, pattern=".xls", full.names=T))

        #open data file
        for (j in seq_along(data_paths)) system2("xdg-open", data_paths[j])

        #check success
        success <- tcltk::tk_select.list(choices = c("Yes", "No"),
                                         title = paste(gse_names[i],
                                                       "formated successfully?"))
        #remove unsuccessful
        if (success == "No") out_names <- setdiff(out_names, gse_names[i])
    }
    return(out_names)
}

#------------------------

# Load and pre-process raw Illum files.
#
# Load raw txt files previously downloaded with \code{get_raw} and checked
# for format with \code{open_raw_illum}. Used by \code{load_raw}.
#
# Data is normalized, SYMBOL and PROBE annotation are added to fData slot, and
# detection p-values are added to pvals slot.
#
# @param gse_names Character vector of Illumina GSE names.
# @param data_dir String specifying directory with GSE folders.
#
# @seealso \code{\link{get_raw}} to obtain raw Illumina data.
#   \code{\link{open_raw_illum}} to ensure their correctness.

# @return List of annotated esets.

load_illum <- function (gse_names, data_dir) {

    esets <- list()
    for (gse_name in gse_names) {

      gse_dir <- paste(data_dir, gse_name, sep="/")

      #get GSEMatrix (for pheno data)
      eset <- GEOquery::getGEO(gse_name, destdir=gse_dir, GSEMatrix=T)[[1]]

      #load non-normalized txt files and normalize
      data_paths <- list.files(gse_dir, pattern="non.norm.*txt", full.names=T)
      data <- limma::read.ilmn(data_paths, probeid="ID_REF")
      data <- tryCatch (
          limma::neqc(data),
        error = function(cond) {
          data <- limma::backgroundCorrect(data, method="normexp",
                                           normexp.method="rma", offset=16)

          return(limma::normalizeBetweenArrays(data, method="quantile"))
        })

      #transfer exprs from data to eset (maintaining eset feature order)
      feature_order <- featureNames(eset)
      pData(eset)$title_GSEMatrix <- pData(eset)$title  #to check if sample order mismatch
      pData(eset)$title <- colnames(data)  #use raw data titles to ensure correct contrasts
      colnames(data) <- sampleNames(eset)
      exprs(eset) <- data$E[feature_order,]

      #transfer pvals from data to eset
      pvals <- data$other$Detection[feature_order, ]
      eset <- add_pvals(eset, pvals)

      #add SYMBOL annotation
      gpl_name <- annotation(eset)
      eset <- symbol_annot(eset, gpl_name)

      esets[[gse_name]] <- eset
      }
    return (esets)
}


# Add detection p-values to Illumina expression set.
#
# Adds detection p-vals to pvals slot of illumina expression set. Used by
# \code{load_illum}.
#
# @param eset Illumina expression set to add pvals slot to.
# @param pvals Detection slot obtained from \link{read.ilmn}
#
# @return Expression set with detection p-values in pvals slot.

add_pvals <- function (eset, pvals) {
    storageMode(eset) = "environment"
    assayData(eset)[["pvals"]] = pvals
    storageMode(eset) = "lockedEnvironment"
    return (eset)
}
