#' Download and unpack Illumina microarray supplementary files from GEO.
#'
#' Downloads Illumina microarray supplementary files from GEO.
#' Files are stored in the supplied data directory under the GSE name.
#'
#'
#' @importFrom GEOquery getGEOSuppFiles gunzip
#'
#' @param gse_names Character vector of GSE names to download.
#' @param data_dir base data directory (a folder for each GSE will be created
#'        here).
#' @export
#' @seealso \code{\link{get_raw_illum}}, \code{\link{get_raw_agil}}.
#' @return NULL (for download/unpack only).
#' @examples \dontrun{
#'
#'}

get_raw_illum <- function(gse_names, data_dir) {

    for (gse_name in gse_names) {
      gse_dir <- paste(data_dir, gse_name, sep="/")

      #get raw data
      if (!file.exists(gse_dir)) {
        getGEOSuppFiles(gse_name, baseDir=data_dir)
      }
      #unzip
      gz_paths <- list.files(gse_dir, pattern=".gz", full.names=T, ignore.case=T)
      sapply(gz_paths, gunzip)
    }
}


#------------------------


#' Opens non-normalized Illumina txt and xls files.
#'
#' Helper utility to open non-normalized Illumina files. User must check and
#' possibly edit the files. To help, I recommend downloading Sublime Text 2
#' (a text editor with regular expression capabilities).
#'
#' @importFrom tcltk tk_select.list
#'
#' @param gse_names Character vector of GSE names to open raw files for.
#' @param data_dir base data directory (contains folder for each GSE).
#'
#' @export
#' @seealso \code{\link{load_illum}}, \code{\link{get_raw_illum}}.
#' @return Character vector of successfully formatted GSEs (only load these).
#' @examples \dontrun{
#'
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
        success <- tk_select.list(choices = c("Yes", "No"),
                                  title = paste(gse_names[i],
                                                "formated successfully?"))
        #remove unsuccessful
        if (success == "No") out_names <- setdiff(out_names, gse_names[i])
    }
    return(out_names)
}

#------------------------

#' Load and pre-process raw Illum files for multiple GSEs.
#'
#' Load raw txt files previously downloaded with \code{get_raw_illum} and checked
#' for format (see \code{\link{open_raw_illum}}). Data is normalized using neqc,
#' detection p-values are added to pvals slot and gene SYMBOLs are added to
#' featureData slot.
#'
#' @importFrom GEOquery getGEO
#' @importFrom limma read.ilmn neqc backgroundCorrect
#' @importFrom Biobase featureNames sampleNames pData exprs
#'
#' @param gse_names Character vector of GSE names to load and pre-process.
#' @param data_dir Character, base data directory (contains a folder with raw
#'        data for each GSE to be loaded).
#' @export
#' @seealso \code{\link{get_raw_illum}} to obtain raw Illumina files and
#'          \code{\link{open_raw_illum}} to ensure their correctness.
#' @return list of processed esets.
#' @examples \dontrun{
#'
#'}

load_illum <- function (gse_names, data_dir) {

    esets <- list()
    for (gse_name in gse_names) {

      gse_dir <- paste(data_dir, gse_name, sep="/")

      #get GSEMatrix (for pheno data)
      eset <- getGEO(gse_name, destdir=gse_dir, GSEMatrix=T)[[1]]

      #load non-normalized txt files and normalize
      data_paths <- list.files(gse_dir, pattern="non.norm.*txt", full.names=T)
      data <- read.ilmn(data_paths, probeid="ID_REF")
      data <- tryCatch (
        neqc(data),
        error = function(cond) {
          return(backgroundCorrect(data,method="normexp"))}
        )

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


#' Add detection p-values to Illumina expression set.
#'
#' Adds detection p-vals to pvals slot of illumina expression set.
#'
#' @importFrom Biobase storageMode assayData
#'
#' @param eset Illumina expression set to add pvals slot to.
#' @param pvals in Detection slot obtained from \link{read.ilmn}
#'
#' @return expression set with pvals slot containing detection p-values.
#' @examples \dontrun{
#'
#'}

add_pvals <- function (eset, pvals) {
    storageMode(eset) = "environment"
    assayData(eset)[["pvals"]] = pvals
    storageMode(eset) = "lockedEnvironment"
    return (eset)
}
