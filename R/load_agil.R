#' Download and unpack Agilent microarray supplementary files from GEO.
#'
#' Downloads Agilent microarray supplementary files from GEO.
#' Files are stored in the supplied data directory under the GSE name.
#'
#' @importFrom GEOquery getGEOSuppFiles gunzip
#' @importFrom utils untar
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
#'
get_raw_agil <- function (gse_names, data_dir) {

    for (gse_name in gse_names) {

        gse_dir <- paste(data_dir, gse_name, sep="/")

        if (!file.exists(gse_dir)) {
            getGEOSuppFiles(gse_name, baseDir=data_dir)
        }
        #untar
        tar_name <- list.files(gse_dir, pattern="tar")
        untar(paste(gse_dir, tar_name, sep="/"), exdir=gse_dir)

        #unzip
        paths <- list.files(gse_dir, pattern=".gz", full.names=T, ignore.case=T)
        sapply(paths, gunzip, overwrite=T)
    }
}



#------------



#' Load and pre-process raw Agilent files for multiple GSEs.
#'
#' Load raw txt files previously downloaded with \code{get_raw_agil}. Data is
#' neqc normalized and gene SYMBOLs are added to featureData slots.
#'
#' @importFrom Biobase fData
#' @importFrom GEOquery getGEO
#' @importFrom limma read.maimages neqc
#' @importFrom stringr str_match
#' @importFrom BiocGenerics annotation
#'
#' @param gse_names Character vector of GSE names to load and pre-process.
#' @param data_dir Character, base data directory (contains a folder with raw
#'        data for each GSE to be loaded).
#' @export
#' @seealso \code{\link{get_raw_agil}} to obtain raw Agilent txt files for
#'          multiple GSEs.
#' @return list of processed esets.
#' @examples \dontrun{
#'
#'}
#'
load_agil <- function (gse_names, data_dir) {

    esets <- list()
    for (gse_name in gse_names) {

      gse_dir <- paste(data_dir, gse_name, sep="/")

      #get GSEMatrix (for pheno data)
      eset <- getGEO(gse_name, destdir=gse_dir, GSEMatrix=T)[[1]]

      #load non-normalized txt files and normalize
      data_paths <- list.files(gse_dir, pattern="GSM.*txt", full.names=T, ignore.case=T)
      data <- read.maimages(data_paths, source="agilent", green.only=TRUE)
      data <- neqc(data, status=data$genes$ControlType, negctrl=-1, regular=0)

      #fix up sample/feature names
      colnames(data) <- str_match(colnames(data), ".*(GSM\\d+).*")[, 2]
      row.names(data$E)     <- data$genes$ProbeName
      data$genes <- data$genes[!duplicated(data$genes$ProbeName),]
      row.names(data$genes) <- data$genes$ProbeName

      #transfer to eset
      exprs(eset) <- data$E
      fData(eset) <- data$genes

      #add SYMBOL annotation
      gpl_name <- annotation(eset)
      eset <- symbol_annot(eset, gpl_name)

      esets[[gse_name]] <- eset
    }
    return(esets)
}
