#' Download and unpack Affymetrix microarray supplementary files from GEO.
#'
#' Downloads Affymetrix microarray supplementary files from GEO.
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

get_raw_affy <- function (gse_names, data_dir) {

    for (gse_name in gse_names) {

      gse_dir <- paste(data_dir, gse_name, sep="/")
      #get raw data
      if (!file.exists(gse_dir)) {
        getGEOSuppFiles(gse_name, baseDir=data_dir)
      }
      #untar
      tar_name <- list.files(gse_dir, pattern="tar")
      untar(paste(gse_dir, tar_name, sep="/"), exdir=gse_dir)

      #unzip
      cel_paths <- list.files(gse_dir, pattern=".CEL.gz",
                              full.names=T, ignore.case=T)
      sapply(cel_paths, gunzip, overwrite=T)
    }
}


#' Extract run date from Affymetrix CEL file.
#'
#' Useful for defining sample batches for \code{\link{ComBat}}. No longer
#' used by crossmeta, which discovers nuissance variables using \code{\link{sva}}.
#'
#' @importFrom affxparser readCelHeader
#'
#' @param cel_paths Charactor vector, full paths to CEL files.
#'
#' @export
#' @seealso \code{\link{ComBat}}
#' @return Factor vector of CEL run dates.
#'
cel_dates <- function(cel_paths) {
  #IN: vector with paths to CEL files
  #OUT: vector of CEL scan dates
  scan_dates <- c()
  for (i in seq_along(cel_paths)) {
    datheader <- readCelHeader(cel_paths[i])$datheader
    scan_date <- gsub(".*([0-9]{2}/[0-9]{2}/[0-9]{2}).*", "\\1", datheader)
    scan_dates[i] <- scan_date
  }
  return (as.factor(scan_dates))
}


#------------------------


#' Load and pre-process raw Affymetrix CEL files for multiple GSEs.
#'
#' Load raw CEL files previously downloaded with \code{get_raw_affy}. Data is
#' RMA normalized, CEL scan dates are extracted and added to phenoData slots,
#' and gene SYMBOLs are added to featureData slots.
#'
#' @importFrom GEOquery getGEO
#'
#' @param gse_names Character vector of GSE names to load and pre-process.
#' @param data_dir Character, base data directory (contains a folder with raw
#'        data for each GSE to be loaded).
#' @export
#' @seealso \code{\link{get_raw_affy}} to obtain raw Affymetrix CEL files for
#'          multiple GSEs.
#' @return list of processed esets (one for each unique GSE/GPL platform).
#' @examples \dontrun{
#'
#'}

load_affy <- function (gse_names, data_dir) {

    esets <- list()
    for (gse_name in gse_names) {
      gse_dir <- paste(data_dir, gse_name, sep="/")

      #get GSEMatrix (for pheno data)
      eset <- getGEO(gse_name, destdir=gse_dir, GSEMatrix=T)

      #load eset for each platform in GSE
      esets[[gse_name]] <- lapply(eset,load_affy_plat, gse_dir)
    }

    eset_names <- get_eset_names(esets, gse_names)
    esets <- unlist(esets)
    names(esets) <- eset_names
    return (esets)
}


#' Get eset names for load_affy.
#'
#' Helper function to get around issue of a single GSE having multiple platforms
#' (and thus \link{getGEO} returns a list of esets). To distinguish these cases,
#'  the GPL platform is appended to the GSE name.
#'
#'
#' @importFrom BiocGenerics annotation
#'
#' @param esets loaded by \code{load_affy}.
#' @param gse_names argument supplied to \code{load_affy}.
#'
#' @seealso \code{\link{load_affy}}
#' @return Character vector of GSE names with GPL appended when multiple
#'         platforms per GSE.
#' @examples \dontrun{
#'
#'}
#'
get_eset_names <- function(esets, gse_names) {
  eset_names <- c()

  for (i in seq_along(esets)) {
    #get gse name
    gse_name <- gse_names[i]

    if (length(esets[[i]]) > 1) {
      #add gpl_name to make gse_name unique
      gpl_name <- sapply(esets[[i]], annotation)
      gse_name <- paste(gse_name, gpl_name, sep=".")
    }
    #add gse_name to eset_names
    eset_names <- c(eset_names, gse_name)
  }
  return(eset_names)
}



#' Helper utility for load_affy.
#'
#' Used by load_affy to load an eset for each GPL platform in a GSE.
#'
#' @importFrom Biobase sampleNames featureNames fData pData
#' @importFrom affy ReadAffy
#' @importFrom oligo read.celfiles
#' @importFrom stringr str_extract
#' @importFrom BiocGenerics annotation
#'
#' @param eset GSEMatrix obtained by load_affy call to getGEO.
#' @param gse_dir directory containing raw data for GSE.
#'
#' @seealso \code{\link{load_affy}}
#' @return eset with scan_date in pData slot and SYMBOL in fData slot.
#' @examples \dontrun{
#'
#'}
#'
load_affy_plat <- function (eset, gse_dir) {
    #used by load_affy to load eset for each platform in GSE

    sample_names <- sampleNames(eset)
    pattern <- paste(".*", sample_names, ".*CEL", collapse="|", sep="")

    cel_paths <- list.files(gse_dir, pattern, full.names=T, ignore.case=T)
    data <- tryCatch (
        {
        raw_data <- ReadAffy (celfile.path=gse_dir)
        affy::rma(raw_data)
        },
        warning = function(cond) {
            raw_data <- read.celfiles(cel_paths)
            return (oligo::rma(raw_data))
        },
        error = function(cond) {
            raw_data <- read.celfiles(cel_paths)
            return (oligo::rma(raw_data))
        }
    )
    #rename samples in data
    sampleNames(data) <- str_extract(sampleNames(data), "GSM[0-9]+")

    #transfer exprs from data to eset (maintaining eset sample/feature order)
    sample_order <- sampleNames(eset)
    feature_order <- featureNames(eset)
    eset <- tryCatch (
        {
        exprs(eset) <- exprs(data)[feature_order, sample_order]
        eset
        },
        #if features don't match: also transfer featureData from data to eset
        error = function(cond) {
            exprs(eset) <- exprs(data)[, sample_order]
            fData(eset) <- fData(data)
            return(eset)
        }
    )

    #add scan dates to pheno data (maintaining eset sample order)
    scan_dates <- cel_dates (cel_paths)
    names(scan_dates) <- sampleNames(data)
    pData(eset)$scan_date <- scan_dates[sample_order]

    #add SYMBOL annotation
    gpl_name <- annotation(eset)
    eset <- symbol_annot(eset, gpl_name)

    return(eset)

}
