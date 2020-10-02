# Load and pre-process raw Agilent files.
#
# Load raw txt files previously downloaded with \code{get_raw_agil}. Used by
# \code{load_raw}.
#
# Data is normalized and SYMBOL and PROBE annotation are added to fData slot.
#
# @param gse_names Character vector of Agilent GSE names.
# @param data_dir String specifying directory with GSE folders.
#
# @seealso \code{\link{get_raw}} to obtain raw data.
# @return List of annotated esets.

load_agil <- function (gse_names, data_dir, gpl_dir, ensql) {
  
  esets  <- list()
  errors <- c()
  for (gse_name in gse_names) {
    
    gse_dir <- file.path(data_dir, gse_name)
    save_name <- paste(gse_name, "eset.rds", sep = "_")
    
    # get GSEMatrix (for pheno dat)
    eset <- NULL
    while (is.null(eset)) {
      eset <- tryCatch(crossmeta:::getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE, getGPL = FALSE),
                       error=function(e) return(NULL))
    }
    
    
    # if _ch2 in pdata => dual channel
    ch2 <- sapply(eset, function(x) any(grepl('ch2', colnames(pData(x)))))
    
    
    # check if have GPL
    gpl_names <- paste0(sapply(eset, annotation), '.soft', collapse = "|")
    gpl_paths <- sapply(gpl_names, function(gpl_name) {
      list.files(gpl_dir, gpl_name, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)[1]
    })
    
    # copy over GPL
    if (length(gpl_paths) > 0)
      file.copy(gpl_paths, gse_dir)
    
    # will use local GPL or download if couldn't copy
    eset <- NULL
    while (is.null(eset)) {
      eset <- try(crossmeta:::getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE))
    }
    
    # name esets
    if (length(eset) > 1) {
      names(eset) <- paste(gse_name, sapply(eset, annotation), sep='.')
    } else {
      names(eset) <- gse_name
    }
    
    # load eset for each platform in GSE
    for (i in seq_along(eset)) {
      eset[[i]] <- tryCatch(
        {load_agil_plat(eset[[i]], ch2[i], gse_dir, gse_name, ensql)},
        error = function(e) {message(e$message, '\n'); return(NA)})
    }
    
    # save to disc
    if (!all(is.na(eset)))
      saveRDS(eset[!is.na(eset)], file.path(gse_dir, save_name))
    
    if (anyNA(eset))
      errors <- c(errors, names(eset[is.na(eset)]))
    
    
    if (!all(is.na(eset)))
      esets[[gse_name]] <- eset[!is.na(eset)]
  }
  
  eset_names <- get_eset_names(esets, gse_names)
  esets <- unlist(esets)
  names(esets) <- eset_names
  return (list(esets = esets, errors = errors))
}


# ------------------------


#' Load Agilent raw data
#'
#' @param eset ExpressionSet from \link{getGEO}
#' @param ch2 Boolean indicating in two-channel array
#' @param gse_dir Direction with Agilent raw data
#' @param gse_name Accession name for \code{eset}.
#' @inheritParams load_raw
#'
#' @return ExpressionSet
#' 
load_agil_plat <- function (eset, ch2, gse_dir, gse_name, ensql) {
  
  try(fData(eset)[fData(eset) == ""] <- NA)
  try(fData(eset)[] <- lapply(fData(eset), as.character))
  
  # get paths to raw files for samples in eset
  pattern <- paste('^', sampleNames(eset), ".*(txt|gpr)$", collapse = "|", sep = "")
  elist_paths <- list.files(gse_dir, pattern, full.names = TRUE, ignore.case = TRUE)
  
  # if multiple with same GSM, take first
  gsm_names  <- stringr::str_extract(elist_paths, "GSM[0-9]+")
  elist_paths <- elist_paths[!duplicated(gsm_names)]
  
  # support for genepix .gpr files
  is.genepix <- grepl('.gpr$', elist_paths[1])
  source <- ifelse(is.genepix, 'genepix', 'agilent')
  
  # load non-normalized txt files and normalize
  elist <- tryCatch(
    {limma::read.maimages(elist_paths, source = source, green.only = !ch2)},
    error = function(e) {
      # determine source of error
      output <- capture.output(tryCatch(
        limma::read.maimages(elist_paths, source = source, green.only = !ch2),
        error = function(e) NULL))
      
      # retry with error excluded
      exclude     <- which(elist_paths == gsub('^Read| ', '', output[length(output)-1])) + 1
      elist_paths <- elist_paths[-exclude]
      limma::read.maimages(elist_paths, source = source, green.only = !ch2)
    })
  
  if (ch2) {
    # follows 'Separate Channel Analysis of Two-Color Data' to make as if single channel
    elist <- limma::backgroundCorrect(elist, method="normexp", offset=50)
    elist <- limma::normalizeWithinArrays(elist, method="loess")
    elist <- limma::normalizeBetweenArrays(elist, method="Aquantile")
    
  } else  {
    elist <- limma::neqc(elist, status = elist$genes$ControlType, negctrl = -1, regular = 0)
  }
  
  # fix up sample names
  colnames(elist) <- stringr::str_match(colnames(elist), ".*(GSM\\d+).*")[, 2]
  eset <- eset[, colnames(elist)]
  
  # merge elist and eset feature data
  elist <- merge_elist(eset, elist)
  
  if (ch2) {
    # not sure 'ID' is general fill-in for 'ProbeName'?
    elist <- elist[!is.na(elist$genes$ID), ]
    row.names(elist$M) <- row.names(elist$A) <- row.names(elist$genes) <- make.unique(elist$genes$ID)
    
    # transfer to eset
    # A: average log-2 expression
    eset <- ExpressionSet(elist$A,
                          phenoData = phenoData(eset),
                          featureData = as(elist$genes, 'AnnotatedDataFrame'),
                          annotation = annotation(eset))
    
    # M: log-2 expression ratios
    # used in fit_ebayes
    Biobase::assayDataElement(eset, 'M') <- elist$M
    
  } else {
    row.names(elist$E) <- row.names(elist$genes) <- make.unique(elist$genes$ProbeName)
    
    # transfer to eset
    eset <- ExpressionSet(elist$E,
                          phenoData = phenoData(eset),
                          featureData = as(elist$genes, 'AnnotatedDataFrame'),
                          annotation = annotation(eset))
  }
  
  
  # add SYMBOL annotation
  eset <- symbol_annot(eset, gse_name, ensql)
  return(eset)
}

