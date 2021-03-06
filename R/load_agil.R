#' Load Agilent raw data
#'
#' @param eset ExpressionSet from \link{getGEO}.
#' @param gse_name Accession name for \code{eset}.
#' @param gse_dir Direction with Agilent raw data.
#' @inheritParams load_raw
#'
#' @return ExpressionSet
#' 
load_agil_plat <- function (eset, gse_name, gse_dir, ensql) {
  
  # if _ch2 in pdata => dual channel
  ch2 <- any(grepl('ch2', colnames(pData(eset))))
  
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
      output <- utils::capture.output(tryCatch(
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
    probe <- 'ID'  # not sure 'ID' is general fill-in for 'ProbeName'?
    
  } else  {
    elist <- limma::neqc(elist, status = elist$genes$ControlType, negctrl = -1, regular = 0)
    probe <- 'ProbeName'
  }
  
  # fix up sample names
  colnames(elist) <- stringr::str_match(colnames(elist), ".*(GSM\\d+).*")[, 2]
  eset <- eset[, colnames(elist)]
  
  # merge elist and eset feature data
  elist <- merge_elist(eset, elist)
  
  # remove empty rows and use probe as row names
  elist <- elist[!is.na(elist$genes[[probe]]), ]
  row.names(elist) <- row.names(elist$genes) <- make.unique(elist$genes[[probe]])
  
  # convert limma object to eset
  eset <- to_eset(elist, eset)
  
  # add SYMBOL annotation
  eset <- symbol_annot(eset, gse_name, ensql)
  return(eset)
}

#' Convert limma object to ExpressionSet
#'
#' @param object an EList of MAList object containing expression data.
#' @param eset ExpressionSet from \link{getGEO}. Used for annotation.
#'
#' @return ExpressionSet using expression data from \code{object} and annotation from \code{eset}.
#' @keywords internal
#'
to_eset <- function(object, eset) {
  
  E <- switch(class(object),
              'EList' = object$E,
              'MAList' = exprs.MA(object))
  
  pdata <- switch(class(object),
              'EList' = phenoData(eset),
              'MAList' = phenoData.ch2(eset))
 
  eset <- ExpressionSet(E,
                        phenoData = pdata,
                        featureData = methods::as(object$genes, 'AnnotatedDataFrame'),
                        annotation = annotation(eset))
  
  return(eset)
}

#' Extract Log-Expression Matrix from MAList
#' 
#' Converts M and A-values to log-expression values. The output matrix will have two columns for each array, in the order all red then all green.
#' Adapted from \link[limma]{plotDensities.MAList} instead of \link[limma]{exprs.MA} so that order is same as \link{phenoData.ch2}.
#'
#' @inheritParams limma::exprs.MA
#'
#' @return A numeric matrix with twice the columns of the input.
#' @export
#'
exprs.MA <- function(MA) {
  y <- cbind(MA$A + MA$M/2, MA$A - MA$M/2)
  colnames(y) <- c(paste0(colnames(MA), '_red'), paste0(colnames(MA), '_green'))
  return(y)
}

#' Covert expression values to MAList
#'
#' @param y Expression values from two-channel agilent array in order all red then all green.
#'
#' @return MAList
#' @keywords internal
#'
#' @examples
#' 
#' A <- matrix(rnorm(100), ncol = 5)
#' M <- matrix(rnorm(100), ncol = 5)
#' MA <- new('MAList', list(M=M, A=A))
#' colnames(MA) <- letters[1:5]
#' 
#' y <- exprs.MA(MA)
#' MA2 <- crossmeta:::to_ma(y)
#' all.equal(MA, MA2)
#' 
to_ma <- function(y) {
  narray <- ncol(y)/2
  R <- y[, 1:narray]
  G <- y[, (narray+1):ncol(y)]
  M <- R - G
  A <- (R + G)/2
  MA <- methods::new("MAList", list(M=M, A=A))
  colnames(MA) <- gsub('_red|_green', '', colnames(MA))
  return(MA)
}

#' Construct AnnotatedDataFrame from Two-Channel ExpressionSet
#' 
#'
#' @param eset ExpressionSet with \code{pData} for two-channel Agilent array.
#'
#' @return AnnotatedDataFrame with twice as many rows as \code{eset}, one for each channel of each array in order all red then all green.
#' @export
#'
phenoData.ch2 <- function(eset) {
  
  pdata <- Biobase::pData(eset)
  cols <- colnames(pdata)
  stopifnot(all(c('label_ch1', 'label_ch2') %in% cols))
  
  # concat rows for non-channel columns
  not.ch <- !grepl('[_:]ch1|[_:]ch2', cols)
  pdata.ch2 <- rbind(pdata[, not.ch], pdata[, not.ch])
  
  # get columns in common between channels
  ch1.col <- grep('[_:]ch1', cols)
  ch2.col <- grep('[_:]ch2', cols)
  
  ch1_names <- gsub('[_:]ch1(.+)?$', '', cols[ch1.col])
  ch2_names <- gsub('[_:]ch2(.+)?$', '', cols[ch2.col])
  
  # exclude duplicates
  # also explicitly characteristics_ch1 and characteristics_ch2 as are often not same category
  dups <- c(ch1_names[duplicated(ch1_names)],
            ch2_names[duplicated(ch2_names)])
  dups <- unique(c('characteristics', dups))
  
  common <- setdiff(intersect(ch1_names, ch2_names), dups)
  
  # add common as one column per channel with corresponding value
  ch1.col.common <- ch1.col[ch1_names %in% common]
  ch2.col.common <- ch2.col[ch2_names %in% common]
  pdata.common.ch1 <- pdata[, ch1.col.common, drop = FALSE]
  pdata.common.ch2 <- pdata[, ch2.col.common, drop = FALSE]
  
  colnames(pdata.common.ch1) <- colnames(pdata.common.ch2) <- common
  pdata.common <- rbind(pdata.common.ch1, pdata.common.ch2)
  
  pdata.ch2 <- cbind(pdata.ch2, pdata.common)
  
  # add columns that are unique to each channel
  ch1.col.unique <- setdiff(ch1.col, ch1.col.common)
  ch2.col.unique <- setdiff(ch2.col, ch2.col.common)
  
  pdata.unique <- rbind.fill(pdata[, ch1.col.unique, drop = FALSE],
                             pdata[, ch2.col.unique, drop = FALSE])
  
  
  if (ncol(pdata.unique)) pdata.ch2 <- cbind(pdata.ch2, pdata.unique)
  
  # remove white-space from label (e.g. GSE10653)
  pdata.ch2$label <- gsub(' ', '', pdata.ch2$label)
  
  # all red then all green
  stopifnot(setequal(c('Cy5', 'Cy3'), pdata.ch2$label))
  is.red <- pdata.ch2$label == 'Cy5'
  is.green <- pdata.ch2$label == 'Cy3'
  pdata.ch2 <- rbind(pdata.ch2[is.red,, drop = FALSE], pdata.ch2[is.green,, drop = FALSE])
  
  # fix row names
  suffix <- rep(c('_red', '_green'), each = length(is.red)/2)
  row.names(pdata.ch2) <- paste0(pdata$geo_accession, suffix)
  
  
  return(AnnotatedDataFrame(pdata.ch2))
}


rbind.fill <- function(..., dfs=list(...)) {
  ns <- unique(unlist(sapply(dfs, names)))
  do.call(rbind, lapply(dfs, function(x) {
    for(n in ns[! ns %in% names(x)]) {x[[n]] <- NA}; x }))
}

