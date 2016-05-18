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

load_agil <- function (gse_names, homologene, data_dir, overwrite) {

    esets <- list()
    for (gse_name in gse_names) {

        gse_dir <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep="_")
        eset_path <- list.files(gse_dir, save_name, full.names=TRUE)

        #check if saved copy
        if (length(eset_path) != 0 & overwrite == FALSE) {
            eset <- readRDS(eset_path)

        } else {
            #get GSEMatrix (for pheno data)
            eset <- GEOquery::getGEO(gse_name, destdir=gse_dir, GSEMatrix=TRUE)[[1]]

            #load non-normalized txt files and normalize
            data_paths <- list.files(gse_dir, pattern="GSM.*txt",
                                     full.names=TRUE, ignore.case=TRUE)
            data <- limma::read.maimages(data_paths,
                                         source="agilent", green.only=TRUE)
            data <- limma::neqc(data, status=data$genes$ControlType,
                                negctrl=-1, regular=0)

            #fix up sample/feature names
            colnames(data) <- stringr::str_match(colnames(data),
                                                 ".*(GSM\\d+).*")[, 2]

            row.names(data$E) <- make.unique(data$genes$ProbeName)
            row.names(data$genes) <- make.unique(data$genes$ProbeName)

            eset <- fix_agil_features(eset, data)

            #transfer to eset
            exprs(eset) <- data$E
            fData(eset) <- merge_fdata(fData(eset), data$genes)

            #add SYMBOL annotation
            eset <- symbol_annot(eset, homologene)

            #save to disc
            saveRDS(eset, file.path(gse_dir, save_name))
        }
        esets[[gse_name]] <- eset
    }
    return(esets)
}


#' Title
#'
#' @param eset
#' @param data
#'
#' @return
#' @export
#'
#' @examples
fix_agil_features <- function(eset, data) {

    #use eset fData column that best matches data probe names
    data_ids <- data$genes$ProbeName
    fData(eset)$rownames <- featureNames(eset)

    cols <- colnames(fData(eset))

    matches <- rep(0, length(cols))
    names(matches) <- cols

    for (col in cols) {
        vals <- fData(eset)[, col]
        matches[col] <- sum(vals %in% data_ids)
    }

    best <- as.character(fData(eset)[, names(which.max(matches))])

    row.names(exprs(eset)) <- make.unique(best)
    row.names(fData(eset)) <- make.unique(best)

    fData(eset)$rownames <- NULL
    return(eset)
}
