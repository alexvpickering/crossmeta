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

load_agil <- function (gse_names, data_dir, overwrite) {

    esets <- list()
    for (gse_name in gse_names) {

        gse_dir <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep = "_")
        eset_path <- list.files(gse_dir, save_name, full.names = TRUE)

        # check if saved copy
        if (length(eset_path) != 0 & overwrite == FALSE) {
            eset <- readRDS(eset_path)

        } else {
            # get GSEMatrices (for pheno/feature data)
            eset <- GEOquery::getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE)

            # load eset for each platform in GSE
            eset <- lapply(eset, load_agil_plat, gse_dir, gse_name)

            # save to disc
            saveRDS(eset, file.path(gse_dir, save_name))
        }
        esets[[gse_name]] <- eset
    }

    eset_names <- get_eset_names(esets, gse_names)
    esets <- unlist(esets)
    names(esets) <- eset_names
    return (esets)
}


# ------------------------


load_agil_plat <- function (eset, gse_dir, gse_name) {

    # get paths to raw files for samples in eset
    pattern <- paste(sampleNames(eset), ".*txt", collapse = "|", sep = "")
    data_paths <- list.files(gse_dir, pattern, full.names = TRUE, ignore.case = TRUE)

    # load non-normalized txt files and normalize
    data <- limma::read.maimages(data_paths, source = "agilent", green.only = TRUE)
    data <- limma::neqc(data, status = data$genes$ControlType, negctrl = -1, regular = 0)

    # fix up sample/feature names
    colnames(data) <- stringr::str_match(colnames(data), ".*(GSM\\d+).*")[, 2]

    row.names(data$genes) <- make.unique(data$genes$ProbeName)
    row.names(data$E)     <- make.unique(data$genes$ProbeName)

    eset <- fix_agil_features(eset, data)

    # transfer to eset
    exprs(eset) <- data$E
    fData(eset) <- merge_fdata(fData(eset), data$genes)

    # add SYMBOL annotation
    eset <- symbol_annot(eset, gse_name)

    return(eset)
}


# ------------------------


# Set eset row names to Agilent probe ids.
#
# Sets eset row names to eset feature data column that best matches raw data
# probe identifiers.
#
# Agilent raw data has probe identifiers (in data$genes$ProbeName). The row
# names of the eset GSEMatrix ('ExpressionSet' object) may not use the same
# identifiers. This is fixed by setting the eset row names to the eset feature
# data column that best matches the raw data probe identifiers.
#
# @param eset Expression set from getGEO with GSEMatrix = TRUE.
# @param data Expression set from raw data (read and processed by limma).
#
# @return eset with row names set to Agilent probe ids.

fix_agil_features <- function(eset, data) {

    # use eset fData column that best matches data probe names
    fData(eset)$rownames <- featureNames(eset)

    cols <- colnames(fData(eset))
    matches <- sapply(cols, function(col){
        sum(fData(eset)[, col] %in% data$genes$ProbeName)
    })

    best <- fData(eset)[, names(which.max(matches))]
    na   <- best %in% c("", NA)
    best <- make.unique(as.character(best))

    eset <- eset[!na, ]
    row.names(eset) <- best[!na]

    fData(eset)$rownames <- NULL
    return(eset)
}
