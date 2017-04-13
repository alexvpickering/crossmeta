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

load_illum <- function (gse_names, data_dir, gpl_dir) {

    esets  <- list()
    errors <- c()
    for (gse_name in gse_names) {

        gse_dir <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep = "_")

        # get GSEMatrix (for pheno data)
        eset <- GEOquery::getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE, getGPL = FALSE)

        # check if have GPL
        gpl_names <- paste0(sapply(eset, annotation), '.soft', collapse = "|")
        gpl_paths <- sapply(gpl_names, function(gpl_name) {
            list.files(gpl_dir, gpl_name, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)[1]
        })

        # copy over GPL
        if (length(gpl_paths) > 0)
            file.copy(gpl_paths, gse_dir)

        # will use local GPL or download if couldn't copy
        eset <- GEOquery::getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE)

        if (length(eset) > 1) {
            warning("Multi-platform Illumina GSEs not supported. ", gse_name)
            errors <- c(errors, gse_name)
            next
        }
        eset <- tryCatch(load_illum_plat(eset[[1]], gse_name, gse_dir),
                         error = function(e) NULL)

        # save to disc
        if (!is.null(eset)) {
            saveRDS(eset, file.path(gse_dir, save_name))
        } else {
            errors <- c(errors, gse_name)
        }

        esets[[gse_name]] <- eset
    }
    eset_names <- get_eset_names(esets, gse_names)
    esets <- unlist(esets)
    names(esets) <- eset_names
    return (list(esets = esets, errors = errors))
}

# -------------------

# Helper utility for load_illum.
#
# Used by load_illum to load an eset.
#
# @param eset Expression set obtained by load_illum call to getGEO.
# @param gse_name String specifying GSE name.
# @param gse_dir String specifying path to GSE folder.
#
# @seealso \code{\link{load_illum}}.
# @return Annotated eset.

load_illum_plat <- function(eset, gse_name, gse_dir) {

    # load header fixed if available
    data_paths <- list.files(gse_dir, pattern = "_fixed\\.txt$",
                             full.names = TRUE, ignore.case = TRUE)

    # otherwise load non-normalized txt files and normalize
    if (!length(data_paths))
        data_paths <- list.files(gse_dir, pattern = "non.norm.*txt$|raw.*txt$|nonorm.*txt$",
                                 full.names = TRUE, ignore.case = TRUE)


    # don't correct if already log transformed
    data <- limma::read.ilmn(data_paths, probeid = "ID_REF")
    logd <- max(data$E) > 1000

    if (!logd) {
        data <- tryCatch (
            limma::neqc(data),
            error = function(cond) {
                data <- limma::backgroundCorrect(data, method = "normexp",
                                                 normexp.method = "rma",
                                                 offset = 16)

                return(limma::normalizeBetweenArrays(data, method = "quantile"))
            })
    }

    # use raw data titles to ensure correct contrasts
    pData(eset)$title <- colnames(data)
    colnames(data) <- sampleNames(eset)

    # add illum colname to warn about pData
    pData(eset)$illum <- NA

    # fix up data feature names
    data <- fix_illum_features(eset, data)

    # reserve '.' for replicates
    row.names(eset) <- gsub(".", "*", row.names(eset), fixed = TRUE)

    # transfer data to eset
    fData(eset) <- merge_fdata(fData(eset), data.frame(row.names = row.names(data)))
    eset <- ExpressionSet(data$E,
                          phenoData = phenoData(eset),
                          featureData = featureData(eset),
                          annotation = annotation(eset))

    # transfer pvals from data to eset
    eset <- add_pvals(eset, data$other$Detection)

    # add SYMBOL annotation
    eset <- symbol_annot(eset, gse_name)
    return(eset)
}


# -------------------


# Convert Illumina data features to eset features.
#
# Maps from raw data row names to eset row names.
#
# Illumina raw data ('EList' object) has no feature information other than
# row names. These rownames may not match those of the eset GSEMatrix
# ('ExpressionSet' object). This is fixed by mapping from raw data row names to
# eset row names through the eset feature data column that best matches the
# raw data row names.
#
# @param eset ExpressionSet from call to getGEO with GSEMatrix = TRUE.
# @param data EList from loading raw Illumina data then limma::neqc.
#
# @return data with rownames matching eset rownames.

fix_illum_features <- function(eset, data) {

    data <- data[row.names(data) != "", ]
    fData(eset)$rownames <- featureNames(eset)

    df <- data.frame(dfn = row.names(data), stringsAsFactors = FALSE)
    ef <- data.frame(efn = featureNames(eset), stringsAsFactors = FALSE)

    # find eset fData column that best matches data features
    cols <- colnames(fData(eset))

    matches <- sapply(cols, function(col) {
        sum(df$dfn %in% fData(eset)[, col]) / dim(eset)[[1]]
    })

    if (max(matches) > 0.5) {
        ef$best <- fData(eset)[, names(which.max(matches))]

        # get map from data features -> best eset match -> eset features
        map <- merge(df, ef, by.x = "dfn", by.y = "best", sort = FALSE)

        # expand 1:many map
        data <- data[map$dfn, ]

        # reserve periods in data row names to indicate replicates
        map$efn <- gsub(".", "*", map$efn, fixed = TRUE)
        row.names(data) <- make.unique(map$efn)
    }

    return(data)
}


# -------------------


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


# -------------------


#' Open raw Illumina microarray files.
#'
#' Helper function to open raw Illumina microarray files in order to check that
#' they are formatted correctly. For details on correct format, please see
#' 'Checking Raw Illumina Data' in vignette.
#'
#' @param gse_names Character vector of Illumina GSE names to open.
#' @param data_dir String specifying directory with GSE folders.
#'
#' @return Character vector of successfully formated Illumina GSE names.
#' @export
#'
#' @examples
#' library(lydata)
#'
#' # Illumina GSE names
#' illum_names <- c("GSE50841", "GSE34817", "GSE29689")
#'
#' # location of raw data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # open raw data files with default text editor
#' # open_raw_illum(illum_names)

open_raw_illum <- function (gse_names, data_dir = getwd()) {

    out_names <- gse_names
    for (i in seq_along(gse_names)) {
        # get data paths
        gse_dir <- paste(data_dir, gse_names[i], sep = "/")
        data_paths <- list.files(gse_dir, pattern = "non.norm.*txt",
                                 full.names = TRUE, ignore.case = TRUE)
        data_paths <- c(data_paths, list.files(gse_dir, pattern = ".xls",
                                               full.names = TRUE))
        # open data file
        for (j in seq_along(data_paths)) pander::openFileInOS(data_paths[j])

        # check success
        success <- tcltk::tk_select.list(choices = c("Yes", "No"),
                                         title = paste(gse_names[i],
                                                       "formated successfully?"))
        # remove unsuccessful
        if (success == "No") out_names <- setdiff(out_names, gse_names[i])
    }
    return(out_names)
}


