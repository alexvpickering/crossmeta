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

load_agil <- function (gse_names, data_dir, gpl_dir) {

    esets  <- list()
    errors <- c()
    for (gse_name in gse_names) {

        gse_dir <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep = "_")

        # get GSEMatrix (for pheno data)
        dl_methods <- c('auto', 'libcurl', 'wget', 'curl')
        eset <- NULL

        for (method in dl_methods) {
            options('download.file.method.GEOquery' = method)

            eset <- tryCatch(getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE, getGPL = FALSE),
                             error = function(e) return(NULL))

            if (inherits(eset, 'list')) break()
            Sys.sleep(5)
            if (method == 'curl') stop("Couldn't get GSEMatrix for: ", gse_names[1])
        }

        # check if have GPL
        gpl_names <- paste0(sapply(eset, annotation), '.soft', collapse = "|")
        gpl_paths <- sapply(gpl_names, function(gpl_name) {
            list.files(gpl_dir, gpl_name, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)[1]
        })

        # copy over GPL
        if (length(gpl_paths) > 0)
            file.copy(gpl_paths, gse_dir)

        # will use local GPL or download if couldn't copy
        eset <- getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE)

        # name esets
        if (length(eset) > 1) {
            names(eset) <- paste(gse_name, sapply(eset, annotation), sep='.')
        } else {
            names(eset) <- gse_name
        }

        # load eset for each platform in GSE
        eset <- lapply(eset, function(eset.gpl) {
            tryCatch(load_agil_plat(eset.gpl, gse_dir, gse_name),
                     error = function(e) NA)
        })


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


load_agil_plat <- function (eset, gse_dir, gse_name) {

    # get paths to raw files for samples in eset
    pattern <- paste(sampleNames(eset), ".*txt", collapse = "|", sep = "")
    data_paths <- list.files(gse_dir, pattern, full.names = TRUE, ignore.case = TRUE)

    # load non-normalized txt files and normalize
    data <- limma::read.maimages(data_paths, source = "agilent", green.only = TRUE)
    data <- limma::neqc(data, status = data$genes$ControlType, negctrl = -1, regular = 0)

    # fix up sample/feature names
    colnames(data) <- stringr::str_match(colnames(data), ".*(GSM\\d+).*")[, 2]
    data <- data[!is.na(data$genes$ProbeName), ]

    row.names(data$genes) <- make.unique(data$genes$ProbeName)
    row.names(data$E)     <- make.unique(data$genes$ProbeName)

    eset <- fix_agil_features(eset, data)
    eset <- eset[, colnames(data)]

    # transfer to eset
    fData(eset) <- merge_fdata(fData(eset), data$genes)
    eset <- ExpressionSet(data$E,
                          phenoData = phenoData(eset),
                          featureData = featureData(eset),
                          annotation = annotation(eset))

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
