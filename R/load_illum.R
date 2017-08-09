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
    logd <- max(data$E, na.rm = TRUE) < 1000

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

    # fix up data feature names
    data <- fix_illum_features(eset, data)

    # determine best sample matches
    match_res <- match_samples(eset, data)
    data <- match_res$data
    warn <- match_res$warn

    # keep gse matrix and raw data title
    pData(eset)$title.gsemat <- pData(eset)$title
    pData(eset)$title.raw    <- colnames(data)


    if (warn) {
        # use raw data titles to ensure correct contrasts
        pData(eset)$title <- colnames(data)

        # add illum colname to warn about pData
        pData(eset)$illum <- NA
    }
    colnames(data) <- sampleNames(eset)


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

match_samples <- function(eset, data) {

    # check if colnames match
    if (all(colnames(eset) %in% colnames(data))) {
        cat('Illumina samples matched by column names.\n')
        return(list(data = data[, colnames(eset)], warn = FALSE))
    }

    # check if any pheno data columns match
    ismatch <- sapply(pData(eset), function(col) {
        all(tolower(colnames(data)) %in% gsub('.+?[;:] ', '', tolower(col)))
        })

    if (any(ismatch)) {
        eset_order <- as.character(pData(eset)[[which(ismatch)[1]]])
        col_order  <- match(tolower(colnames(data)), gsub('.+?: ', '', tolower(eset_order)))

        cat('Illumina samples matched by pData column.\n')
        return(list(data = data[, col_order], warn = FALSE))
    }

    # make sure eset is log2 transformed
    logd <- max(exprs(eset), na.rm = TRUE) < 1000
    if (!logd) {
        exprs(eset) <- log2(exprs(eset) + abs(min(exprs(eset), na.rm = TRUE)) + 16)
    }

    # determine most similar data sample for each sample in eset
    qres <- list()

    for (i in 1:ncol(eset)) {

        # query sample
        qsamp <- exprs(eset)[, i]
        qres[[colnames(eset)[i]]] <- ccmap::query_drugs(qsamp, data$E, sorted = FALSE, ngenes = nrow(eset))

    }
    # eset sample to most similar data sample
    qres <- as.data.frame(qres)
    best <- sapply(qres, which.max)

    if (length(best) == length(unique(best))) {
        cat('Illumina samples matched by similarity.\n')
        return(list(data = data[, best], warn = FALSE))

    } else {
        # cat('checking non-first query results.\n')
        # look for misses in non-first query results
        dups   <- unique(best[duplicated(best)])
        misses <- setdiff(1:nrow(qres), unique(best))

        # cat('Duplicated:', paste0(row.names(qres)[dups], collapse = ', '), '\n')
        # cat('Missing:', paste0(row.names(qres)[misses], collapse = ', '), '\n\n')

        n <- nrow(qres)

        for (dup in dups) {

            # cat('Checking duplicate:', row.names(qres)[dup], '\n')

            # query results for duplicate
            i <- 1
            qres_dup  <- qres[, best == dup]

            while (dup %in% dups & i < n) {

                # cat('Checking rank:', i+1, '\n')
                ibest_dup <- sapply(qres_dup, function(col) which(col == sort(col, partial=n-i)[n-i]))

                # for each miss
                for (miss in misses) {

                    # cat('Checking if', row.names(qres)[miss], 'is in rank', i+1, 'exactly once.\n')

                    # check if one ibest is miss
                    imiss <- ibest_dup == miss

                    if (sum(imiss) == 1){

                        # cat('yes it is!\n')

                        # if so, replace best with ibest
                        ibest_repl <- ibest_dup[imiss]
                        best[names(ibest_repl)] <- ibest_repl

                        # also update duplicates and misses
                        dups   <- best[duplicated(best)]
                        misses <- setdiff(1:nrow(qres), unique(best))

                        # if no more misses, break
                        if (!length(misses)) {
                            # cat('no more missing!\n')
                            break()
                        }
                    }
                }
                i <- i + 1
            }
        }

        if (!length(dups)) {
            cat('Illumina samples matched by similarity using non-first ranks.\n')
            return(list(data = data[, best], warn = TRUE))
        } else {
            cat('Illumina samples not matched.\n')
            return(list(data = data, warn = TRUE))
        }
    }
}



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


