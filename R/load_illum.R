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

load_illum <- function (gse_names, data_dir, gpl_dir, entrez_dir) {

    esets  <- list()
    errors <- c()
    for (gse_name in gse_names) {

        gse_dir <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep = "_")


        # get GSEMatrix (for pheno dat)
        eset <- crossmeta:::getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE, getGPL = FALSE)


        # check if have GPL
        gpl_names <- paste0(sapply(eset, annotation), '.soft', collapse = "|")
        gpl_paths <- sapply(gpl_names, function(gpl_name) {
            list.files(gpl_dir, gpl_name, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)[1]
        })

        # copy over GPL
        if (length(gpl_paths) > 0)
            file.copy(gpl_paths, gse_dir)

        # will use local GPL or download if couldn't copy
        eset <- crossmeta:::getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE)

        if (length(eset) > 1) {
            warning("Multi-platform Illumina GSEs not supported. ", gse_name)
            errors <- c(errors, gse_name)
            next
        }
        eset <- tryCatch(load_illum_plat(eset[[1]], gse_name, gse_dir, entrez_dir),
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

load_illum_plat <- function(eset, gse_name, gse_dir, entrez_dir) {

    # fix header issues
    data_paths <- list.files(gse_dir, pattern = "non.norm.*txt$|raw.*txt$|nonorm.*txt$", full.names = TRUE, ignore.case = TRUE)
    data_paths <- data_paths[!grepl('fixed[.]txt$', data_paths)]
    anncols    <- crossmeta:::fix_illum_headers(data_paths, eset)

    # load fixed data paths
    data_paths <- gsub(".txt", "_fixed.txt", data_paths, fixed = TRUE)

    # don't correct if already log transformed (already corrected?)
    data <- limma::read.ilmn(data_paths, probeid = "ID_REF", annotation = anncols)
    logd <- max(data$E, na.rm = TRUE) < 1000

    if (!logd) {
        data <- tryCatch (
            limma::neqc(data),
            error = function(c) {
                # PMID:19068485 recommends mle and offset 50
                data <- limma::backgroundCorrect(data, method = "normexp",
                                                 normexp.method = "mle",
                                                 offset = 50)

                return(limma::normalizeBetweenArrays(data, method = "quantile"))
            })
    }

    # fix up data feature names
    data <- crossmeta:::fix_illum_features(eset, data)

    # determine best sample matches
    res  <- crossmeta:::match_samples(eset, data)
    data <- data[, res$data_order]
    eset <- eset[, res$eset_order]
    warn <- res$warn

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

    # transfer data to eset
    if (!is.null(data$genes)) {
        dfdat <- data.frame(data$genes, row.names = row.names(data))
    } else {
        dfdat <- data.frame(row.names = row.names(data))
    }

    fData(eset) <- crossmeta:::merge_fdata(fData(eset), dfdat)
    eset <- ExpressionSet(data$E,
                          phenoData = phenoData(eset),
                          featureData = featureData(eset),
                          annotation = annotation(eset))

    # transfer pvals from data to eset
    eset <- crossmeta:::add_pvals(eset, data$other$Detection)

    # add SYMBOL annotation
    eset <- crossmeta:::symbol_annot(eset, gse_name, entrez_dir)
    return(eset)
}


# -------------------

# like base pmatch
# partial match occurs if the whole of the element of x matches any part of the element of table
fuzzy_pmatch <- function(x, table) {
    x <- tolower(x)
    table <- tolower(table)

    # first look for perfect matches
    perfect <- match(x, table)

    # is every x has a perfect match in table, return
    if (sum(is.na(perfect)) == 0) return(perfect)

    # otherwise first grep
    tomatch <- x[is.na(perfect)]
    gmatch  <- sapply(tomatch, function(val) {
        res <- grep(val, table, fixed = TRUE)[1]
        if (!length(res)) return(NA_integer_)
        return(res)
    })

    # fill in grep result where NA in perfect
    perfect[is.na(perfect)] <- gmatch[is.na(perfect)]
    return(perfect)
}

match_samples <- function(eset, data) {

    # determine if data has fewer samples
    data_fewer <- ncol(data) < ncol(eset)

    # check if colnames match ----
    if (data_fewer) {
        # check if all data colnames in eset colnames
        if (all(colnames(data) %in% colnames(eset))) {
            cat('Illumina samples matched by column names.\n')
            return(list(data_order = colnames(data), eset_order = colnames(data), warn = FALSE))
        }

    } else {
        # check if all eset colnames in data colnames
        if (all(colnames(eset) %in% colnames(data))) {
            cat('Illumina samples matched by column names.\n')
            return(list(data_order = colnames(eset), eset_order = colnames(eset), warn = FALSE))
        }
    }

    # check if eset pdata col matches data colnames ----

    if (!is.null(colnames(data))) {

        # matrix of positions of matches for data colnames among those for each pdata column
        matches <- sapply(pData(eset), function(col) {
            fuzzy_pmatch(colnames(data), col)
        })

        # number of unique non NA matches for each pdata column
        nunique <- apply(matches, 2, function(match) length(unique(match[!is.na(match)])))

        # number of unique non NA matches should be the min of number of eset or pdata samples
        nmin <- min(ncol(eset), ncol(data))
        if (any(nunique == nmin)) {
            cat('Illumina samples matched by pdata column.\n')

            # matches where satisfied
            bestcol <- names(which(nunique == nmin))[1]
            matches <- matches[, bestcol]

            # data_order is positions where matches are not NA
            data_order <- which(!is.na(matches))

            # eset_order is non NA matches
            eset_order <- matches[!is.na(matches)]

            return(list(data_order = data_order, eset_order = eset_order, warn = FALSE))
        }
    }

    # check if similarity offers unique match ----


    # make sure eset is log2 transformed
    logd <- max(exprs(eset), na.rm = TRUE) < 1000
    if (!logd) {
        exprs(eset) <- log2(exprs(eset) + abs(min(exprs(eset), na.rm = TRUE)) + 16)
    }

    qres  <- list()
    eset  <- eset[complete.cases(exprs(eset)), ]
    data  <- data[complete.cases(data$E), ]
    ngenes <- min(nrow(eset), nrow(data))

    if (data_fewer) {
        # determine most similar eset sample for each sample in data
        for (i in 1:ncol(data)) {
            qsamp <- data$E[, i]
            qres[[colnames(data)[i]]] <- ccmap::query_drugs(qsamp, exprs(eset), sorted = FALSE, ngenes = ngenes)
        }

    } else {
        # determine most similar data sample for each sample in eset
        for (i in 1:ncol(eset)) {
            qsamp <- exprs(eset)[, i]
            qres[[colnames(eset)[i]]] <- ccmap::query_drugs(qsamp, data$E, sorted = FALSE, ngenes = ngenes)
        }
    }

    # eset sample to most similar data sample
    qres <- as.data.frame(qres)
    best <- sapply(qres, which.max)

    if (length(best) == length(unique(best))) {
        cat('Illumina samples matched by similarity.\n')

        if (data_fewer) {
            data_order <- colnames(data)
            eset_order <- best
        } else {
            data_order <- best
            eset_order <- colnames(eset)
        }

        return(list(data_order = data_order, eset_order = eset_order, warn = FALSE))

    } else {
        # look for misses in non-first query results
        dups   <- unique(best[duplicated(best)])
        misses <- setdiff(1:nrow(qres), unique(best))

        n <- nrow(qres)
        for (dup in dups) {
            # query results for duplicate
            i <- 1
            qres_dup  <- qres[, best == dup]

            while (dup %in% dups & i < n) {
                ibest_dup <- sapply(qres_dup, function(col) which(col == sort(col, partial=n-i)[n-i]))

                # for each miss
                for (miss in misses) {
                    # check if one ibest is miss
                    imiss <- ibest_dup == miss

                    if (sum(imiss) == 1){
                        # if so, replace best with ibest
                        ibest_repl <- ibest_dup[imiss]
                        best[names(ibest_repl)] <- ibest_repl

                        # also update duplicates and misses
                        dups   <- best[duplicated(best)]
                        misses <- setdiff(1:nrow(qres), unique(best))

                        # if no more misses, break
                        if (!length(misses)) {
                            break()
                        }
                    }
                }
                i <- i + 1
            }
        }

        if (!length(dups)) {
            cat('Illumina samples matched by similarity using non-first ranks.\n')
            if (data_fewer) {
                data_order <- colnames(data)
                eset_order <- best
            } else {
                data_order <- best
                eset_order <- colnames(eset)
            }
            return(list(data_order = data_order, eset_order = eset_order, warn = FALSE))

        } else {
            cat('Illumina samples not matched.\n')
            return(list(data_order = colnames(data), eset_order = colnames(eset), warn = TRUE))
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

    genes_null <- is.null(data$genes)
    names_null <- is.null(row.names(data))

    if (genes_null & names_null) stop('Raw data lacks feature names.')

    datacols <- c(data$genes, list(row.names(data)))
    datacols <- datacols[!sapply(datacols, is.null)]

    fData(eset)$rn <- gsub('[.]\\d+$', '', row.names(eset))


    # find eset fData column that best matches data features
    best  <- c(esetcol=NA, datacol=NA)
    bestf <- 0

    for (i in seq_along(datacols)) {

        datacol <- datacols[[i]]

        # get fraction of fdata column that has a match
        matches <- sapply(fvarLabels(eset), function(fdatacol) {
            sum(datacol %in% fData(eset)[, fdatacol]) / length(datacol)
        })

        # update best
        if (max(matches) > bestf) {
            bestf <- max(matches)
            best['datacol'] <- names(datacols)[i]
            best['esetcol'] <- names(matches[which.max(matches)])
        }
    }

    if (bestf > 0.5) {

        # map from best datacol to eset row names
        df <- data.frame(datacols)[best['datacol']]
        ef <- fData(eset)[, best['esetcol'], drop=FALSE]
        ef$rn <- fData(eset)$rn

        map <- merge(df, ef, all.x = TRUE, by.x = best['datacol'], by.y = best['esetcol'], sort = FALSE)
        map <- unique(map)


        # positions where best datacol is in map
        rowind <- match(data$genes[, best['datacol']], map[, best['datacol']])
        row.names(data) <- make.unique(map[rowind, 'rn'])
        data <- data[!is.na(row.names(data)), ]

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


