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

load_illum <- function (gse_names, data_dir, gpl_dir, ensql) {

    esets  <- list()
    errors <- c()
    for (gse_name in gse_names) {

        gse_dir <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep = "_")


        # get GSEMatrix (for pheno dat)
        eset <- NULL
        while (is.null(eset)) {
            eset <- try(crossmeta:::getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE, getGPL = FALSE))
        }

        if (length(eset) > 1) {
            warning("Multi-platform Illumina GSEs not supported. ", gse_name)
            errors <- c(errors, gse_name)
            next
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
        eset <- NULL
        while (is.null(eset)) {
            eset <- try(crossmeta:::getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE))
        }
        eset <- tryCatch(load_illum_plat(eset[[1]], gse_name, gse_dir, ensql),
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

load_illum_plat <- function(eset, gse_name, gse_dir, ensql) {

    try(fData(eset)[fData(eset) == ""] <- NA)
    try(fData(eset)[] <- lapply(fData(eset), as.character))

    # fix header issues
    elist_paths <- list.files(gse_dir, pattern = "non.*norm.*txt$|raw.*txt$|nonorm.*txt$", full.names = TRUE, ignore.case = TRUE)
    elist_paths <- elist_paths[!grepl('fixed[.]txt$', elist_paths)]
    annotation  <- fix_illum_headers(elist_paths, eset)

    # load fixed elist paths
    elist_paths <- gsub(".txt", "_fixed.txt", elist_paths, fixed = TRUE)
    elist <- limma::read.ilmn(elist_paths, probeid = "ID_REF", annotation = annotation)


    # don't correct if already log transformed (already corrected?)
    logd <- max(elist$E, na.rm = TRUE) < 100
    if (!logd) {
        elist <- tryCatch (
            limma::neqc(elist),
            error = function(c) {
                # PMID:19068485 recommends mle and offset 50
                elist <- limma::backgroundCorrect(elist, method = "normexp",
                                                 normexp.method = "mle",
                                                 offset = 50)

                return(limma::normalizeBetweenArrays(elist, method = "quantile"))
            })
    }

    # merge eset and elist fdata
    elist <- merge_elist(eset, elist)

    if ('ID_REF' %in% colnames(elist$genes)) {
        row.names(elist$E) <- row.names(elist$genes) <- make.unique(elist$genes$ID_REF)
    } else {
        row.names(elist$E) <- row.names(elist$genes) <- NULL
    }

    # determine best sample matches
    res   <- match_samples(eset, elist)
    elist <- elist[, res$elist_order]
    eset  <- eset[, res$eset_order]
    warn  <- res$warn

    # keep gse matrix and raw elist title
    pData(eset)$title.gsemat <- pData(eset)$title
    pData(eset)$title.raw    <- colnames(elist)


    if (warn) {
        # use raw elist titles to ensure correct contrasts
        pData(eset)$title <- colnames(elist)

        # add illum colname to warn about pData
        pData(eset)$illum <- NA
    }
    colnames(elist) <- sampleNames(eset)

    # transfer elist to eset
    eset <- ExpressionSet(elist$E,
                          phenoData = phenoData(eset),
                          featureData = as(elist$genes, 'AnnotatedDataFrame'),
                          annotation = annotation(eset))

    # add SYMBOL annotation
    eset <- symbol_annot(eset, gse_name, ensql)
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

match_samples <- function(eset, elist) {

    # determine if elist has fewer samples
    data_fewer <- ncol(elist) < ncol(eset)

    # check if colnames match ----
    if (data_fewer) {
        # check if all elist colnames in eset colnames
        if (all(colnames(elist) %in% colnames(eset))) {
            cat('Illumina samples matched by column names.\n')
            return(list(elist_order = colnames(elist), eset_order = colnames(elist), warn = FALSE))
        }

    } else {
        # check if all eset colnames in elist colnames
        if (all(colnames(eset) %in% colnames(elist))) {
            cat('Illumina samples matched by column names.\n')
            return(list(elist_order = colnames(eset), eset_order = colnames(eset), warn = FALSE))
        }
    }

    # check if eset pdata col matches elist colnames ----

    if (!is.null(colnames(elist))) {

        # matrix of positions of matches for elist colnames among those for each pdata column
        matches <- sapply(pData(eset), function(col) {
            fuzzy_pmatch(colnames(elist), col)
        })

        # number of unique non NA matches for each pdata column
        nunique <- apply(matches, 2, function(match) length(unique(match[!is.na(match)])))

        # number of unique non NA matches should be the min of number of eset or pdata samples
        nmin <- min(ncol(eset), ncol(elist))
        if (any(nunique == nmin)) {
            cat('Illumina samples matched by pdata column.\n')

            # matches where satisfied
            bestcol <- names(which(nunique == nmin))[1]
            matches <- matches[, bestcol]

            # elist_order is positions where matches are not NA
            elist_order <- which(!is.na(matches))

            # eset_order is non NA matches
            eset_order <- matches[!is.na(matches)]

            return(list(elist_order = elist_order, eset_order = eset_order, warn = FALSE))
        }
    }

    # check if similarity offers unique match ----


    # make sure eset is log2 transformed
    logd <- max(exprs(eset), na.rm = TRUE) < 1000
    if (!logd) {
        exprs(eset) <- log2(exprs(eset) + abs(min(exprs(eset), na.rm = TRUE)) + 16)
    }

    # row names are the best match columns for elist and eset
    elist <- elist[!is.na(elist$genes[[elist$elistcol]]), ]
    row.names(elist) <- make.unique(elist$genes[[elist$elistcol]])

    eset <- eset[!is.na(fData(eset)[[elist$esetcol]]), ]
    row.names(eset)  <- make.unique(fData(eset)[[elist$esetcol]])

    # only include rows without missing values
    eset   <- eset[complete.cases(exprs(eset)), ]
    elist  <- elist[complete.cases(elist$E), ]

    qres   <- list()
    ngenes <- min(nrow(eset), nrow(elist))

    if (data_fewer) {
        # determine most similar eset sample for each sample in elist
        for (i in 1:ncol(elist)) {
            qsamp <- elist$E[, i]
            qres[[colnames(elist)[i]]] <- ccmap::query_drugs(qsamp, exprs(eset), sorted = FALSE, ngenes = ngenes)
        }

    } else {
        # determine most similar elist sample for each sample in eset
        for (i in 1:ncol(eset)) {
            qsamp <- exprs(eset)[, i]
            qres[[colnames(eset)[i]]] <- ccmap::query_drugs(qsamp, elist$E, sorted = FALSE, ngenes = ngenes)
        }
    }

    # eset sample to most similar elist sample
    qres <- as.data.frame(qres)
    best <- sapply(qres, which.max)

    if (length(best) == length(unique(best))) {
        cat('Illumina samples matched by similarity.\n')

        if (data_fewer) {
            elist_order <- colnames(elist)
            eset_order <- best
        } else {
            elist_order <- best
            eset_order <- colnames(eset)
        }

        return(list(elist_order = elist_order, eset_order = eset_order, warn = FALSE))

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
                elist_order <- colnames(elist)
                eset_order <- best
            } else {
                elist_order <- best
                eset_order <- colnames(eset)
            }
            return(list(elist_order = elist_order, eset_order = eset_order, warn = FALSE))

        } else {
            cat('Illumina samples not matched.\n')
            return(list(elist_order = colnames(elist), eset_order = colnames(eset), warn = TRUE))
        }
    }
}




merge_elist <- function(eset, elist) {

    if (is.null(elist$genes)) stop('Raw elist lacks feature names.')

    # get eset and elist fdata columns
    esetcols  <- fData(eset)
    elistcols <- elist$genes


    # find eset fData column that best matches elist features
    best  <- c(esetcol=NA, elistcol=NA)
    bestf <- 0

    for (i in seq_along(elistcols)) {

        elistcol <- elistcols[[i]]

        # get fraction of fdata column that has a match
        matches <- sapply(names(esetcols), function(esetcol) {
            sum(elistcol %in% esetcols[, esetcol]) / length(elistcol)
        })

        # update best
        if (max(matches) >= bestf) {
            bestf <- max(matches)
            best['elistcol'] <- names(elistcols)[i]
            best['esetcol']  <- names(matches[which.max(matches)])
        }
    }

    if (bestf > 0.3) {
        # merge eset and elist fdata columns
        esetcols  <- esetcols[!duplicated(esetcols[best['esetcol']]),, drop = FALSE]
        elistcols <- merge(elistcols, esetcols, all.x = TRUE, by.x = best['elistcol'], by.y = best['esetcol'], sort = FALSE)
        elistcols[elistcols == ""] <- NA

        elist$genes <- elistcols

        # add best info for illumina sample matching
        elist$elistcol <- best[['elistcol']]
        elist$esetcol  <- best[['esetcol']]
    }
    elist$genes[] <- lapply(elist$genes, as.character)
    return(elist)
}

