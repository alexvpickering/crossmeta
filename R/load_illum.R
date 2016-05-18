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

load_illum <- function (gse_names, homologene, data_dir, overwrite) {

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
            data_paths <- list.files(gse_dir, pattern="non.norm.*txt",
                                     full.names=TRUE, ignore.case=TRUE)
            data <- limma::read.ilmn(data_paths, probeid="ID_REF")
            data <- tryCatch (
                limma::neqc(data),
                error = function(cond) {
                    data <- limma::backgroundCorrect(data, method="normexp",
                                                     normexp.method="rma",
                                                     offset=16)

                    return(limma::normalizeBetweenArrays(data, method="quantile"))
                })

            #to check if sample order mismatch
            pData(eset)$title_GSEMatrix <- pData(eset)$title

            #use raw data titles to ensure correct contrasts
            pData(eset)$title <- colnames(data)
            colnames(data) <- sampleNames(eset)

            #fix up feature names and transfer data
            data <- fix_illum_features(eset, data)

            exprs(eset) <- data$E[row.names(data$E) != "", ]
            fData(eset) <- merge_fdata(fData(eset),
                                       data.frame(row.names = make.unique(row.names(data))))
            fData(eset) <- fData(eset)[featureNames(eset), ]

            #transfer pvals from data to eset
            pvals <- data$other$Detection
            eset <- add_pvals(eset, pvals)

            #add SYMBOL annotation
            eset <- symbol_annot(eset, homologene)

            #save to disc
            saveRDS(eset, file.path(gse_dir, save_name))
        }
        esets[[gse_name]] <- eset
    }
    return (esets)
}


# used to fix raw data feature names
fix_illum_features <- function(eset, data) {

    fData(eset)$rownames <- featureNames(eset)

    dataf <- data.frame(dfn = row.names(data), stringsAsFactors = FALSE)
    esetf <- data.frame(efn = featureNames(eset), stringsAsFactors = FALSE)

    # find eset fData column that best matches data features
    cols <- colnames(fData(eset))

    matches <- rep(0, length(cols))
    names(matches) <- cols

    for (col in cols) {
        vals <- fData(eset)[, col]
        matches[col] <- sum(dataf$dfn %in% vals)
    }
    best <- names(which.max(matches))
    esetf$best <- fData(eset)[, best]

    # get map from data features -> best eset match -> eset features
    map <- merge(dataf, esetf, by.x="dfn", by.y="best", sort = FALSE)

    # expand 1:many map and set data row names to eset features
    data <- data[map$dfn, ]

    row.names(data) <- make.unique(map$efn)
    return(data)
}


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




open_raw_illum <- function (gse_names, data_dir=getwd()) {

    #OUT: names of successfully formated (probeid = "ID_REF",
    #                                     exprs = "AVG_Signal-sample_name",
    #                                     pvals = "Detection-sample_name",
    #                                     sep = "\t")
    out_names <- gse_names
    for (i in seq_along(gse_names)) {
        #get data paths
        gse_dir <- paste(data_dir, gse_names[i], sep="/")
        data_paths <- list.files(gse_dir, pattern="non.norm.*txt",
                                 full.names=TRUE, ignore.case=TRUE)
        data_paths <- c(data_paths, list.files(gse_dir, pattern=".xls",
                                               full.names=TRUE))
        #open data file
        for (j in seq_along(data_paths)) system2("xdg-open", data_paths[j])

        #check success
        success <- tcltk::tk_select.list(choices = c("Yes", "No"),
                                         title = paste(gse_names[i],
                                                       "formated successfully?"))
        #remove unsuccessful
        if (success == "No") out_names <- setdiff(out_names, gse_names[i])
    }
    return(out_names)
}
