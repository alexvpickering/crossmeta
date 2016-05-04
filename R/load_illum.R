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

load_illum <- function (gse_names, data_dir) {

    esets <- list()
    for (gse_name in gse_names) {

        gse_dir <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep="_")
        eset_path <- list.files(gse_dir, save_name, full.names=TRUE)

        #check if saved copy
        if (length(eset_path) != 0) {
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

            #transfer exprs from data to eset (maintaining eset feature order)
            feature_order <- featureNames(eset)

            #to check if sample order mismatch
            pData(eset)$title_GSEMatrix <- pData(eset)$title

            #use raw data titles to ensure correct contrasts
            pData(eset)$title <- colnames(data)
            colnames(data) <- sampleNames(eset)
            exprs(eset) <- data$E[feature_order,]

            #transfer pvals from data to eset
            pvals <- data$other$Detection[feature_order, ]
            eset <- add_pvals(eset, pvals)

            #add SYMBOL annotation
            gpl_name <- annotation(eset)
            eset <- symbol_annot(eset, gpl_name)

            #save to disc
            saveRDS(eset, file.path(gse_dir, save_name))
        }
        esets[[gse_name]] <- eset
    }
    return (esets)
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




open_raw_illum <- function (gse_names, data_dir) {

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
