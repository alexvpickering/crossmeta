library(GEOquery)
library(stringr)

source("shared_utils.R")

#------------------------

get_raw <- function(gse_name, data_dir) {
    #wrapper for get_raw_one
    lapply(gse_name, get_raw_one, data_dir)
}


get_raw_one <- function (gse_name, data_dir) {
    #IN
    #OUT
    gse_dir <- paste(data_dir, gse_name, sep="/")
    
    #get raw data
    if (!file.exists(gse_dir)) {
        getGEOSuppFiles(gse_name, baseDir=data_dir)
    }   
    #unzip
    gz_paths <- list.files(gse_dir, pattern=".gz", full.names=T, ignore.case=T)
    sapply(gz_paths, gunzip)  
}

#------------------------

open_raw <- function (gse_names, data_dir) {
    #IN: opens non-normalized illumina data (change to below format)
    #OUT: names of successfully formated (probeid = "ID_REF",
    #                                     exprs = "AVG_Signal-sample_name",
    #                                     pvals = "Detection-sample_name",
    #                                     sep = ",") 
    out_names <- gse_names
    for (i in seq_along(gse_names)) {
        #open data
        gse_dir <- paste(data_dir, gse_names[i], sep="/")
        data_paths <- list.files(gse_dir, pattern="non.norm.*txt", full.names=T)
        data_paths <- c(data_paths, list.files(gse_dir, pattern=".xls", full.names=T))
        for (i in seq_along(data_paths)) {
            system2("xdg-open", data_paths[i])
        }
        #check success
        success <- tk_select.list(choices = c("Yes", "No"),
                                  title = paste(gse_names[i], "formated successfully?"))
        #remove unsuccessful
        if (success == "No") {
            out_names <- setdiff(out_names, gse_names[i])  
        }
    }
    return(out_names)
}

#------------------------

load_eset <- function (gse_name, data_dir) {
    #wrapper for load_eset_one
    eset <- lapply(gse_name, load_eset_one, data_dir)
    names(eset) <- gse_name
    return (eset)
}


load_eset_one <- function (gse_name, data_dir) {
    #IN: non-normalized illumina data (format: probeid = "ID_REF",
    #                                          exprs = "AVG_Signal-sample_name",
    #                                          pvals = "Detection-sample_name",
    #                                          sep = ",")
    #OUT
    gse_dir <- paste(data_dir, gse_name, sep="/")
    
    #get GSEMatrix (for pheno data)
    eset <- getGEO(gse_name, destdir=gse_dir, GSEMatrix=T)[[1]]
    
    #load non-normalized txt files and normalize
    data_paths <- list.files(gse_dir, pattern="non.norm.*txt", full.names=T)
    data <- read.ilmn(data_paths, probeid="ID_REF")
    data <- neqc(data)
    
    #transfer exprs from data to eset (maintaining eset feature order)
    feature_order <- featureNames(eset)
    pData(eset)$title_GSEMatrix <- pData(eset)$title  #to check if sample order mismatch
    pData(eset)$title <- colnames(data)  #use raw data titles to ensure correct contrasts
    colnames(data) <- sampleNames(eset)
    exprs(eset) <- data$E[feature_order,]

    #transfer pvals from data to eset
    pvals <- data$other$Detection[feature_order, ]
    eset <- add_pvals(eset, pvals)

    #add SYMBOL annotation
    eset <- symbol_annot(eset, gse_name)

    return (eset)
}


add_pvals <- function (eset, pvals) {
    storageMode(eset) = "environment"
    assayData(eset)[["pvals"]] = pvals
    storageMode(eset) = "lockedEnvironment"
    return (eset)
}

#------------------------

diff_expr <- function (eset, data_dir) {
    #wrapper for diff_expr_one
    diff_expr <- mapply(diff_expr_one, eset, names(eset), data_dir, SIMPLIFY=F)
    names(diff_expr) <- names(eset)
    return (diff_expr)
}


diff_expr_one <- function (eset, name, data_dir) {
    #INPUT:
    #OUTPUT:
    gse_dir <- paste(data_dir, name, sep="/")

    #add blocking variables
    choices <- paste(sampleNames(eset), pData(eset)$title)
    block <- add_blocking(eset, choices)

    #select contrasts
    cons <- add_contrasts(block$eset, name)

    #differential expression
    setup <- diff_setup(cons$eset, cons$levels, block$names)

    eset_filt <- pval_filt(cons$eset, cons$samples, cons$levels)

    anal <- diff_anal(eset_filt, cons$contrasts, cons$levels,
                      setup$mod, setup$modsv, setup$svobj, gse_dir, name)

    return (anal)
}


pval_filt <- function (eset, selected_samples, group_levels) {

    filt_num <- round (length(selected_samples) / length(group_levels))
    pval_cnt <- rowSums (assayData(eset)$pvals < 0.05)
    eset <- eset[pval_cnt >= filt_num,]

    return (eset)
}