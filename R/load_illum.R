library(GEOquery)
library(stringr)
library(limma)

source("load_utils.R")

#------------------------

get_raw_illum <- function(gse_name, data_dir) {
    #wrapper for get_raw_illum_one
    lapply(gse_name, get_raw_illum_one, data_dir)
}


get_raw_illum_one <- function (gse_name, data_dir) {
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

open_raw_illum <- function (gse_names, data_dir) {
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

load_illum <- function (gse_name, data_dir) {
    #wrapper for load_illum_one
    eset <- lapply(gse_name, load_illum_one, data_dir)
    names(eset) <- gse_name
    return (eset)
}


load_illum_one <- function (gse_name, data_dir) {
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
    gpl_name <- annotation(eset)
    eset <- symbol_annot(eset, gpl_name)

    return (eset)
}


add_pvals <- function (eset, pvals) {
    storageMode(eset) = "environment"
    assayData(eset)[["pvals"]] = pvals
    storageMode(eset) = "lockedEnvironment"
    return (eset)
}

#------------------------