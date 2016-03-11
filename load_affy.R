library(GEOquery)
library(affy)
library(oligo)
library(stringr)
library(affxparser)

source("load_utils.R")

#------------------------

cel_dates <- function(cel_paths) {
    #IN: vector with paths to CEL files
    #OUT: vector of CEL scan dates
    scan_dates <- c()
    for (i in seq_along(cel_paths)) {
        datheader <- readCelHeader(cel_paths[i])$datheader
        scan_date <- gsub(".*([0-9]{2}/[0-9]{2}/[0-9]{2}).*", "\\1", datheader)
        scan_dates[i] <- scan_date       
    }
    return (as.factor(scan_dates))
}


#------------------------

get_raw_affy <- function(gse_name, data_dir) {
    #wrapper for get_raw_affy_one
    lapply(gse_name, get_raw_affy_one, data_dir)
}

get_raw_affy_one <- function (gse_name, data_dir) {
    #IN
    #OUT
    gse_dir <- paste(data_dir, gse_name, sep="/")
    #get raw data
    if (!file.exists(gse_dir)) {
        getGEOSuppFiles(gse_name, baseDir=data_dir)
    }
    #untar
    tar_name <- list.files(gse_dir, pattern="tar")
    untar(paste(gse_dir, tar_name, sep="/"), exdir=gse_dir)

    #unzip
    cel_paths <- list.files(gse_dir, pattern=".CEL.gz", full.names=T, ignore.case=T)
    sapply(cel_paths, gunzip, overwrite=T)   
}

#------------------------

load_affy <- function (gse_name, data_dir) {
    #wrapper for load_affy_one
    eset <- lapply(gse_name, load_affy_one, data_dir)
    names(eset) <- gse_name
    return (eset)
}

load_affy_one <- function (gse_name, data_dir) {
    #IN
    #OUT
    gse_dir <- paste(data_dir, gse_name, sep="/")

    #get GSEMatrix (for pheno data)
    eset <- getGEO(gse_name, destdir=gse_dir, GSEMatrix=T)[[1]]

    #load celfiles and normalize
    cel_paths <- list.files(gse_dir, pattern=".CEL", full.names=T, ignore.case=T)
    data <- tryCatch (
        {
        raw_data <- ReadAffy (celfile.path=gse_dir)
        affy::rma(raw_data)
        },
        warning = function(cond) {
            raw_data <- read.celfiles(cel_paths)
            return (oligo::rma(raw_data))  
        } 
    )
    #rename samples in data
    sampleNames(data) <- str_extract(sampleNames(data), "GSM[0-9]+")

    #transfer exprs from data to eset (maintaining eset sample/feature order)
    sample_order <- sampleNames(eset)
    feature_order <- featureNames(eset)
    exprs(eset) <- exprs(data)[feature_order, sample_order]

    #add scan dates to pheno data (maintaining eset sample order)
    scan_dates <- cel_dates (cel_paths)
    names(scan_dates) <- sampleNames(data)
    pData(eset)$scan_date <- scan_dates[sample_order]

    #add SYMBOL annotation
    eset <- symbol_annot(eset, gse_name)

    return (eset)
}

#------------------------

