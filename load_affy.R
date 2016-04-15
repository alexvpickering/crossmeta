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
get_eset_names <- function(esets, gse_names) {
    #used by load_affy to name esets

    eset_names <- c()

    for (i in seq_along(esets)) {
        #get gse name
        gse_name <- gse_names[i]
        
        if (length(esets[[i]]) > 1) {
            #add gpl_name to make gse_name unique
            gpl_name <- sapply(esets[[i]], annotation)
            gse_name <- paste(gse_name, gpl_name, sep=".")
        }
        #add gse_name to eset_names
        eset_names <- c(eset_names, gse_name)
    }
    return(eset_names)
}


load_affy <- function (gse_names, data_dir) {
    #wrapper for load_affy_one
    esets <- lapply(gse_names, load_affy_one, data_dir)
    eset_names <- get_eset_names(esets, gse_names)

    esets <- unlist(esets)
    names(esets) <- eset_names
    return (esets)
}

load_affy_one <- function (gse_name, data_dir) {
    #IN
    #OUT
    gse_dir <- paste(data_dir, gse_name, sep="/")

    #get GSEMatrix (for pheno data)
    esets <- getGEO(gse_name, destdir=gse_dir, GSEMatrix=T)

    #load eset for each platform in GSE
    esets <- lapply(esets,load_affy_plat, gse_dir)

    return(esets)
}


load_affy_plat <- function (eset, gse_dir) {
    #used by load_affy_one to load eset for each platform in GSE
    
    sample_names <- sampleNames(eset)
    pattern <- paste(".*", sample_names, ".*CEL", collapse="|", sep="")
    
    cel_paths <- list.files(gse_dir, pattern, full.names=T, ignore.case=T)
    data <- tryCatch (
        {
        raw_data <- ReadAffy (celfile.path=gse_dir)
        affy::rma(raw_data)
        },
        warning = function(cond) {
            raw_data <- read.celfiles(cel_paths)
            return (oligo::rma(raw_data))  
        },
        error = function(cond) {
            raw_data <- read.celfiles(cel_paths)
            return (oligo::rma(raw_data))  
        } 
    )
    #rename samples in data
    sampleNames(data) <- str_extract(sampleNames(data), "GSM[0-9]+")

    #transfer exprs from data to eset (maintaining eset sample/feature order)
    sample_order <- sampleNames(eset)
    feature_order <- featureNames(eset)
    eset <- tryCatch (
        {
        exprs(eset) <- exprs(data)[feature_order, sample_order]
        eset
        },
        #if features don't match: also transfer featureData from data to eset
        error = function(cond) {
            exprs(eset) <- exprs(data)[, sample_order]
            fData(eset) <- fData(data)
            return(eset)
        } 
    )

    #add scan dates to pheno data (maintaining eset sample order)
    scan_dates <- cel_dates (cel_paths)
    names(scan_dates) <- sampleNames(data)
    pData(eset)$scan_date <- scan_dates[sample_order]

    #add SYMBOL annotation
    gpl_name <- annotation(eset)
    eset <- symbol_annot(eset, gpl_name)
    
    return(eset)
    
}


#------------------------

