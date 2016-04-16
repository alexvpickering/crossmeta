library(GEOquery)
library(stringr)

source("load_utils.R")

#------------

get_raw_agil <- function (gse_names, data_dir) {
    
    for (gse_name in gse_names) {
    
        gse_dir <- paste(data_dir, gse_name, sep="/")

        if (!file.exists(gse_dir)) {
            getGEOSuppFiles(gse_name, baseDir=data_dir)
        }
        #untar
        tar_name <- list.files(gse_dir, pattern="tar")
        untar(paste(gse_dir, tar_name, sep="/"), exdir=gse_dir)

        #unzip
        paths <- list.files(gse_dir, pattern=".gz", full.names=T, ignore.case=T)
        sapply(paths, gunzip, overwrite=T) 
    }        
}

#------------

load_agil <- function (gse_name, data_dir) {
    #wrapper for load_illum_one
    eset <- lapply(gse_name, load_agil_one, data_dir)
    names(eset) <- gse_name
    return (eset)
}

#-----------

load_agil_one <- function (gse_name, data_dir) {
    
    gse_dir <- paste(data_dir, gse_name, sep="/")

    #get GSEMatrix (for pheno data)
    eset <- getGEO(gse_name, destdir=gse_dir, GSEMatrix=T)[[1]]

    #load non-normalized txt files and normalize
    data_paths <- list.files(gse_dir, pattern="GSM.*txt", full.names=T, ignore.case=T)
    data <- read.maimages(data_paths, source="agilent", green.only=TRUE)
    data <- neqc(data, status=data$genes$ControlType, negctrl=-1, regular=0)
    
    #fix up sample/feature names
    colnames(data) <- str_match(colnames(data), ".*(GSM\\d+).*")[, 2]
    row.names(data$E)     <- data$genes$ProbeName
    data$genes <- data$genes[!duplicated(data$genes$ProbeName),]
    row.names(data$genes) <- data$genes$ProbeName
    
    #transfer to eset
    exprs(eset) <- data$E
    fData(eset) <- data$genes

    #add SYMBOL annotation
    gpl_name <- annotation(eset)
    eset <- symbol_annot(eset, gpl_name)
    
    return(eset)
}