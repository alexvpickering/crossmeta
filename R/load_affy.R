# Extract scan date from Affymetrix CEL file.
#
# Useful for defining sample batches for \code{ComBat}. No longer
# used by crossmeta, which discovers nuissance variables using \code{sva}.
#
# @param cel_paths Charactor vector specifying full paths to CEL files.
#
# @seealso \code{\link{ComBat}}
# @return Factor vector of CEL scan dates.

cel_dates <- function(cel_paths) {
    
    scan_dates <- c()
    for (i in seq_along(cel_paths)) {
        datheader <- affxparser::readCelHeader(cel_paths[i])$datheader
        scan_date <- gsub(".*([0-9]{2}/[0-9]{2}/[0-9]{2}).*", "\\1", datheader)
        scan_dates[i] <- scan_date
    }
    return (as.factor(scan_dates))
}



# Affymetrix loader for load_plat.
#
# Used by load_plat to load an eset for each GPL platform in a GSE.
#
# @param eset Expression set obtained by load_plat call to getGEO.
# @param gse_dir String specifying path to GSE folder.
#
# @seealso \code{\link{load_plat}}.
# @return Annotated eset with scan_date in pData slot.

load_affy_plat <- function (eset, gse_name, gse_dir, ensql) {
    
    try(Biobase::fData(eset)[Biobase::fData(eset) == ""] <- NA)
    try(Biobase::fData(eset)[] <- lapply(Biobase::fData(eset), as.character))
    
    sample_names <- Biobase::sampleNames(eset)
    pattern <- paste(sample_names, ".*CEL$", collapse = "|", sep = "")
    
    cel_paths <- tryCatch (
        list.files(gse_dir, pattern, full.names = TRUE, ignore.case = TRUE),
        
        # list.files fails if too many files
        error = function(c) {
            n <- length(sample_names)
            p1 <- paste(sample_names[1:(n/2)], ".*CEL$", collapse = "|", sep = "")
            p2 <- paste(sample_names[(n/2+1):n], ".*CEL$", collapse = "|", sep = "")
            
            pth1 <- list.files(gse_dir, p1, full.names = TRUE, ignore.case = TRUE)
            pth2 <- list.files(gse_dir, p2, full.names = TRUE, ignore.case = TRUE)
            
            return(c(pth1, pth2))
        }
    )
    
    # if multiple with same GSM, take first
    gsm_names  <- stringr::str_extract(cel_paths, "GSM[0-9]+")
    cel_paths <- cel_paths[!duplicated(gsm_names)]
    
    abatch <- tryCatch (
        {
            raw_abatch <- affy::ReadAffy(filenames = cel_paths)
            affy::rma(raw_abatch)
        },
        warning = function(c) {
            # is the warning to use oligo/xps?
            if (grepl("oligo", c$message)) {
                raw_abatch <- oligo::read.celfiles(cel_paths)
                return(oligo::rma(raw_abatch))
                # if not, use affy
            } else {
                raw_abatch <- affy::ReadAffy(filenames = cel_paths)
                return(affy::rma(raw_abatch))
            }
        },
        error = function(c) {
            # is the error a corrupted CEL?
            if (grepl('corrupted', c$message)) {
                # exclude corrupted and try again
                corrupted <- stringr::str_extract(c$message, 'GSM\\d+')
                cel_paths <- cel_paths[!grepl(corrupted, cel_paths)]
                raw_abatch  <- affy::ReadAffy(filenames = cel_paths)
                return(affy::rma(raw_abatch))
                
            } else {
                raw_abatch <- tryCatch(oligo::read.celfiles(cel_paths),
                                       error = function(d) {
                                           if (grepl('pd.huex.1.0.st.v1', d$message))
                                               return(oligo::read.celfiles(cel_paths, pkgname = 'pd.huex.1.0.st.v2'))
                                           if (grepl('pd.hugene.2.0.st.v1', d$message))
                                               return(oligo::read.celfiles(cel_paths, pkgname = 'pd.hugene.2.0.st'))
                                           if (grepl('pd.mogene.2.0.st.v1', d$message))
                                               return(oligo::read.celfiles(cel_paths, pkgname = 'pd.mogene.2.0.st'))
                                       })
                return (oligo::rma(raw_abatch))
            }
        }
    )
    # rename samples in abatch
    Biobase::sampleNames(abatch) <- stringr::str_extract(Biobase::sampleNames(abatch), "GSM[0-9]+")
    
    # transfer exprs from abatch to eset (maintaining eset sample order)
    sample_order <- Biobase::sampleNames(eset)[Biobase::sampleNames(eset) %in% Biobase::sampleNames(abatch)]
    
    eset <- eset[, sample_order]
    abatch <- abatch[, sample_order]
    Biobase::assayData(eset) <- Biobase::assayData(abatch)
    
    # transfer merged fdata
    Biobase::fData(eset) <- merge_fdata(Biobase::fData(eset), Biobase::fData(abatch))
    
    # add SYMBOL annotation
    eset <- symbol_annot(eset, gse_name, ensql)
    
    return(eset)
}



# Merge feature data from eset and raw data.
#
# Merges feature data from eset GSEMatrix and raw data.
#
# Data.frames are merged on feature names. Result has same row names as raw
# feature data. NAs are added where eset feature data is missing a feature
# in raw data.
#
# @param efdat data.frame with eset feature data (fData).
# @param dfdat data.frame with raw feature data (varies).
#
# @return Data.frame with all columns present in efdat and dfdat. Row names
#    are same as dfdat.

merge_fdata <- function(eset_fdata, abatch_fdata) {
    
    # merge feature data from raw data and eset
    abatch_fdata$ID <- row.names(abatch_fdata)
    abatch_fdata <- merge(abatch_fdata, eset_fdata, by = "ID", all.x = TRUE, sort = FALSE)
    row.names(abatch_fdata) <- make.unique(abatch_fdata$ID)
    abatch_fdata[] <- lapply(abatch_fdata, as.character)
    return(abatch_fdata)
}
