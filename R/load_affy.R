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


# ------------------------


# Load and pre-process raw Affymetrix CEL files for multiple GSEs.
#
# Load raw CEL files previously downloaded with \code{get_raw_affy}. Used by
# \code{load_raw}.

# Data is normalized, SYMBOL and PROBE annotation are added to fData slot, and
# scan dates are added to pData slot.
#
# @param gse_names Character vector of Affymetrix GSE names.
# @param data_dir String specifying directory with GSE folders.
#
# @seealso \code{\link{get_raw}} to obtain raw data.
# @return List of annotated esets (one for each unique GSE/GPL platform).

load_affy <- function (gse_names, data_dir, gpl_dir, ensql) {
    
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
        
        # check if have GPL
        gpl_names <- paste0(sapply(eset, Biobase::annotation), '.soft', collapse = "|")
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
        
        # name esets
        if (length(eset) > 1) {
            names(eset) <- paste(gse_name, sapply(eset, Biobase::annotation), sep='.')
        } else {
            names(eset) <- gse_name
        }
        
        # load eset for each platform in GSE
        eset <- lapply(eset, function(eset.gpl) {
            tryCatch(load_affy_plat(eset.gpl, gse_dir, gse_name, ensql),
                     error = function(e) {message(e$message, '\n'); return(NA)}
                     )
        })
        
        
        # save to disc
        if (!all(is.na(eset)))
            saveRDS(eset[!is.na(eset)], file.path(gse_dir, save_name))
        
        if (anyNA(eset))
            errors <- c(errors, names(eset[is.na(eset)]))
        
        
        if (!all(is.na(eset)))
            esets[[gse_name]] <- eset[!is.na(eset)]
    }
    eset_names <- get_eset_names(esets, gse_names)
    esets <- unlist(esets)
    names(esets) <- eset_names
    return (list(esets = esets, errors = errors))
}


# ------------------------


# Helper utility for load_affy.
#
# Used by load_affy to load an eset for each GPL platform in a GSE.
#
# @param eset Expression set obtained by load_affy call to getGEO.
# @param gse_dir String specifying path to GSE folder.
#
# @seealso \code{\link{load_affy}}.
# @return Annotated eset with scan_date in pData slot.

load_affy_plat <- function (eset, gse_dir, gse_name, ensql) {
    
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

# ----

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
