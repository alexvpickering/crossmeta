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

load_affy <- function (gse_names, data_dir, gpl_dir) {


    esets  <- list()
    errors <- c()
    for (gse_name in gse_names) {

        gse_dir <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep = "_")

        # get GSEMatrix (for pheno dat)
        dl_methods <- c('auto', 'libcurl', 'wget', 'curl')
        eset <- NULL

        for (method in dl_methods) {
            options('download.file.method.GEOquery' = method)

            eset <- tryCatch(getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE, getGPL = FALSE, limit_gpls = TRUE),
                             error = function(e) {
                                 return(NULL)
                                 })

            if (inherits(eset, 'list')) break()
            Sys.sleep(5)
            if (method == 'curl') stop("Couldn't get GSEMatrix for: ", gse_names[1])
        }

        # check if have GPL
        gpl_names <- paste0(sapply(eset, annotation), '.soft')
        gpl_paths <- sapply(gpl_names, function(gpl_name) {
            list.files(gpl_dir, gpl_name, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)[1]
        })

        # copy over GPL
        if (length(gpl_paths) > 0)
            file.copy(gpl_paths, gse_dir)

        # will use local GPL or download if couldn't copy
        eset <- getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE, limit_gpls = TRUE)

        # name esets
        if (length(eset) > 1) {
            names(eset) <- paste(gse_name, sapply(eset, annotation), sep='.')
        } else {
            names(eset) <- gse_name
        }

        # load eset for each platform in GSE
        eset <- lapply(eset, function(eset.gpl) {
            tryCatch(load_affy_plat(eset.gpl, gse_dir, gse_name),
                     error = function(e) NA)
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

load_affy_plat <- function (eset, gse_dir, gse_name) {

    sample_names <- sampleNames(eset)
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

    data <- tryCatch (
        {
            raw_data <- affy::ReadAffy(filenames = cel_paths)
            affy::rma(raw_data)
        },
        warning = function(c) {
            # is the warning to use oligo/xps?
            if (grepl("oligo", c$message)) {
                raw_data <- oligo::read.celfiles(cel_paths)
                return (oligo::rma(raw_data))
                # if not, use affy
            } else {
                raw_data <- affy::ReadAffy(filenames = cel_paths)
                return(affy::rma(raw_data))
            }
        },
        error = function(c) {
            raw_data <- oligo::read.celfiles(cel_paths)
            return (oligo::rma(raw_data))
        }
    )
    # rename samples in data
    sampleNames(data) <- stringr::str_extract(sampleNames(data), "GSM[0-9]+")

    # reserve '.' for duplicate features
    row.names(eset) <- gsub(".", "*", row.names(eset), fixed = TRUE)
    row.names(data) <- gsub(".", "*", row.names(data), fixed = TRUE)

    # transfer exprs from data to eset (maintaining eset sample order)
    sample_order <- sampleNames(eset)[sampleNames(eset) %in% sampleNames(data)]

    eset <- eset[, sample_order]
    data <- data[, sample_order]
    assayData(eset) <- assayData(data)

    # transfer merged fdata
    fData(eset) <- merge_fdata(fData(eset), fData(data))

    # add scan dates to pheno data (maintaining eset sample order)
    scan_dates <- cel_dates(cel_paths)
    names(scan_dates) <- sampleNames(data)
    pData(eset)$scan_date <- scan_dates[sample_order]

    # add SYMBOL annotation
    eset <- symbol_annot(eset, gse_name)

    return(eset)
}
