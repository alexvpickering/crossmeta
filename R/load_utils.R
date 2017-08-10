#' Download and unpack microarray supplementary files from GEO.
#'
#' Downloads and unpacks microarray supplementary files from GEO. Files are
#' stored in the supplied data directory under the GSE name.
#'
#' @param gse_names Character vector of GSE names to download.
#' @param data_dir  String specifying directory for GSE folders.
#'
#' @seealso \code{\link{load_raw}}.
#' @return NULL (for download/unpack only).
#' @export
#' @examples
#' get_raw("GSE41845")

get_raw <- function (gse_names, data_dir = getwd()) {

    for (gse_name in gse_names) {

        gse_dir <- paste(data_dir, gse_name, sep = "/")
        work_dir  <- getwd()

        # get raw data
        if (!file.exists(gse_dir)) {
            crossmeta:::getGEOSuppFiles(gse_name, baseDir = data_dir)
        }

        # untar
        tar_names <- list.files(gse_dir, pattern = "\\.tar$")
        if (length(tar_names) > 0) {
            setwd(gse_dir)
            tryCatch(utils::untar(tar_names),
                     error = function(e) setwd(work_dir))
            setwd(work_dir)
        }
        # unzip
        paths <- list.files(gse_dir, pattern = "\\.gz$",
                            full.names = TRUE, ignore.case = TRUE)
        sapply(paths, GEOquery::gunzip, overwrite = TRUE)
    }
}

# -----------

#' Load and annotate raw data downloaded from GEO.
#'
#' Loads and annotates raw data previously downloaded with \code{\link{get_raw}}.
#' Supported platforms include Affymetrix, Agilent, and Illumina.
#'
#' @import data.table
#'
#' @param gse_names Character vector of GSE names.
#' @param data_dir  String specifying directory with GSE folders.
#' @param gpl_dir   String specifying parent directory to search for previously downloaded GPL.soft files.
#' @param overwrite Do you want to overwrite saved esets from previous \code{load_raw}?
#'
#' @return List of annotated esets.
#' @export
#' @examples
#' library(lydata)
#' data_dir <- system.file("extdata", package = "lydata")
#' eset <- load_raw("GSE9601", data_dir = data_dir)


load_raw <- function(gse_names, data_dir = getwd(), gpl_dir = '..', overwrite = FALSE) {

    # no duplicates allowed (causes mismatched names/esets)
    gse_names <- unique(gse_names)

    affy_names  <- c()
    agil_names  <- c()
    illum_names <- c()

    errors <- c()
    saved <- list()
    for (gse_name in gse_names) {

        gse_dir   <- paste(data_dir, gse_name, sep = "/")
        save_name <- paste(gse_name, "eset.rds", sep = "_")
        eset_path <- list.files(gse_dir, save_name, full.names = TRUE)

        # check if saved copy
        if (length(eset_path) > 0 & overwrite == FALSE) {
            saved[[gse_name]] <- readRDS(eset_path)
            next()
        }

        # determine platform (based on filenames)
        affy  <- list.files(gse_dir, ".CEL$", ignore.case = TRUE)
        agil  <- list.files(gse_dir, "^GSM.*txt$", ignore.case = TRUE)
        illum <- list.files(gse_dir, "non.norm.*txt$|raw.*txt$|nonorm.*txt$", ignore.case = TRUE)

        # add to appropriate names vector
        if (length(affy) != 0) {
            affy_names  <- c(affy_names, gse_name)

        } else if  (length(agil) != 0) {
            agil_names  <- c(agil_names, gse_name)

        } else if (length(illum) != 0) {
            illum_names <- c(illum_names, gse_name)
        } else {
            errors <- c(errors, gse_name)
        }
    }

    # fix up saved names
    eset_names <- get_eset_names(saved, names(saved))
    saved <- unlist(saved)
    names(saved) <- eset_names

    # load non-saved esets
    affy  <- load_affy(affy_names, data_dir, gpl_dir)
    agil  <- load_agil(agil_names, data_dir, gpl_dir)
    illum <- load_illum(illum_names, data_dir, gpl_dir)


    # no raw data found
    if (length(errors) > 0) {
        message(paste0("Couldn't find raw data for: ",
                       paste(errors, collapse=", ")))
    }

    # couldn't load raw data
    if (length(affy$errors) > 0) {
        message(paste0("Couldn't load raw Affymetrix data for: ",
                       paste(affy$errors, collapse=", ")))
    }

    if (length(agil$errors) > 0) {
        message(paste0("Couldn't load raw Agilent data for: ",
                       paste(agil$errors, collapse=", ")))
    }

    if (length(illum$errors) > 0) {

        # try to fix headers
        fixed   <- fix_illum_headers(illum$errors, data_dir)
        unfixed <- setdiff(illum$errors, fixed)

        # try to load fixed
        if (!is.null(fixed)) {
            illum_fixed <- load_illum(fixed, data_dir, gpl_dir)

            if (!is.null(illum_fixed$esets))
                cat('\nFixed headers for Illumina data:', paste(names(illum_fixed$esets), collapse=", "), '\n')

            illum$esets <- c(illum$esets, illum_fixed$esets)
            unfixed     <- c(unfixed, illum_fixed$errors)
        }

        # message any unfixed
        if (length(unfixed) > 0) {

            # remove previous fix attempt
            for (gse_name in unfixed) {
                gse_dir <- file.path(data_dir, gse_name)
                prev    <- list.files(gse_dir, "_fixed\\.txt$", full.names = TRUE, ignore.case = TRUE)
                unlink(prev)
            }

            message(paste0("Couldn't load raw Illumina data for: ",
                           paste(unfixed, collapse=", "),
                           ".\nsee 'Checking Raw Illumina Data' in vignette."))
        }
    }


    return (c(saved, affy$esets, agil$esets, illum$esets))
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

merge_fdata <- function(efdat, dfdat) {

    # merge feature data from raw data and eset
    efdat$rn <- as.character(sapply(strsplit(row.names(efdat), "\\."), `[[`, 1))
    dfdat$rn <- as.character(sapply(strsplit(row.names(dfdat), "\\."), `[[`, 1))

    dt1 <- data.table(efdat, key = "rn")
    dt2 <- data.table(dfdat, key = "rn")

    fdat <- merge(unique(dt1), dt2, by = "rn", all.y = TRUE)
    fdat <- as.data.frame(fdat)
    row.names(fdat) <- make.unique(fdat$rn)
    fdat$rn <- NULL

    return(fdat[row.names(dfdat), ])
}


# Downloads bioconductor package.
#
# Used by symbol_annot to download annotation data packages from bioconductor.
#
# @param biocpack_name String specifying bioconductor package to download.
#
# @seealso \link{get_biocpack_name}, \link{symbol_annot}.
# @return NULL (downloads and loads requested package).

get_biocpack <- function(biocpack_name) {

    if (!requireNamespace(biocpack_name, quietly = TRUE)) {
        BiocInstaller::biocLite(biocpack_name, suppressUpdates = TRUE)
    }
    db <- get(biocpack_name, getNamespace(biocpack_name))
    return (db)
}


# Uses sysdata to obtain bioconductor annotation data package name.
#
# Function checks \code{gpl_bioc data.frame} from sysdata to obtain
# bioconductor annotation data package name for a given platform (GPL).

# If package name not available, user must look up manually on bioconductor
# website and then enter. Reference data from gpl table in GEOmetadb.sqlite.
#
# @param gpl_name String specifying platform for GSE.
#
# @seealso \link{symbol_annot}.
# @return String specifying name of bioconductor annotation data package.

get_biocpack_name <- function (gpl_name) {

    # get from gpl_bioc
    biocpack_name <- gpl_bioc[gpl_name, "bioc_package"]

    # manual entry if needed
    if (is.na(biocpack_name)) {
        biocpack_name <- readline(prompt = paste("Enter biocpack_name for",
                                                 paste0(gpl_name, ":")))
    }
    return (paste(biocpack_name, ".db", sep = ""))
}


# ------------------------


#' Add hgnc symbol to expression set.
#'
#' Function first maps entrez gene ids to homologous human entrez gene ids and
#' then to hgnc symbols.
#'
#' Initial entrez gene ids are obtained from bioconductor annotation data
#' packages or from feature data of supplied expression set. Homologous human
#' entrez ids are obtained from homologene and then mapped to hgnc symbols
#' using org.Hs.eg.db. Expression set is expanded if 1:many mappings occur.
#'
#' @param eset Expression set to annotate.
#' @param gse_name GSE name for eset.
#'
#' @export
#' @seealso \code{\link{load_raw}}.
#'
#' @return Expression set with hgnc symbols ("SYMBOL") and row names ("PROBE")
#'    added to fData slot.
#'
#' @examples
#' library(lydata)
#'
#' # location of raw data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # load eset
#' eset <- load_raw("GSE9601", data_dir)[[1]]
#'
#' # annotate eset (need if load_raw failed to annotate)
#' eset <- symbol_annot(eset)

symbol_annot <- function (eset, gse_name = "") {

    cat("Annotating")

    # --------------------- Map Features to Homologous Human Entrez ID


    # get map from features to entrez ids
    map <- entrez_map(eset)

    # get map from entrez ids to species gene name
    map <- gene_map(eset, map)

    # get map from entrez ids to homologous human entrez ids
    map <- merge(map, homologene, by = "ENTREZID", all.x = TRUE, sort = FALSE)

    exist_homologues <- sum(!is.na(map$ENTREZID_HS)) != 0

    if (exist_homologues) {

        # get organism name
        org_col <- grep('organism', colnames(pData(eset)))
        org     <- unique(as.character(pData(eset)[, org_col]))

        # If human, use original entrez id (useful if human platform):
        if (org == 'Homo sapiens') {
            filt <- is.na(map$ENTREZID_HS)
            map[filt, "ENTREZID_HS"] <- map$ENTREZID[filt]
        }


        # expand one-to-many mappings
        eset <- eset[map$PROBEID, ]

        # --------------------- Map Human Entrez ID to Gene Symbol

        hs <- get_biocpack("org.Hs.eg.db")

        SYMBOL <- tryCatch (
            {
                suppressMessages(
                    SYMBOL <- AnnotationDbi::select(hs, as.character(map$ENTREZID_HS),
                                                    "SYMBOL", "ENTREZID")$SYMBOL)
            },
            error = function(c) {
                if (grepl("valid keys", c$message)) {
                    message("Annotation failed.")
                    c$message <- paste0(gse_name, ": Annotation failed. Please see ",
                                        "'Loading and Annotating Data' in vignette")
                    warning(c)
                    return(NULL)
                }
            }
        )

        # --------------------- Add 'PROBE', 'Org_SYMBOL' and 'SYMBOL' to Feature Data

        # PROBE is feature names
        # added '.' to identify replicates (remove)
        PROBE <- sapply(strsplit(featureNames(eset), "\\."), `[[`, 1)

        # replaced '.' with '*' in features (reverse)
        PROBE <- gsub("*", ".", PROBE, fixed = TRUE)

        # add PROBE to fData and use for unique row names
        fData(eset)$PROBE <- PROBE
        row.names(eset) <- make.unique(PROBE)

        # add Org_SYMBOL to fData
        org_col <- colnames(map)[grep('_SYMBOL$', colnames(map))]
        if (length(org_col)) fData(eset)[[org_col]] <- map[[org_col]]

        if (!is.null(SYMBOL)) {
            # add uppercase gene symbols to fData
            fData(eset)$SYMBOL <- toupper(SYMBOL)

            # add entrez ids to fData
            fData(eset)$ENTREZID    <- map$ENTREZID
            fData(eset)$ENTREZID_HS <- map$ENTREZID_HS
        }

    } else {

        # expand one-to-many mappings
        eset <- eset[map$PROBEID, ]

        # --------------------- Add 'PROBE' and 'SYMBOL' to Feature Data

        # PROBE is feature names
        # added '.' to identify replicates (remove)
        PROBE <- sapply(strsplit(featureNames(eset), "\\."), `[[`, 1)

        # replaced '.' with '*' in features (reverse)
        PROBE <- gsub("*", ".", PROBE, fixed = TRUE)

        # add PROBE to fData and use for unique row names
        fData(eset)$PROBE <- PROBE
        row.names(eset) <- make.unique(PROBE)

        # add entrez ids to fData
        fData(eset)$ENTREZID <- map$ENTREZID

        # add Org_SYMBOL to fData
        org_col <- colnames(map)[grep('_SYMBOL$', colnames(map))]
        if (length(org_col)) fData(eset)[[org_col]] <- map[[org_col]]

    }
    return (eset)
}



# ------------------------


# Get map from eset features to entrez id.
#
#
# @param eset Expression set.
#
# @return Data.frame with columns 'PROBEID' and 'ENTREZID' which maps from eset
#    feature names to corresponding entrez gene ids.

entrez_map <- function(eset) {

    # default
    map <- data.frame(PROBEID = 1:nrow(eset), ENTREZID = NA)

    # check internal data or UI for biocpack_name
    biocpack_name <- get_biocpack_name(annotation(eset))

    # if biocpack_name not empty, use to get entrez id
    if (!biocpack_name %in% c("", ".db")) {
        biocpack <- get_biocpack(biocpack_name)

        # added '.' to identify replicates (remove)
        ID <- sapply(strsplit(featureNames(eset), "\\."), `[[`, 1)

        # replaced '.' with '*' in features (reverse)
        ID <- gsub("*", ".", ID, fixed = TRUE)
        suppressMessages(map <- AnnotationDbi::select(biocpack, ID, "ENTREZID"))

        # replace '.' with '*' (to match eset row names)
        map$PROBEID <- gsub(".", "*", map$PROBEID, fixed = TRUE)

    } else {
        # if not, check fData column for entrez id
        cols <- grep("gene_id|^gene$|entrez",
                     fvarLabels(eset), ignore.case = TRUE, value = TRUE)

        if (length(cols) != 0) {
            # pick col with most homologene matches (min 1/4)
            matches <- sapply(cols, function(col) {
                sum(fData(eset)[, col] %in% homologene$ENTREZID)
            })
            best <- names(which.max(matches))

            # expand one-to-many (row-to-entrez)
            if (matches[best] >= 0.25 * nrow(eset)) {
                entrez <- as.character(fData(eset)[, best])
                entrez <- strsplit(entrez, "\\D+")
                rn <- sapply(seq_along(entrez),
                             function(x) rep(x, length(entrez[[x]])))
                map <- data.frame(PROBEID = unlist(rn),
                                  ENTREZID = as.integer(unlist(entrez)))
            }
        }
    }
    return(map)
}


gene_map <- function(eset, map) {

    # get organism name
    org_col <- grep('organism', colnames(pData(eset)))
    org     <- unique(as.character(pData(eset)[, org_col]))

    if (!org %in% names(org_pkg)) return(map)

    # get organism package
    pkg_name <- org_pkg[[org]]
    pkg      <- get_biocpack(pkg_name)

    # organism shorthand
    org_short <- gsub('^org[.](.+?)[.].+?db$', '\\1', pkg_name)

    # map from entrezid to organism gene name
    gene_map <- suppressMessages(AnnotationDbi::select(pkg, unique(map$ENTREZID), 'SYMBOL', 'ENTREZID'))
    colnames(gene_map) <- c('ENTREZID', paste0(org_short, '_SYMBOL'))

    # merge with original map
    gene_map <- merge(gene_map, map, by='ENTREZID', all.x = TRUE)

    return(gene_map)
}

# ------------------------


# Get eset names for load_affy and load_agil.
#
# Helper function to get around issue of a single GSE having multiple platforms
# (and thus \code{getGEO} returns a list of esets). To distinguish these cases,
# the GPL platform is appended to the GSE name.
#
# @param List of annotated esets. Created by \code{load_affy} or \code{load_agil}.
# @param gse_names Character vector of GSE names for each eset.
#
# @seealso \code{\link{load_affy}}, \code{\link{load_agil}}
# @return Character vector of GSE names with GPL appended when multiple
#   platforms per GSE.

get_eset_names <- function(esets, gse_names) {
    eset_names <- c()

    for (i in seq_along(esets)) {
        # get gse name
        gse_name <- gse_names[i]

        if (length(esets[[i]]) > 1) {
            # add gpl_name to make gse_name unique
            gpl_name <- sapply(esets[[i]], annotation)
            gse_name <- paste(gse_name, gpl_name, sep = ".")
        }
        # add gse_name to eset_names
        eset_names <- c(eset_names, gse_name)
    }
    return(eset_names)
}



# GEOquery functions ----

getGEO <- function(GEO=NULL,
                   filename=NULL,
                   destdir=tempdir(),
                   GSElimits=NULL,GSEMatrix=TRUE,
                   AnnotGPL=FALSE,
                   getGPL=TRUE,
                   limit_gpls=FALSE) {
    con <- NULL
    if(!is.null(GSElimits)) {
        if(length(GSElimits)!=2) {
            stop('GSElimits should be an integer vector of length 2, like (1,10) to include GSMs 1 through 10')
        }
    }
    if(is.null(GEO) & is.null(filename)) {
        stop("You must supply either a filename of a GEO file or a GEO accession")
    }
    if(is.null(filename)) {
        GEO <- toupper(GEO)
        geotype <- toupper(substr(GEO,1,3))
        if(GSEMatrix & geotype=='GSE') {
            return(getAndParseGSEMatrices(GEO,destdir,AnnotGPL=AnnotGPL,getGPL=getGPL,limit_gpls=limit_gpls))
        }
        filename <- GEOquery::getGEOfile(GEO,destdir=destdir,AnnotGPL=AnnotGPL)
    }
    ret <- GEOquery:::parseGEO(filename,GSElimits,destdir,AnnotGPL=AnnotGPL,getGPL=getGPL)
    return(ret)
}


getAndParseGSEMatrices <- function(GEO, destdir, AnnotGPL, getGPL=TRUE, limit_gpls = FALSE) {
    GEO <- toupper(GEO)
    ## This stuff functions to get the listing of available files
    ## for a given GSE given that there may be many GSEMatrix
    ## files for a given GSE.
    stub = gsub('\\d{1,3}$','nnn',GEO,perl=TRUE)
    gdsurl <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/'
    b = crossmeta:::getDirListing(sprintf(gdsurl,stub,GEO), limit_gpls)
    message(sprintf('Found %d file(s)',length(b)))
    ret <- list()
    ## Loop over the files, returning a list, one element
    ## for each file
    for(i in 1:length(b)) {
        message(b[i])
        destfile=file.path(destdir,b[i])
        if(file.exists(destfile)) {
            message(sprintf('Using locally cached version: %s',destfile))
        } else {
            download.file(sprintf('https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s',
                                  stub,GEO,b[i]),destfile=destfile,mode='wb',
                          method=getOption('download.file.method.GEOquery'))
        }
        ret[[b[i]]] <- GEOquery:::parseGSEMatrix(destfile,destdir=destdir,AnnotGPL=AnnotGPL,getGPL=getGPL)$eset
    }
    return(ret)
}


getGEOSuppFiles <- function(GEO,makeDirectory=TRUE,baseDir=getwd()) {
    geotype <- toupper(substr(GEO,1,3))
    storedir <- baseDir
    fileinfo <- list()
    stub = gsub('\\d{1,3}$','nnn',GEO,perl=TRUE)
    if(geotype=='GSM') {
        url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/",stub,GEO)
    }
    if(geotype=='GSE') {
        url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/",stub,GEO)
    }
    if(geotype=='GPL') {
        url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platform/%s/%s/suppl/",stub,GEO)
    }
    fnames <- try(crossmeta:::getDirListing(url),silent=TRUE)
    if(inherits(fnames,'try-error')) {
        message('No supplemental files found.')
        message('Check URL manually if in doubt')
        message(url)
        return(NULL)
    }
    if(makeDirectory) {
        suppressWarnings(dir.create(storedir <- file.path(baseDir,GEO)))
    }
    for(i in fnames) {
        download.file(paste(file.path(url,i),'tool=geoquery',sep="?"),
                      destfile=file.path(storedir,i),
                      mode='wb',
                      method=getOption('download.file.method.GEOquery'))
        fileinfo[[file.path(storedir,i)]] <- file.info(file.path(storedir,i))
    }
    invisible(do.call(rbind,fileinfo))
}

getDirListing <- function(url, limit_gpls = FALSE) {
    message(url)
    # Takes a URL and returns a character vector of filenames
    a <- RCurl::getURL(url)
    ## Renaud Gaujoux reported problems behind firewall
    ## where the ftp index was converted to html content
    ## The IF statement here is his fix--harmless for the rest
    ## of us.
    if( grepl("<HTML", a, ignore.case=T) ){ # process HTML content
        doc <- XML::htmlParse(a)
        links <- XML::xpathSApply(doc, "//a/@href")
        XML::free(doc)

        # make sure link doesn't end with '/'
        links <- links[!grepl('/$', links)]

        if (limit_gpls) {
            gpls  <- stringr::str_extract(links, 'GPL\\d+')
            if (!all(is.na(gpls))) links <- links[gpls %in% row.names(gpl_bioc)]
        }

        b <- as.matrix(links)
        message('OK')
    } else { # standard processing of txt content
        tmpcon <- textConnection(a, "r")
        b <- read.table(tmpcon)
        close(tmpcon)
    }
    b <- as.character(b[,ncol(b)])
    return(b)
}
