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
        tar_names <- list.files(gse_dir, pattern = "\\.tar$", full.names = TRUE)
        if (length(tar_names) > 0) {
            res <- 1; try <- 0
            while (res != 0 & try != 3) {
                res <- try(utils::untar(tar_names,  exdir=gse_dir))
                if (res != 0) Sys.sleep(10); try <- try + 1
            }
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
#' @param ensql For development. Path to sqlite file with ENTREZID and SYMBOL columns created in data-raw/entrezdt.
#'
#' @return List of annotated esets.
#' @export
#' @examples
#' library(lydata)
#' data_dir <- system.file("extdata", package = "lydata")
#' eset <- load_raw("GSE9601", data_dir = data_dir)


load_raw <- function(gse_names, data_dir = getwd(), gpl_dir = '..', overwrite = FALSE, ensql = NULL) {

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
        illum <- list.files(gse_dir, "non.*norm.*txt$|raw.*txt$|nonorm.*txt$", ignore.case = TRUE)

        # add to appropriate names vector
        if (length(affy) != 0) {
            affy_names  <- c(affy_names, gse_name)
        } else if (length(illum) != 0) {
            illum_names <- c(illum_names, gse_name)
        } else if  (length(agil) != 0) {
            agil_names  <- c(agil_names, gse_name)

        }  else {
            errors <- c(errors, gse_name)
        }
    }

    # fix up saved names
    eset_names <- get_eset_names(saved, names(saved))
    saved <- unlist(saved)
    names(saved) <- eset_names

    # load non-saved esets
    affy  <- load_affy(affy_names, data_dir, gpl_dir, ensql)
    agil  <- load_agil(agil_names, data_dir, gpl_dir, ensql)
    illum <- load_illum(illum_names, data_dir, gpl_dir, ensql)


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
        message(paste0("Couldn't load raw Illumina data for: ",
                       paste(illum$errors, collapse=", ")))
    }


    return (c(saved, affy$esets, agil$esets, illum$esets))
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
        BiocManager::install(biocpack_name, update=FALSE)
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

    if (is.na(biocpack_name)) biocpack_name <- ""
    if (grepl('probeset$', biocpack_name))
        biocpack_name <- c(biocpack_name, gsub('probeset$', 'transcriptcluster', biocpack_name))

    return (paste0(biocpack_name, ".db"))
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

symbol_annot <- function (eset, gse_name = "", ensql = NULL) {
    cat("Annotating")

    # get map from features to organism entrez ids and  symbols
    map <- entrez_map(eset, ensql)

    # get map from entrez ids to homologous human entrez ids
    map <- merge(map, homologene, by = "ENTREZID", all.x = TRUE, sort = FALSE)

    exist_homologues <- sum(!is.na(map$ENTREZID_HS)) != 0

    if (exist_homologues) {

        # use original entrez id and symbols if human (taxid 9606)
        if ('SYMBOL_9606' %in% colnames(map)) {
            map$ENTREZID_HS <- map$ENTREZID
            map$SYMBOL_9606 <- map$SYMBOL <- toupper(map$SYMBOL_9606)
        } else {
            # map human entrez to gene symbol
            map$SYMBOL <- toupper(hs[map$ENTREZID_HS, SYMBOL_9606])
        }
    }

    # merge map and exprs
    PROBE <- featureNames(eset)
    dt <- data.table(exprs(eset), PROBE, key='PROBE')
    dt <- merge(unique(dt), map, by = 'PROBE', all.x=TRUE, sort=FALSE)
    dt <- data.frame(dt, row.names = make.unique(dt$PROBE))

    # transfer to eset
    eset <- ExpressionSet(as.matrix(dt[, sampleNames(eset)]),
                          phenoData(eset),
                          AnnotatedDataFrame(dt[, colnames(map)]),
                          annotation = annotation(eset))
    return(eset)
}



# ------------------------


# Get map from eset features to entrez id.
#
#
# @param eset Expression set.
#
# @return Data.frame with columns 'PROBEID' and 'ENTREZID' which maps from eset
#    feature names to corresponding entrez gene ids.

entrez_map <- function(eset, ensql) {

    # default map
    map <- data.frame(PROBE = featureNames(eset), ENTREZID = NA)

    # try to get ENTREZ from biocpack ----

    biocpack_names <- get_biocpack_name(annotation(eset))
    PROBE <-  featureNames(eset)

    # if biocpack_name not empty, use to get entrez id
    if (!biocpack_names[1] %in% c("", ".db")) {

        for(biocpack_name in biocpack_names) {
            # get biocpack and keys
            biocpack <- get_biocpack(biocpack_name)
            biockeys <- AnnotationDbi::keys(biocpack)

            # find eset column that best matches biocpack keys
            matches <- sapply(fvarLabels(eset), function(fdatacol) {
                sum(fData(eset)[, fdatacol] %in% biockeys) / length(PROBE)
            })
            best <- names(which.max(matches))

            if (max(matches) > 0.5) {
                # map from feature names to best to entrezid
                idmap  <- data.table(PROBE, fData(eset)[, best])
                colnames(idmap) <- c('PROBE', 'PROBEID')

                # get entrezid map and merge
                map <- suppressMessages(AnnotationDbi::select(biocpack, idmap$PROBEID, "ENTREZID"))
                map <- merge(idmap, map, by='PROBEID', all.x=TRUE, sort=FALSE, allow.cartesian=TRUE)
                map <- unique(map[, c('PROBE', 'ENTREZID')])
            }
        }

        # use biocpack species
        org   <- AnnotationDbi::species(biocpack)
        taxid <- org_taxid[org]

    } else {
        # use pdata species
        # TODO: there is usually taxid column
        org_col <- grep('organism', colnames(pData(eset)))[1]
        org     <- unique(as.character(pData(eset)[, org_col]))
        taxid   <- org_taxid[org]
    }

    # fraction of probes with entrez ids
    fentrez <- sum(!is.na(map$ENTREZID)) / nrow(map)

    orgpack_name   <- org_pkg[org]
    orgpack_exists <- !is.na(orgpack_name)
    if (orgpack_exists) orgpack <- get_biocpack(orgpack_name)


    # try to get ENTREZ from fdata ----
    if (fentrez < 0.2) {

        # check fData column for organism entrez id
        entrezcols <- grep("gene_id|^gene$|entrez",
                           fvarLabels(eset), ignore.case = TRUE, value = TRUE)



        if (length(entrezcols) != 0) {
            # entrez ids
            if (is.null(ensql)) {
                enids <- AnnotationDbi::keys(orgpack, 'ENTREZID')
            } else {
                db <- DBI::dbConnect(RSQLite::SQLite(), ensql)
            }

            # pick col with most organism entrez id matches (min 1/4)
            matches <- sapply(entrezcols, function(entrezcol) {

                colvals <- as.character(fData(eset)[, entrezcol])
                colvals <- unlist(strsplit(colvals, "\\D+"))

                if(is.null(ensql)) {
                    sum(colvals %in% enids) / length(PROBE)

                } else {
                    enids <- paste(shQuote(colvals), collapse=', ')
                    statement <- paste0('select count(ENTREZID) from ensql where ENTREZID in (', enids, ')')
                    sum(DBI::dbGetQuery(db, statement)[,1] / length(PROBE))
                }
            })
            best <- names(which.max(matches))

            # expand one-to-many (row-to-entrez)
            if (matches[best] >= 0.25) {
                entrez <- as.character(fData(eset)[, best])
                entrez <- strsplit(entrez, "\\D+")
                rn <- sapply(seq_along(entrez),
                             function(x) rep(x, length(entrez[[x]])))

                map <- data.frame(PROBE = PROBE[unlist(rn)],
                                  ENTREZID = unlist(entrez), stringsAsFactors = FALSE)
            }
            if (exists('db')) DBI::dbDisconnect(db)
        }
    }

    # update fraction of probes with entrez ids
    fentrez <- sum(!is.na(map$ENTREZID)) / nrow(map)

    # try to get ENTREZ from orgpack ----
    if (fentrez < 0.2) {

        # find fdata column that best matches organism column
        fdatacols <- fvarLabels(eset)

        keytypes <- AnnotationDbi::keytypes(orgpack)
        orgcols  <- keytypes[keytypes %in% c('ACCNUM', 'ALIAS', 'ENSEMBL', 'ENSEMBLPROT', 'ENSEMBLTRANS', 'REFSEQ', 'SYMBOL')]

        best  <- c(orgcol=NA, fdatacol=NA)
        bestf <- 0

        for (orgcol in orgcols) {
            # get keys of type orgcol
            orgkeys <- AnnotationDbi::keys(orgpack, orgcol)

            # get fraction of fdata column that has a match
            matches <- sapply(fdatacols, function(fdatacol) {
                sum(fData(eset)[, fdatacol] %in% orgkeys) / length(PROBE)
            })

            # update best
            if (max(matches) > bestf) {
                bestf <- max(matches)
                best['orgcol']   <- orgcol
                best['fdatacol'] <- names(matches[which.max(matches)])
            }
        }

        # min 1/5 match to use
        if (bestf > 0.2) {

            # map from best fdatacol to entrez id using best orgcol
            orgkeys <- AnnotationDbi::keys(orgpack, best['orgcol'])
            suppressMessages(map <- AnnotationDbi::select(orgpack, orgkeys, "ENTREZID", best['orgcol']))

            # map from probe id to best fdatacol
            idmap <- data.frame(PROBE, fData(eset)[, best['fdatacol']], stringsAsFactors = FALSE)
            names(idmap) <- c('PROBEID', best['orgcol'])

            # merge
            map <- merge(idmap, map, by = best['orgcol'], all.x = TRUE, sort = FALSE)
            map <- map[, c('PROBEID', 'ENTREZID')]
        }
    }

    colnames(map) <- c('PROBE', 'ENTREZID')

    # get SYMBOL_taxid ----

    if (is.null(ensql)) {
        # map from organism entrez id to organism symbol
        entrezdt <- suppressMessages(AnnotationDbi::select(orgpack, unique(map$ENTREZID), 'SYMBOL', 'ENTREZID'))

    } else {
        db <- DBI::dbConnect(RSQLite::SQLite(), ensql)
        map_enids <- paste(shQuote(map$ENTREZID), collapse=', ')
        statement <- paste0('select * from ensql where ENTREZID in (', map_enids, ')')
        entrezdt  <- DBI::dbGetQuery(db, statement)
        DBI::dbDisconnect(db)
    }

    # merge with original map
    colnames(entrezdt) <- c('ENTREZID', paste0('SYMBOL_', taxid))
    entrezdt <- data.table(entrezdt, key = 'ENTREZID')

    map <- data.table(map, key = 'ENTREZID')
    map <- merge(map, entrezdt, by='ENTREZID', all.x = TRUE, sort = FALSE)

    # remove duplicated rows
    return(unique(map))
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
                   getGPL=TRUE) {
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
            return(getAndParseGSEMatrices(GEO,destdir,AnnotGPL=AnnotGPL,getGPL=getGPL))
        }
        filename <- GEOquery::getGEOfile(GEO,destdir=destdir,AnnotGPL=AnnotGPL)
    }
    ret <- GEOquery:::parseGEO(filename,GSElimits,destdir,AnnotGPL=AnnotGPL,getGPL=getGPL)
    return(ret)
}


getAndParseGSEMatrices <- function(GEO, destdir, AnnotGPL, getGPL=TRUE) {
    GEO <- toupper(GEO)
    ## This stuff functions to get the listing of available files
    ## for a given GSE given that there may be many GSEMatrix
    ## files for a given GSE.
    stub = gsub('\\d{1,3}$','nnn',GEO,perl=TRUE)
    gdsurl <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/'
    b = crossmeta:::getDirListing(sprintf(gdsurl,stub,GEO))
    message(sprintf('Found %d file(s)',length(b)))
    ret <- list()
    ## Loop over the files, returning a list, one element
    ## for each file
    for(i in 1:length(b)) {
        message(b[i])
        destfile=list.files(destdir, gsub('.gz$', '', b[i]), full.names = TRUE)[1]
        if(file.exists(destfile)) {
            message(sprintf('Using locally cached version: %s',destfile))
        } else {
            destfile=file.path(destdir,b[i])
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

getDirListing <- function(url) {
    message(url)
    # Takes a URL and returns a character vector of filenames
    a <- RCurl::getURL(url)
    ## Renaud Gaujoux reported problems behind firewall
    ## where the ftp index was converted to html content
    ## The IF statement here is his fix--harmless for the rest
    ## of us.
    if( grepl("<HTML", a, ignore.case=TRUE) ){ # process HTML content
        doc <- XML::htmlParse(a)
        links <- XML::xpathSApply(doc, "//a/@href")
        XML::free(doc)

        # make sure link doesn't end with '/'
        links <- links[!grepl('/$', links)]

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
