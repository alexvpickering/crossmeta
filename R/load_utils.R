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

        gse_dir <- file.path(data_dir, gse_name)
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
#' 
load_raw <- function(gse_names, data_dir = getwd(), gpl_dir = '..', overwrite = FALSE, ensql = NULL) {

    # no duplicates allowed
    gse_names <- unique(gse_names)
    esets <- list()
    for (gse_name in gse_names) {

        gse_dir   <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep = "_")
        eset_path <- list.files(gse_dir, save_name, full.names = TRUE)

        # check if saved copy
        if (length(eset_path) > 0 & overwrite == FALSE) {
            eset <- readRDS(eset_path)
            
            # ensure in correct format
            if (methods::is(eset, 'ExpressionSet')) {
              eset <- list(eset)
              names(eset) <- gse_name
            }
            
            esets <- c(esets, eset)
            
        } else {
            esets <- c(esets, load_plat(gse_name, data_dir, gpl_dir, ensql))
        }
    }
    
    esets <- unlist(esets[!is.na(esets)])
    return (esets)
}

#' Load and pre-process raw Affymetrix, Illumina, and Agilent microarrays.
#'
#' Load raw files previously downloaded with \code{get_raw}. Used by \code{load_raw}.
#'
#' Data is normalized, SYMBOL and PROBE annotation are added to fData slot.
#'
#' @param gse_name GSE names.
#' @param data_dir String specifying directory with GSE folder.
#' @inheritParams load_raw
#' @keywords internal
#'
#' @seealso \code{\link{get_raw}} to obtain raw data.
#'
#' @return List of annotated esets, one for each platform in \code{gse_name}.

load_plat <- function(gse_name, data_dir, gpl_dir, ensql) {
  
  gse_dir <- file.path(data_dir, gse_name)
  
  # get GSEMatrix for phenotype data
  eset <- NULL
  while (is.null(eset)) 
    eset <- try(crossmeta:::getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE, getGPL = FALSE))
  
  # check if have GPLs
  gpl_names <- paste0(sapply(eset, Biobase::annotation), '.soft', collapse = "|")
  gpl_paths <- sapply(gpl_names, function(gpl_name) {
    list.files(gpl_dir, gpl_name, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)[1]
  })
  
  # copy over GPLs
  if (length(gpl_paths) > 0) file.copy(gpl_paths, gse_dir)
  
  # will use local GPL or download if couldn't copy
  eset <- NULL
  while (is.null(eset)) 
    eset <- try(crossmeta:::getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE))
  
  # make sure eset uses GSM accession for sample names (e.g. GSE10653 doesn't)
  eset <- lapply(eset, function(x) {colnames(x) <- x$geo_accession; x})
  
  # name esets
  if (length(eset) > 1) {
    names(eset) <- paste(gse_name, sapply(eset, Biobase::annotation), sep='.')
  } else {
    names(eset) <- gse_name
  }
  
  # determine manufacturer
  supported <- c('Affymetrix', 'Illumina', 'Agilent')
  gpl_names <- sapply(eset, Biobase::annotation)
  mft_names <- sapply(gpl_names, get_mft)
  
  # load eset for each platform in GSE
  esets <- list()
  nplat <- length(eset)
  for (i in 1:nplat) {
    pre_msg <- paste0('Error: ', gse_name, ' (', i, ' of ', nplat, '): ')
    
    mft_name <- mft_names[i]
    if (!mft_name %in% supported) {
      message(pre_msg, mft_name, ' platforms (', gpl_names[i], ') not supported.\n')
      next()
    }
    
    load_fun <- switch(mft_name,
                       'Affymetrix' = load_affy_plat, 
                       'Illumina' = load_illum_plat, 
                       'Agilent' = load_agil_plat)
    
    gse.gpl <- names(eset)[i]
    esets[[gse.gpl]] <- tryCatch(
      {load_fun(eset[[i]], gse_name, gse_dir, ensql)},
      error = function(e) {message(pre_msg, e$message, '\n'); return(NA)})
  }
  
  # save to disc
  save_path <- file.path(gse_dir, paste(gse_name, "eset.rds", sep = "_"))
  na.esets <- is.na(esets)
  if (!all(na.esets)) saveRDS(esets[!na.esets], save_path)
  
  return(esets)
}

get_mft <- function(gpl_name) {
  
  # check gpl_bioc first
  if (gpl_name %in% row.names(gpl_bioc)) {
    mft <- gpl_bioc[gpl_name, 'manufacturer']
    
  } else {
    mft <- scrape_mft(gpl_name)
  }
  return(mft)
}

scrape_mft <- function(gpl_name) {
  gpl_text <- crawl_acc(gpl_name)
  mft <- grep('^!Platform_manufacturer = ', gpl_text, value = TRUE)
  mft <- gsub('^!Platform_manufacturer = ', '', mft)
  
  if (grepl('illumina', mft, ignore.case = TRUE)) mft <- 'Illumina'
  else if (grepl('affymetrix', mft, ignore.case = TRUE)) mft <- 'Affymetrix'
  else if (grepl('agilent', mft, ignore.case = TRUE)) mft <- 'Agilent'
  
  return(mft)
}

crawl_acc <- function(acc_name) {
  
  # get html text for GPL page
  acc_url  <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", acc_name, "&targ=self&form=text&view=full")
  
  acc_text <- NULL
  attempt <- 1
  while(is.null(acc_text) && attempt <= 3) {
    con <- url(acc_url)
    try(acc_text <- readLines(con))
    if(is.null(acc_text)) Sys.sleep(15)
    close(con)
  }
  return(acc_text)
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
#' @inheritParams load_raw
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
    SYMBOL_9606 <- NULL
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
    
    # restore M slot if two channel array
    M <- assayDataElement(eset, 'M')

    # transfer to eset
    eset <- ExpressionSet(as.matrix(dt[, sampleNames(eset)]),
                          phenoData(eset),
                          AnnotatedDataFrame(dt[, colnames(map)]),
                          annotation = annotation(eset))
    
    suppressWarnings(Biobase::assayDataElement(eset, 'M') <- M)
    return(eset)
}


# Get map from eset features to entrez id.
#
#
# @param eset Expression set.
#
# @return Data.frame with columns 'PROBEID' and 'ENTREZID' which maps from eset
#    feature names to corresponding entrez gene ids.

entrez_map <- function(eset, ensql) {
  
    # default map
    map <- data.frame(PROBE = Biobase::featureNames(eset), ENTREZID = NA)

    # try to get ENTREZ from biocpack ----

    biocpack_names <- get_biocpack_name(Biobase::annotation(eset))
    PROBE <-  Biobase::featureNames(eset)

    # if biocpack_name not empty, use to get entrez id
    if (!biocpack_names[1] %in% c("", ".db")) {

        for(biocpack_name in biocpack_names) {
            # get biocpack and keys
            biocpack <- get_biocpack(biocpack_name)
            biockeys <- AnnotationDbi::keys(biocpack)

            # find eset column that best matches biocpack keys
            matches <- sapply(Biobase::fvarLabels(eset), function(fdatacol) {
                sum(Biobase::fData(eset)[, fdatacol] %in% biockeys) / length(PROBE)
            })
            best <- names(which.max(matches))

            if (max(matches) > 0.5) {
                # map from feature names to best to entrezid
                idmap  <- data.table::data.table(PROBE, Biobase::fData(eset)[, best])
                idmap[] <- lapply(idmap, as.character)
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
    entrezdt <- data.table::data.table(entrezdt, key = 'ENTREZID')

    map <- data.table::data.table(map, key = 'ENTREZID')
    map <- merge(map, entrezdt, by='ENTREZID', all.x = TRUE, sort = FALSE)

    # remove duplicated rows
    return(unique(map))
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
    ret <- GEOquery::parseGEO(filename,GSElimits,destdir,AnnotGPL=AnnotGPL,getGPL=getGPL)
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
            utils::download.file(sprintf('https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s',
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
        utils::download.file(paste(file.path(url,i),'tool=geoquery',sep="?"),
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
        b <- utils::read.table(tmpcon)
        close(tmpcon)
    }
    b <- as.character(b[,ncol(b)])
    return(b)
}
