#' Download and unpack microarray supplementary files from GEO.
#'
#' Downloads and unpacks microarray supplementary files from GEO.
#' Files are stored in the supplied data directory under the GSE name.
#'
#' @importFrom GEOquery getGEOSuppFiles gunzip
#' @importFrom utils untar
#'
#' @param gse_names Character vector of GSE names to download.
#' @param data_dir base data directory (a folder for each GSE will be created
#'        here).
#' @export
#' @seealso
#' @return NULL (for download/unpack only).
#' @examples \dontrun{
#'
#'}
#'

get_raw <- function (gse_names, data_dir) {

  for (gse_name in gse_names) {

    gse_dir <- paste(data_dir, gse_name, sep="/")
    #get raw data
    if (!file.exists(gse_dir)) {
      getGEOSuppFiles(gse_name, baseDir=data_dir)
    }
    #untar
    tar_names <- list.files(gse_dir, pattern="tar")
    if (length(tar_names) > 0) {
      untar(paste(gse_dir, tar_names, sep="/"), exdir=gse_dir)
    }
    #unzip
    paths <- list.files(gse_dir, pattern=".gz", full.names=T, ignore.case=T)
    sapply(paths, gunzip, overwrite=T)
  }
}


#' Downloads bioconductor package.
#'
#' Used by symbol_annot to download annotation data packages from bioconductor.
#'
#' @importFrom BiocInstaller biocLite
#'
#' @param biocpack_name name of bioconductor package to download.
#'
#' @seealso \link{get_biocpack_name}, \link{symbol_annot}.
#' @return NULL (downloads and loads requested package).
#' @examples \dontrun{
#'
#'}

get_biocpack <- function(biocpack_name) {
    #IN:
    #OUT:

    if (!require(biocpack_name , character.only=T)) {
        source("https://bioconductor.org/biocLite.R")
        biocLite(biocpack_name)
        require(biocpack_name , character.only=T)
    }
}


#' Queries GEOmetadb sqlite file to obtain bioconductor annotation data package name.
#'
#' Function queries GEOmetadb.sqlite in order to obtain bioconductor annotation
#' data package name for a given platform (GPL). If bioconductor package name
#' not available, user must look up manually on bioconductor website and then
#' enter.
#'
#' @importFrom DBI dbConnect dbGetQuery dbDisconnect
#' @importFrom RSQLite SQLite
#' @importFrom GEOmetadb getSQLiteFile
#'
#' @param gpl_name platform name of GSE.
#'
#' @seealso \link{symbol_annot}.
#' @return name of bioconductor package.
#' @examples \dontrun{
#'
#'}

get_biocpack_name <- function (gpl_name) {

    #connect to GEOmetadb database
    meta_path <- file.path(getwd(), "GEOmetadb.sqlite")

    if (!file.exists(meta_path)) {
        stop("GEOmetadb.sqlite not in working directory. See
              ?GEOmetadb::getSQLiteFile")
    }
    con <- dbConnect(SQLite(), meta_path)

    #query metadata database
    query <- paste(
        "SELECT gpl.bioc_package, gpl.title",
        "FROM gpl",
        "WHERE gpl.gpl=", shQuote(gpl_name))

    biocpack_name <- dbGetQuery(con, query)[, "bioc_package"]
    if (length(biocpack_name) == 0) biocpack_name <- NA

    #manual entry if needed
    if (is.na(biocpack_name)) {
        title <- dbGetQuery(con, query)[, "title"]
        biocpack_name <- inputs(box1="Enter biocpack_name", def1=title)
    }

    dbDisconnect(con)
    return (paste(biocpack_name, ".db", sep=""))
}


#------------------------


#' Adds gene name (SYMBOL) column to fData slot of expression set.
#'
#' Function uses platform (GPL) to identify and download corresponding
#' bioconductor annotation data package. Gene name ("SYMBOL") is then added to
#' featureData slot of supplied eset.
#'
#' @importFrom Biobase featureNames fvarLabels fData
#' @importFrom tcltk tk_select.list
#'
#' @param eset expression set object to add annotation data to.
#' @param gpl_name platform name used to look up bioconductor
#'        annotation data package.
#'
#' @seealso \link{load_affy}, \link{load_illum}, \link{load_agil}.
#' @return eset with gene names (in "SYMBOL" column of fData slot).
#' @examples \dontrun{
#'
#'}

symbol_annot <- function (eset, gpl_name) {
    biocpack_name <- get_biocpack_name(gpl_name)
    SYMBOL <- NA  #value if no biocpack/selection

    if (biocpack_name != "") {
        #get map for SYMBOL from biocpack
        get_biocpack(biocpack_name)
        ID <- featureNames(eset)
        map <- AnnotationDbi::select(get(biocpack_name), ID, "SYMBOL")
        eset <- eset[map$PROBEID,]  #expands one-to-many mappings
        SYMBOL <- map$SYMBOL
    } else {
        #TODO: remove NAs and subset eset
        #try fData column
        choices <- setdiff(fvarLabels(eset), "SYMBOL")
        column <- tk_select.list(choices, title="select SYMBOL column")
        if (column != "") {
            SYMBOL <- fData(eset)[, column]
        }
    }
    fData(eset)$PROBE <- featureNames(eset)
    fData(eset)$SYMBOL <- SYMBOL
    return (eset)
}


#' Keep only common features for a list of esets.
#'
#' Used prior to differential expression analysis (diff_expr) to remove
#' non-common features (either gene or probe name). The subsequent meta-analysis
#' uses common features only, so eliminating non-common features reduces the number
#' of comparisons made during differential expression analysis (increases power).
#'
#' @importFrom Biobase fData
#'
#' @param esets list of annotated expression sets to commonize.
#' @param annot feature names to commonize by. Either "SYMBOL" or "PROBE".
#'
#' @export
#' @seealso \link{load_affy}, \link{load_illum}, and \link{load_agil} to obtain
#'          list of annotated esets.
#'          \link{diff_expr} to run differential expression subsequent to
#'          \code{commonize}.
#' @return list of esets with where all features (probe or gene) are common.
#' @examples \dontrun{
#'
#'}

commonize <- function(esets, annot="SYMBOL") {

    esets <- lapply(esets, function(eset) {

      #make gene symbols uppercase
      fData(eset)[, "SYMBOL"] <- toupper(fData(eset)[, "SYMBOL"])
      eset
      })

    #get common genes
    all_genes <- lapply(esets, function(x) unique(fData(x)[, annot]))
    common_genes <- Reduce(intersect, all_genes)

    for (i in seq_along(esets)) {
        #keep rows with annot (ID) in common_genes
        ID <- fData(esets[[i]])[, annot]
        filter <- ID %in% common_genes
        esets[[i]] <- esets[[i]][filter, ]

    }
    return (esets)
}


#-------------------


#' Query user to provide description.
#'
#' Uses tcltk to request input from user for group or tissue names.
#'
#' @import tcltk
#'
#' @param msg Message to display in title bar.
#' @param box1 Description beside box1.
#' @param def1 Default value for box1.
#' @param box2 Description beside box2.
#' @param def2 Default value for box2.
#' @param two Do you want two input boxes (one if FALSE)?
#'
#' @seealso \link{add_contrasts}
#' @return Character vector with inputs typed into box1 and/or box2.

inputs <- function(msg="", box1="eg. AL.obob", def1="AL", box2="eg. CR.obob", def2="CR", two=F) {
    #IN:
    #OUT:

    xvar <- tclVar(def1)
    if (two) {
    yvar <- tclVar(def2)
    }

    tt <- tktoplevel()
    tkwm.title(tt,"")
    x.entry <- tkentry(tt, textvariable=xvar)
    if (two) {
    y.entry <- tkentry(tt, textvariable=yvar)
    }

    submit <- function() {
        e <- parent.env(environment())
        x <- as.character(tclvalue(xvar))
        e$x <- x
        if (two) {
        y <- as.character(tclvalue(yvar))
        e$y <- y
        }
        tkdestroy(tt)
    }
    submit.but <- tkbutton(tt, text="Submit", command=submit)

    reset <- function() {
        tclvalue(xvar)<-""
        if (two) {
        tclvalue(yvar)<-""
        }
    }
    reset.but <- tkbutton(tt, text="Reset", command=reset)

    tkgrid(tklabel(tt,text=msg),columnspan=2)
    tkgrid(tklabel(tt,text=box1), x.entry, pady = 10, padx =10)
    if (two) {
    tkgrid(tklabel(tt,text=box2), y.entry, pady = 10, padx =10)
    }
    tkgrid(submit.but, reset.but)

    tkwait.window(tt)
    if (two) {
    return (c(x,y))
    } else {
    return (x)
    }
}
