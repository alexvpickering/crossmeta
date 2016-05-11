#' Download and unpack microarray supplementary files from GEO.
#'
#' Downloads and unpacks microarray supplementary files from GEO. Files are
#' stored in the supplied data directory under the GSE name.
#'
#' @param gse_names Character vector of GSE names to download.
#' @param data_dir  String specifying directory for GSE folders.
#'
#' @seealso \code{\link{load_raw}}, \code{\link{open_raw_illum}}.
#' @return NULL (for download/unpack only).
#' @export
#' @examples
#' get_raw("GSE41845")

get_raw <- function (gse_names, data_dir=getwd()) {

    for (gse_name in gse_names) {

        gse_dir <- paste(data_dir, gse_name, sep="/")
        #get raw data
        if (!file.exists(gse_dir)) {
            GEOquery::getGEOSuppFiles(gse_name, baseDir=data_dir)
        }
        #untar
        tar_names <- list.files(gse_dir, pattern="tar")
        if (length(tar_names) > 0) {
            utils::untar(paste(gse_dir, tar_names, sep="/"), exdir=gse_dir)
        }
        #unzip
        paths <- list.files(gse_dir, pattern=".gz",
                            full.names=TRUE, ignore.case=TRUE)
        sapply(paths, GEOquery::gunzip, overwrite=TRUE)
    }
}

#-----------

#' Load and annotate raw data downloaded from GEO.
#'
#' Loads and annotates raw data previously downloaded with \code{get_raw}.
#' Supported platforms include Affymetrix, Agilent, and Illumina.
#'
#' @import tcltk data.table
#'
#' @param gse_names Character vector of GSE names.
#' @param data_dir  String specifying directory with GSE folders.
#'
#' @seealso \code{\link{get_raw}} to obtain raw data.
#'
#' @return List of annotated esets.
#' @export
#' @examples
#' library(lydata)
#' data_dir <- system.file("extdata", package = "lydata")
#' eset <- load_raw("GSE9601", data_dir)


load_raw <- function(gse_names, data_dir=getwd()) {

    affy_names  <- c()
    agil_names  <- c()
    illum_names <- c()

    for (gse_name in gse_names) {

        #determine platform (based on filenames)
        gse_dir <- paste(data_dir, gse_name, sep="/")

        affy  <- list.files(gse_dir, ".CEL", ignore.case=TRUE)
        agil  <- list.files(gse_dir, "^GSM.*txt", ignore.case=TRUE)
        illum <- list.files(gse_dir, "non.norm.*txt", ignore.case=TRUE)

        #add to appropriate names vector
        if (length(affy) != 0) {
            affy_names  <- c(affy_names, gse_name)

        } else if  (length(agil) != 0) {
            agil_names  <- c(agil_names, gse_name)

        } else if (length(illum) != 0) {
            illum_names <- c(illum_names, gse_name)
        }
    }

    #load esets
    affy_esets  <- load_affy(affy_names, data_dir)
    agil_esets  <- load_agil(agil_names, data_dir)
    illum_esets <- load_illum(illum_names, data_dir)

    return (c(affy_esets, agil_esets, illum_esets))
}

# Title
#
# @import data.table
# @param eset
# @param data
#
# @return
#
# @examples

# !
# USE FOR load_illum and load_agil if possible
# !

merge_fdata <- function(eset, data) {

    #merge feature data from raw data and eset
    fData(eset)$rn <- row.names(fData(eset))
    fData(data)$rn <- row.names(fData(data))

    dt1 <- data.table(fData(eset), key="rn")
    dt2 <- data.table(fData(data), key="rn")

    feature_data <- dt1[dt2]
    feature_data <- as.data.frame(feature_data)
    row.names(feature_data) <- feature_data$rn
    feature_data$rn <- NULL

    return(feature_data[featureNames(data), ])
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

    if (!require(biocpack_name, character.only=TRUE)) {
        source("https://bioconductor.org/biocLite.R")
        biocLite(biocpack_name)
        require(biocpack_name, character.only=TRUE)
    }
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

    #get from gpl_bioc
    biocpack_name <- gpl_bioc[gpl_name, "bioc_package"]

    #manual entry if needed
    if (is.na(biocpack_name)) {
        biocpack_name <- inputs(box1="Enter biocpack_name", def1=gpl_name)
    }
    return (paste(biocpack_name, ".db", sep=""))
}


#------------------------


# Add gene symbol to expression set.
#
# Function uses platform (GPL) to identify and download corresponding
# bioconductor annotation data package. Gene ("SYMBOL") and probe ("PROBE")
# names are then added to featureData slot of supplied eset.
#
# @param eset Expression set to annotate.
# @param gpl_name Platform name used to look up annotation data package.
#
# @seealso \code{\link{load_raw}}, \code{\link{get_biocpack_name}},
#   \code{\link{get_biocpack}}.
# @return Annotated eset.

symbol_annot <- function (eset, gpl_name) {
    biocpack_name <- get_biocpack_name(gpl_name)
    SYMBOL <- NA  #value if no biocpack/selection

    if (!biocpack_name %in% c("", ".db")) {
        #get map for SYMBOL from biocpack
        get_biocpack(biocpack_name)
        ID <- featureNames(eset)
        map <- AnnotationDbi::select(get(biocpack_name), ID, "SYMBOL")
        eset <- eset[map$PROBEID, ]  #expands one-to-many mappings
        SYMBOL <- map$SYMBOL
    } else {
        #try fData column
        choices <- setdiff(fvarLabels(eset), "SYMBOL")
        column <- tcltk::tk_select.list(choices, title="select SYMBOL column")
        if (column != "") {
            SYMBOL <- fData(eset)[, column]
        }
    }
    fData(eset)$PROBE <- featureNames(eset)
    fData(eset)$SYMBOL <- SYMBOL
    return (eset)
}


#-------------------


# Query user to provide description.
#
# Uses tcltk to request input from user.
#
# @param msg Message to display in title bar.
# @param box1 Description beside box1.
# @param def1 Default value for box1.
# @param box2 Description beside box2.
# @param def2 Default value for box2.
# @param two Do you want two input boxes? If FALSE, one box.
#
# @seealso \code{\link{add_contrasts}}, \code{\link{get_biocpack_name}}.
# @return Character vector with inputs typed into box1 and/or box2.

inputs <- function(msg="", box1="eg. CTRL",
                   def1="CTRL", box2="eg. TEST", def2="TEST", two=FALSE) {

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
