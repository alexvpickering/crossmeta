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


load_raw <- function(gse_names, homologene_path, data_dir=getwd(),
                     overwrite=FALSE) {

    #no duplicates allowed (somehow causes mismatched names/esets)
    gse_names <- unique(gse_names)

    #get homologene
    homologene <- get_homologene(homologene_path)

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
    affy_esets  <- load_affy(affy_names, homologene, data_dir, overwrite)
    agil_esets  <- load_agil(agil_names, homologene, data_dir, overwrite)
    illum_esets <- load_illum(illum_names, homologene, data_dir, overwrite)

    return (c(affy_esets, agil_esets, illum_esets))
}



#' Get homologene data frame.
#'
#' Function loads and, if necessary, sets up homologene data frame.
#'
#' Setup results in a dataframe with all entrez ids in one column ('entrez') and
#' the homologous human entrez ids in another column ('entrez_HS').
#'
#' @param homologene_path
#'
#' @return
#' @export
#'
#' @examples
get_homologene <- function(homologene_path) {

    homologene <- read.delim(homologene_path, header = FALSE)

    if (ncol(homologene) == 6) {

        message("Setting up homologene data. Will take a while (one time only).")

        #get homologous human (9606) entrez ids for all entrez ids (V3)
        entrez_HS <- annotationTools::getHOMOLOG(homologene$V3,
                                                 9606,
                                                 homologene)

        #remove entrez ids with homologous human entrez id
        homologene <- homologene[!is.na(entrez_HS), ]
        entrez_HS  <- entrez_HS[!is.na(entrez_HS)]

        #expand homologene where multiple homologous human entrez ids
        homologene$entrez_HS <- sapply(entrez_HS,
                                       function(x) paste(x, collapse=","))

        homologene <- data.table::data.table(homologene)
        homologene <- homologene[,
                                 list(entrez_HS = unlist(strsplit(entrez_HS, ","))),
                                 by = V3]

        homologene <- as.data.frame(homologene)

        #save to disc
        write.table(homologene, file = homologene_path,
                    sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    colnames(homologene) <- c("entrez", "entrez_HS")
    return(homologene)
}

# Title
#
# Returns data.frame, same order/nrow as data_fdat with feature info from eset
#
# @import data.table
# @param eset
# @param data
#
# @return
#
# @examples


merge_fdata <- function(eset_fdat, data_fdat) {

    #merge feature data from raw data and eset
    eset_fdat$rn <- sapply(strsplit(row.names(eset_fdat), "\\."), `[[`, 1)
    data_fdat$rn <- sapply(strsplit(row.names(data_fdat), "\\."), `[[`, 1)

    dt1 <- data.table(eset_fdat, key="rn")
    dt2 <- data.table(data_fdat, key="rn")

    feature_data <- dt1[dt2]
    feature_data <- as.data.frame(feature_data)
    row.names(feature_data) <- make.unique(feature_data$rn)
    feature_data$rn <- NULL

    return(feature_data[row.names(data_fdat), ])
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
        source("https://bioconductor.org/biocLite.R")
        biocLite(biocpack_name)
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

symbol_annot <- function (eset, homologene, gse_name = "") {

    cat("Annotating")


    #--------------------- Get Entrez ID


    # check internal data or UI for biocpack_name
    biocpack_name <- get_biocpack_name(annotation(eset))
    hs <- get_biocpack("org.Hs.eg.db")

    # if biocpack_name not empty, use to get entrez id
    if (!biocpack_name %in% c("", ".db")) {

        biocpack <- get_biocpack(biocpack_name)
        ID <- sapply(strsplit(featureNames(eset), "\\."), `[[`, 1)

        suppressMessages(map <- AnnotationDbi::select(biocpack, ID, "ENTREZID"))
        eset <- eset[map$PROBEID, ] #expands one-to-many mappings
        entrez <- map$ENTREZID

    # if not, ask user for fData column with entrez id
    } else {

        # default if no selection
        entrez <- rep(NA, nrow(eset))

        # ask for fData column
        while (TRUE) {
            choices <- fvarLabels(eset)
            column <- tcltk::tk_select.list(choices, title="select ENTREZID column")
            if (column == "") break

            entrez <- as.character(fData(eset)[, column])

            #test if column is mostly digits
            chk <- grepl("^[[:digit:]]*$", unique(entrez))
            if (sum(chk) >= 0.5 * length(unique(entrez))) break

            message(column, " not mostly digits - unlikely to be entrez ids.")
        }

        if (column != "") {
            #expand one-to-many
            entrez <- strsplit(entrez, "\\D+")
            rn <- sapply(seq_along(entrez),
                         function(x) rep(x, length(entrez[[x]])))

            entrez <- as.integer(unlist(entrez))
            eset <- eset[unlist(rn), ]
        }
    }


    #--------------------- Map Entrez ID to Homologous Human Entrez ID


    # where entrez id in homologene:
    endf <- data.frame(entrez, rn = 1:length(entrez))
    filt <- entrez %in% homologene$entrez
    map <- merge(endf[filt, ], homologene, by="entrez")

    # where no homology, use original entrez id (useful if human platform):
    endf$entrez_HS <- entrez
    map <- rbind(endf[!filt, ], map)
    eset <- eset[map$rn, ] #expands one-to-many mappings



    #--------------------- Map Human Entrez ID to Gene Symbol


    SYMBOL <- tryCatch (
        {
            suppressMessages(
                SYMBOL <- AnnotationDbi::select(hs, as.character(map$entrez_HS),
                                      "SYMBOL", "ENTREZID")$SYMBOL)
        },
        error = function(c) {
            if (grepl("valid keys", c$message)) {
                warning(gse_name, ": Annotation failed. ",
                        "Add human gene symbols to 'SYMBOL' column in fData.")
                return(NULL)
            }
        }
    )
    # PROBE is feature names (remove '.' then restore '*' to '.')
    PROBE <- sapply(strsplit(featureNames(eset), "\\."), `[[`, 1)
    PROBE <- gsub("*", ".", PROBE, fixed = TRUE)

    # add PROBE to fData and use for unique row names
    fData(eset)$PROBE <- PROBE
    row.names(eset) <- make.unique(PROBE)

    # add uppercase gene symbols to fData
    if (!is.null(SYMBOL)) {
        fData(eset)$SYMBOL <- toupper(SYMBOL)
    }
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
