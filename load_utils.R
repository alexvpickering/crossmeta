library(GEOquery)
library(GEOmetadb)
library(annotate)
library(tcltk)

#-------------------
#  Load Utilities
#-------------------

# GPL platform maps are supplied by investigator and not updated.
# Bioc platform maps are updated.

get_biocpack <- function(biocpack_name) {
    #IN:
    #OUT:

    if (!require(biocpack_name , character.only=T)) {
        source("https://bioconductor.org/biocLite.R")
        biocLite(biocpack_name)
        require(biocpack_name , character.only=T)
    }
}


get_biocpack_name <- function (gse_name) {
    #IN:
    #OUT:

    #connect to GEOmetadb database
    if (!file.exists("~/Documents/Batcave/GEO/GEOmetadb.sqlite")) {
        getSQLiteFile("~/Documents/Batcave/GEO/")
    }
    con <- dbConnect(SQLite(), "~/Documents/Batcave/GEO/GEOmetadb.sqlite")

    #query metadata database
    query <- paste(
        "SELECT gpl.bioc_package, gpl.title",
        "FROM gpl JOIN gse_gpl ON gpl.gpl=gse_gpl.gpl",
        "WHERE gse_gpl.gse=", shQuote(gse_name))

    biocpack_name <- dbGetQuery(con, query)[, "bioc_package"]

    #manual entry if needed
    if (is.na (biocpack_name)) {
        title <- dbGetQuery(con, query)[, "title"]
        biocpack_name <- inputs(box1="Enter biocpack_name", def1=title)
    }

    dbDisconnect(con)
    return (paste(biocpack_name, ".db", sep=""))
}

#------------------------

symbol_annot <- function (eset, gse_name) {
    #IN:
    #OUT:
    biocpack_name <- get_biocpack_name(gse_name) 
    SYMBOL <- NA  #value if no biocpack/selection

    if (biocpack_name != "") {
        #get map for SYMBOL from biocpack
        get_biocpack(biocpack_name)
        ID <- featureNames(eset)
        map <- AnnotationDbi::select(get(biocpack_name), ID, "SYMBOL")
        map <- map[!is.na(map$SYMBOL),]
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
    fData(eset)$SYMBOL <- SYMBOL
    return (eset)
}


commonize <- function(esets) {
    #IN: esets
    #OUT: esets with common genes
    all_genes <- lapply(esets, function(x) unique(fData(x)$SYMBOL))   
    common_genes <- Reduce(intersect, all_genes)

    for (i in seq_along(esets)) {
        #keep rows with SYMBOL in common_genes
        SYMBOL <- fData(esets[[i]])$SYMBOL
        filter <- SYMBOL %in% common_genes
        esets[[i]] <- esets[[i]][filter, ]

    }
    return (esets)
}

#-------------------

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