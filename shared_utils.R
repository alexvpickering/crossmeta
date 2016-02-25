library(GEOquery)
library(GEOmetadb)
library(annotate)
library(dplyr)
library(RColorBrewer)
library(tcltk)
library(sva)
library(limma)

#------------------------

get_biocpack <- function(biocpack_name) {
    #IN:
    #OUT:
    biocpack_db <- paste(biocpack_name, ".db", sep="")

    if (!require(biocpack_db , character.only=T)) {
        source("https://bioconductor.org/biocLite.R")
        biocLite(biocpack_db)
        require(biocpack_db , character.only=T)
    }
}

#------------------------

#GPL platform maps are supplied by investigator and not updated.
#Bioc platform maps are updated.

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
    return (biocpack_name)
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
        SYMBOL <- getSYMBOL(ID, biocpack_name)
        SYMBOL <- SYMBOL[!is.na(SYMBOL)]  #removes IDs where SYMBOL=NA
    } else { 
        #try fData column
        choices <- setdiff(fvarLabels(eset), "SYMBOL")
        column <- tk_select.list(choices, title="select SYMBOL column")
        if (column != "") {
            SYMBOL <- fData(eset)[, column]
        }
    }
    eset <- eset[names(SYMBOL),]  #subset: only rows with SYMBOL
    fData(eset)$SYMBOL <- SYMBOL
    return (eset)
}

#------------------------

cleanY <- function(y, mod, svs) {
    #IN:  y   = exprs(eset), 
    #     mod = full model matrix supplied to sva, 
    #     svs = surrogate variables returned by sva (svobj$sv)
    #OUT: exprs(eset) adjusted for surrogate variables

    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    rm(Hat)
    gc()
    P = ncol(mod)
    return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

#------------------------

add_blocking <- function (eset, choices) {

    #ask if want to add blocking variable?
    i <- 1  #block count
    block_names <- c()
    while (TRUE) {
        block <- tk_select.list(c(choices, "YES", "NO"), title="Add blocking variable?")
        if (block != "YES") {break}

        #if YES, setup
        block_name <- paste("block", i, sep="_")
        block_names <- c(block_names, block_name)
        pData(eset)[, block_name] <- "level_0"  #baseline level
        i <- i + 1

        #select samples in each level of blocking variable (except 1)
        j <- 1  #level count 
        while (TRUE) {
            level <- tk_select.list(choices, multiple=T, title="Select samples in each level (except 1)")    
            level <- str_extract(level, "GSM[0-9]+")
            if (length(level) == 0) {break}  

            #add level to pheno
            level_name <- paste("level", j, sep="_")
            pData(eset)[level, block_name] <- level_name
            j <- j + 1
        }
    }
    return (list("eset"=eset, "names"=block_names))
}


add_contrasts <- function (eset) {

    choices <- paste(sampleNames(eset), pData(eset)$title)
    contrasts <- c()
    group_levels <- c()
    selected_samples <- c()

    #repeat until all contrasts selected
    while (TRUE) {
        #select AL samples
        AL <- tk_select.list(choices, multiple=T, title="select AL (control) samples for contrast")
        AL <- str_extract(AL, "GSM[0-9]+")
        if (length(AL) == 0) {break}

        #select CR samples
        CR <- tk_select.list(choices, multiple=T, title="select CR (test) samples for contrast")
        CR <- str_extract(CR, "GSM[0-9]+")

        #add group names to pheno
        group_names <- inputs("Enter group names", two=T)
        pData(eset)[AL, "group"] <- group_names[1]
        pData(eset)[CR, "group"] <- group_names[2]

        #add to contrasts
        contrasts <- c(contrasts, paste(group_names[2], group_names[1], sep="-"))

        #add to group_levels
        group_levels <- unique(c(group_levels, group_names[1], group_names[2]))

        #add to selected_samples
        selected_samples <- unique(c(selected_samples, AL, CR))
    }
    #retain selected samples only
    eset <- eset[, selected_samples]

    return (list("eset"=eset, "contrasts"=contrasts, 
                 "levels"=group_levels, "samples"=selected_samples))
}

diff_setup <- function(eset, group_levels, block_names){

    #make full model matrix   
    vars <- c("~0", "group", block_names)
    fmla <- as.formula(paste(vars, collapse="+"))
    mod <- model.matrix(fmla, data=pData(eset))
    colnames(mod)[seq_along(group_levels)] <- group_levels

    #make null model matrix (sva)
    vars0 <- c("~1", block_names)
    fmla0 <- as.formula(paste(vars0, collapse="+"))
    mod0 <- model.matrix(fmla0, data=pData(eset))

    #surrogate variable analysis
    svobj <- sva(exprs(eset), mod, mod0)
    modsv <- cbind(mod, svobj$sv)
    colnames(modsv) <- c(colnames(mod), paste("SV", 1:svobj$n.sv, sep=""))

    return (list("mod"=mod, "modsv"=modsv, "svobj"=svobj))
}


unique_genes <- function(top_genes) {
    #IN: result from limma topTable containing gene_symbols column
    #OUT: result with unique gene_symbols
    #           - min adj.P.Val
    #           - if tie, max logFC
    #           - sorted by descending absolute logFC
    top_genes %>% 
        group_by(gene_symbols) %>%
        arrange(adj.P.Val, desc(logFC)) %>%
        slice(which.min(adj.P.Val)) %>%
        ungroup() %>%
        arrange(desc(abs(logFC)))


}


diff_anal <- function(eset, contrasts, group_levels, mod, modsv, svobj, gse_dir, name){

    #differential expression
    contrast_matrix <- makeContrasts(contrasts=contrasts, levels=modsv)
    fit <- contrasts.fit (lmFit(exprs(eset),modsv), contrast_matrix)
    ebayes <- eBayes(fit)

    #annotate/store results
    top_tables <- list()
    cat("\n\n")
    for (i in seq_along(contrasts)){
        top_genes <- topTable(ebayes, coef=i, n=Inf)  #return all (for MetaDE)
        gene_symbols <- fData(eset)[row.names(top_genes), "SYMBOL"]
        top_genes <- unique_genes(cbind(gene_symbols, top_genes))
        num_sig <- dim(top_genes[top_genes$adj.P.Val < 0.05 & top_genes$logFC > 1,])[1]
        top_tables[[contrasts[i]]] <- top_genes
        cat (contrasts[i], "(n significant):", num_sig, "\n")
    }
    cat("\n\n")

    #save to disk
    diff_expr <- list(eset, modsv, contrast_matrix, top_tables)
    names(diff_expr) <- c("eset", "modsv", "contrast_matrix", "top_tables")
    save_name <- paste (name, "diff_expr.rds", sep="_")
    saveRDS(diff_expr, file = paste(gse_dir, save_name, sep="/"))

    #plot MDS 
    group <- factor(pData(eset)$group, levels=group_levels)
    palette <- brewer.pal(12, "Paired")
    colours <- palette[group]

    sva_exprs <- cleanY(exprs(eset), mod, svobj$sv)

    plotMDS(sva_exprs, pch=19, main = name, col = colours)
    legend("center", legend=group_levels, fill=unique(colours))

    return (diff_expr)
}


#-------------------

inputs <- function(msg="", box1="AL group name", def1="AL", box2="CR group name", def2="CR", two=F) {
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