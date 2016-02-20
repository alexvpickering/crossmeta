#load libraries

library(GEOquery)
library(stringr)
library(tcltk)
library(RColorBrewer)
library(sva)
library(limma)

#------------------------

inputs <- function(instructions="", two=F, box1="AL group name", box2="CR group name"){
    #IN:
    #OUT:

    xvar <- tclVar("")
    if (two) {
    yvar <- tclVar("")
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

    tkgrid(tklabel(tt,text=instructions),columnspan=2)
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

#------------------------

get_raw <- function(gse_name, data_dir) {
    #wrapper for get_raw_one
    lapply(gse_name, get_raw_one, data_dir)
}


get_raw_one <- function (gse_name, data_dir) {
    #IN
    #OUT
    gse_dir <- paste(data_dir, gse_name, sep="/")
    
    #get raw data
    if (!file.exists(gse_dir)) {
        getGEOSuppFiles(gse_name, baseDir=data_dir)
    }   
    #unzip
    rawtxt_paths <- list.files(gse_dir, pattern="non-norm.*txt.gz", full.names=T)
    sapply(rawtxt_paths, gunzip)  
}

#------------------------

open_raw <- function (gse_names, data_dir) {
    #IN: opens non-normalized illumina data (change to below format)
    #OUT: names of successfully formated (probeid = "ID_REF",
    #                                     exprs = "AVG_Signal-sample_name",
    #                                     pvals = "Detection-sample_name",
    #                                     sep = ",") 
    out_names <- gse_names
    for (i in seq_along(gse_names)) {
        #open data
        gse_dir <- paste(data_dir, gse_names[i], sep="/")
        data_paths <- list.files(gse_dir, pattern="non-norm.*txt", full.names=T)
        system2("xdg-open", data_paths)
        #check success
        success <- tk_select.list(choices = c("Yes", "No"),
                                  title = paste(gse_name, "formated successfully?"))
        #remove unsuccessful
        if (success == "No") {
            out_names <- setdiff(out_names, gse_names[i])  
        }
    }
    return(out_names)
}

#------------------------

load_eset <- function (gse_name, data_dir) {
    #wrapper for load_eset_one
    esets <- lapply(gse_name, load_eset_one, data_dir)
    names(eset) <- gse_name
    return (eset)
}


load_eset_one <- function (gse_name, data_dir) {
    #IN: non-normalized illumina data (format: probeid = "ID_REF",
    #                                          exprs = "AVG_Signal-sample_name",
    #                                          pvals = "Detection-sample_name",
    #                                          sep = ",")
    #OUT
    gse_dir <- paste(data_dir, gse_name, sep="/")
    
    #get GSEMatrix (for pheno data)
    eset <- getGEO(gse_name, destdir=gse_dir, GSEMatrix=T)[[1]]
    
    #load non-normalized txt files and normalize
    data_paths <- list.files(gse_dir, pattern="non-norm.*txt", full.names=T)
    data <- read.ilmn(data_paths, probeid="ID_REF", sep=",")
    data <- neqc(data)$E
    
    #transfer exprs from data to eset (maintaining eset feature order)
    feature_order <- featureNames(eset)
    pData(eset)$title_GSEMatrix <- pData(eset)$title  #to check if sample order mismatch
    pData(eset)$title <- colnames(data)  #use raw data titles to ensure correct contrasts
    colnames(data) <- sampleNames(eset)
    exprs(eset) <- data[feature_order,]
    
    return (eset)
}

#------------------------

diff_expr <- function (eset) {
    #wrapper for diff_expr_one
    top_tables <- mapply(diff_expr_one, eset, names(eset))
    names(top_tables) <- names(eset)
    return (top_tables)
}


diff_expr_one <- function (eset, name) {
    #INPUT:
    #OUTPUT:

    #BLOCKING VARIABLES:
    #-------------------
    block_names <- c()

    #ask if want to add blocking variable?
    choices <- paste(sampleNames(eset), pData(eset)$title)
    i <- 1  #block count
    while (TRUE){
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
  
  
    #CONTRASTS:        
    #----------

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
  
    #MDS PLOT:
    #---------

    group <- factor(pData(eset)$group, levels=group_levels)
    palette <- brewer.pal(8, "Dark2")
    colours <- palette[group]

    plotMDS(exprs(eset), labels = sampleNames(eset), 
            main = name, col = colours)


    #DIFFERENTIAL EXPRESSION (limma):
    #-------------------------------

    #make full model matrix   
    vars <- c("~0", "group", block_names)
    fmla <- as.formula(paste(vars, collapse="+"))
    mod <- model.matrix(fmla, data=pData(eset))
    colnames(mod)[seq_along(group_levels)] <- group_levels

    #make null model matrix (sva)
    vars0 <- c("~1", block_names)
    fmla0 <- as.formula(paste(vars0, collapse="+"))
    mod0 <- model.matrix(fmla0, data=pData(eset))

    svobj <- sva(exprs(eset), mod, mod0)
    modSv <- cbind(mod, SV1=svobj$sv)  #TODO: if #SV > 1
    contrast_matrix <- makeContrasts(contrasts=contrasts, levels=modSv)
    print (modSv)
    
    fit <- contrasts.fit (lmFit(exprs(eset),modSv), contrast_matrix)
    ebayes <- eBayes(fit)

    top_tables <- list()
    for (i in seq_along(contrasts)){
        top_genes <- topTable(ebayes, coef=i, n=Inf, resort.by="logFC", p.value=0.05, lfc=1)
        top_tables[[contrasts[i]]] <- top_genes
    }
    return (top_tables)
}