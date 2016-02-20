#load libraries

library(GEOquery)
library(affy)
library(oligo)
library(stringr)
library(affxparser)
library(tcltk)
library(RColorBrewer)
library(limma)

#------------------------

cel_dates <- function(cel_paths) {
  #IN: vector with paths to CEL files
  #OUT: vector of CEL scan dates
  scan_dates <- c()
  for (i in seq_along(cel_paths)) {
    datheader <- readCelHeader(cel_paths[i])$datheader
    scan_date <- gsub(".*([0-9]{2}/[0-9]{2}/[0-9]{2}).*", "\\1", datheader)
    scan_dates[i] <- scan_date       
  }
  return (as.factor(scan_dates))
}

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
  #untar
  tar_name <- list.files(gse_dir, pattern="tar")
  untar(paste(gse_dir, tar_name, sep="/"), exdir=gse_dir)
  
  #unzip
  cel_paths <- list.files(gse_dir, pattern=".CEL.gz", full.names=T)
  sapply(cel_paths, gunzip, overwrite=T)   
}

#------------------------

load_eset <- function (gse_name, data_dir) {
  #wrapper for load_eset_one
  eset <- lapply(gse_name, load_eset_one, data_dir)
  names(eset) <- gse_name
  return (eset)
}

load_eset_one <- function (gse_name, data_dir) {
  #IN
  #OUT
  gse_dir <- paste(data_dir, gse_name, sep="/")
  
  #get GSEMatrix (for pheno data)
  eset <- getGEO(gse_name, destdir=gse_dir, GSEMatrix=T)[[1]]
  
  #load celfiles and normalize
  cel_paths <- list.files(gse_dir, pattern=".CEL", full.names=T)
  data <- tryCatch (
    {
      raw_data <- ReadAffy (celfile.path=gse_dir)
      affy::rma(raw_data)
    },
    warning = function(cond) {
      raw_data <- read.celfiles(cel_paths)
      return (oligo::rma(raw_data))  
    } 
  )
  #rename samples in data
  sampleNames(data) <- str_extract(sampleNames(data), "GSM[0-9]+")
  
  #transfer exprs from data to eset (maintaining eset sample/feature order)
  sample_order <- sampleNames(eset)
  feature_order <- featureNames(eset)
  exprs(eset) <- exprs(data)[feature_order, sample_order]
  
  #add scan dates to pheno data (maintaining eset sample order)
  scan_dates <- cel_dates (cel_paths)
  names(scan_dates) <- sampleNames(data)
  pData(eset)$scan_date <- scan_dates[sample_order]
  
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
  choices <- paste(pData(eset)$scan_date, sampleNames(eset), pData(eset)$title)
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
  
  #make model/contrast matrix   
  variables <- c("~0", "group", block_names)
  formula <- as.formula(paste(variables, collapse="+"))
  model <- model.matrix(formula, data=pData(eset))
  colnames(model)[seq_along(group_levels)] <- group_levels
  contrast_matrix <- makeContrasts(contrasts=contrasts, levels=model)
  
  print (model)
  print (contrast_matrix)
  
  
  #MDS PLOT:
  #---------
  
  group <- factor(pData(eset)$group, levels=group_levels)
  palette <- brewer.pal(8, "Dark2")
  colours <- palette[group]
    
  plotMDS(exprs(eset), labels = group, 
          main = name, col = colours)
          
  
  
  #DIFFERENTIAL EXPRESSION (limma):
  #-------------------------------
  
  fit <- contrasts.fit (lmFit(exprs(eset), model), contrast_matrix)
  ebayes <- eBayes(fit)
  
  top_tables <- list()
  for (i in seq_along(contrasts)){
    top_genes <- topTable(ebayes, coef=i, n=Inf, resort.by="logFC", p.value=0.05, lfc=1)
    top_tables[[contrasts[i]]] <- top_genes
  }
  return (top_tables)
}
