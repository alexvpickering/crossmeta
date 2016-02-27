library(dplyr)
library(magrittr)
library(MetaDE)
options(warn=-1)

gse_names <- c("GSE60596")#, "GSE68646","GSE50789",
               #"GSE16790","GSE51885","GSE34773","GSE14202")

data_dir <- paste(getwd(), "data", sep="/")

load_anal <- function(gse_name, data_dir) {
    anals <- lapply(gse_name, load_anal_one, data_dir)
    names(anals) <- gse_name
    return (anals)
}

load_anal_one <- function(gse_name, data_dir) {
    gse_dir <- file.path (data_dir, gse_name)
    diff_name <- paste (gse_name, "_diff_expr.rds", sep="")
    diff_path <- file.path (gse_dir, diff_name)
    return (readRDS (diff_path))    
}

anals <- load_anal(gse_names, data_dir)


#-----------------------------------
#             MetaDE
#-----------------------------------

merge_intra <- function(anal, val_name) {
    intra_vals <- lapply(anal, merge_intra_one, val_name)
    names(intra_vals) <- names(anal)
    return(intra_vals)
}

merge_intra_one <- function(anal, val_name) {
    #IN: result of diff_expr (contains top_tables for each contrast)
    #OUT: dataframe with values of val_name for each contrast
    tts <- anal$top_tables
    tts %>% 
        Reduce(function(x, y) full_join(x, y, by=c("gene_symbols")), .) %>%
        select(gene_symbols, contains(val_name)) ->
        intra_vals

    names(intra_vals)[-1] <- names(tts)
    return (intra_vals)
}
      
               
               
merge_inter <- function(intra_vals) {
    #IN: result of merge_intra
    #OUT: dataframe 
    intra_vals %>% 
        Reduce(function(x, y) full_join(x, y, by=c("gene_symbols")), .) ->
        inter_vals
    
    rownames(inter_vals) <- inter_vals$gene_symbols
    inter_vals %<>% data.matrix(rownames.force=T)
    inter_vals <- inter_vals[,-1]
    return (inter_vals)
}


run_MetaDE <- function(anals) {
    #merge p/t-values from within study contrasts
    intra_pvals <- merge_intra(anals, "adj.P.Val")
    intra_tstats <- merge_intra(anals, "t")
    
    #merge above results between studies
    inter_pvals <- merge_inter(intra_pvals)
    inter_tstats <- merge_inter(intra_tstats)
    
    MetaDE_in <- list(stat=inter_tstats, p=inter_pvals)
    MetaDE_out <- MetaDE.pvalue(MetaDE_in, meta.method='AW.OC')
    return (MetaDE_out)
}