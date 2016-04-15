library(MAMA)
library(VennDiagram)
library(GeneExpressionSignature)
source("~/Documents/Batcave/GEO/1-meta/MAMA_edits.R")
#-------------------

load_diff <- function(gse_names, data_dir, probe=F) {
    #wrapper for load_diff_one
    anals <- lapply(gse_names, load_diff_one, data_dir, probe)
    anals <- unlist(anals, recursive=F)
    return (anals)
}


load_diff_one <- function(gse_name, data_dir, probe=F) {
    #loads saved diff_data (rds file)
    gse_dir <- file.path (data_dir, gse_name)

    if (probe) {
        anal_paths <- list.files(gse_dir, pattern=".*_diff_expr_probe.rds", full=T)
    } else {
        anal_paths <- list.files(gse_dir, pattern=".*_diff_expr.rds", full=T)
    }
    #load each diff path 
    #multiple if more than one platform per GSE)
    anals <- list()
    for (path in anal_paths) {
        anal <- readRDS(path)
        anal_name <- strsplit(names(anal$top_tables), "_")[[1]][1]
        anals[[anal_name]] <- anal
    }
    return (anals)    
}

#-------------------


get_ebayes <- function(diff_exprs, sva) {
    # wrapper for ebayes_info_one
    ebayes <- mapply (get_ebayes_one, diff_exprs, sva, SIMPLIFY=F)
    ebayes <- unlist(ebayes, recursive=F, use.names=F)

    names(ebayes) <- unlist(lapply(diff_exprs, function(x) names(x$top_tables)))
    return(ebayes)
}


get_ebayes_one <- function(diff_exprs, sva) {
    # IN: diff_exprs
    # OUT: ebayes_info for each contrast

    ebayes_list <- list()
    contrast_names <- names(diff_exprs$top_tables)

    if (sva) {
        ebayes <- diff_exprs$ebayes_sv
    } else {
        ebayes <- diff_exprs$ebayes
    }
    for (i in seq_along(contrast_names)) {

        eb <- new("MArrayLM")
        eb$coefficients <- ebayes$coefficients[,i]
        eb$p.value <- ebayes$p.value[,i]
        eb$df.residual <- ebayes$df.residual
        eb$df.prior <- ebayes$df.prior
        eb$t <- ebayes$t[,i]

        ebayes_list[[contrast_names[i]]] <- eb    
    }
    return (ebayes_list)
}

#---------------------------

make_ma <- function(diff_exprs, sva) {
    #IN: diff_exprs
    #OUT: MetaArray object
    
    GEDM <- list()  #gene expression data matrices
    all_clinicals <- list()  #sample descriptions

    mama_datas <- lapply(diff_exprs, function(x) x$mama_data)
    
    for (i in seq_along(mama_datas)) {
        #add exprs to GEDM 
        if (sva) {
            exprs <- lapply(mama_datas[[i]]$esets_sva, exprs)
        } else {
            exprs <- lapply(mama_datas[[i]]$esets, exprs)
        }
        GEDM <- c(GEDM, exprs)
            
        #add gse_clinicals to all_clinicals
        gse_clinicals <- mama_datas[[i]]$clinicals
        all_clinicals <- c(all_clinicals, gse_clinicals)        
    }
    ma <- new("MetaArray", GEDM=GEDM, clinical=all_clinicals, datanames=names(GEDM))
    return(ma)
}


#-------------------------

get_tts <- function(diff_exprs) {
    #returns complete list of topTables for each contrast (SVs modeled)
    #needed for merge_ranks meta-analysis
    top_tables <- list()

    for (diff in diff_exprs) {
        
        contrast_names <- names(diff$top_tables)
        ebayes_sv <- diff$ebayes_sv
        
        for (j in seq_along(contrast_names)) {

            top_genes <- topTable(ebayes_sv, coef=j, n=Inf)
            top_tables[[contrast_names[j]]] <- top_genes   
        }
    }
    return (top_tables)
}


merge_ranks <- function(diff_exprs, n=NULL) {
    #employs RankMerging from GeneExpressionSignature

    #get complete topTables
    top_tables <- get_tts(diff_exprs)

    #setup rank matrix
    gene_names <- featureNames(diff_exprs[[1]]$eset)
    rank_matrix <- matrix(nrow = length(gene_names), 
                          ncol = length(top_tables), 
                          dimnames = list(gene_names, seq_along(top_tables)))

    #generate featureData
    featureData <- as (as.data.frame(rank_matrix), "AnnotatedDataFrame")

    for (i in seq_along(top_tables)) {
        
        #get order from most positive to most negative modt-statistic
        order.t <- order(top_tables[[i]]$t, decreasing=T)
        
        #add rank to top_table
        top_tables[[i]][order.t, "rank"] <- 1:nrow(top_tables[[i]])
        
        #add rank to rank_matrix (using consistent gene_names order)
        rank_matrix[, i] <- top_tables[[i]][gene_names, "rank"]
                
    }

    #generate ExpressionSet needed for RankMerging
    rank_eset <- ExpressionSet(rank_matrix, featureData=featureData)
    pData(rank_eset)$state <- "test"

    #merge ranks to obtain prototype ranked list (PRL)
    prl <- exprs(RankMerging(rank_eset,"Spearman",weighted=TRUE))
    prl_sorted <- sort(prl[,])

    if (!is.null(n)) {
        #get top n up/down-regulated genes from prl_sorted
        prl_up <- names(head(prl_sorted, n))
        prl_dn <- names(tail(prl_sorted, n))
        prl_genes <- c(prl_up, prl_dn)
        return(prl_genes)
    } else {
        #return complete list
        return (prl_sorted)
    }    
}