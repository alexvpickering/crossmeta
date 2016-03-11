library(MAMA)
library(VennDiagram)
library(GeneExpressionSignature)
source("MAMA_edits.R")
#-------------------

load_diff <- function(gse_name, data_dir) {
    #wrapper for load_diff_one
    anals <- lapply(gse_name, load_diff_one, data_dir)
    names(anals) <- gse_name
    return (anals)
}


load_diff_one <- function(gse_name, data_dir) {
    #loads saved diff_data (rds file)
    gse_dir <- file.path (data_dir, gse_name)
    diff_name <- paste (gse_name, "_diff_expr.rds", sep="")
    diff_path <- file.path (gse_dir, diff_name)
    return (readRDS (diff_path))    
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
        order.t <- order(top_tables[[i]]$t)
        
        #add rank to top_table
        top_tables[[i]][order.t, "rank"] <- 1:nrow(top_tables[[i]])
        
        #add rank to rank_matrix (using consistent gene_names order)
        rank_matrix[, i] <- top_tables[[i]][gene_names, "rank"]
                
    }

    #generate ExpressionSet needed for RankMerging
    rank_eset <- ExpressionSet(rank_matrix, featureData=featureData)
    pData(rank_eset)$state <- "CR"

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

#-------------------------


random_TAB <- function(MAT) {

    #number of genes in each method
    gene_names <- row.names(MAT)
    methods_n <- colSums(MAT)
    method_names <- colnames(MAT)
    TABS <- list()

    for (i in 1:1000) {
        #get random samples of methods_n genes
        lists <- list()
        for (j in seq_along(methods_n)) {

            rs_genes <- sample(gene_names, methods_n[j])
            lists[[j]] <- rs_genes
        }
        #turn it into a conting table
        names(lists)<-c("PvalCom", "ESCom","ESCom2", "RankProduct", "MergeRanks")
        TABS[[i]] <- conting.tab(lists)
    }
    return(TABS)      
}

#-------------------------
#adapted from DOI: 10.1186/1471-2105-12-35

draw_venn <- function(genes_list, filename) {

    n_lists <- length(genes_list)

    
    if (n_lists == 2) {
        venn.diagram(
        x = genes_list,
        filename = filename,
        lwd = 4,
        fill = c("cornflowerblue", "darkorchid1"),
        alpha = 0.75,
        label.col = "white",
        cex = 4,
        fontfamily = "serif",
        fontface = "bold",
        cat.col = c("cornflowerblue", "darkorchid1"),
        cat.cex = 3,
        cat.fontfamily = "serif",
        cat.fontface = "bold",
        cat.dist = c(0.03, 0.03),
        cat.pos = c(-20, 14)
        )
    }


    if (n_lists == 3) {
        venn.diagram(
        x = genes_list,
        filename = filename,
        col = "transparent",
        fill = c("red", "blue", "green"),
        alpha = 0.5,
        label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),
        cex = 2.5,
        fontfamily = "serif",
        fontface = "bold",
        cat.default.pos = "text",
        cat.col = c("darkred", "darkblue", "darkgreen"),
        cat.cex = 2.5,
        cat.fontfamily = "serif",
        cat.dist = c(0.06, 0.06, 0.03),
        cat.pos = 0
        )
    }

    if (n_lists == 4) {
        venn.diagram(
        x = genes_list,
        filename = filename,
        col = "black",
        lty = "dotted",
        lwd = 4,
        fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
        alpha = 0.50,
        label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
        cex = 2.5,
        fontfamily = "serif",
        fontface = "bold",
        cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
        cat.cex = 2.5,
        cat.fontfamily = "serif"
        )
    }
}