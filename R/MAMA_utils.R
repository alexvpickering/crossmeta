#' Load previous differential expression analysis.
#'
#' Previous runs of diff_expr are loaded.
#'
#' @param gse_names Character vector of GSE names to be loaded.
#' @param data_dir base data directory containing a folder for each GSE name.
#' @param probe Load probe level analysis? If FALSE, loads gene level analysis.
#'
#' @export
#' @seealso \link{diff_expr}
#' @return saved results from previous call to \code{diff_expr}.
#' @examples
#'
load_diff <- function(gse_names, data_dir, probe=F) {

  anals <- list()
  for (gse_name in gse_names) {
    gse_dir <- file.path (data_dir, gse_name)

    if (probe) {
      anal_paths <- list.files(gse_dir, pattern=".*_diff_expr_probe.rds", full=T)
    } else {
      anal_paths <- list.files(gse_dir, pattern=".*_diff_expr.rds", full=T)
    }
    #load each diff path
    #multiple if more than one platform per GSE)
    for (path in anal_paths) {
      anal <- readRDS(path)
      anal_name <- strsplit(names(anal$top_tables), "_")[[1]][1]
      anals[[anal_name]] <- anal
    }
  }
  return (anals)
}

#-------------------


#' Get ebayes info for each contrast.
#'
#' Helper utility in order to obtain ebayes info for each contrast.
#'
#' @param diff_exprs result of previous call to \code{diff_expr} of \code{load_diff}.
#' @param sva Do you want ebayes info from model with surrogate variables?
#'
#' @export
#' @seealso \link{diff_expr}, or \link{load_diff} to obtain diff_exprs.
#'
#'          \link{pvalcombination}, and \link{EScombination} for functions that
#'          require result of \code{get_ebayes}.
#' @return list of MArrayLM objects (one per contrast - contain ebayes info).
#' @examples \dontrun{
#'
#'}

get_ebayes <- function(diff_exprs, sva) {

    ebayes <- list()
    for (anal in diff_exprs) {
      contrast_names <- names(anal$top_tables)

      if (sva) {
        ebayes <- anal$ebayes_sv
      } else {
        ebayes <- anal$ebayes
      }
      for (i in seq_along(contrast_names)) {

        eb <- new("MArrayLM")
        eb$coefficients <- ebayes$coefficients[,i]
        eb$p.value <- ebayes$p.value[,i]
        eb$df.residual <- ebayes$df.residual
        eb$df.prior <- ebayes$df.prior
        eb$t <- ebayes$t[,i]

        ebayes[[contrast_names[i]]] <- eb
      }
    }
    return(ebayes)
}

#---------------------------

#' Make MetaArray object from results of differential expression.
#'
#' Function is used to construct a MetaArray object from results of call to
#' \code{diff_expr} or \code{load_diff}. MetaArray object allows for use of
#' all meta-analysis methods present in MAMA package.
#'
#' @importFrom Biobase exprs
#'
#' @param diff_exprs result of call to \code{diff_expr} or \code{load_diff}.
#' @param sva Use esets with effect of surrogate variables removed? If FALSE,
#'            none-adjusted esets will be used.
#' @export
#' @seealso \link[MAMA]{MetaArray-class} for possible meta-analysis methods.
#'
#'          \link{diff_expr} or \link{load_diff} to obtain \code{diff_expr}.
#' @return MetaArray object.

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


#' Get topTables for each contrast.
#'
#' Helper function for \code{merge_ranks} meta-analysis method.
#'
#' @param diff_exprs result of call to \code{diff_expr} or \code{load_diff}.
#'
#' @seealso \link{merge_ranks}.
#' @return list of topTables (one for each contrast in diff_exprs).
#' @examples \dontrun{
#'
#'}

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


#' Employ RankMerging meta-analysis from GeneExpressionSignature package.
#'
#' Wrapper for RankMerging function in GeneExpressionSignature.
#'
#' @importFrom Biobase ExpressionSet featureNames exprs
#' @importFrom GeneExpressionSignature RankMerging
#'
#' @param diff_exprs result of call to \code{diff_expr} or \code{load_diff}.
#' @param n number of top ranked (up-regulated) and bottom ranked (down-regulated)
#'        genes to return.
#' @export
#' @seealso \link{RankMerging}, \link{diff_expr}, \link{load_diff}.
#' @return Character vector of gene identifiers sorted by rank.
#' @examples \dontrun{
#'
#'}

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
