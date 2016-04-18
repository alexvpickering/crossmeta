#' Edit of \link[metaMA]{pvalcombination} for multiple contrasts per study.
#'
#' Results from seperate differential expression can be supplied. This allows for
#' multiple contrasts to be analysed within a study.
#'
#' @importFrom metaMA directpvalcombi IDDIRR
#'
#' @param esets list of expression sets (one per contrast).
#' @param classes list of factor vectors indicating class labels for each contrast.
#' @param ebayes result of call to limma eBayes function.
#' @param BHth cutoff for false discovery rate.
#'
#' @export
#' @seealso \link{make_ma} to make MetaArray object from result of call to diff_expr.
#'
#'          \link{GEDM} and \link{selectClass} to get esets and classes respectively.
#'
#'          \link{get_ebayes} to extract ebayes info for each contrast.
#'
#'          \link[metaMA]{pvalcombination} for original method.
#' @return List
#'         \item{Study1}{Vector of indices of differentially expressed genes in
#'          study 1. Similar names are given for the other individual studies.}
#'         \item{AllIndStudies}{Vector of indices of differentially expressed
#'         genes found by at least one of the individual studies.}
#'         \item{Meta}{Vector of indices of differentially expressed genes in
#'         the meta-analysis.}
#'         \item{TestStatistic}{Vector with test statistics for differential
#'         expression in the meta-analysis.}
#' @examples \dontrun{
#'
#'}

pvalcombination <- function (esets, classes, ebayes, BHth = 0.05) {

    nbstudies = length(esets)


    for (i in 1:nbstudies) {
        Ref <- levels(classes[[i]])[1]
        classes[[i]] <- sapply(classes[[i]], function(x) ifelse(x == Ref, 0, 1))
    }

    listgd = vector("list", (nbstudies + 3))

    for (i in 1:nbstudies) {

        listgd[[i]] = which(p.adjust(ebayes[[i]]$p.value, method = "BH") <= BHth)
        p1sidedLimma = pt(ebayes[[i]]$t, df = (ebayes[[i]]$df.prior + ebayes[[i]]$df.residual))
        assign(paste("p1sidedLimma", i, sep = ""), p1sidedLimma)
    }
    tempvec = paste("p1sidedLimma", 1:nbstudies, sep = "")

    lsinglep = lapply(tempvec, FUN = function(x) get(x, inherits = TRUE))
    nrep = unlist(lapply(classes, FUN = function(x) length(x)))
    listgd[[(nbstudies + 1)]] = unique(unlist(listgd[1:nbstudies]))
    restempdirect = directpvalcombi(lsinglep, nrep, BHth)
    listgd[[(nbstudies + 2)]] = restempdirect$DEindices
    listgd[[(nbstudies + 3)]] = restempdirect$TestStatistic
    names(listgd) = c(paste("study", 1:nbstudies, sep = ""),
        "AllIndStudies", "Meta", "TestStatistic")
    restemp = IDDIRR(listgd$Meta, listgd$AllIndStudies)
    print(restemp)
    invisible(listgd)
}


#--------------------------


#' Edit of metaMA::\link{EScombination} for multiple contrasts per study.
#'
#' Results from seperate differential expression can be supplied. This allows for
#' multiple contrasts to be analysed within a study.
#'
#' @importFrom metaMA effectsize directEScombi IDDIRR
#'
#' @param esets list of expression sets (one per contrast).
#' @param classes list of factor vectors indicating class labels for each contrast.
#' @param ebayes result of call to \code{get_ebayes}.
#' @param BHth cutoff for false discovery rate.
#'
#' @export
#' @seealso \link{make_ma} to make MetaArray object from result of call to diff_expr.
#'
#'          \link{GEDM} and \link{selectClass} to get esets and classes respectively.
#'
#'          \link{get_ebayes} to extract ebayes info for each contrast.
#'
#'          \link[metaMA]{pvalcombination} for original method.
#' @return List
#'         \item{Study1}{Vector of indices of differentially expressed genes in
#'          study 1. Similar names are given for the other individual studies.}
#'         \item{AllIndStudies}{Vector of indices of differentially expressed
#'         genes found by at least one of the individual studies.}
#'         \item{Meta}{Vector of indices of differentially expressed genes in
#'         the meta-analysis.}
#'         \item{TestStatistic}{Vector with test statistics for differential
#'         expression in the meta-analysis.}
#' @examples \dontrun{
#'
#'}

EScombination <- function (esets, classes, ebayes, BHth = 0.05) {

    nbstudies = length(esets)

    for (i in 1:nbstudies) {
        Ref <- levels(classes[[i]])[1]
        classes[[i]] <- sapply(classes[[i]], function(x) ifelse(x == Ref, 0, 1))
    }

    listgd = vector("list", (nbstudies + 3))
    ES = array(dim = c(dim(esets[[1]])[1], 4, nbstudies))


    for (i in 1:nbstudies) {

        listgd[[i]] = which(p.adjust(ebayes[[i]]$p.value, method = "BH") <= BHth)
        n1i = length(which(classes[[i]] == 1))
        n2i = length(which(classes[[i]] == 0))
        ES[, , i] = effectsize(ebayes[[i]]$t, ((n1i * n2i)/(n1i +
            n2i)), (ebayes[[i]]$df.prior + ebayes[[i]]$df.residual))
    }

    listgd[[(nbstudies + 1)]] = unique(unlist(listgd[1:nbstudies]))
    restempdirect = directEScombi(ES[, 3, ], ES[, 4, ], BHth)
    listgd[[(nbstudies + 2)]] = restempdirect$DEindices
    listgd[[(nbstudies + 3)]] = restempdirect$TestStatistic
    names(listgd) = c(paste("study", 1:nbstudies, sep = ""),
        "AllIndStudies", "Meta", "TestStatistic")
    restemp = IDDIRR(listgd$Meta, listgd$AllIndStudies)
    print(restemp)
    invisible(listgd)
}

#---------------------------

#' Fixes bug in MAMA join.DEG.
#'
#' Bug was present under condition \code{if (type[i] == 3)}. Column was selected
#' by number, but is not in general constant. Fix refers to column by name
#' ("FDR").
#'
#' @param ... Outputs from different function for methods of meta-analysis of
#'        microarray
#' @param genenames a character vector - names of all genes (or probe ID)
#'        included in meta-analysis. It can be NULL if the wrapper functions were
#'        used for the analysis.
#' @param type a numeric vector idicating from which function the output is,
#'        kth element in type corresponds to kth element of .... It is not
#'        needed when wrapper functions where used.
#' @param cutoff a numeric value - a cutoff level for p-value to select
#'        significant genes
#' @seealso \link[MAMA]{join.DEG} for original method.
#' @return A list in which each slot refers to one meta-analytical method and
#'         contains names of differentially expressed genes found by the method.
#' @examples \dontrun{
#'
#'}

join.DEG <- function (..., genenames = NULL, type = NULL, cutoff)
{
    args <- list(...)
    N <- length(args)
    if (!(is.null(type)) & N != length(type))
        stop("Vector type has not correct length")
    genelist <- list()
    if (is.null(type)) {
        for (i in 1:N) {
            if ("metaMA.res" %in% class(args[[i]]))
                genelist[[i]] <- args[[i]]$gene.names[args[[i]]$Meta]
            if ("ES.GeneMeta.res" %in% class(args[[i]]))
                genelist[[i]] <- rownames(args[[i]]$ScoresFDR$two.sided)[args[[i]]$ScoresFDR$two.sided[,
                  "FDR"] < cutoff]
            if ("RankProduct.res" %in% class(args[[i]]))
                genelist[[i]] <- unique(c(rownames(args[[i]]$Table1),
                  rownames(args[[i]]$Table2)))
            if ("SOGLresult" %in% class(args[[i]]))
                genelist[[i]] <- args[[i]]$genes
            if ("posterior.mean" %in% class(args[[i]]))
                genelist[[i]] <- rownames(args[[i]])[args[[i]]$Pvalue <
                  cutoff & !(is.nan(args[[i]]$Pvalue))]
            if ("MAP.Matches.res" %in% class(args[[i]]))
                genelist[[i]] <- unique(unlist(args[[i]]$genes))
        }
    }
    else {
        if (is.null(genenames))
            stop("The 'genenames' must be provided")
        for (i in 1:N) {
            if (type[i] == 1) {
                genelist[[i]] <- genenames[args[[i]]$Meta]
            }
            if (type[i] == 2) {
            }
            if (type[i] == 3) {
                genelist[[i]] <- rownames(args[[i]]$two.sided)[args[[i]]$two.sided[,
                  "FDR"] < cutoff]
            }
            if (type[i] == 4) {
                genelist[[i]] <- args[[i]]$genes
            }
            if (type[i] == 5) {
                genelist[[i]] <- unique(c(rownames(args[[i]]$Table1),
                  rownames(args[[i]]$Table2)))
            }
            if (type[i] == 6) {
                genelist[[i]] <- rownames(args[[i]])[args[[i]]$Pvalue <
                  cutoff]
            }
            if (type[i] == 7) {
                genelist[[i]] <- unique(unlist(args[[i]]))
            }
        }
    }
    return(genelist)
}
