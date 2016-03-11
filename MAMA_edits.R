
library(MAMA)

#---------------------
#   Edits of MAMA
#---------------------

pvalcombination_edit <- function (esets, classes, ebayes, BHth = 0.05) {

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

EScombination_edit <- function (esets, classes, ebayes, BHth = 0.05) {
    
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

join.DEG_edit <- function (..., genenames = NULL, type = NULL, cutoff) 
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