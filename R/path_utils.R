#' Differential expression of KEGG pathways.
#'
#' Performs PADOG pathway analysis using KEGG database (downloaded Feb 2017).
#'
#' If you wish to perform source-specific pathway meta-analyses,
#' \code{\link{add_sources}} must be used before \code{diff_paths}.
#'
#' For each GSE, analysis results are saved in the corresponding GSE
#' folder in \code{data_dir} that was created by \code{\link{get_raw}}. PADOG outperforms
#' other pathway analysis algorithms at prioritizing expected pathways (see references).
#'
#' @importFrom doRNG %dorng%
#'
#' @param esets List of annotated esets. Created by \code{\link{load_raw}}.
#' @param prev_anals Previous result of \code{\link{diff_expr}}, which can
#'    be reloaded using \code{\link{load_diff}}.
#' @param data_dir String specifying directory for GSE folders.
#'
#'
#' @return List of named lists, one for each GSE. Each named list contains:
#'    \item{padog_tables}{data.frames containing \code{\link[PADOG]{padog}}
#'       pathway analysis results for each contrast.}
#'
#'    If \code{\link{add_sources}} is used first:
#'    \item{sources}{Named vector specifying selected sample source for each contrast.
#'       Vector names identify the contrast.}
#'    \item{pairs}{List of character vectors indicating tissue sources that should be
#'       treated as the same source for subsequent pathway meta-analysis.}
#'
#' @references Tarca AL, Bhatti G, Romero R. A Comparison of Gene Set Analysis Methods
#'    in Terms of Sensitivity, Prioritization and Specificity. Chen L, ed. PLoS ONE.
#'    2013;8(11):e79217. doi:10.1371/journal.pone.0079217.
#'
#'    Dong X, Hao Y, Wang X, Tian W. LEGO: a novel method for gene set over-representation
#'    analysis by incorporating network-based gene weights. Scientific Reports.
#'    2016;6:18871. doi:10.1038/srep18871.
#'
#' @export
#'
#' @examples
#'
#' library(lydata)
#'
#' # location of data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # gather GSE names
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load esets
#' esets <- load_raw(gse_names, data_dir)
#'
#' # load previous differential expression analysis
#' anals <- load_diff(gse_names, data_dir)
#'
#' # add tissue sources to perform seperate meta-analyses for each source (recommended)
#' # anals <- add_sources(anals)
#'
#' # perform pathway analysis for each contrast
#' # path_anals <- diff_path(esets, anals, data_dir)

diff_path <- function(esets, prev_anals, data_dir = getwd()) {

    prev_anals <- prev_anals[names(prev_anals) %in% names(esets)]
    esets      <- esets[names(prev_anals)]
    anals      <- list()

    for (i in seq_along(esets)) {

        eset <- esets[[i]]
        gse_name <- names(esets)[i]
        prev <- prev_anals[[i]]

        cat('Working on', paste0(gse_name, ':'), '\n')

        gse_folder <- strsplit(gse_name, "\\.")[[1]][1]  # name can be "GSE.GPL"
        gse_dir <- file.path(data_dir, gse_folder)
        
        # select groups/contrasts
        if (is.null(prev)) prev <- select_contrasts(eset, gse_name)
        if (is.null(prev)) next
        
        # add groups from selection
        eset <- match_prev_eset(eset, prev)
        
        # possibly subset two-channel to be like one-channel
        eset <- ch2_subset(eset, prev)
        
        # group/contrast info from previous analysis
        cons <- colnames(prev$ebayes_sv$contrasts)
        group_levels <- unique(eset$group)
        
        # run surrogate variable analysis
        sva_mods <- get_mods(eset@phenoData)
        svobj <- run_sva(sva_mods, eset, svanal)
        
        # add surrogate variable/pair adjusted ("clean") expression matrix for iqr_replicates
        eset <- add_adjusted(eset, svobj)
        
        # remove rows with duplicated/NA annot (SYMBOL or ENTREZID)
        eset <- iqr_replicates(eset, annot = 'ENTREZID_HS', rm.dup = TRUE)
        groups <- pData(eset)$group

        pts <- list()
        # padog for each contrast
        for (con in cons) {

            cat('    ', paste0(con, '\n'))

            # contrast levels
            ctrl <- gsub('^.+?-', '', con)
            test <- gsub('-.+?$', '', con)

            # group
            incon <- groups %in% c(ctrl, test)
            group <- ifelse(groups[incon] == ctrl, 'c', 'd')

            # other padog inputs
            esetm  <- exprs(eset[, incon])
            pair   <- pData(eset[, incon])$pair
            paired <- !anyNA(pair)
            block  <- NULL
            if (paired)
                block <- pair

            anal_name <- paste(gse_name, con, sep='_')

            # run padog
            pts[[anal_name]] <- padog(esetm, group, paired, block, crossmeta::gslist, gs.names = crossmeta::gs.names,
                                      parallel = TRUE, ncr = parallel::detectCores(), verbose = FALSE)

        }

        if (length(pts) == 0) next()

        anal <- list(padog_tables = pts, sources = prev$sources, pair = prev$pair)

        # save to disk
        save_name <- paste(gse_name, "diff_path", sep = "_")
        save_name <- paste0(save_name, ".rds")
        saveRDS(anal, file.path(gse_dir, save_name))

        anals[[gse_name]] <- anal
    }
    return (anals)
}



# Pathway Analysis with Down-weighting of Overlapping Genes (PADOG)
#
# This is a general purpose gene set analysis method that downplays the
# importance of genes that apear often accross the sets of genes analyzed.
#
# The original implementation of \code{\link[PADOG]{padog}} was modified to allow two-sample
# groups and eliminate defunct KEGG.db warning.
#
# @export
#
# @importFrom foreach foreach
# @importFrom doRNG %dorng%
#
padog <- function (esetm = NULL, group = NULL, paired = FALSE, block = NULL,
                   gslist = NULL, organism = "hsa", annotation = NULL,
                   gs.names = NULL, NI = 1000, plots = FALSE, targetgs = NULL,
                   Nmin = 3, verbose = TRUE, parallel = FALSE, dseed = NULL,
                   ncr = NULL)
{



    getFDR = function(p0, p1) {
        # p1: observed p value vector; p0: permutation p value matrix.
        fdr0 = sapply(p1, function(z) {
            ifelse(is.na(z), NA, min(sum(p0 <= z, na.rm=TRUE) * sum(!is.na(p1))
                                     / sum(p1 <= z, na.rm=TRUE) / sum(!is.na(p0)), 1))
        })
        nna = !is.na(fdr0)
        fdr = fdr0[nna]
        ord = order(p1[nna], decreasing=TRUE)
        fdr[ord] = cummin(fdr[ord])
        fdr0[nna] = fdr
        fdr0
    }



    stopifnot(is(esetm, "matrix"))
    # stopifnot(all(dim(esetm) > 4))
    stopifnot(is(group, "factor") | is(group, "character"))
    stopifnot(length(group) == dim(esetm)[2])
    stopifnot(all(group %in% c("c", "d")))
    # stopifnot(all(table(group) > 2))
    if (paired) {
        stopifnot(length(block) == length(group))
        stopifnot(all(table(block) == 2))
    }
    stopifnot(is(gslist, "list"))
    stopifnot(length(gslist) >= 3)
    if (!is.null(gs.names)) {
        stopifnot(length(gslist) == length(gs.names))
    }
    stopifnot(is(NI, "numeric"))
    stopifnot(NI > 5)


    stopifnot(sum(rownames(esetm) %in% as.character(unlist(gslist))) >
                  10 & !any(duplicated(rownames(esetm))))

    Block = block
    gf = table(unlist(gslist))
    if (!all(gf == 1)) {
        if (stats::quantile(gf, 0.99) > mean(gf) + 3 * stats::sd(gf)) {
            gf[gf > stats::quantile(gf, 0.99)] <- stats::quantile(gf, 0.99)
        }
        gff <- function(x) {
            1 + ((max(x) - x)/(max(x) - min(x)))^0.5
        }
        gf = gff(gf)
    }
    else {
        fdfd = unique(unlist(gslist))
        gf = rep(1, length(fdfd))
        names(gf) <- fdfd
    }
    allGallP = unique(unlist(gslist))

    restg = setdiff(rownames(esetm), names(gf))
    appendd = rep(1, length(restg))
    names(appendd) <- restg
    gf = c(gf, appendd)
    stopifnot(all(!duplicated(rownames(esetm))))
    stopifnot(sum(rownames(esetm) %in% allGallP) > 10)
    if (verbose) {
        cat(paste("Starting with ", length(gslist), " gene sets!",
                  sep = ""))
        cat("\n")
    }
    gslist = gslist[unlist(lapply(gslist, function(x) {
        length(intersect(rownames(esetm), x)) >= Nmin
    }))]
    gs.names = gs.names[names(gslist)]
    stopifnot(length(gslist) >= 3)
    if (verbose) {
        cat(paste("Analyzing ", length(gslist), " gene sets with ",
                  Nmin, " or more genes!", sep = ""))
        cat("\n")
    }
    # if (!is.null(dseed))
    #     set.seed(dseed)
    G = factor(group)
    Glen = length(G)
    tab = table(G)
    idx = which.min(tab)
    minG = names(tab)[idx]
    minGSZ = tab[idx]
    bigG = rep(setdiff(levels(G), minG), length(G))
    block = factor(Block)
    topSigNum = dim(esetm)[1]
    combFun = function(gi, countn = TRUE) {
        g = G[gi]
        tab = table(g)
        if (countn) {
            minsz = min(tab)
            ifelse(minsz > 10, -1, choose(length(g), minsz))
        }
        else {
            dup = which(g == minG)
            cms = utils::combn(length(g), tab[minG])
            del = apply(cms, 2, setequal, dup)
            if (paired) {
                cms = cms[, order(del, decreasing = TRUE), drop = FALSE]
                cms[] = gi[c(cms)]
                cms
            }
            else {
                cms[, !del, drop = FALSE]
            }
        }
    }
    if (paired) {
        bct = tapply(seq_along(G), block, combFun, simplify = TRUE)
        nperm = ifelse(any(bct < 0), -1, prod(bct))
        if (nperm < 0 || nperm > NI) {
            btab = tapply(seq_along(G), block, `[`, simplify = FALSE)
            bSamp = function(gi) {
                g = G[gi]
                tab = table(g)
                bsz = length(g)
                minsz = tab[minG]
                cms = do.call(cbind, replicate(NI, sample.int(bsz,
                                                              minsz), simplify = FALSE))
                cms[] = gi[c(cms)]
                cms
            }
            combidx = do.call(rbind, lapply(btab, bSamp))
        }
        else {
            bcomb = tapply(seq_along(G), block, combFun, countn = FALSE,
                           simplify = FALSE)
            colb = expand.grid(lapply(bcomb, function(x) 1:ncol(x)))[-1,
                                                                     , drop = FALSE]
            combidx = mapply(function(x, y) x[, y, drop = FALSE],
                             bcomb, colb, SIMPLIFY = FALSE)
            combidx = do.call(rbind, combidx)
        }
    }
    else {
        nperm = combFun(seq_along(G))
        if (nperm < 0 || nperm > NI) {
            combidx = do.call(cbind, replicate(NI, sample.int(Glen,
                                                              minGSZ), simplify = FALSE))
        }
        else {
            combidx = combFun(seq_along(G), countn = FALSE)
        }
    }
    NI = ncol(combidx)

    deINgs = intersect(rownames(esetm), unlist(gslist))
    gslistINesetm = lapply(gslist, match, table = deINgs, nomatch = 0)
    MSabsT <- MSTop <- matrix(NA, length(gslistINesetm), NI +
                                  1)
    gsScoreFun <- function(G, block) {
        force(G)
        force(block)
        if (ite > 1) {
            G = bigG
            G[combidx[, ite - 1]] = minG
            G = factor(G)
        }
        if (paired) {
            design <- stats::model.matrix(~0 + G + block)
            colnames(design) <- substr(colnames(design), 2, 100)
        }
        else {
            design <- stats::model.matrix(~0 + G)
            colnames(design) <- levels(G)
        }
        fit <- limma::lmFit(esetm, design)
        cont.matrix <- limma::makeContrasts(contrasts = "d-c", levels = design)
        fit2 <- limma::contrasts.fit(fit, cont.matrix)
        fit2 <- limma::eBayes(fit2)
        aT1 <- limma::topTable(fit2, coef = 1, number = topSigNum)
        aT1$ID = rownames(aT1)
        de = abs(aT1$t)
        names(de) <- aT1$ID
        degf = scale(cbind(de, de * gf[names(de)]))
        rownames(degf) = names(de)
        degf = degf[deINgs, , drop = FALSE]
        sapply(gslistINesetm, function(z) {
            X = stats::na.omit(degf[z, , drop = FALSE])
            colMeans(X, na.rm = TRUE) * sqrt(nrow(X))
        })
    }
    if (parallel && requireNamespace("doParallel", quietly = TRUE) &&
        requireNamespace("parallel", quietly = TRUE)) {
        ncores = parallel::detectCores()
        if (!is.null(ncr))
            ncores = min(ncores, ncr)

        clust = parallel::makeCluster(ncores)

        doParallel::registerDoParallel(clust)
        tryCatch({
            parRes = foreach::foreach(ite = 1:(NI + 1), .combine = "c",
                             .packages = "limma") %dorng% {
                                 Sres <- gsScoreFun(G, block)
                                 tmp <- list(t(Sres))
                                 names(tmp) <- ite
                                 if (verbose && (ite%%10 == 0)) {
                                     cat(ite, "/", NI, "\n")
                                 }
                                 tmp
                             }
            parRes = do.call(cbind, parRes[order(as.numeric(names(parRes)))])
            evenCol = (1:ncol(parRes))%%2 == 0
            MSabsT[] = parRes[, !evenCol]
            MSTop[] = parRes[, evenCol]
            rm(parRes)
        }, finally = parallel::stopCluster(clust))
    }
    else {
        if (parallel)
            message("Execute in serial! Packages 'doParallel' and 'parallel' \n                       needed for parallelization!")
        for (ite in 1:(NI + 1)) {
            Sres <- gsScoreFun(G, block)
            MSabsT[, ite] <- Sres[1, ]
            MSTop[, ite] <- Sres[2, ]
            if (verbose && (ite%%10 == 0)) {
                cat(ite, "/", NI, "\n")
            }
        }
    }
    meanAbsT0 = MSabsT[, 1]
    padog0 = MSTop[, 1]
    plotIte = min(NI, 21)
    MSabsT_raw = MSabsT
    MSTop_raw = MSTop
    MSabsT = scale(MSabsT)
    MSTop = scale(MSTop)
    mff = function(x) {
        mean(x[-1] > x[1], na.rm = TRUE)
    }
    PSabsT = apply(MSabsT, 1, mff)
    PSTop = apply(MSTop, 1, mff)
    PSabsT[PSabsT == 0] <- 1/NI/100
    PSTop[PSTop == 0] <- 1/NI/100

    if (!is.null(gs.names)) {
        myn = gs.names
    }
    else {
        myn = names(gslist)
    }
    SIZE = unlist(lapply(gslist, function(x) {
        length(intersect(rownames(esetm), x))
    }))
    res = data.frame(Name = myn, ID = names(gslist), Size = SIZE,
                     meanAbsT0, padog0, PmeanAbsT = PSabsT, Ppadog = PSTop,
                     stringsAsFactors = FALSE)
    ord = order(res$Ppadog, -res$padog0)
    res = res[ord, ]
    res
}



#' Load previous pathway analyses.
#'
#' @param gse_names Character vector of GSE names.
#' @param data_dir String specifying directory for GSE folders.
#'
#' @return Result of previous call to \code{\link{diff_path}}.
#' @export
#'
#' @examples
#'
#' library(lydata)
#'
#' # location of data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # gather GSE names
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load previous pathway analyses
#' # path_anals <- load_path(gse_names, data_dir)

load_path <- function(gse_names, data_dir = getwd()) {

    anals <- list()
    for (gse_name in gse_names) {
        gse_dir <- file.path(data_dir, gse_name)

        # get paths
        pattern <- paste0(".*_diff_path", ".rds")
        paths   <- list.files(gse_dir, pattern, full.names = TRUE)

        for (path in paths) {
            anal_name <- stringr::str_extract(path, "GSE[0-9]+.GPL[0-9]+")

            if (is.na(anal_name))
                anal_name <- gse_name

            anals[[anal_name]] <- readRDS(path)
        }
    }
    return (anals)
}




#' Pathway p-value meta analysis.
#'
#' Uses Fisher's method to combine p-values from PADOG pathway analyses.
#'
#' Permutation p-values are determined by shuffling pathway names associated with PADOG p-values prior
#' to meta-analysis. Permutation p-values are then adjusted using the Benjamini & Hochberg method to  obtain
#' false discovery rates.
#'
#' @importFrom foreach foreach %dopar%
#'
#' @param path_anals Previous result of \code{\link{diff_path}}, which can
#'    be reloaded using \code{\link{load_path}}.
#' @param ncores Number of cores to use. Default is all available.
#' @param nperm Number of permutation to perform to calculate p-values.
#' @param by_source Should seperate meta-analyses be performed for each tissue
#'    source added with \code{\link{add_sources}}?
#'
#' @return A list of matrices, one for each tissue source. Each matrix contains
#'    a column of PADOG p-values for each contrast and permutation p- and fdr-values
#'    for the meta analysis.
#' @export
#' @seealso \code{\link[metap]{sumlog}}, \code{\link[PADOG]{padog}}.
#'
#' @examples
#'
#' library(lydata)
#'
#' # location of data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # gather GSE names
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load previous pathway analyses
#' # path_anals <- load_path(gse_names, data_dir)
#'
#' # perform pathway meta analysis
#' # path_res <- path_meta(path_anals, ncores = 1, nperm = 100)
#'
path_meta <- function(path_anals, ncores = parallel::detectCores(), nperm = ncores*10000, by_source = FALSE) {


    path_meta_src <- function(path_anals, src_name, ncores, nperm) {
        # bindings to pass check
        i = NULL

        # get padog tables
        path_anals <- unlist(unname(lapply(path_anals, `[[`, 'padog_tables')), recursive = FALSE)

        # same order
        ord <- row.names(path_anals[[1]])
        names(ord) <- path_anals[[1]]$Name
        path_anals <- lapply(path_anals, function(res) res[ord, ])

        # matrix of padog pvals
        pval_mat <- lapply(path_anals, `[[`, 'Ppadog')
        pval_mat <- do.call(cbind, pval_mat)
        row.names(pval_mat) <- names(ord)

        if (ncol(pval_mat) == 1) {
            # fdr values are Ppadog
            fdr  <- pval_mat[,1]
            return(cbind(pval_mat, fdr))

        } else {
            # get metap sumlog pvalues
            mpval_fun <- function(row) metap::sumlog(row[!is.na(row)])$p
            mpval     <- apply(pval_mat, 1, mpval_fun)

            # permute results
            cat(paste0("Calculating permutation p-values (", nperm, " permutations, sources: ", src_name, ") ... \n"))

            cl <- parallel::makeCluster(ncores)
            doParallel::registerDoParallel(cl)

            bins <- split(1:nperm, sort(1:nperm%%ncores))

            resl <- foreach::foreach(i=1:min(ncores, length(bins))) %dopar% {

                bin <- bins[[i]]
                B   <- rep(0, length(mpval))

                for (j in bin) {
                    # scramble pvals for each column
                    perm_mat <- apply(pval_mat, 2, sample)

                    # get perm metap sumlog pvalues
                    mpval_perm <- apply(perm_mat, 1, mpval_fun)

                    # accumulate counts where
                    B <- B + (mpval_perm <= mpval)
                }
                B
            }
            parallel::stopCluster(cl)

            # add results from each core
            B <- Reduce(`+`, resl)

            # permutation p and fdr values
            pval <- sort((B+1)/(nperm+1))
            fdr  <- stats::p.adjust(pval, 'BH')

            pval_mat <- pval_mat[names(fdr), ]
            return(cbind(pval_mat, pval, fdr))
        }
    }

    if (by_source) {

        # check for sources
        null_sources <- sapply(path_anals, function(anal) is.null(anal$sources))
        if (any(null_sources)) {
            message("Sources missing from path_anals (to add, use add_sources then re-run diff_path).\nContinuing with by_source = FALSE.")
            path_res <- list(all = path_meta_src(path_anals, 'all', ncores, nperm))

        } else {
            anals_src <- list(all = path_anals)
            anals_src <- c(anals_src, setup_src(path_anals, "padog_tables"))
            path_res <- mapply(path_meta_src, anals_src, names(anals_src),
                               MoreArgs = list(ncores = ncores, nperm = nperm))
        }

    } else {
        path_res <- list(all = path_meta_src(path_anals, 'all', ncores, nperm))
    }

    return(path_res)
}
