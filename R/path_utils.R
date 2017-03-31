#' Differential expression of KEGG pathways.
#'
#' Performs PADOG pathway analysis using KEGG database (downloaded Feb 2017).
#'
#' Only contrasts with more than 2 samples per group can be analysed by PADOG.
#'
#' @import ggplot2
#'
#' @param esets List of annotated esets. Created by \code{load_raw}.
#' @param prev_anals Previous result of \code{diff_expr}.
#' @param data_dir String specifying directory of GSE folders. Defaul is working directory.
#'
#' @seealso \code{\link[PADOG]{padog}} for details on PADOG pathway analysis.
#'
#' @return List of data.frames containing PADOG pathway analysis results for each contrast.
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
#' # perform pathway analysis for each contrast
#' # path_anals <- diff_path(esets, anals, data_dir)

diff_path <- function(esets, prev_anals, data_dir = getwd()) {
    # bindings to pass check
    gslist = gs.names = NULL

    prev_anals <- prev_anals[names(prev_anals) %in% names(esets)]
    esets      <- esets[names(prev_anals)]
    anals      <- list()

    # load KEGG data
    utils::data("gslist", "gs.names", package = "crossmeta", envir = environment())

    for (i in seq_along(esets)) {

        eset <- esets[[i]]
        gse_name <- names(esets)[i]
        prev_anal <- prev_anals[[i]]

        cat('Working on', paste0(gse_name, ':'), '\n')

        gse_folder <- strsplit(gse_name, "\\.")[[1]][1]  # name can be "GSE.GPL"
        gse_dir <- paste(data_dir, gse_folder, sep = "/")

        # select contrasts
        cons <- add_contrasts(eset, gse_name, prev_anal)
        contrasts <- cons$contrasts

        # remove replicates/duplicates and annotate with human ENTREZID
        dups <- iqr_replicates(cons$eset, annot = 'ENTREZID_HS', rm.dup = TRUE)
        eset <- dups$eset

        groups <- pData(eset)$group

        pts <- list()
        # padog for each contrast
        for (con in contrasts) {

            cat('    ', paste0(con, '\n'))

            # contrast levels
            ctrl <- gsub('^.+?-', '', con)
            test <- gsub('-.+?$', '', con)

            # group
            incon <- groups %in% c(ctrl, test)
            group <- ifelse(groups[incon] == ctrl, 'c', 'd')

            # other padog inputs
            esetm  <- exprs(eset[, incon])
            pairs  <- pData(eset[, incon])$pairs
            paired <- !anyNA(pairs)
            block  <- NULL
            if (paired)
                block <- pairs

            anal_name <- paste(gse_name, con, sep='_')

            # run padog
            pts[[anal_name]] <- suppressPackageStartupMessages(tryCatch(PADOG::padog(esetm, group, paired, block, gslist, gs.names = gs.names,
                                                                                     parallel = TRUE, ncr = parallel::detectCores(),
                                                                                     verbose = FALSE),
                                                                        error = function(e) NULL))

        }

        if (length(pts) == 0) next()

        anal <- list(padog_tables = pts, sources = prev_anal$sources, pairs = prev_anal$pairs)

        # save to disk
        save_name <- paste(gse_name, "diff_path", sep = "_")
        save_name <- paste0(save_name, ".rds")
        saveRDS(anal, file = paste(gse_dir, save_name, sep = "/"))

        anals[[gse_name]] <- anal
    }
    return (anals)
}


# -------------------


#' Load previous pathway analyses.
#'
#' @param gse_names Character vector of GSE names.
#' @param data_dir String specifying directory for GSE folders.
#'
#' @return Result of previous call to \code{diff_path}.
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
#' path_anals <- load_path(gse_names, data_dir)

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


# -------------------


#' Fisher's p-value combination meta analysis of PADOG pathway analyses.
#'
#' \code{\link[metap]{sumlog}} is used to perform meta-analysis of \code{\link[PADOG]{padog}} p-values.
#' Permutation p-values are determined by shuffling pathway names associated with PADOG p-values prior
#' to meta-analysis. Permutation p-values are adjusted using the Benjamini & Hochberg method to  obtain
#' false discovery rates.
#'
#' @importFrom foreach foreach %dopar%
#'
#' @param path_anals Result of previous call to \code{diff_path} or \code{load_path}.
#' @param ncores Number of cores to use. Default is all available.
#' @param nperm Number of permutation to perform to calculate p-values.
#'
#' @return Matrix with PADOG p-values for each contrast and permutation p- and fdr-values for meta analysis.
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
#' path_anals <- load_path(gse_names, data_dir)
#'
#' # perform pathway meta analysis
#' path_res <- path_meta(path_anals, ncores = 1, nperm = 100)
#'
path_meta <- function(path_anals, ncores = parallel::detectCores(), nperm = ncores*10000, by_source = FALSE) {

    if (by_source) {
        anals_src <- setup_src(path_anals, "padog_tables")
        path_res <- mapply(path_meta_src, anals_src, names(anals_src), MoreArgs = list(ncores = ncores, nperm = nperm))

    } else {
        path_res <- list(all = path_meta_src(path_anals, 'all', ncores, nperm))
    }

    return(path_res)
}

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


# -------------------



#' Plot meta analysis of pathway.
#'
#' For all transcripts in a given pathway, dprime values are plotted as grey circles for each contrast. Black
#' error bars are also plotted to indicate +/- one standard deviation from the overall mean effect sizes. The
#' parameters \code{drugs} and \code{drug_info} are used to plot
#'
#' @param path Character vector giving pathway name. Valid values are row.names of path_meta result.
#' @param es Result of call to \code{es_meta}.
#' @param drugs Character vector giving drug names to plot alongside meta analysis result.
#' @param drug_info Matrix of expression values that include columns named for each drug.
#' @param xlim Numeric vector of length 2 specifying lower and upper bounds for dprime values.
#'
#' @return plot
#' @export
#'
#' @examples
#' library(ccmap)
#' library(ccdata)
#' library(lydata)
#'
#' data_dir <- system.file("extdata", package = "lydata")
#' data(cmap_es)
#'
#' # gather GSE names
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load previous differential expression analysis
#' anals <- load_diff(gse_names, data_dir)
#'
#' # run meta-analysis
#' es <- es_meta(anals)
#'
#' # highest ranking pathway from path_meta
#' path <- 'Amino sugar and nucleotide sugar metabolism'
#'
#' # plot meta-analysis of pathway
#' plot_path(path, es, cmap_es)
#'
#' # get meta-analysis effect size values
#' dprimes <- get_dprimes(es)
#'
#' # drug with most similar transcriptional profile
#' topd <- names(query_drugs(dprimes$all$meta, cmap_es))[1]
#'
#' # drug with most similar transcriptional profile for top pathway
#' topd_path <- names(query_drugs(dprimes$all$meta, cmap_es, path=path))[1]
#'
#' # plot meta-analysis of pathway alongside top overall drug and top drug for pathway
#' plot_path(path, es, cmap_es, c(topd, topd_path))

plot_path <- function(path, es, drug_info = NULL, drugs = NULL, xlim = NULL) {
    # bindings to pass check
    gslist = gs.names = NULL

    utils::data("gslist", "gs.names", package = "crossmeta", envir = environment())

    # dprime columns
    isdp <- grepl('^dprime', colnames(es$filt))
    ndp  <- sum(isdp)

    # path symbols
    path_num <- names(gs.names)[gs.names == path]
    path_sym <- unique(names(gslist[[path_num]]))
    path_sym <- path_sym[path_sym %in% row.names(es$filt)]
    nsym <- length(path_sym)

    # default drug effect sizes and names for plot
    des <- dnm <- rep(NA, nsym*ndp)

    if (!is.null(drugs)) {
        # path symbols also must be in drug_info
        path_sym <- path_sym[path_sym %in% row.names(drug_info)]
        nsym <- length(path_sym)

        des <- dnm <- rep(NA, nsym*ndp)
        for (i in seq_along(drugs)) {
            des[(1:nsym)*ndp-i] <- drug_info[path_sym, drugs[i]]
            dnm[(1:nsym)*ndp-i] <- drugs[i]
        }
    } else {
        path <- paste0(path, '\n')
    }



    # construct data.frame
    qes <- es$filt[path_sym, ]
    mus <- sds <- rep(NA, nsym*ndp)
    mus[(1:nsym)*ndp]  <- qes$mu
    sds[(1:nsym)*ndp] <- sqrt(qes$var)

    if (is.null(xlim))
        xlim <- c(floor(min(qes[, isdp], na.rm = TRUE)), ceiling(max(qes[, isdp], na.rm = TRUE)))

    df <- data.frame(qes  = c(t(qes[, isdp])),
                     gene = as.factor(rep(path_sym, each=ndp)),
                     mus  = mus,
                     sds  = sds,
                     ymin = mus-sds,
                     ymax = mus+sds,
                     des  = des,
                     drug = dnm)

    cols <- RColorBrewer::brewer.pal(9, 'Set1')[1:length(drugs)]

    pl <-  ggplot(df, aes_string(y = 'qes', x = 'gene')) +
        geom_point(shape=1, colour = '#999999', na.rm=TRUE) +
        geom_errorbar(aes_string(x = 'gene', ymin = 'ymin', ymax = 'ymax'),
                      colour = 'black', width = 0.2, size=.5, na.rm=TRUE)

    if (!is.null(drugs))
        pl <- pl + geom_point(aes_string(y = 'des', x = 'gene', colour = 'drug'), na.rm=TRUE)

    pl <- pl +
        xlab("Transcript") +
        ylab("Dprime") +
        geom_hline(yintercept = 0, colour = '#999999') +
        scale_y_continuous(breaks=xlim[1]:xlim[2], limits = c(xlim[1], xlim[2])) +
        scale_color_manual(breaks=drugs, values = cols) +
        ggtitle(path)+
        theme(panel.background = element_rect(fill = "#f8f8f8"),
              plot.background  = element_rect(fill = "#f8f8f8", colour = "#f8f8f8"),
              legend.background = element_rect(fill = "#f8f8f8"),
              panel.grid.major = element_line(colour = "#dddddd"),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(margin=margin(5,5,10,5,"pt"), angle = 90, vjust=0.5),
              axis.text.y = element_text(margin=margin(5,5,10,5,"pt")),
              legend.position = "top",
              legend.title=element_blank())

    pl
}



