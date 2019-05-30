#' Differential expression analysis of esets.
#'
#' After selecting control and test samples for each contrast, surrogate variable
#' analysis (\code{\link[sva]{sva}}) and differential expression analysis is performed.
#'
#' The \strong{Samples} tab is used to select control and test samples for
#' each contrast. To do so: select rows for control samples, type a group name
#' in the \emph{Control group name} text input box and click the \emph{Add Group}
#' button. Repeat for test samples. While adding additional contrasts, a previous
#' control group can be quickly reselected from the \emph{Previous selections}
#' dropdown box. After control and test samples have been added for all contrasts
#' that you wish to include, click the \emph{Done} button. Repeat for all GSEs.
#'
#' Paired samples (e.g. the same subject before and after treatment) can be
#' specified by selecting sample rows to pair and then clicking \emph{Pair Samples}.
#' The author does not usually specify paired samples and instead allows surrogate
#' variable analysis to discover these inter-sample relationships from the data itself.
#'
#' The \strong{Contrasts} tab is used to view and delete contrasts that have
#' already been added.
#'
#' For each GSE, analysis results are saved in the corresponding GSE
#' folder in \code{data_dir} that was created by \code{\link{get_raw}}. If analyses
#' needs to be repeated, previous results can be reloaded with \code{\link{load_diff}}
#' and supplied to the \code{prev_anals} parameter. In this case, previous
#' selections, names, and pairs will be reused.
#'
#' @import Biobase shiny miniUI
#' @importFrom BiocGenerics annotation
#'
#' @param esets List of annotated esets. Created by \code{\link{load_raw}}.
#' @param data_dir String specifying directory of GSE folders.
#' @param annot String, column name in fData common to all esets. For duplicated
#'   values in this column, the row with the highest interquartile range
#'   across selected samples will be kept. If meta-analysis will follow, appropriate
#'   values are "SYMBOL" (default - for gene level analysis) or, if all esets are
#'   from the same platform, "PROBE" (for probe level analysis).
#' @param prev_anals Previous result of \code{\link{diff_expr}}, which can
#'    be reloaded using \code{\link{load_diff}}. If present, previous
#'   selections, names, and pairs will be reused.
#'
#' @export
#'
#' @return List of named lists, one for each GSE. Each named list contains:
#'   \item{pdata}{data.frame with phenotype data for selected samples.
#'      Columns \code{treatment} ('ctrl' or 'test'), \code{group}, and \code{pairs} are
#'      added based on user selections.}
#'   \item{top_tables}{List with results of \code{\link[limma]{topTable}} call (one per
#'      contrast). These results account for the effects of nuissance variables
#'      discovered by surrogate variable analysis.}
#'   \item{ebayes_sv}{Results of call to \code{\link[limma]{eBayes}} with surrogate
#'      variables included in the model matrix.}
#'   \item{annot}{Value of \code{annot} variable.}
#'
#' @examples
#' library(lydata)
#'
#' # location of raw data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # gather GSE names
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load first eset
#' esets <- load_raw(gse_names[1], data_dir)
#'
#' # run analysis
#' # anals <- diff_expr(esets, data_dir)
#'
#' # re-run analysis on first eset
#' prev <- load_diff(gse_names[1], data_dir)
#' # anals <- diff_expr(esets[1], data_dir, prev_anals = prev)


diff_expr <- function (esets, data_dir = getwd(),
                       annot = "SYMBOL", prev_anals = list(NULL)) {

    # within organism symbol
    if (annot == 'SPECIES') {

        # set annot to Org_SYMBOL of first eset
        eset <- esets[[1]]
        annot <- grep('^\\d+_SYMBOL$', colnames(Biobase::fData(eset)), value = TRUE)
    }

    # check for annot column
    chk <- sapply(esets, function(x) annot %in% colnames(Biobase::fData(x)))

    if (FALSE %in% chk) {
        stop(annot, " column in fData missing for esets: ",
             paste(names(which(!chk)), collapse = ", "))
    }

    prev_anals <- prev_anals[names(esets)]
    anals <- list()
    for (i in seq_along(esets)) {

        eset <- esets[[i]]
        gse_name <- names(esets)[i]
        prev_anal <- prev_anals[[i]]

        gse_folder <- strsplit(gse_name, "\\.")[[1]][1]  # name can be "GSE.GPL"
        gse_dir <- file.path(data_dir, gse_folder)

        # select contrasts
        cons <- add_contrasts(eset, gse_name, prev_anal)
        if (is.null(cons)) next

        # setup for differential expression
        setup <- diff_setup(cons$eset, cons$levels, gse_name)


        # remove rows with duplicated/NA annot (SYMBOL or ENTREZID)
        dups <- tryCatch (
            {iqr_replicates(cons$eset, setup$mod, setup$svobj, annot)},

            error = function(c) {
                message(gse_name, ": couldn't fit model - skipping GSE.",
                        " Could try non-GUI selection (see ?setup_prev).")
                return(c)
            })
        if(inherits(dups, "error")) next

        # option to not account for surrogate variables
        # to restore add 'use_sva = TRUE' as function parameter

        # if (!use_sva) setup$modsv <- setup$mod

        # differential expression
        anal <- diff_anal(dups$eset, dups$exprs_sva, cons$contrasts, cons$levels,
                          setup$modsv, gse_dir, gse_name, annot)

        anals[[gse_name]] <- anal
    }
    return (anals)
}


# ---------------------


# Reuse contrast selections from previous analysis.
#
# Transfers user-supplied selections from previous call of diff_expr.
#
# @param eset Annotated eset. Created by \code{load_raw}.
# @param prev_anal One item (for eset) from previous result of \code{diff_expr}.
#    If present, previous selections and names will be reused.
#
# @seealso \code{\link{diff_expr}}
# @return Expression set with samples and pData as in prev_anal.

match_prev_eset <- function(eset, prev_anal) {

    # retain previously selected samples only
    # TODO: keep all samples: https://support.bioconductor.org/p/73107/
    selected_samples <- row.names(prev_anal$pdata)
    eset <- eset[, selected_samples]

    # transfer previous treatment, group, and pairs to eset
    pData(eset)$treatment <- prev_anal$pdata$treatment
    pData(eset)$group     <- prev_anal$pdata$group
    pData(eset)$pairs     <- NA

    if ("pairs" %in% colnames(prev_anal$pdata)) {
        pData(eset)$pairs <- prev_anal$pdata$pairs
    }

    return (eset)
}


# ------------------------


# Select contrasts for each GSE.
#
# Function is used by \code{diff_expr} to get sample selections for each
# contrast from user.
#
# @inheritParams match_prev_eset
# @param gse_name String specifying GSE name for eset.
#
# @seealso \code{\link{diff_expr}}
# @return List with
#    \item{eset}{Expression set with selected samples only.}
#    \item{contrasts}{Character vector of contrast names.}
#    \item{levels}{Character vector of group variable names.


add_contrasts <- function (eset, gse_name, prev_anal) {

    if (!is.null(prev_anal)) {
        # re-use selections/sample labels from previous analysis
        eset <- match_prev_eset(eset, prev_anal)

        # get contrast info from previous analysis
        contrasts    <- colnames(prev_anal$ebayes_sv$contrasts)
        group_levels <- unique(prev_anal$pdata$group)


    } else {
        # get contrast info from user input
        sels <- select_contrasts(gse_name, eset)

        # add pairs info to pheno data
        pData(eset)$pairs <- sels$pairs

        if (length(sels$cons$Control) ==  0) return(NULL)

        # setup data for each contrast
        for (i in seq_along(sels$cons$Control)) {
            # get group names
            cgrp <- sels$cons[i, "Control"]
            tgrp <- sels$cons[i, "Test"]

            # get sample names
            ctrl <- sampleNames(eset)[ sels$rows[[cgrp]] ]
            test <- sampleNames(eset)[ sels$rows[[tgrp]] ]

            # add treatment/group labels to pheno
            pData(eset)[ctrl, "treatment"] <- "ctrl"
            pData(eset)[test, "treatment"] <- "test"

            pData(eset)[ctrl, "group"] <- cgrp
            pData(eset)[test, "group"] <- tgrp
        }
        # create contrast names
        contrasts <- paste(sels$cons$Test, sels$cons$Control, sep = "-")

        # store levels for group name variable
        group_levels <- names(sels$rows)

        # retain selected samples only
        # TODO: keep all samples: https://support.bioconductor.org/p/73107/
        eset <- eset[, unique(unlist(sels$rows))]
    }

    # put data together
    contrast_data <- list(eset = eset, contrasts = contrasts, levels = group_levels)

    return (contrast_data)
}



# ------------------------


# Generate model matrix with surrogate variables.
#
# Used by \code{diff_expr} to create model matrix with surrogate variables
# in order to run \code{diff_anal}.
#
# @param eset Annotated eset with samples selected during \code{add_contrasts}.
# @param group_levels Character vector of unique group names created by
#    \code{add_contrasts}.
#
# @seealso \code{\link{add_contrasts}}, \code{\link{diff_expr}}.
# @return List with model matrix(mod), model matrix with surrogate
#         variables(modsv), and result of \code{sva} function.

diff_setup <- function(eset, group_levels, gse_name, svanal = TRUE){

    # incase svanal FALSE
    svobj <- list("sv" = NULL)

    # make full and null model matrix
    group <- factor(pData(eset)$group, levels = group_levels)
    pairs <- factor(pData(eset)$pairs)

    if (length(levels(pairs)) > 1) {
        mod <- stats::model.matrix(~0 + group + pairs)
        mod0 <- stats::model.matrix(~1 + pairs)
    } else {
        mod <- stats::model.matrix(~0 + group)
        mod0 <- stats::model.matrix(~1, data = group)
    }
    colnames(mod)[1:length(group_levels)] <- group_levels

    if (svanal) {
        # surrogate variable analysis
        # remove duplicated rows (from 1:many PROBE:SYMBOL) as affect sva
        PROBE <- fData(eset)$PROBE
        expr  <- unique(data.table(exprs(eset), PROBE))[, PROBE:=NULL]
        svobj <- tryCatch (
            {utils::capture.output(svobj <- sva::sva(as.matrix(expr), mod, mod0)); svobj},

            error = function(cond) {
                message(gse_name, ": sva failed - continuing without.")
                return(list("sv" = NULL))
            })
    }

    if (is.null(svobj$sv) || svobj$n.sv ==  0) {
        svobj$sv <- NULL
        modsv <- mod
    } else {
        modsv <- cbind(mod, svobj$sv)
        colnames(modsv) <- c(colnames(mod), paste("SV", 1:svobj$n.sv, sep = ""))
    }
    return (list("mod" = mod, "modsv" = modsv, "svobj" = svobj))
}


# ------------------------



# Removes features with replicated annotation.
#
# For rows with duplicated annot, highested IQR retained.
#
# @inheritParams diff_expr
# @inheritParams diff_setup
# @param mod Model matrix without surrogate variables. generated by \code{diff_setup}.
# @param svobj Result from \code{sva} function called during \code{diff_setup}.
# @param annot feature to use to remove replicates.
# @param rm.dup remove duplicates (same measure, multiple ids)?
#
# @return List with:
#    \item{eset}{Expression set with unique features at probe or gene level.}
#    \item{exprs_sva}{Expression data from eset with effect of surrogate
#       variable removed.}

iqr_replicates <- function (eset, mod = NULL, svobj = NULL, annot = "SYMBOL", rm.dup = FALSE) {

    # for R CMD check
    iqrange = SYMBOL = NULL

    if (length(svobj) > 0) {
        # get eset with surrogate variables modeled out
        exprs_sva <- clean_y(exprs(eset), mod, svobj$sv)
    } else {
        exprs_sva <- exprs(eset)
    }
    
    iqr_rows <- which_max_iqr(eset, annot, exprs_sva)

    # use row number to keep selected features
    eset <- eset[iqr_rows, ]
    exprs_sva <- exprs_sva[iqr_rows, ]

    # use annot for feature names
    featureNames(eset) <- fData(eset)[, annot]
    row.names(exprs_sva)   <- fData(eset)[, annot]

    if (rm.dup) {
        not.dup <- !duplicated(exprs_sva)
        eset <- eset[not.dup, ]
        exprs_sva <- exprs_sva[not.dup, ]
    }

    return (list(eset = eset, exprs_sva = exprs_sva))
}

#' Get row indices of maximum IQR within annotation groups
#' 
#' Groups by \code{group_by} and determines row with maximum IQR.
#'
#' @param eset \code{ExpressionSet}
#' @param groub_by Column in \code{fData(eset)} to group by
#' @param x \code{matrix} of expression values to use for IQR
#'
#' @return Integer vector of row numbers representing rows with the maximum IQR after grouping by \code{group_by}
#' @export
#'
#' @examples
which_max_iqr <- function(eset, groub_by, x = exprs(eset)) {
    
    # add inter-quartile ranges, row, and feature data to exprs data
    data <- as.data.frame(x)
    data$iqrange <- matrixStats::rowIQRs(x)
    data$row <- 1:nrow(data)
    data[, colnames(Biobase::fData(eset))] <- Biobase::fData(eset)
    
    # remove rows with NA groub_by (occurs if groub_by is SYMBOL)
    data <- data[!is.na(data[, groub_by]), ]
    
    # for rows with same groub_by, keep highest IQR
    data <- data.table(data)
    data <- data[data[, .I[which.max(iqrange)], by = eval(groub_by)]$V1]
    
    return(data$row)
}




# ------------------------


# Run limma analysis.
#
# Runs limma differential expression analysis on all contrasts selected by
# \code{add_contrasts}. Analysis performed with and without surrogate
# variables discovered by \code{diff_setup}. Also prints MDS plot and saves
# results.
#
# @param eset Annotated eset created by \code{load_raw}. Replicate features and
#   non-selected samples removed by \code{iqr_replicates}.
# @param exprs_sva Expression data with surrogate variables removed. Created by
#    \code{iqr_replicates}
# @param contrasts Character vector generated by \code{add_contrasts}.
# @param group_levels Character vector of group names generated by
#    \code{add_contrasts}.
# @param mod, modsv Model matrix generated by \code{diff_setup}. With
#   and without surrogate variables.
# @param svobj Result from \code{sva} function called during \code{diff_setup}.
# @param gse_dir String, path to directory with GSE folders.
# @param gse_name String, name of GSE.
# @param annot String, either "ENTREZID" or "SYMBOL" for probe or gene level
#   analysis respectively. If "ENTREZID", appends "_entrezid.rds" to save name.
#
# @seealso \code{\link{diff_expr}}.
# @return List, final result of \code{diff_expr}. Used for subsequent
#   meta-analysis.

diff_anal <- function(eset, exprs_sva, contrasts, group_levels,
                      modsv, gse_dir, gse_name, annot = "SYMBOL"){

    # differential expression (surrogate variables modeled and not)
    ebayes_sv <- fit_ebayes(eset, contrasts, modsv)

    # annotate/store results
    top_tables <- list()
    contrast_names <- paste(gse_name, contrasts, sep = "_")

    for (i in seq_along(contrast_names)) {
        top_genes <- limma::topTable(ebayes_sv, coef = i, n = Inf)
        num_sig <- sum(top_genes$adj.P.Val < 0.05)
        top_tables[[contrast_names[i]]] <- top_genes
        cat (contrast_names[i], "(# p < 0.05):", num_sig, "\n")
    }
    cat("\n")

    # setup plot items
    group <- factor(pData(eset)$group, levels = group_levels)
    palette <- RColorBrewer::brewer.pal(12, "Paired")
    colours <- palette[group]

    # Add extra space to right of plot area
    graphics::par(mai = c(1, 1, 1, 1.4))

    # plot MDS
    limma::plotMDS(exprs_sva, pch = 19, main = gse_name, col = colours)
    graphics::legend("topright", inset = c(-0.18, 0), legend = group_levels,
                     fill = unique(colours), xpd = TRUE, bty = "n", cex = 0.65)

    # only store phenoData (exprs and fData large)
    pdata <- pData(eset)

    # save to disk
    diff_expr <- list(pdata = pdata, top_tables = top_tables, ebayes_sv = ebayes_sv, annot = annot)
    save_name <- paste(gse_name, "diff_expr", tolower(annot), sep = "_")
    save_name <- paste0(save_name, ".rds")

    saveRDS(diff_expr, file.path(gse_dir, save_name))
    return (diff_expr)
}


# ------------------------


# Perform eBayes analysis from limma.
#
# Generates contrast matrix then runs eBayes analysis from limma.
#
# @param eset Annotated eset created by \code{load_raw}. Non-selected samples
#    and duplicate features removed by \code{add_contrasts} and
#    \code{iqr_replicates}.
# @param contrasts Character vector of contrast names generated by
#    \code{add_contrasts}.
# @param mod Model matrix generated by \code{diff_setup}. With
#   or without surrogate variables.
#
# @return result from call to limma \code{eBayes}.

fit_ebayes <- function(eset, contrasts, mod) {
    contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = mod)
    fit <- limma::contrasts.fit(limma::lmFit(exprs(eset), mod), contrast_matrix)
    return (limma::eBayes(fit))
}



# ------------------------


# Adjusts expression data for surrogate variables.
#
# Factors out effect of surrogate variables discovered during surrogate variable
# analysis.
#
# @param y Expression data of eset.
# @param mod Full model matrix supplied to \code{sva}.
# @param svs Surrogate variables returned by \code{sva} (svobj$sv).
#
# @seealso \code{\link{get_contrast_esets}}.
# @return Expression data with effects of svs removed.

clean_y <- function(y, mod, svs) {

    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    rm(Hat)
    gc()
    P = ncol(mod)
    return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}


# ------------------------



#' Load previous differential expression analyses.
#'
#' Loads previous differential expression analyses.
#'
#' @param gse_names Character vector specifying GSE names to be loaded.
#' @param data_dir String specifying directory of GSE folders.
#' @param annot Level of previous analysis (e.g. "SYMBOL" or "PROBE").
#'
#' @export
#' @return Result of previous call to \code{\link{diff_expr}}.
#' @examples
#' library(lydata)
#'
#' data_dir <- system.file("extdata", package = "lydata")
#' gse_names <- c("GSE9601", "GSE34817")
#' prev <- load_diff(gse_names, data_dir)

load_diff <- function(gse_names, data_dir = getwd(), annot = "SYMBOL") {

    anals <- list()
    for (gse_name in gse_names) {
        gse_dir <- file.path (data_dir, gse_name)

        # need until rename all previous
        if (annot == "SYMBOL") {
            paths <- list.files(gse_dir, pattern = ".*_diff_expr.rds", full.names = TRUE)
            to    <- gsub("expr", "expr_symbol", paths, fixed = TRUE)
            file.rename(paths, to)
        }

        # get paths
        pattern <- paste0(paste(".*_diff_expr", tolower(annot), sep="_"), ".rds")
        anal_paths <- list.files(gse_dir, pattern, full.names = TRUE)

        # load each diff path
        # multiple if more than one platform per GSE)
        for (path in anal_paths) {
            anal <- readRDS(path)

            # only need pdata from prev_anal (previously saved eset)
            if ("eset" %in% names(anal)) {
                anal$pdata <- pData(anal$eset)
                anal$eset  <- NULL
                saveRDS(anal, path)
            }

            anal_name <- strsplit(names(anal$top_tables), "_")[[1]][1]
            anals[[anal_name]] <- anal
        }
    }
    return (anals)
}
