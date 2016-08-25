#' Differential expression of esets.
#'
#' User selects contrasts, then surrogate variable analysis (sva) and
#' differential expression analysis (limma) is performed.
#'
#' For each GSE, analysis results are saved in the corresponding GSE
#' folder (in \code{data_dir}) that was created by \code{get_raw}. If analysis
#' needs to be repeated, previous results can be reloaded with \code{load_diff}
#' and supplied to the \code{prev_anals} parameter. In this case, previous
#' selections/names will be reused.
#'
#' @import Biobase shiny miniUI
#' @importFrom BiocGenerics annotation
#'
#' @param esets List of annotated esets. Created by \code{load_raw}.
#' @param data_dir String specifying directory of GSE folders.
#' @param annot String, column name in fData common to all esets. For duplicated
#'   values in this column, the row with the highest interquartile range
#'   across selected samples will be kept. If meta-analysis will follow, appropriate
#'   values are "SYMBOL" (default - for gene level analysis) or, if all esets are
#'   from the same platform, "PROBE" (for probe level analysis).
#' @param prev_anals Previous result of \code{diff_expr}. If Present, previous
#'   selections and names will be reused.
#'
#' @export
#' @seealso \code{\link{get_raw}}, \code{\link{load_raw}}, and
#'   \code{\link{load_diff}}.
#'
#' @return List of lists (one per GSE), each containing:
#'   \item{eset}{Expression set without expression or feature data.
#'      \code{Treatment} ('ctl' or 'test') and \code{group} columns have been
#'      added to the \code{pData} slot. Only selected samples kept.}
#'   \item{top_tables}{List with results of \code{\link{topTable}} call (one per
#'      contrast). These results account for the effects of nuissance variables
#'      discovered by surrogate variable analysis.}
#'   \item{ebayes_sv}{Results of call to \code{\link{eBayes}} with surrogate
#'      variables included in the model matrix.}
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
#' anals <- diff_expr(esets[1], data_dir, prev_anals = prev)


diff_expr <- function (esets, data_dir = getwd(),
                       annot = "SYMBOL", prev_anals = list(NULL)) {

    # check for annot column
    chk <- sapply(esets, function(x) annot %in% colnames(fData(x)))

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
        gse_dir <- paste(data_dir, gse_folder, sep = "/")

        # select contrasts
        cons <- add_contrasts(eset, gse_name, prev_anal)
        if (is.null(cons)) next

        # setup for differential expression
        setup <- diff_setup(cons$eset, cons$levels, gse_name)


        # remove rows with duplicated/NA annot (SYMBOL or ENTREZID)
        dups <- tryCatch (
            {iqr_duplicates(cons$eset, setup$mod, setup$svobj, annot)},

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

diff_setup <- function(eset, group_levels, gse_name){

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


    # surrogate variable analysis
    # remove duplicated rows (from 1:many PROBE:SYMBOL) as affect sva
    expr <- unique(data.table(exprs(eset)))
    svobj <- tryCatch (
        {utils::capture.output(svobj <- sva::sva(as.matrix(expr), mod, mod0)); svobj},

        error = function(cond) {
            message(gse_name, ": sva failed - continuing without.")
            return(list("sv" = NULL))
        })

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



# Removes features with duplicated annotation.
#
# For rows with duplicated annot, highested IQR retained.
#
# @inheritParams diff_expr
# @inheritParams diff_setup
# @param mod Model matrix without surrogate variables. generated by \code{diff_setup}.
# @param svobj Result from \code{sva} function called during \code{diff_setup}.
#
# @return List with:
#    \item{eset}{Expression set with unique features at probe or gene level.}
#    \item{exprs_sva}{Expression data from eset with effect of surrogate
#       variable removed.}

iqr_duplicates <- function (eset, mod, svobj, annot = "SYMBOL") {

    # for R CMD check
    iqrange = SYMBOL = NULL

    # get eset with surrogate variables modeled out
    exprs_sva <- clean_y(exprs(eset), mod, svobj$sv)

    # add inter-quartile ranges, row, and feature data to exprs data
    data <- as.data.frame(exprs_sva)
    data$iqrange <- matrixStats::rowIQRs(exprs_sva)
    data$row <- 1:nrow(data)
    data[, colnames(fData(eset))] <- fData(eset)

    # remove rows with NA annot (occurs if annot is SYMBOL)
    data <- data[!is.na(data[, annot]), ]

    # for rows with same annot, keep highest IQR
    data <- data.table(data)
    data <- data[, .SD[which.max(iqrange)], by = eval(annot)]

    # use row number to keep selected features
    eset <- eset[data$row, ]
    exprs_sva <- exprs_sva[data$row, ]

    # use annot for feature names
    featureNames(eset) <- fData(eset)[, annot]
    row.names(exprs_sva)   <- fData(eset)[, annot]

    return (list(eset = eset, exprs_sva = exprs_sva))
}


# ------------------------


# Run limma analysis.
#
# Runs limma differential expression analysis on all contrasts selected by
# \code{add_contrasts}. Analysis performed with and without surrogate
# variables discovered by \code{diff_setup}. Also prints MDS plot and saves
# results.
#
# @param eset Annotated eset created by \code{load_raw}. Duplicate features and
#   non-selected samples removed by \code{iqr_duplicates}.
# @param exprs_sva Expression data with surrogate variables removed. Created by
#    \code{iqr_duplicates}
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
    graphics::par(mar = c(5, 4, 2, 6))

    # plot MDS
    limma::plotMDS(exprs_sva, pch = 19, main = gse_name, col = colours)
    graphics::legend("topright", inset = c(-0.4, 0), legend = group_levels,
                     fill = unique(colours), xpd = TRUE, bty = "n", cex = 0.65)

    # only store phenoData (exprs and fData large)
    pdata <- pData(eset)

    # save to disk
    diff_expr <- list(pdata = pdata, top_tables = top_tables, ebayes_sv = ebayes_sv)
    save_name <- paste(gse_name, "diff_expr", tolower(annot), sep = "_")
    save_name <- paste0(save_name, ".rds")

    saveRDS(diff_expr, file = paste(gse_dir, save_name, sep = "/"))
    return (diff_expr)
}


# ------------------------


# Perform eBayes analysis from limma.
#
# Generates contrast matrix then runs eBayes analysis from limma.
#
# @param eset Annotated eset created by \code{load_raw}. Non-selected samples
#    and duplicate features removed by \code{add_contrasts} and
#    \code{iqr_duplicates}.
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



#' Load previous differential expression analysis.
#'
#' Previous runs of \code{diff_expr} are loaded.
#'
#' @param gse_names Character vector specifying GSE names to be loaded.
#' @param data_dir String specifying directory of GSE folders.
#' @param annot Level of previous analysis (e.g. "SYMBOL" or "PROBE").
#'
#' @export
#' @seealso \code{\link{diff_expr}}.
#' @return Result of previous call to \code{diff_expr}.
#' @examples
#' library(lydata)
#'
#' data_dir <- system.file("extdata", package = "lydata")
#' gse_names<- c("GSE9601", "GSE34817")
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
