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
#' @param svanal Use surrogate variable analysis? Default is \code{TRUE}.
#'
#' @export
#'
#' @return List of named lists, one for each GSE. Each named list contains:
#'   \item{pdata}{data.frame with phenotype data for selected samples.
#'      Columns \code{treatment} ('ctrl' or 'test'), \code{group}, and \code{pair} are
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
                       annot = "SYMBOL", prev_anals = list(NULL), svanal = TRUE) {

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
        prev <- prev_anals[[i]]

        gse_folder <- strsplit(gse_name, "\\.")[[1]][1]  # name can be "GSE.GPL"
        gse_dir <- file.path(data_dir, gse_folder)
      
        # select groups/contrasts
        if (is.null(prev)) prev <- select_contrast(eset)
        if (is.null(prev)) next
        
        # add groups from selection
        eset <- match_prev_eset(eset, prev)
        
        # possibly subset two-channel to be like one-channel
        eset <- ch2_subset(eset, prev)
        
        # group/contrast info from previous analysis
        contrasts <- colnames(prev$ebayes_sv$contrasts)
        group_levels <- unique(eset$group)
        
        # run surrogate variable analysis
        sva_mods <- get_sva_mods(eset@phenoData)
        svobj <- run_sva(sva_mods, eset, svanal)

        # add surrogate variable/pair adjusted ("clean") expression matrix for iqr_replicates
        eset <- add_adjusted(eset, svobj)

        # remove rows with duplicated/NA annot (SYMBOL or ENTREZID)
        eset <- iqr_replicates(eset, annot)

        # differential expression
        anal <- diff_anal(eset, svobj, contrasts, group_levels, gse_dir, gse_name, annot)

        anals[[gse_name]] <- anal
    }
    return (anals)
}


# ---------------------
#' Add expression data adjusted for pairs/surrogate variables
#'
#' @param eset ExpressionSet
#' @param svobj surrogate variable object
#'
#' @return eset with \code{adjusted} element added
#' 
add_adjusted <- function(eset, svobj = list(sv = NULL)) {
  
  # get mods with group and pair effects
  mods <- get_sva_mods(eset@phenoData)
  mod <- mods$mod
  mod0 <- mods$mod0
  
  # remove pairs from alternative model so that get cleaned
  pair_cols <- colnames(mod0)[-1]
  
  svs <- svobj$sv
  mod <- mod[, !colnames(mod) %in% pair_cols]
  mod.clean <- cbind(mod0[, pair_cols], svs)
  
  y <- Biobase::exprs(eset)
  adj <- clean_y(y, mod, mod.clean)
  Biobase::assayDataElement(eset, 'adjusted') <- adj
  return(eset)
}

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
    # to keep all samples (https://support.bioconductor.org/p/73107/) specify groups without contrasts
    selected_samples <- row.names(prev_anal$pdata)
    eset <- eset[, selected_samples]

    # transfer previous treatment, group, and pairs to eset
    eset$treatment <- prev$treatment
    eset$group     <- prev$group
    eset$pair      <- NA

    if ("pair" %in% colnames(prev)) {
        eset$pair <- prev$pair
    }

    return (eset)
}


#' Subset for Paired Two-Channel ExpressionSet
#' 
#' Two-channel esets use intraspotCorrelation and lmscFit so can't use duplicateCorrelation.
#' If not using one channel in contrasts (e.g. because all reference RNA) and have paired design,
#' better to treat as single channel so that can use duplicateCorrelation.
#'
#' @param eset ExpressionSet
#' @inheritParams diff_expr
#'
#' @return ExpressionSet. If two-channel, paired and one channel not used will subset to used channel.
#' @export
#'
ch2_subset <- function(eset, prev_anal) {
  cons <- colnames(prev_anal$ebayes_sv$contrasts)
  gsm_names <- colnames(eset)
  
  # treat as single-channel if two-channel, paired, and not using both channels
  # so that can use duplicateCorrelation
  ch2 <- any(grepl('_red|_green', gsm_names))
  paired <- length(unique(eset$pair)) > 1
  
  if (ch2 & paired) {
    
    sel <- strsplit(cons, '-')[[1]]
    sel.gsm <- gsm_names[eset$group %in% sel]
    sel.colrs <- gsub('^.+_(red|green)$', '\\1', sel.gsm)
    sel.colrs <- unique(sel.colrs)
    one.colr <- length(sel.colrs == 1)
    
    if (!one.colr) stop('two-color paired designs not implemented')
    
    eset <- eset[, sel.gsm]
    colnames(eset) <- gsub('^(.+)_(red|green)$', '\\1', sel.gsm)
  }
  
  return(eset)
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
        group_levels <- unique(eset$group)


    } else {
        # get contrast info from user input
        sels <- select_contrasts(gse_name, eset)

        # add pairs info to pheno data
        pData(eset)$pair <- sels$pair

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



#' Get model matrices for surrogate variable analysis
#'
#' Used by \code{add_adjusted} to create model matrix with surrogate variables.
#'
#' @param eset Annotated eset with samples selected during \code{add_contrasts}.
#'
#' @return List with model matrix(mod) and null model matrix (mod0) used for \code{sva}.
#'
get_sva_mods <- function(pdata) {
  
  # make full and null model matrix
  group_levels <- setdiff(unique(pdata$group), NA)
  group <- factor(pdata$group, levels = group_levels)
  pair <- factor(pdata$pair)
  
  contrasts.fun <- function(l)lapply(l, stats::contrasts, contrasts = FALSE)
  
  # if pairs include in alternative and null model
  if (length(pair)) {
    mod <- stats::model.matrix(~0 + group + pair, contrasts.arg = contrasts.fun(list(group=group, pair=pair)))
    mod0 <- stats::model.matrix(~1 + pair, data = group, contrasts.arg = contrasts.fun(list(pair=pair)))
    
  } else {
    mod <- stats::model.matrix(~0 + group)
    mod0 <- stats::model.matrix(~1, data = group)
  }
  
  # rename group columns
  colnames(mod)[1:length(group_levels)] <- group_levels
  
  if (length(pair)) {
    # remove non matched pairs
    pair_cols <- colnames(mod0)[-1]
    has.pair <- colSums(mod0[, pair_cols, drop = FALSE]) >= 2
    has.pair <- names(which(has.pair))
    
    # if multiple pair columns then remove first
    if (length(has.pair) > 1) has.pair <- has.pair[-1]
    
    mod <- mod[, c(group_levels, has.pair), drop = FALSE]
    mod0 <- mod0[, c("(Intercept)", has.pair), drop = FALSE]
  }
  
  return(list("mod" = mod, "mod0" = mod0))
}

#' Run surrogate variable analysis
#'
#' @param mods result of \code{get_mods}
#' @param eset ExpressionSet
#' @param svanal Should surrogate variable analysis be run? If \code{FALSE}, returns dummy result.
#'
run_sva <- function(mods, eset, svanal) {
  if (!svanal) return(list("sv" = NULL))
  
  # remove duplicated rows (from 1:many PROBE:SYMBOL) as affect sva
  PROBE <- Biobase::fData(eset)$PROBE
  
  expr <- unique(data.table::data.table(Biobase::exprs(eset), PROBE))[, PROBE := NULL]
  expr <- as.matrix(expr)
  
  set.seed(100)
  svobj <- tryCatch (
    {utils::capture.output(svobj <- sva::sva(expr, mods$mod, mods$mod0)); svobj},
    
    error = function(e) {
      message(gse_name, ": sva failed - continuing without.")
      return(list("sv" = NULL))
    })
  
  return(svobj)
}

# ------------------------
#' Removes features with replicated annotation.
#'
#' For rows with duplicated annot, highested IQR retained.
#'
#' @param mod Model matrix without surrogate variables. generated by \code{diff_setup}.
#' @param svobj Result from \code{sva} function called during \code{diff_setup}.
#' @param annot feature to use to remove replicates.
#' @param rm.dup remove duplicates (same measure, multiple ids)? Used for Pathway analysis so that doesn't treat
#'  probes that map to multiple genes as distinct measures.
#'
#' @return Expression set with unique features at probe or gene level.
#' @export
iqr_replicates <- function(eset, annot = "SYMBOL", rm.dup = FALSE) {
  
  # for R CMD check
  iqrange = SYMBOL = NULL
  
  # do less work if possible as can take seconds
  fdata <- fData(eset)
  annot.all <- fData(eset)[, annot]
  annot.na  <- is.na(annot.all)
  annot.dup <- duplicated(annot.all[!annot.na])
  
  if (!any(annot.dup)) {
    eset <- eset[!annot.na, ]
    featureNames(eset) <- fdata[!annot.na, annot]
    
  } else {
    # use sva adjusted to compute IQRs
    y <- assayDataElement(eset, 'adjusted')
    
    iqr_rows <- which_max_iqr(eset, annot, y)
    eset <- eset[iqr_rows, ]
    
    # use annot for feature names
    featureNames(eset) <- fData(eset)[, annot]
  }
  
  # remove probes that map to multiple genes (for pathway analysis)
  if (rm.dup) {
    not.dup <- !duplicated(assayDataElement(eset, 'adjusted'))
    eset <- eset[not.dup, ]
  }
  
  return(eset)
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
# @param svobj Surrogate variable analysis object returned by \code{run_sva}.
# @param contrasts Character vector generated by \code{add_contrasts}.
# @param group_levels Character vector of group names generated by
#    \code{add_contrasts}.
# @param svobj Result from \code{sva} function called during \code{diff_setup}.
# @param gse_dir String, path to directory with GSE folders.
# @param gse_name String, name of GSE.
# @param annot String, either "ENTREZID" or "SYMBOL" for probe or gene level
#   analysis respectively. If "ENTREZID", appends "_entrezid.rds" to save name.
#
# @seealso \code{\link{diff_expr}}.
# @return List, final result of \code{diff_expr}. Used for subsequent
#   meta-analysis.

diff_anal <- function(eset, svobj, contrasts, group_levels, gse_dir, gse_name, annot = "SYMBOL"){
  
    # setup model matrix with surrogate variables
    group <- eset$group
    mod <- model.matrix(~0 + group)
    colnames(mod) <- gsub('^group', '', colnames(mod))
    svmod <- svobj$sv
    svind <- seq_len(ncol(svmod))
    if (length(svind)) colnames(svmod) <- paste0('SV', svind)
    mod <- cbind(mod, svmod)

    # run lmFit and eBayes
    ebayes_sv <- fit_ebayes(eset, contrasts, mod)

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

    # plot MDS using sv/pair cleaned data
    adj <- Biobase::assayDataElement(eset, 'adjusted')
    limma::plotMDS(adj, pch = 19, main = gse_name, col = colours)
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


#' Get design matrix for two-channel array
#'
#' @param eset ExpressionSet with \code{colnames} that end in '_red' and '_green'
#'   indicating channel and \code{eset$group} indicating group membership.
#'
#' @return model matrix for use by \link[limma]{intraspotCorrelation} and \link[limma]{lmscFit}
#'
get_ch2_mod <- function(eset) {
  
  gsm_name <- gsub('_red|_green', '', colnames(eset))
  color <- gsub('^.+?_(red|green)$', '\\1', colnames(eset))
  group <- eset$group
  
  is.green <- color == 'green'
  targets <- data.frame(row.names = gsm_name[is.green],
                        Cy3 = group[is.green],
                        Cy5 = group[!is.green])
  
  targets2 <- limma::targetsA2C(targets)
  levels <- unique(targets2$Target)
  group <- factor(targets2$Target, levels=levels)
  mod <- stats::model.matrix(~0+group)
  colnames(mod) <- levels
  return(mod)
}

# Perform lmFit and eBayes analysis from limma.
#
# If paired samples, runs \code{\link[limma]{duplicateCorrelation}} to estimate intra-patient variance.
# If non-paired two-channel Agilent, runs \code{\link[limma]{intraspotCorrelation}} and \code{\link[limma]{lmscFit}}.
#
# @param eset Annotated eset created by \code{load_raw}. Non-selected samples
#    and duplicate features removed by \code{add_contrasts} and
#    \code{iqr_replicates}.
# @param contrasts Character vector of contrast names generated by
#    \code{add_contrasts}.
# @param mod Model matrix generated by \code{diff_anal}. With
#   or without surrogate variables.
#
# @return result from call to limma \code{eBayes}.

fit_ebayes <- function(eset, contrasts, mod) {
  
    pair <- eset$pair
    y <- Biobase::exprs(eset)
    
    # check for two-channel Agilent array
    ch2 <- any(grepl('_red', colnames(eset)))
    
    if(ch2) {
      # unpaired two-channel agilent
      # covert to MAList
      MA <- to_ma(y)
      
      # get two-channel design matrix
      mod <- get_ch2_mod(eset)
      
      # run fit using intra-spot correlation
      corfit <- limma::intraspotCorrelation(MA, mod)
      fit <- limma::lmscFit(MA, mod, correlation=corfit$consensus)
      
    } else if (length(pair)) {
      # paired
      corfit <- limma::duplicateCorrelation(y, mod, block = pair)
      
      # if couldn't estimate within-block correlation, model pair as fixed effect
      if (is.nan(corfit$consensus.correlation)) {
        fit <- limma::lmFit(y, mod)
        
        # if no dof, drop pairs and retry
        if (fit$df.residual == 0) {
          eset$pair <- NULL
          mod <- get_sva_mods(eset@phenoData)$mod
          fit <- limma::lmFit(y, mod)
        }
        
      } else {
        fit <- limma::lmFit(y, mod, correlation = corfit$consensus.correlation, block = pair)
      }
      
    } else {
      # not paired
      fit <- limma::lmFit(y, mod)
    }
  
    # fit contrast
    contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = mod)
    
    fit <- limma::contrasts.fit(fit, contrast_matrix)
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
