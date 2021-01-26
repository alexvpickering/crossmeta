
#' Reuse contrast selections from previous analysis.
#'
#' Transfers user-supplied selections from previous call of diff_expr.
#'
#' @param eset Annotated eset. Created by \code{load_raw}.
#' @param prev_anal One item (for eset) from previous result of \code{diff_expr}.
#'    If present, previous selections and names will be reused.
#'
#' @seealso \code{\link{diff_expr}}
#' @return Expression set with samples and pData as in prev_anal.
#' @keywords internal
#
match_prev_eset <- function(eset, prev_anal) {
  
    # retain selected samples only
    # to keep all samples (https://support.bioconductor.org/p/73107/) specify groups without contrasts
    prev <- prev_anal$pdata
    sel <- row.names(prev)[!is.na(prev$group)]
    
    if (length(sel) < ncol(eset))
      warning('Not all samples grouped: see https://support.bioconductor.org/p/73107/#73109',
              call. = FALSE)
    
    
    # for two-channel, keep both channels for lmscFit
    # order: all R, all G
    ch2 <- any(grepl('_red|_green', colnames(eset)))
    if (ch2) {
      sel <- gsub('_red|_green', '', sel)
      sel <- unique(sel)
      sel <- paste0(sel, rep(c('_red', '_green'), each = length(sel)))
    }
    
  
    eset <- eset[, sel]
    prev <- prev[sel, ]

    # transfer previous treatment, group, and pairs to eset
    eset$treatment <- prev$treatment
    eset$group     <- prev$group

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
#' @param eset Annotated ExpressionSet. Created by \code{load_raw}.
#' @param prev_anal One item (for eset) from previous result of \code{diff_expr}.
#'
#' @return ExpressionSet. If two-channel, paired and one channel not used will subset to used channel.
#' @keywords internal
#'
ch2_subset <- function(eset, prev_anal) {
  cons <- colnames(prev_anal$ebayes_sv$contrasts)
  gsm_names <- colnames(eset)
  
  # treat as single-channel if two-channel, paired, and not using both channels
  # so that can use duplicateCorrelation and not intraspotCorrelation
  ch2 <- any(grepl('_red|_green', gsm_names))
  paired <- length(unique(eset$pair)) > 1
  
  if (ch2 & paired) {
    
    sel <- unique(unlist(strsplit(cons, '-')))
    sel.gsm <- gsm_names[eset$group %in% sel]
    sel.colrs <- gsub('^.+_(red|green)$', '\\1', sel.gsm)
    sel.colrs <- unique(sel.colrs)
    one.colr <- length(sel.colrs) == 1
    
    if (!one.colr) stop('paired two-color designs that use samples from both colors not yet implemented')
    
    # keep samples with selected color and group assignment
    keep <- grepl(sel.colrs, gsm_names) & !is.na(eset$group)
    eset <- eset[, keep]
    colnames(eset) <- gsub('^(.+)_(red|green)$', '\\1', colnames(eset))
  }
  
  return(eset)
}



#' Get model matrices for surrogate variable analysis
#'
#' Used by \code{add_adjusted} to create model matrix with surrogate variables.
#'
#' @param pdata \code{data.frame} of phenotype data with column \code{'group'}
#'   and \code{'pair'} (optional).
#'
#' @return List with model matrix(mod) and null model matrix (mod0) used for \code{sva}.
#' @export
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


#' Get row indices of maximum IQR within annotation groups
#' 
#' Groups by \code{group_by} and determines row with maximum IQR.
#'
#' @param eset \code{ExpressionSet}
#' @param groub_by Column in \code{fData(eset)} to group by
#' @param x \code{matrix} of expression values to use for IQR
#'
#' @return Integer vector of row numbers representing rows with the maximum IQR after grouping by \code{group_by}
#' @keywords internal
#'
which_max_iqr <- function(eset, groub_by, x = exprs(eset)) {
    iqrange <- NULL
    
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


#' Get design matrix for two-channel array
#'
#' @param eset ExpressionSet with \code{colnames} that end in '_red' and '_green'
#'   indicating channel and \code{eset$group} indicating group membership.
#'
#' @return model matrix for use by \link[limma]{intraspotCorrelation} and \link[limma]{lmscFit}
#' @keywords internal
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



#' Adjusts expression data for surrogate variables.
#'
#' Factors out effect of surrogate variables discovered during surrogate variable
#' analysis.
#'
#' @param y Expression data of eset.
#' @param mod Full model matrix supplied to \code{sva}.
#' @param mod.clean Model matrix with factors to clean.
#'
#' @return Expression data with effects of svs removed.
#' 
#' @keywords internal
#' 
clean_y <- function(y, mod, mod.clean) {
  
  # if no factors to clean return original
  if (!ncol(mod.clean)) return(y)
  
  X = cbind(mod, mod.clean)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}


#' Load previous differential expression analyses.
#'
#' Loads previous differential expression analyses.
#'
#' @param gse_names Character vector specifying GSE names to be loaded.
#' @param data_dir String specifying directory of GSE folders.
#' @param annot Level of previous analysis (e.g. "SYMBOL" or "PROBE").
#' @inheritParams diff_expr
#'
#' @export
#' @return Result of previous call to \code{\link{diff_expr}}.
#' @examples
#' library(lydata)
#'
#' data_dir <- system.file("extdata", package = "lydata")
#' gse_names <- c("GSE9601", "GSE34817")
#' prev <- load_diff(gse_names, data_dir)

load_diff <- function(gse_names, data_dir = getwd(), annot = "SYMBOL", postfix = NULL) {

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
        pattern <- paste(".*_diff_expr", tolower(annot), sep = "_")
        if (!is.null(postfix)) pattern <- paste(pattern, postfix, sep = "_")
        pattern <- paste0(pattern, '.rds')
        
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
