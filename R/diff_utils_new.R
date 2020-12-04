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
#' @param recheck Would you like to recheck previous group/contrast annotations? Requires 
#'   \code{prev_anals}. Default is FALSE.
#' @param postfix Optional string to append to saved results. Useful if need to run multiple
#'   meta-analyses on the same series but with different contrasts.
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
#' anals_old <- diff_expr(esets, data_dir)
#'
#' # re-run analysis on first eset
#' prev <- load_diff(gse_names[1], data_dir)
#' # anals <- diff_expr(esets[1], data_dir, prev_anals = prev)
#' 
diff_expr_new <- function(esets, data_dir = getwd(),
                          annot = "SYMBOL", prev_anals = list(NULL), svanal = TRUE, recheck = FALSE, postfix = NULL) {
  
  # within organism symbol
  if (annot == 'SPECIES') {
    
    # set annot to Org_SYMBOL of first eset
    eset <- esets[[1]]
    annot <- grep('^\\d+_SYMBOL$', colnames(Biobase::fData(eset)), value = TRUE)
  }
  
  # check for annot column
  chk <- sapply(esets, function(x) annot %in% colnames(Biobase::fData(x)))
  if (any(!chk)) {
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
    
    # setup contrasts
    # select groups/contrasts
    if (is.null(prev) | recheck) prev <- run_select_contrasts(eset, gse_name, prev)
    
    # add groups from selection
    eset <- match_prev_eset(eset, prev)
    
    # possibly subset two-channel to be like one-channel
    eset <- ch2_subset(eset, prev)
    
    # run sva
    sva_mods <- get_sva_mods(eset@phenoData)
    svobj <- run_sva_new(sva_mods, eset, svanal)
    numsv <- svobj$n.sv
    
    # run limma
    lm_fit <- run_limma(eset, annot, svobj, numsv)
    contrasts <- colnames(prev$ebayes_sv$contrasts)
    
    top_tables <- list()
    for (con in contrasts) {
      groups <- strsplit(con, '-')[[1]]
      tt <- get_top_table(lm_fit, groups)
      num_sig <- sum(tt$adj.P.Val < 0.05)
      con_name <- paste0(gse_name, '_', con)
      cat (con_name, "(# p < 0.05):", num_sig, "\n")
      top_tables[[paste0(gse_name, '_', con)]] <- tt
    }

    anals[[gse_name]] <- c(prev, list(top_tables = top_tables, annot = annot))
    
    # save to disk
    # uses prev so that saved will have _red|_green even if treated like single-channel
    diff_expr <- c(prev, list(top_tables = top_tables, annot = annot))
    
    save_name <- paste(gse_name, "diff_expr", tolower(annot), sep = "_")
    if (!is.null(postfix)) save_name <- paste(save_name, postfix, sep = "_")
    save_name <- paste0(save_name, ".rds")
    
    saveRDS(diff_expr, file.path(gse_dir, save_name))
  }
  return(anals)
}

#' Get top table
#'
#' @param lm_fit Result of \link{run_limma}
#' @param groups Test and Control group as strings.
#'
#' @return result of \link[limma]{toptable}
#' @export
get_top_table <- function(lm_fit, groups = c('test', 'ctrl'), with.es = TRUE, robust = FALSE) {
  contrast <- paste(make.names(groups[1]), make.names(groups[2]), sep = '-')
  
  ebfit <- fit_ebayes_new(lm_fit, contrast, robust = robust)
  tt <- limma::topTable(ebfit, coef = contrast, n = Inf, sort.by = 'p')
  if (with.es) tt <- add_es_new(tt, ebfit, groups = groups)
  
  return(tt)
}


#' Linear model fitting of eset with limma.
#'
#' After selecting control and test samples for a contrast, surrogate variable
#' analysis (\code{\link[sva]{sva}}) and linear model fitting with \link[limma]{lmFit} is performed.
#'
#'
#' If analyses need to be repeated, previous results can be reloaded with \code{\link[base]{readRDS}}
#' and supplied to the \code{prev_anal} parameter. In this case, previous selections will be reused.
#'
#' @param eset Annotated eset created by \code{\link[GEOkallisto]{load_seq}}.
#' @param data_dir String specifying directory of GSE folders.
#' @param annot String, column name in fData. For duplicated
#'   values in this column, the row with the highest interquartile range
#'   across selected samples will be kept. Appropriate values are \code{"SYMBOL"} (default - for gene level analysis)
#'   or \code{"ENTREZID_HS"} (for probe level analysis).
#' @param svobj Surrogate variable analysis results. Returned from \link{run_sva}.
#' @param numsv Number of surrogate variables to model.
#' @param prev_anal Previous result of \code{run_limma}. If present, previous group
#'   selections will be reused.
#' @param filter For RNA-seq. Should genes with low counts be filtered? drugseqr shiny app performs this step
#'   separately. Should be \code{TRUE} (default) if used outside of drugseqr shiny app.
#'
#' @export
#'
#' @return List with:
#'   \item{fit}{result of \code{\link[limma]{lmFit}}.}
#'   \item{mod}{\code{model.matrix} used for \code{fit}}
#'
#'
run_limma <- function (eset, annot = "SYMBOL", svobj = list('sv' = NULL), numsv = 0, filter = TRUE) {
  
  # determine if this is rna seq data
  rna_seq <- 'norm.factors' %in% colnames(Biobase::pData(eset))
  
  # filtering low counts (as in tximport vignette)
  if (filter & rna_seq) eset <- GEOkallisto::filter_genes(eset)
  
  # add vsd element for cleaning
  eset <- add_vsd(eset, rna_seq = rna_seq)
  
  # add surrogate variable/pair adjusted ("clean") expression matrix for iqr_replicates
  eset <- add_adjusted_new(eset, svobj, numsv = numsv)
  
  # remove rows with duplicated/NA annot (SYMBOL or ENTREZID)
  eset <- iqr_replicates_new(eset, annot)
  
  # lmFit
  lm_fit <- fit_lm(eset = eset,
                   svobj = svobj,
                   numsv = numsv,
                   rna_seq = rna_seq)
  return (lm_fit)
}


#' Add VST normalized assay data element to expression set
#'
#' For microarray datasets duplicates exprs slot into vsd slot.
#'
#' @param eset ExpressionSet with group column in \code{pData(eset)}
#' @param vsd_path Path to save result to. Allows skipping running transform
#'   on each load.
#'
#' @return \code{eset} with \code{'vsd'} \code{assayDataElement} added.
#' @export
add_vsd <- function(eset, rna_seq = TRUE, pbulk = FALSE, vsd_path = NULL) {
  
  # for cases where manually added (e.g. nanostring dataset)
  els <- Biobase::assayDataElementNames(eset)
  if ('vsd' %in% els) return(eset)
  
  if (!rna_seq) {
    # for microarray use exprs
    vsd <- Biobase::assayDataElement(eset, 'exprs')
    
  } else if (!is.null(vsd_path) && file.exists(vsd_path)) {
    vsd <- readRDS(vsd_path)
    
  } else if (pbulk) {
    pdata <- Biobase::pData(eset)
    dds <- DESeq2::DESeqDataSetFromMatrix(Biobase::exprs(eset), pdata, design = ~group)
    dds <- DESeq2::estimateSizeFactors(dds)
    vsd <- tryCatch(DESeq2::rlog(dds, blind = FALSE),
                    warning = function(e) {
                      if (grepl('varianceStabilizingTransformation', e$message))
                        DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
                    })
    
    vsd <- SummarizedExperiment::assay(vsd)
    if (!is.null(vsd_path)) saveRDS(vsd, vsd_path)
    
  } else {
    # rlog fast enough for most rna-seq
    vsd <- GEOkallisto::get_vsd(eset)
    vsd <- SummarizedExperiment::assay(vsd)
    if (!is.null(vsd_path)) saveRDS(vsd, vsd_path)
  }
  
  Biobase::assayDataElement(eset, 'vsd') <- vsd
  return(eset)
}


#' Add expression data adjusted for pairs/surrogate variables
#'
#' @param eset ExpressionSet
#' @param svobj surrogate variable object
#' @param numsv Number of surrogate variables to adjust for
#' @param adj_path Path to file to restore/save adjuted data from/to for future speed
#'
#' @return eset with \code{adjusted} element added
#' @export
add_adjusted_new <- function(eset, svobj = list(sv = NULL), numsv = 0, adj_path = NULL) {
  
  
  if (!is.null(adj_path) && file.exists(adj_path)) {
    Biobase::assayDataElement(eset, 'adjusted') <- readRDS(adj_path)
    return(eset)
  }
  
  # get mods with group and pair effects
  mods <- get_sva_mods(eset@phenoData)
  mod <- mods$mod
  mod0 <- mods$mod0
  
  # remove pairs from alternative model so that get cleaned
  pair_cols <- colnames(mod0)[-1]
  
  svs <- svobj$sv[, seq_len(numsv), drop = FALSE]
  
  mod <- mod[, !colnames(mod) %in% pair_cols]
  mod.clean <- cbind(mod0[, pair_cols], svs)
  
  # used DESeq::rlog transformed counts for RNA-Seq
  # for microarray this will be exprs slot
  y <- Biobase::assayDataElement(eset, 'vsd')
  
  adj <- tryCatch(
    clean_y(y, mod, mod.clean),
    error = function(e) {
      message("adjusting for sva/pairs failed - using non-adjusted.") 
      return(y)
    })
  
  Biobase::assayDataElement(eset, 'adjusted') <- adj
  if (!is.null(adj_path)) saveRDS(adj, adj_path)
  return(eset)
}


#' @param mod Model matrix without surrogate variables. generated by \code{diff_setup}.
#' @param svobj Result from \code{sva} function called during \code{diff_setup}.
#' @param annot feature to use to remove replicates.
#' @param rm.dup remove duplicates (same measure, multiple ids)?
#'
#' @return Expression set with unique features at probe or gene level.
#' @export
iqr_replicates_new <- function(eset, annot = "SYMBOL", rm.dup = FALSE, keep_path = NULL) {
  
  # for R CMD check
  iqrange = SYMBOL = NULL
  
  # do less work if possible as can take seconds
  fdata <- Biobase::fData(eset)
  annot.all <- Biobase::fData(eset)[, annot]
  annot.na  <- is.na(annot.all)
  annot.dup <- duplicated(annot.all[!annot.na])
  
  if (!is.null(keep_path) && file.exists(keep_path)) {
    keep <- readRDS(keep_path)
    eset <- eset[keep, ]
    Biobase::featureNames(eset) <- fdata[keep, annot]
    
  } else if (!any(annot.dup)) {
    eset <- eset[!annot.na, ]
    Biobase::featureNames(eset) <- fdata[!annot.na, annot]
    
  } else {
    # use sva adjusted to compute IQRs
    # use rlog transformed or RMA-normalized to compute IQRs
    y <- Biobase::assayDataElement(eset, 'adjusted')
    
    # use row number to keep selected features
    iqr_rows <- which_max_iqr(eset, annot, y)
    eset <- eset[iqr_rows, ]
    
    # use annot for feature names
    Biobase::featureNames(eset) <-  Biobase::fData(eset)[, annot]
    
    # possibly save
    if (!is.null(keep_path)) saveRDS(iqr_rows, keep_path)
  }
  
  if (rm.dup) {
    not.dup <- !duplicated(Biobase::assayDataElement(eset, 'adjusted'))
    eset <- eset[not.dup, ]
  }
  
  return(eset)
}

#' Run surrogate variable analysis
#'
#' @param mods result of \code{get_sva_mods}
#' @param eset Expression set.
#' @param rna_seq Boolean is it RNA-seq?
#'
#' @export
run_sva_new <- function(mods, eset, svanal = TRUE) {
  if (!svanal) return(list("sv" = NULL))
  
  # determine if this is rna seq data
  rna_seq <- 'norm.factors' %in% colnames(Biobase::pData(eset))
  
  # remove duplicated rows (from 1:many PROBE:SYMBOL) as affect sva
  if (rna_seq) {
    PROBE <- Biobase::fData(eset)[,1]
  } else {
    PROBE <- Biobase::fData(eset)$PROBE
  }
  
  expr <- data.frame(Biobase::exprs(eset), PROBE)
  expr <- unique(expr)
  expr$PROBE <- NULL
  expr <- as.matrix(expr)
  
  # sva or svaseq
  sva_fun <-ifelse(rna_seq, sva::svaseq, sva::sva)
  
  set.seed(100)
  svobj <- tryCatch (
    {utils::capture.output(svobj <- sva_fun(expr, mods$mod, mods$mod0))
      colnames(svobj$sv) <- paste0('SV', seq_len(svobj$n.sv))
      return(svobj)
    },
    
    error = function(e) {
      message("sva failed - continuing without.")
      return(list("sv" = NULL))
    })
  return(svobj)
}


#' Run limma analysis.
#'
#' Runs limma differential expression analysis on all contrasts selected by
#' \code{add_contrast}. Analysis performed with and without surrogate
#' variables discovered by \code{diff_setup}. Also prints MDS plot and saves
#' results.
#'
#' @param eset Annotated eset created by \code{load_raw}. Replicate features and
#'   non-selected samples removed by \code{iqr_replicates}.
#'
#' @return list with slots:
#'   * \code{fit} Result of \link[limma]{lmFit}.
#'   * \code{mod} model matrix used for fit.
#'
fit_lm <- function(eset, svobj = list(sv = NULL), numsv = 0, rna_seq = TRUE){
  
  # setup model matrix with surrogate variables
  group <- Biobase::pData(eset)$group
  mod <- model.matrix(~0 + group)
  colnames(mod) <- gsub('^group', '', colnames(mod))
  svind <- seq_len(numsv)
  svmod <- svobj$sv[, svind, drop = FALSE]
  if (length(svind)) colnames(svmod) <- paste0('SV', svind)
  mod <- cbind(mod, svmod)
  
  lm_fit <- run_lmfit(eset, mod, rna_seq)
  
  # add enids for go/kegg pathway analyses
  # add PROBE for microarray dataset
  genes <- Biobase::fData(eset)
  genes <- genes[, colnames(genes) %in% c('ENTREZID', 'PROBE'), drop = FALSE]
  lm_fit$fit$genes <- genes
  
  return(lm_fit)
}


#' Perform lmFit analysis from limma.
#'
#' If paired samples, runs \code{\link{duplicateCorrelation}} to estimate intra-patient variance.
#'
#' @param eset Annotated eset created by \code{load_raw}. Non-selected samples
#'    and duplicate features removed by \code{add_contrasts} and
#'    \code{iqr_replicates}.
#' @param mod Model matrix generated by \code{diff_setup}. With
#'   or without surrogate variables.
#' @param rna_seq Is this an RNA-seq \code{eset}? Default is \code{TRUE}.
#'
#' @return result from call to limma \code{lmFit}.
run_lmfit <- function(eset, mod, rna_seq = TRUE) {
  
  pdata <- Biobase::pData(eset)
  pair <- pdata$pair
  y <- Biobase::exprs(eset)
  if (rna_seq) lib.size <- pdata$lib.size * pdata$norm.factors
  
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
    
  } else if (length(pair) & rna_seq) {
    # rna-seq paired
    # see https://support.bioconductor.org/p/110780/ for similar
    
    # first round
    v <- limma::voomWithQualityWeights(y, mod, lib.size = lib.size)
    corfit <- limma::duplicateCorrelation(v, mod, block = pair)
    
    # second round
    fit <- NULL
    tryCatch({
      v <- limma::voomWithQualityWeights(y, mod, lib.size = lib.size, block = pair, correlation = corfit$consensus.correlation, plot = TRUE)
      corfit <- limma::duplicateCorrelation(v, mod, block = pair)
      fit <- limma::lmFit(v, mod, correlation = corfit$consensus.correlation, block = pair)
    }, error = function(e) NULL)
    
    
    # if couldn't estimate within-block correlation, model pair as fixed effect
    if (is.null(fit)) {
      fit_fun <- function() {
        v <- limma::voomWithQualityWeights(y, mod, lib.size = lib.size)
        limma::lmFit(v, design = mod)
      }
      
      mod <- get_sva_mods(eset@phenoData)$mod
      fit <- fit_fun()
      
      # if no dof, drop pairs and retry
      if (fit$df.residual[1] == 0) {
        eset$pair <- NULL
        mod <- get_sva_mods(eset@phenoData)$mod
        fit <- fit_fun()
      }
    }
    
  } else if (rna_seq) {
    # rna-seq not paired
    # get normalized lib size and voom
    v <- limma::voomWithQualityWeights(y, mod, lib.size = lib.size)
    fit  <- limma::lmFit(v, design = mod)
    
  } else if (length(pair) & !rna_seq) {
    # microarray paired
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
    # microarray not paired
    fit <- limma::lmFit(y, mod)
  }
  
  return(list(fit = fit, mod = mod))
}


#' Fit ebayes model
#'
#' @param lm_fit Result of call to \link{run_limma}
#' @param contrasts Character vector of contrasts to fit.
#'
#' @return result of \link[limma]{eBayes}
#' @export
fit_ebayes_new <- function(lm_fit, contrasts, robust = TRUE) {
  colnames(lm_fit$fit$coefficients) <- colnames(lm_fit$mod) <- make.names(colnames(lm_fit$mod))
  contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = lm_fit$mod)
  
  eb_fit <- limma::contrasts.fit(lm_fit$fit, contrast_matrix)
  eb_fit <- limma::eBayes(eb_fit, robust = robust)
  return (eb_fit)
}

#' Add metaMA effectsize values to top table.
#'
#' Adds moderated unbiased standardised effect sizes (dprimes) to top table
#' from differential expression analysis.
#'
#' @param diff_exprs Result from call.
#' @param cols Columns from \code{\link[metaMA]{effectsize}} result to add to
#'    top table.
#'
#' @export
#'
#' @return diff_exprs with specified columns added to top_table.
#'
#' @examples
#' # location of previously saved analysis
#' data_dir <- system.file('extdata', 'IBD', package = 'drugseqr')
#'
#' # load previous analysis for eset
#' anal <- readRDS(file.path(data_dir, 'diff_exprs_symbol.rds'))
#'
#' # add dprime and vardprime to top tables
#' anal <- add_es(anal)
add_es_new <- function(tt, ebfit, cols = c("dprime", "vardprime"), groups = c('test', 'ctrl')) {
  
  
  # get study degrees of freedom and group classes
  df <- ebfit$df.total
  
  # get sample sizes for groups
  ni <- sum(ebfit$design[, groups[2]])
  nj <- sum(ebfit$design[, groups[1]])
  
  # bind effect size values with top table
  es <- metaMA::effectsize(tt$t, ((ni * nj)/(ni + nj)), df)[, cols, drop = FALSE]
  tt <- cbind(tt, es)
  
  return(tt)
}
