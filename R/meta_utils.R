#' Effect size combination meta anaylsis.
#'
#' Builds on GeneMeta implementation by allowing for genes that were not
#' measured in all studies.
#'
#' In addition to allowing for genes that were not measured in all studies, this
#' method uses moderated unbiased effect sizes calculated by \code{metaMA} and
#' determines false discovery rates using \code{fdrtool}.
#'
#' @param diff_exprs Result of previous call to \code{diff_expr}.
#' @param cutoff Minimum fraction of contrasts that must have measured each gene.
#'    Between 0 and 1.
#'
#' @return A list with two named data.frames. The first ('filt') has all the
#'    columns below for genes present in cutoff or more fraction of contrasts.
#'    The second ('raw') has only \code{dprime} and \code{vardprime} columns but
#'    for all genes (NAs for genes not measured by a given contrast).
#'    \item{dprime}{Unbiased effect sizes (one column per contrast).}
#'    \item{vardprime}{Variances of unbiased effect sizes (one column per contrast).}
#'    \item{mu}{Overall mean effect sizes.}
#'    \item{var}{Variances of overall mean effect sizes.}
#'    \item{z}{Overall z score = \code{mu / sqrt(var)}.}
#'    \item{fdr}{False discovery rates calculated by \code{fdrtool}.}
#'
#' @export
#' @seealso \code{\link[GeneMeta]{zScores}}, \code{\link[metaMA]{effectsize}},
#'    \code{\link{fdrtool}}.
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
#' # load previous analysis
#' anals <- load_diff(gse_names, data_dir)
#'
#' # perform meta-analysis
#' es <- es_meta(anals)

es_meta <- function(diff_exprs, cutoff = 0.3) {

    # get dp and vardp
    es <- get_es(diff_exprs, cutoff)
    df <- es$filt

    dp  <- df[, seq(1, ncol(df), 2)]
    var <- df[, seq(2, ncol(df), 2)]

    # get Cochran Q statistic
    Q <- f.Q(dp, var)

    # get tau (between study variance)
    tau <- tau2.DL(Q,
                   num.studies = apply(var, 1, function(x) sum(!is.na(x))),
                   my.weights  = 1 / var)

    # add tau to vardp then calculate mean effect sizes and variance
    var <- var + tau
    df$mu  <- mu.tau2(dp, var)
    df$var <- var.tau2(var)

    # get z-score and fdr
    df$z   <- df$mu/sqrt(df$var)
    df$fdr <- fdrtool::fdrtool(df$z, plot = FALSE, verbose = FALSE)$qval

    es$filt <- df[order(df$fdr), ]

    return(es)
}


# ---------------------


# Get dprimes and vardprimes values for each contrast.
#
# @inheritParams es_meta
# @return data.frame with dprime and vardprime values.

get_es <- function(diff_exprs, cutoff = 0.3) {

    # add dprimes and vardprimes to top tables
    diff_exprs <- add_es(diff_exprs)

    # get top tables
    es <- lapply(diff_exprs, function(study) study$top_tables)
    es <- unlist(es, recursive = FALSE)

    # get desired top table columns
    es <- lapply(es, function(top) {
        top$SYMBOL <- row.names(top)
        top[, c("SYMBOL", "dprime", "vardprime")]
    })

    # merge dataframes
    es <- merge_dataframes(es)

    # only keep genes where more than cutoff fraction of studies have data
    filt <- apply(es, 1, function(x) sum(!is.na(x))) >= (ncol(es) * cutoff)

    return(list(filt = es[filt, ], raw = es))
}


# ---------------------


# Add metaMA effectsize values to top tables.
#
# Used internally by \code{setup_combo_data} and \code{\link[crossmeta]{es_meta}}
# to add moderated unbiased standardised effect sizes (dprimes) to top tables
# from differential expression analysis.
#
# @param diff_exprs Result from call to \code{\link[crossmeta]{diff_expr}}.
# @param cols Columns from \code{\link[metaMA]{effectsize}} result to add to
#    top tables.
#
# @export
# @seealso \link[crossmeta]{diff_expr}, \link[crossmeta]{es_meta}.
#
# @return diff_exprs with specified columns added to top_tables for each contrast.
#
# @examples
# library(crossmeta)
# library(lydata)
#
# # location of raw data
# data_dir <- system.file("extdata", package = "lydata")
#
# # load previous analysis for eset
# anal <- load_diff("GSE9601", data_dir)
#
# # add dprime and vardprime to top tables
# anal <- add_es(anal)

add_es <- function(diff_exprs, cols = c("dprime", "vardprime")) {

    for (i in seq_along(diff_exprs)) {

        # get study degrees of freedom and group classes
        study <- diff_exprs[[i]]

        df <- study$ebayes_sv$df.residual + study$ebayes_sv$df.prior
        classes <- Biobase::pData(study$eset)$group

        for (con in names(study$top_tables)) {
            # get group names for contrast
            groups <- gsub("GSE.+?_", "", con)
            groups <- strsplit(groups, "-")[[1]]

            # get sample sizes for groups
            ni <- sum(classes == groups[2])
            nj <- sum(classes == groups[1])

            # bind effect size values with top table
            tt <- study$top_tables[[con]]
            es <- metaMA::effectsize(tt$t, ((ni * nj)/(ni + nj)), df)[, cols, drop = FALSE]
            tt <- cbind(tt, es)

            # store results
            study$top_tables[[con]] <- tt
        }
        diff_exprs[[i]] <- study
    }
    return(diff_exprs)
}


# ---------------------


# Merge a list of data.frames.
#
# @param ls List of data.frames.
# @param key Column to merge data.frames on.
#
# @return A merged data.frame with \code{key} set to row names.


merge_dataframes <- function(ls, key = "SYMBOL") {

    # ensure non 'by' names are not duplicated
    ls = Map(function(x, i)
        stats::setNames(x, ifelse(names(x) %in% key,
                                  names(x),
                                  sprintf('%s.%d', names(x), i))),
        ls, seq_along(ls))

    # merge list
    res <- Reduce(function(...) merge(..., by = key, all = TRUE), ls)

    # format result
    row.names(res) <- res[, key]
    res[, key] <- NULL
    return(res)
    print(class(res))
}

# ---------------------

# Modifed f.Q from GeneMeta (allows NAs)
#
# Compute Cochran's Q statistic. Allows genes that were not measured in all
# studies.
#
# @param dadj Dataframe of unbiased effect sizes (dprimes) for each contrast.
# @param varadj Dataframe of variances for unbiased effect sizes (vardprimes)
#    for each contrast.
#
# @return A vector of length equal to the number of rows of dadj with the Q
#    statistics.

f.Q <- function (dadj, varadj) {
    w <- 1/varadj
    tmp1 <- w * dadj
    mu <- rowSums(tmp1, na.rm = TRUE)/rowSums(w, na.rm = TRUE)
    Q <- rowSums(w * (dadj - mu)^2, na.rm = TRUE)
}

# ---------------------

# Modifed tau2.DL from GeneMeta (allows NAs)

# tau2.DL is an estimation of tau in a random effects model (REM) using
# Cochran's Q statistic. Allows genes that were not measured in all studies.
#
# @param Q A vector of Cochran's Q statistics.
# @param num.studies A vector specifying the number of experiments in which each
#    gene was measured.
# @param my.weights A matrix with one column for each experiment containing the
#    variances of the effects that should be combined.
#
# @return A vector of tau values.

tau2.DL <- function (Q, num.studies, my.weights) {
    tmp1 <- rowSums(my.weights, na.rm = TRUE)
    tmp2 <- rowSums(my.weights^2, na.rm = TRUE)
    value <- cbind((Q - (num.studies - 1))/(tmp1 - (tmp2/tmp1)), 0)
    apply(value, 1, max)
}

# ---------------------

# Modifed mu.tau2 from GeneMeta (allows NAs)
#
# Estimate overall mean effect sizes. Allows genes that were not measured in all
# studies.
#
# @param my.d A matrix, with one column for each experiment, containing the
#    effects that should be combined.
# @param my.vars.new A matrix, with one column for each experiment, containing
#    the variances of the effects that should be combined.
#
# @return A vector with the estimates of the overall mean effect sizes.

mu.tau2 <- function (my.d, my.vars.new) {
    w <- 1/my.vars.new
    tmp1 <- w * my.d
    mu <- rowSums(tmp1, na.rm = TRUE)/rowSums(w, na.rm = TRUE)
}

# ---------------------

# Modifed var.tau2 from GeneMeta (allows NAs)
#
# Estimate variances of overall mean effect sizes. Allows genes that were not
# measured in all studies.
#
# @inheritParams mu.tau2
#
# @return

var.tau2 <- function (my.vars.new) {
    w <- 1/my.vars.new
    my.var <- 1/rowSums(w, na.rm = TRUE)
}


# ---------------------


#' Contribute results of meta-analysis to public database.
#'
#' Contributed results will be used to build a freely searchable database of
#' gene expression meta-analyses.
#'
#' Performs meta-analysis on \code{diff_exprs} using \code{es_meta}. Sends
#' overall mean effect size values and minimal information needed to reproduce
#' meta-analysis.
#'
#'
#' @param diff_exprs Result of call to \code{diff_expr}.
#' @param subject String identifying meta-analysis subject (e.g. "rapamycin" or
#'    "prostate_cancer").
#'
#' @export
#'
#' @return NULL (used to contribute meta-analysis).
#'
#' @examples
#' library(lydata)
#'
#' # location of data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # gather GSE names
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load differential expression analyses
#' anals <- load_diff(gse_names, data_dir)
#'
#' # contribute results of meta-analysis
#' # contribute(anals, subject = "LY294002")

contribute <- function(diff_exprs, subject) {

    # get pdata
    pcols <- c("treatment", "group", "pairs")
    pdata <- lapply(diff_exprs, function(x) pData(x$eset)[, pcols])

    # get contrasts
    cons  <- lapply(diff_exprs, function(x) colnames(x$ebayes_sv$contrasts))

    # get effect size values
    es <- es_meta(diff_exprs)
    mu <- es$filt$mu
    names(mu) <- row.names(es$filt)

    # put it all together
    meta_info <- list(pdata = pdata, contrasts = cons, effectsize = mu)

    # upload to dropbox
    tstamp    <- format(Sys.time(), "%Y%m%d_%H%M%S_")
    save_name <- paste0(tstamp, subject, ".rds")

    saveRDS(meta_info, save_name)
    rdrop2::drop_upload(dtoken = token, save_name)
    message("Thank you for your contribution!")
    file.remove(save_name)
    return(NULL)
}
