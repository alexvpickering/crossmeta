#' Effect size combination meta anaylsis.
#'
#' Builds on GeneMeta implementation by allowing for genes that were not
#' measured in all studies.
#'
#' In addition to allowing for genes that were not measured in all studies, this
#' method uses moderated unbiased effect sizes calculated by metaMA and
#' determines false discovery rates using fdrtool.
#'
#' @param diff_exprs Result of previous call to \code{diff_expr}.
#' @param cutoff Minimum fraction of contrasts that must have measured each gene.
#'    Between 0 and 1.
#'
#' @return A data.frame with columns:
#'    \item{dprime}{Unbiased effect sizes(one column per contrast).}
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

es_meta <- function(diff_exprs, cutoff = 0.3) {

    # get dp and vardp
    scores <- get_scores(diff_exprs, cutoff)
    dp  <- scores[, seq(1, ncol(scores), 2)]
    var <- scores[, seq(2, ncol(scores), 2)]

    # get Cochran Q statistic
    Q <- f.Q(dp, var)

    # get tau (between study variance)
    tau <- tau2.DL(Q,
                   num.studies = apply(var, 1, function(x) sum(!is.na(x))),
                   my.weights  = 1 / var)

    # add tau to vardp then calculate mean effect sizes and variance
    var <- var + tau
    scores <- as.data.frame(scores)
    scores$mu  <- mu.tau2(dp, var)
    scores$var <- var.tau2(var)

    # get z-score and fdr
    scores$z   <- scores$mu/sqrt(scores$var)
    scores$fdr <- fdrtool::fdrtool(scores$z, plot = FALSE, verbose = FALSE)$qval

    return(scores)
}

#---------------------

# Get dprimes and vardprimes for each contrast.
#
# @inheritParams es_meta
#
# @return data.frame with dprime and vardprime values.

get_scores <- function(diff_exprs, cutoff = 0.3) {

    scores <- list()

    for (study in names(diff_exprs)) {
        #get study degrees of freedom
        diff <- diff_exprs[[study]]
        df <- diff$ebayes_sv$df.residual + diff$ebayes_sv$df.prior

        scores_cons <- list()

        for (con in names(diff$top_tables)) {
            #get sample sizes and top table for contrast
            classes <- pData(diff$eset)$treatment
            ni <- length(classes[classes == "ctrl"])
            nj <- length(classes[classes == "test"])

            tt <- diff$top_tables[[con]]

            #get dprime and vardprime
            res <-  metaMA::effectsize(tt$t, ((ni * nj)/(ni + nj)), df)
            res <- as.data.frame(res)
            res$SYMBOL <- toupper(row.names(tt))

            #store result
            scores_cons[[con]] <- res[, c("SYMBOL", "dprime", "vardprime")]
        }
        scores[[study]] <- scores_cons


    }
    #merge dataframes
    scores <- merge_dataframes(scores)

    #only keep genes where more than cutoff fraction of studies have data
    filt <- apply(scores, 1, function(x) sum(!is.na(x))) >= (ncol(scores) * cutoff)

    return(scores[filt, ])
}

#---------------------

# Merge a list of data.frames.
#
# @param ls List of data.frames.
# @param key Column to merge data.frames on.
#
# @return A merged data.frame with \code{key} set to row names.


merge_dataframes <- function(ls, key = "SYMBOL") {

    ls <- unlist(ls, recursive = FALSE)

    #ensure non 'by' names are not duplicated
    ls = Map(function(x, i)
        setNames(x, ifelse(names(x) %in% key,
                           names(x),
                           sprintf('%s.%d', names(x), i))),
        ls, seq_along(ls))

    #merge list
    res <- Reduce(function(...) merge(..., by=key, all=TRUE), ls)

    #format result
    row.names(res) <- res[, key]
    res[, key] <- NULL
    return(res)
}

#---------------------

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

#---------------------

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

#---------------------

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

#---------------------

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
