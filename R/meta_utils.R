

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

    # add z-score and fdr
    scores$z   <- scores$mu/sqrt(scores$var)
    scores$fdr <- fdrtool::fdrtool(scores$z, plot = FALSE, verbose = FALSE)$qval

    return(scores)
}


get_scores <- function(diff_exprs, cutoff = 0.3) {

    scores <- list()

    for (study in names(diff_exprs)) {
        #get study degrees of freedom
        diff <- diff_exprs[[study]]
        df <- diff$ebayes_sv$df.residual + diff$ebayes_sv$df.prior

        scores_cons    <- list()

        for (con in names(diff$top_tables)) {
            #get sample sizes and top table for contrast
            classes <- diff$mama_data$clinicals[[con]]$treatment
            ni <- length(classes[classes == "ctrl"])
            nj <- length(classes[classes == "test"])

            tt <- diff$top_tables[[con]]

            #get dprime and vardprime
            res <-  metaMA::effectsize(tt$t, ((ni * nj)/(ni + nj)), df)
            res <- as.data.frame(res)
            res$SYMBOL <- row.names(tt)

            #store result
            scores_cons[[con]] <- res[, c("SYMBOL", "dprime", "vardprime")]
        }
        scores[[study]]    <- scores_cons


    }
    #merge dataframes
    scores <- merge_dataframes(scores)

    #only keep genes where more than cutoff fraction of studies have data
    filt <- apply(scores, 1, function(x) sum(!is.na(x))) > (ncol(scores) * 0.3)

    return(scores[filt, ])
}


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



f.Q <- function (dadj, varadj) {
    w <- 1/varadj
    tmp1 <- w * dadj
    mu <- rowSums(tmp1, na.rm = TRUE)/rowSums(w, na.rm = TRUE)
    Q <- rowSums(w * (dadj - mu)^2, na.rm = TRUE)
}


tau2.DL <- function (Q, num.studies, my.weights) {
    tmp1 <- rowSums(my.weights, na.rm = TRUE)
    tmp2 <- rowSums(my.weights^2, na.rm = TRUE)
    value <- cbind((Q - (num.studies - 1))/(tmp1 - (tmp2/tmp1)), 0)
    apply(value, 1, max)
}


mu.tau2 <- function (my.d, my.vars.new) {
    w <- 1/my.vars.new
    tmp1 <- w * my.d
    mu <- rowSums(tmp1, na.rm = TRUE)/rowSums(w, na.rm = TRUE)
}


var.tau2 <- function (my.vars.new) {
    w <- 1/my.vars.new
    my.var <- 1/rowSums(w, na.rm = TRUE)
}











