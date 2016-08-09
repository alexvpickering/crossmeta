#' Setup selections when many samples.
#'
#' Function useful when number of samples makes manual (GUI) selection  error
#' prone and time-consuming. Particularly useful for large clinical data sets.
#'
#' @param eset List containing one expression set with pData 'group' and 'pairs'
#'    (optional) columns. Name of \code{eset} should be the GSE name.
#' @param contrasts Character vector specifying contrasts to analyse. Each
#'    contrast must take the form "B-A" where both "B" and "A" are present in
#'    \code{eset} pData 'group' column. "B" is the treatment group and "A" is
#'    the control group.
#'
#' @return List containing necessary information for \code{prev_anal} parameter
#'    of \code{\link{diff_expr}}.
#' @export
#'
#' @examples
#'
#' library(lydata)
#' library(Biobase)
#'
#' # location of raw data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # load eset
#' gse_name  <- c("GSE34817")
#' eset <- load_raw(gse_name, data_dir)
#'
#' # inspect pData of eset
#' # View(pData(eset$GSE34817))  # if using RStudio
#' head(pData(eset$GSE34817))    # otherwise
#'
#' # get group info from pData (differs based on eset)
#' group <- pData(eset$GSE34817)$characteristics_ch1.1
#'
#' # make group names concise and valid
#' group <- gsub("treatment: ", "", group)
#' group <- make.names(group)
#'
#' # add group to eset pData
#' pData(eset$GSE34817)$group <- group
#'
#' # setup selections
#' sel <- setup_prev(eset, contrasts = "LY-DMSO")
#'
#' # run differential expression analysis
#' anal <- diff_expr(eset, data_dir, prev_anal = sel)

setup_prev <- function(eset, contrasts) {

    gse_name <- names(eset)
    eset <- eset[[1]]

    # check for group info
    if (!"group" %in% colnames(pData(eset)))
        stop("'group' column missing from pData(eset)")

    # remove expression and feature data
    exprs(eset) <- matrix()
    fData(eset) <- data.frame()

    # setup treatment and keep samples used in contrasts
    used <- strsplit(contrasts, "-")
    ctrl <- sapply(used, `[`, 2)
    test <- sapply(used, `[`, 1)
    used <- c(ctrl, test)

    group <- pData(eset)$group
    pData(eset)$treatment <- ifelse(group %in% ctrl, "ctrl", "test")
    pData(eset) <- pData(eset)[group %in% used, ]

    # setup fake ebayes_sv
    eb <- list(contrasts = matrix(ncol=length(contrasts),
                                  dimnames=list("", contrasts)))

    # setup prev
    prev <- list(list("eset" = eset, "ebayes_sv" = eb))
    names(prev) <- gse_name

    return(prev)
}
