#' Load previous differential expression analysis.
#'
#' Previous runs of \code{diff_expr} are loaded.
#'
#' @param gse_names Character vector specifying GSE names to be loaded.
#' @param data_dir String specifying directory of GSE folders.
#' @param probe Load probe level analysis? Default loads gene level analysis.
#'
#' @export
#' @seealso \code{\link{diff_expr}}.
#' @examples
#' library(lydata)
#' library(crossmeta)
#'
#' data_dir <- system.file("extdata", package = "lydata")
#' gse_names<- c("GSE9601", "GSE34817")
#' prev <- load_diff(gse_names, data_dir)

load_diff <- function(gse_names, data_dir, probe=F) {

  anals <- list()
  for (gse_name in gse_names) {
    gse_dir <- file.path (data_dir, gse_name)

    if (probe) {
      anal_paths <- list.files(gse_dir, pattern=".*_diff_expr_probe.rds", full=T)
    } else {
      anal_paths <- list.files(gse_dir, pattern=".*_diff_expr.rds", full=T)
    }
    #load each diff path
    #multiple if more than one platform per GSE)
    for (path in anal_paths) {
      anal <- readRDS(path)
      anal_name <- strsplit(names(anal$top_tables), "_")[[1]][1]
      anals[[anal_name]] <- anal
    }
  }
  return (anals)
}


#---------------------------


#' Make MetaArray object from results of differential expression.
#'
#' Function is used to construct a MetaArray object from results of
#' \code{diff_expr}. MetaArray object allows for use of meta-analysis methods
#' present in MAMA package.
#'
#' @param diff_exprs Result from \code{diff_expr} or \code{load_diff}.
#' @param sva Use esets with effect of surrogate variables removed? If FALSE,
#'            non-adjusted esets will be used.
#' @export
#' @seealso \link[MAMA]{MetaArray-class}, \link{diff_expr}, and \link{load_diff}.
#'
#' @return MetaArray object.

make_ma <- function(diff_exprs, sva) {

    GEDM <- list()  #gene expression data matrices
    all_clinicals <- list()  #sample descriptions

    mama_datas <- lapply(diff_exprs, function(x) x$mama_data)

    for (i in seq_along(mama_datas)) {
        #add exprs to GEDM
        if (sva) {
            exprs <- lapply(mama_datas[[i]]$esets_sva, exprs)
        } else {
            exprs <- lapply(mama_datas[[i]]$esets, exprs)
        }
        GEDM <- c(GEDM, exprs)

        #add gse_clinicals to all_clinicals
        gse_clinicals <- mama_datas[[i]]$clinicals
        all_clinicals <- c(all_clinicals, gse_clinicals)
    }
    ma <- new("MetaArray", GEDM=GEDM, clinical=all_clinicals, datanames=names(GEDM))
    return(ma)
}
