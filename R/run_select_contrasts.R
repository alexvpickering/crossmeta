#' Shiny gadget to upload groups and select contrasts
#'
#' @param eset ExpressionSet
#' @param gse_name GEO accession for the series.
#' @param prev Previous result of \code{diff_expr}. Used to allow rechecking previous selections.
#' @param app_dir Directory to shiny app. For local development use 'inst/select_contrasts'. Default is
#'   in 'select_contrasts' sub directory of crossmeta package.
#'
#' @return result of \link{setup_prev}. Used to specify sample groups and contrasts for differential expression analysis.
#' 
#' library(lydata)
#' # location of data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' # gather GSE names
#' gse_name  <- "GSE9601"
#'
#' # load previous analysis
#' eset <- load_raw(gse_name, data_dir)[[1]]
#' run_select_contrasts(eset, gse_name)
#' 
run_select_contrasts <- function(eset, 
                                 gse_name,
                                 prev,
                                 app_dir = system.file('select_contrasts', package = 'crossmeta', mustWork = TRUE),
                                 port = 3838) {
  
  # pass arguments to app through options then run
  shiny::shinyOptions(eset = eset, gse_name = gse_name, prev = prev)
  
  # auto-reload if update app files
  options(shiny.autoreload = TRUE)
  shiny::runGadget(shinyAppDir(app_dir),port = port,  viewer = browserViewer(), stopOnCancel = FALSE)
}
