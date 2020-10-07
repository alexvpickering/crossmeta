#' Shiny gadget to upload groups and select contrasts
#'
#' @param eset ExpressionSet
#' @param app_dir Directory to shiny app. For local development use 'inst/select_contrasts'. Default is
#'   in 'select_contrasts' sub directory of crossmeta package.
#'
#' @return result of \link{setup_prev}. Used to specify sample groups and contrasts for differential expression analysis.
#'
run_select_contrasts <- function(eset, 
                                 app_dir = system.file('select_contrasts', package = 'crossmeta', mustWork = TRUE),
                                 port = 3838) {
  
  # pass arguments to app through options then run
  shiny::shinyOptions(eset = eset)
  
  # auto-reload if update app files
  options(shiny.autoreload = TRUE)
  shiny::runGadget(shinyAppDir(app_dir),port = port,  viewer = browserViewer(), stopOnCancel = FALSE)
}
