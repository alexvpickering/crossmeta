run_select_contrasts <- function(eset, 
                                 app_dir = system.file('select_contrasts', package = 'crossmeta', mustWork = TRUE),
                                 port = 3838) {
  
  # pass arguments to app through options then run
  shiny::shinyOptions(eset = eset)
  
  # auto-reload if update app files
  options(shiny.autoreload = TRUE)
  shiny::runGadget(shinyAppDir(app_dir),port = port,  viewer = browserViewer(), stopOnCancel = FALSE)
}


gse_name <- 'GSE11975'
eset <- crossmeta::load_raw(gse_name)[[1]]


app_dir <- 'inst/select_contrasts'
run_select_contrasts(eset, app_dir)
