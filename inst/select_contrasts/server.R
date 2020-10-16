server <- function(input, output, session) {
  
  # get arguments from calling function
  eset <- getShinyOption('eset')
  gse_name <- getShinyOption('gse_name')
  prev <- getShinyOption('prev')
  
  # title for widget
  output$title <- shiny::renderText(paste('Select Contrasts:', gse_name))
  
  observeEvent(input$done, {
    js$closeWindow()
    stopApp(res())
  })
  
  observeEvent(input$goto_geo, {
    # click genecards
      geo_link <- paste0('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', gse_name)
      runjs(paste0("window.open('", geo_link, "')"))
  })
  
  # return value in format of previously saved analysis
  res <- reactive({
    pdata <- bulkPage$pdata()
    contrasts <- bulkPage$contrasts()
    if (!length(contrasts)) return(NULL)
    
    eset$group <- pdata$`Group name`
    if (!all(is.na(pdata$Pair))) eset$pair  <- pdata$Pair
    setup_prev(list(eset), contrasts)[[1]]
  })
  
  bulkPage <- callModule(bulkPage, 'bulk', 
                         eset = eset,
                         gse_name = gse_name,
                         prev = prev)
  
}
