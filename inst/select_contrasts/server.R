server <- function(input, output, session) {
  
  # get arguments from calling function
  eset <- getShinyOption('eset')
  
  observeEvent(input$done, {
    js$closeWindow()
    stopApp(prev())
  })
  
  observeEvent(input$cancel, {
    js$closeWindow()
    stopApp(stop("User cancel", call. = FALSE))
  })
  
  # return value in format of previously saved analysis
  prev <- reactive({
    pdata <- bulkPage$pdata()
    contrasts <- bulkPage$contrasts()
    if (!length(contrasts)) return(NULL)
    
    eset$group <- pdata$`Group name`
    if (!all(is.na(pdata$Pair))) eset$pair  <- pdata$Pair
    setup_prev(list(eset), contrasts)[[1]]
  })
  
  bulkPage <- callModule(bulkPage, 'bulk', 
                         eset = eset,
                         gse_name = gse_name)
  
}
