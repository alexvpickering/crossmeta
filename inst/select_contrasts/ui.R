jscode <- "shinyjs.closeWindow = function() { window.close(); }"

bootstrapPage(
  shinyjs::useShinyjs(),
  shinyjs::extendShinyjs(text = jscode, functions = c("closeWindow")),
  includeScript(path = 'www/select_contrasts.js'),
  includeCSS(path = 'www/select_contrasts.css'),
  tags$div(
    miniUI::gadgetTitleBar('Select Contrasts', left = miniUI::miniTitleBarButton('goto_geo', 'GEO')),
    shiny::fluidPage(
      tags$div(bulkPageUI('bulk'), style='padding-top: 15px;'
               
      )
      
    )
  )
)