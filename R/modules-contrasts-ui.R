
#' UI for Bulk Data page
#' @export
#' @keywords internal
bulkPageUI <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = 'row',
        div(class = 'col-lg-4', bulkFormInput(ns('form')))
    ),
    hr(),
    div(id = ns('anal_table_container'), style = '',
        bulkAnnotInput(ns('anal')),
        bulkTableOuput(ns('explore'))
    )
  )
}

#' UI for Bulk Data annotation upload/download
#' @export
#' @keywords internal
bulkAnnotInput <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(id=ns('validate-up'), class='validate-wrapper',
             shinyWidgets::actionGroupButtons(
               inputIds = c(ns('click_dl'), ns('click_up')),
               labels = list(icon('download', 'fa-fw'), icon('upload', 'fa-fw'))
             ),
             tags$div(class='help-block', id = ns('error_msg'))
             
    ),
    # hidden dl/upload buttons
    div(style = 'display: none',
        fileInput(ns('up_annot'), '', accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
    ),
    downloadLink(ns('dl_annot'), ''),
    
    shinyBS::bsTooltip(id = ns('click_dl'), title = 'Download metadata to fill: <b>Group name</b> and <b>Pair</b> (optional)', options = list(container = 'body')),
    shinyBS::bsTooltip(id = ns('click_up'), title = 'Upload filled metadata', options = list(container = 'body'))
  )
}


#' Input form for Bulk Data page
#' @export
#' @keywords internal
bulkFormInput <- function(id) {
  ns <- NS(id)
  
  withTags({
    div(class = "well-form well-bg",
        div(id = ns('anal_dataset_panel'), style = '',
            bulkFormAnalInput(ns('anal_form'))
        )
    )
  })
}



#' Differential expression analysis inputs for bulkFormInput
#' @export
#' @keywords internal
bulkFormAnalInput <- function(id) {
  ns <- NS(id)
  
  tagList(
    bulkAnalInput(ns('ds'), label = 'Add two-group comparison:')
  )
}


#' Tables for datasets page
#' @export
#' @keywords internal
bulkTableOuput <- function(id) {
  ns <- NS(id)
  withTags({
    div(class = 'dt-container',
        DT::dataTableOutput(ns("pdata"))
    )
  })
}



#' Bulk Differential expression analysis input
#' @export
#' @keywords internal
bulkAnalInput <- function(id, with_dl = TRUE, label = 'Select groups to compare:') {
  ns <- NS(id)
  
  options <- list(maxItems = 2, placeholder = 'Select test then control group')
  if (with_dl) {
    input <- tags$div(id = 'bulk-intro-comparison',
                      shinypanel::selectizeInputWithButtons(
                        ns('contrast_groups'),
                        label,
                        actionButton(ns('plus'), label = NULL, icon = icon('plus', 'fa-fw'), title = 'Add contrast'),
                        options = options,
                        container_id = ns('run_anal_container')
                      )
    )
    
  } else {
    input <- tags$div(
      id = ns('run_anal_container'),
      shinypanel::selectizeInputWithValidation(
        ns('contrast_groups'),
        label,
        options = options
      )
      
    )
  }
  
  return(input)
}


select_contrasts <- function(eset, gse_name) {
  
  server <- function(input, output, session) {
    
    bulkPage <- callModule(bulkPage, 'bulk', 
                           eset = eset,
                           gse_name = gse_name)
    
  }
  
  
  ui <- bootstrapPage(
    shiny::tags$head(
      shiny::tags$style(".dt-container {white-space: nowrap;}"), # table text on 1 line
      shiny::tags$style(".dt-fake-height {height: 1px;}"), # to make 100% height div work
      shiny::tags$style("td.dt-nopad {padding: 0px !important; height: 100%;}"), # td for bg color group column
      shiny::tags$style("td.dt-nopad div {height: 100%; width: 100%; text-align: center;}"), # div inside td for bg color group column
      shiny::tags$style("td.dt-nopad span {display: inline-block; padding: 8px 10px; color: white;}") # span inside div inside td for bg color group column
    ),
    shinyjs::useShinyjs(),
    fluidPage(
      tags$div(bulkPageUI('bulk'), style='padding-top: 15px;'
      )
    )
  )
  
  runGadget(shinyApp(ui, server), viewer = paneViewer())
  
}

select_contrasts(eset, gse_name)
