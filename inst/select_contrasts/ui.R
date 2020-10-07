
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
    div(style = 'display: none;',
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
            addContrastInput(ns('add_contrast')),
            delContrastsInput(ns('del_contrasts'))
        )
    )
  })
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



#' Add contrast input
#' @export
#' @keywords internal
addContrastInput <- function(id) {
  ns <- NS(id)
  
  options_add <- list(maxItems = 2, placeholder = 'Select test then control group')
  tags$div(shinypanel::selectizeInputWithButtons(
    ns('select_groups'),
    'Select groups to compare:',
    actionButton(ns('add_contrast'), label = NULL, icon = icon('plus', 'fa-fw'), title = 'Add contrast'),
    options = options_add,
    container_id = ns('add_contrast_container')
  )
  )
}

#' Delete contrasts input
#' @export
#' @keywords internal
delContrastsInput <- function(id) {
  ns <- NS(id)
  
  options_del <- list(multiple = TRUE)
  tags$div(shinypanel::selectizeInputWithButtons(
    ns('select_contrasts'),
    'Select comparisons to remove:',
    actionButton(ns('del_contrasts'), label = NULL, icon = icon('minus', 'fa-fw'), title = 'Remove contrasts'),
    options = options_del,
    container_id = ns('del_contrast_container')
  )
  )
  
}

jscode <- "shinyjs.closeWindow = function() { window.close(); }"

bootstrapPage(
  shinyjs::useShinyjs(),
  shinyjs::extendShinyjs(text = jscode, functions = c("closeWindow")),
  includeScript(path = 'www/select_contrasts.js'),
  includeCSS(path = 'www/select_contrasts.css'),
  miniUI::miniPage(
    miniUI::gadgetTitleBar('Select Contrasts'),
    tags$div(bulkPageUI('bulk'), style='padding: 15px;'
    )
  )
)