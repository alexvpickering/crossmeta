# Title
#
# @param title
# @param pdata
# @param contrasts
# @param previous
#
# @return
# @export
#
# @examples
select_samples <- function(gse_name, title, pdata, contrasts, previous) {

    gse_link <- paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
                      gse_name, sep="")
    ui <- miniPage(
        gadgetTitleBar(title, left = miniTitleBarButton("geo",
                                                        a("GEO",
                                                          href=gse_link))),
        miniTabstripPanel(
            miniTabPanel("Samples", icon = icon("table"),
                         miniContentPanel(
                             fillCol(flex = c(1, 1, 4),
                                     textInput(
                                         "group",
                                         "Group name",
                                         width = "400px"
                                     ),
                                     selectInput(
                                         "prev",
                                         NULL,
                                         choices = c("", names(previous)),
                                         width = "400px"
                                     ),
                                     DT::dataTableOutput("pdata")
                             )
                         )
            ),

            miniTabPanel("Contrasts", icon = icon("list-ol"),
                         miniContentPanel(
                             DT::dataTableOutput("contrasts")
                         )
            )
        )
    )

    server <- function(input, output, session) {



        output$contrasts <- DT::renderDataTable(
            DT::datatable(
                contrasts,
                selection = 'none',
                options = list(
                    paging = FALSE,
                    searching = FALSE
                )
            )
        )


        output$pdata <- DT::renderDataTable(
            DT::datatable(
                pdata,
                options = list(
                    paging = FALSE,
                    searching = FALSE
                )
            )
        )


        proxy = DT::dataTableProxy('pdata')


        observeEvent(input$prev, {

            if (input$prev == "") {
                DT::selectRows(proxy, NULL)
                updateTextInput(session, "group", value = "")
            }
            if (input$prev != "") {
                DT::selectRows(proxy, previous[[input$prev]])
                updateTextInput(session, "group", value = input$prev)
            }
        })


        observeEvent(input$done, {
            selected_rows <- input$pdata_rows_selected

            if (length(selected_rows) == 0 & input$group == "") {
                stopApp(NULL)

            } else if (input$group == "") {
                message("Enter group name or clear selections.")

            } else if (length(selected_rows) == 0) {
                message("Select rows or clear group name.")
            } else {
                res <- list(selected_rows)
                names(res) <- input$group
                stopApp(res)
            }
        })

    }

    runGadget(shinyApp(ui, server), viewer = paneViewer())
}
