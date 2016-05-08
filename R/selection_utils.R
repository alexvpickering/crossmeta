# Shiny gadget to select samples.
#
# @param gse_name String, used to create link to Gene Expression Omnibus.
# @param title String, displayed in gadget title bar.
# @param pdata Dataframe with phenotype data to help with selection.
# @param contrasts Dataframe with columns "Test" and "Control" specifying
#    group names for previously selected contrasts.
# @param previous Named lists of numeric vectors. List names are previous groups
#    and numeric vectors are previously selected rows. Created by concatenating
#    results of previous calls to select_samples. Used for reselection of
#    previous groups.
#
# @return Named list with a numeric vector. Name is typed group name,
#    numeric vector is selected rows.

select_contrasts <- function(gse_name, pdata) {

    #------------------- setup


    # objects we will update
    previous <- list()
    groups <- c()
    contrasts <- data.frame(Control = character(0),
                            Test = character(0), stringsAsFactors = FALSE)

    # link for GSE
    gse_link <- paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
                      gse_name, sep="")



    #------------------- user interface


    ui <- miniPage(
        # title bar
        gadgetTitleBar("Select Contrasts", left = miniTitleBarButton("geo",
                                                                     a("GEO",
                                                                       href=gse_link))),
        miniTabstripPanel(
            miniTabPanel("Samples", icon = icon("table"),
                         miniContentPanel(
                             fillCol(flex = c(NA, NA, NA, 1),
                                     textInput(
                                         "group",
                                         "Control group name:",
                                         width = "200px"
                                     ),
                                     selectInput(
                                         "prev",
                                         NULL,
                                         choices = c("", names(previous)),
                                         width = "200px"
                                     ),
                                     hr(),
                                     DT::dataTableOutput("pdata")
                             )
                         ),
                         miniButtonBlock(
                             actionButton("add", "Add Group")
                         )
            ),

            miniTabPanel("Contrasts", icon = icon("list-ol"),
                         miniContentPanel(
                             DT::dataTableOutput("contrasts")
                         ),
                         miniButtonBlock(
                             actionButton("delete", "Delete Contrast(s)")
                         )
            )
        )
    )



    #------------------------- server


    server <- function(input, output, session) {


        # show selected contrasts
        output$contrasts <- DT::renderDataTable({

            # invalidate when add/delete click
            input$add
            input$delete

            DT::datatable(
                contrasts,
                options = list(
                    paging = FALSE,
                    searching = FALSE,
                    bInfo = 0
                )
            )},
            server = FALSE
        )

        # show phenotype data
        output$pdata <- DT::renderDataTable(
            DT::datatable(
                pdata,
                options = list(
                    paging = FALSE,
                    bInfo = 0
                )
            )
        )

        # make reactive state value to keep track of ctrl vs test group
        state <- reactiveValues(ctrl = 1)

        # need pdata proxy to reselect rows of pheno data
        proxy = DT::dataTableProxy('pdata')



        #------------------- click 'Add'


        observeEvent(input$add, {

            #construct group data
            rows  <- input$pdata_rows_selected
            group <- input$group
            group_data <- list()
            group_data[[group]] <- rows


            #check for incomplete input
            if (group == "") {
                message("Enter group name.")

            } else if (length(rows) == 0) {
                message("Select rows.")


            #check for wrong input
            } else if (make.names(group) != group) {
                message("Group name invalid.")

            } else if (group %in% names(previous) &&
                       !setequal(previous[[group]], rows)) {
                message("Group name in use with different samples.")

            } else if (Position(function(x) setequal(x, rows), previous, nomatch = 0) > 0 &&
                        !group %in% names(previous)) {
                message("Selection in use with different group name.")


            #add ctrl group data to previous and groups
            } else if (state$ctrl == 1) {
                if (!group %in% names(previous))
                    previous <<- c(previous, group_data)
                groups <<- c(groups, names(group_data))

                #update inputs
                updateTextInput(session, "group",
                                label = "Test group name:", value = "")
                updateSelectInput(session, "prev",
                                  choices = c("", names(previous)))
                DT::selectRows(proxy, NULL)
                state$ctrl <- 0


            #add test group data to previous and groups
            } else {
                if (!group %in% names(previous))
                    previous <<- c(previous, group_data)
                groups <<- c(groups, names(group_data))

                #update inputs
                updateTextInput(session, "group",
                                label = "Control group name:", value = "")
                updateSelectInput(session, "prev",
                                  choices = c("", names(previous)))
                DT::selectRows(proxy, NULL)

                #add groups to contrasts table then blank groups
                contrasts[nrow(contrasts) + 1, ] <<- groups
                groups <<- c()
                state$ctrl <- 1
            }
        })



        #------------------- click 'Delete'


        observeEvent(input$delete, {
            rows  <- input$contrasts_rows_selected

            # check if no selection
            if (length(rows) == 0) {
                message("Select contrasts.")

            } else {
                #for groups in rows to delete
                for (group in unlist(contrasts[rows, ])){

                    # remove from previos if not in another row
                    if (!group %in% unlist(contrasts[-rows, ]))
                        previous[[group]] <<- NULL

                    updateSelectInput(session, "prev",
                                      choices = c("", names(previous)))
                }
                #remove rows from contrasts
                contrasts <<- contrasts[-rows, ]
                row.names(contrasts) <<- NULL
            }
        })



        #------------------- select previous


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



        #------------------- click 'Done'


        observeEvent(input$done, {
            if (state$ctrl == 0) {
                message("Need to add test group.")
            } else if (nrow(contrasts) == 0) {
                message("No contrasts selected.")
                stopApp(NULL)
            } else {
                stopApp(list(rows=previous, cons=contrasts))
            }
        })

    }

    runGadget(shinyApp(ui, server), viewer = paneViewer())
}
