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

select_contrasts <- function(gse_name, eset) {

    # ------------------- Setup

    # objects we will update
    previous <- list()

    pairs <- 0

    pdata <- data.frame(Accession = sampleNames(eset),
                        Title = pData(eset)$title,
                        Pair = NA,
                        row.names = 1:ncol(eset))

    contrasts <- data.frame(Control = character(0),
                            Test = character(0), stringsAsFactors = FALSE)



    # link for GSE
    gse_name <- strsplit(gse_name, ".", fixed = TRUE)[[1]][1]
    gse_link <- paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
                      gse_name, sep = "")



    # ------------------- user interface


    ui <- miniPage(
        # title bar
        gadgetTitleBar("Select Contrasts", left = miniTitleBarButton("geo", "GEO")),
        miniTabstripPanel(
            miniTabPanel("Samples", icon = icon("table"),
                         miniContentPanel(
                             fillCol(flex = c(NA, NA, 1),
                                     fillRow(flex = c(NA, 0.025, NA),
                                             textInput(
                                                 "group",
                                                 "Control group name:",
                                                 width = "200px"
                                             ),
                                             br(),
                                             selectInput(
                                                 "prev",
                                                 br(),
                                                 choices = c("", names(previous)),
                                                 width = "200px"
                                             )
                                     ),
                                     hr(),
                                     DT::dataTableOutput("pdata")
                             )
                         ),
                         miniButtonBlock(
                             actionButton("add", "Add Group"),
                             actionButton("pair", "Pair Samples")
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



    # ------------------------- server


    server <- function(input, output, session) {


        # show selected contrasts
        output$contrasts <- DT::renderDataTable({

            # invalidate when add/delete click
            state$contrast

            DT::datatable(
                contrasts,
                options = list(
                    scrollY = FALSE,
                    paging = FALSE,
                    searching = FALSE,
                    bInfo = 0
                )
            )},
            server = FALSE
        )

        # show phenotype data
        output$pdata <- DT::renderDataTable({

            DT::datatable(
                isolate(loopData()),
                #rownames = FALSE,
                options = list(
                    scrollY = FALSE,
                    paging = FALSE,
                    bInfo = 0
                )
            )}
        )


        loopData = reactive({

            # invalidate if click 'Pair Samples'
            input$pair

            # initial load: set pairs to 1
            if (pairs == 0) {
                pairs <<- 1

            } else {
                rows <- input$pdata_rows_selected

                # check for valid selection
                if (length(rows) <= 1){
                    message("Select at least two samples to pair.")

                    # valid selection: update pdata and increment pairs
                } else  {
                    pdata[rows, 'Pair'] <<- pairs
                    pairs <<- pairs + 1
                }
            }
            pdata
        })


        # make reactive state value to keep track of ctrl vs test group
        state <- reactiveValues(ctrl = 1, contrast = 0)

        # need pdata proxy to reselect rows of pheno data
        proxy = DT::dataTableProxy('pdata')



        # ------------------- click 'Add'


        observeEvent(input$add, {

            # construct group data
            rows  <- input$pdata_rows_selected
            group <- input$group
            group_data <- list()
            group_data[[group]] <- rows


            # check for incomplete/wrong input
            if (group == "") {
                message("Enter group name.")

            } else if (length(rows) == 0) {
                message("Select rows.")

            } else if (make.names(group) != group) {
                message("Group name invalid.")

            } else if (group %in% names(previous) &&
                       !setequal(previous[[group]], rows)) {
                message("Group name in use with different samples.")

            } else if (Position(function(x) setequal(x, rows), previous, nomatch = 0) > 0 &&
                       !group %in% names(previous)) {
                message("Selection in use with different group name.")


            } else if (state$ctrl == 1) {
                # add ctrl group data to previous and contrasts
                if (!group %in% names(previous))
                    previous <<- c(previous, group_data)
                contrasts[nrow(contrasts) + 1, ] <<- c(group, NA)

                # update inputs
                updateTextInput(session, "group",
                                label = "Test group name:", value = "")
                updateSelectInput(session, "prev",
                                  choices = c("", names(previous)))
                DT::selectRows(proxy, NULL)

                # update states
                state$ctrl <- 0
                state$contrast <- state$contrast + 1


            } else {
                # add test group data to previous and contrasts
                if (!group %in% names(previous))
                    previous <<- c(previous, group_data)
                contrasts[nrow(contrasts), "Test"] <<- group

                # update inputs
                updateTextInput(session, "group",
                                label = "Control group name:", value = "")
                updateSelectInput(session, "prev",
                                  choices = c("", names(previous)))
                DT::selectRows(proxy, NULL)

                #update states
                state$ctrl <- 1
                state$contrast <- state$contrast + 1
            }
        })

        #------------------- click 'Pair Samples'

        observeEvent(input$pair, {
            DT::replaceData(proxy, loopData(), resetPaging = FALSE)
        })


        #------------------- click 'Delete'


        observeEvent(input$delete, {
            rows    <- input$contrasts_rows_selected
            test.na <- is.na(contrasts[rows, 'Test'])

            # check if no selection
            if (length(rows) == 0) {
                message("Select contrasts.")

            } else {
                # for groups in rows to delete
                for (group in unlist(contrasts[rows, ])){

                    # remove from previos if not in another row
                    if (!group %in% unlist(contrasts[-rows, ]))
                        previous[[group]] <<- NULL

                    updateSelectInput(session, "prev",
                                      choices = c("", names(previous)))
                }
                # remove rows from contrasts
                contrasts <<- contrasts[-rows, ]
                row.names(contrasts) <<- NULL

                # update contrast state
                state$contrast <- state$contrast + 1

                if (test.na) {
                    # put back to control group
                    updateTextInput(session, "group",
                                    label = "Control group name:", value = "")
                    state$ctrl <- 1
                }
            }
        })



        # ------------------- select previous


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



        # ------------------- click 'Done'


        observeEvent(input$done, {
            if (state$ctrl == 0) {
                message("Need to add test group.")
            } else if (nrow(contrasts) == 0) {
                message("No contrasts selected.")
                stopApp(NULL)
            } else {
                stopApp(list(rows = previous, cons = contrasts, pairs = pdata$Pair))
            }
        })


        # ------------------- click 'GEO'


        observeEvent(input$geo, {
            utils::browseURL(gse_link)
        })

    }

    runGadget(shinyApp(ui, server), viewer = paneViewer())
}
