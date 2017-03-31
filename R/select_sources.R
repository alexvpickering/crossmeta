


# Select sources for contrasts.
#
# @param srcdf
# @param added_prs
# @param sources
#
# @return
# @export
#
# @examples
select_sources <- function(srcdf, added_prs, sources) {

    # ------------------- Setup

    # setup
    added_src <- setdiff(srcdf$Source, NA)

    prsht  <- paste0(nrow(srcdf) * 21.32, 'px')
    srcht  <- paste0(nrow(srcdf) * 42.63, 'px')

    if (is.null(added_prs)) {
        prsdf  <- data.frame(Pairs = character(0), stringsAsFactors = FALSE)

    } else {
        prsvec <- sapply(added_prs, paste, collapse = ', ')
        prsdf  <- data.frame(Pairs = prsvec, stringsAsFactors = FALSE)
    }


    # link for GSE
    orig_names <- srcdf$GSE
    gse_names  <- sapply(strsplit(srcdf$GSE, ".", fixed = TRUE), `[`, 1)
    srcdf$GSE  <- paste0('<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',
                         gse_names, '">', gse_names, '</a>')



    # ------------------- user interface


    ui <- miniPage(
        # title bar
        gadgetTitleBar("Select Sources"),
        miniTabstripPanel(
            miniTabPanel("Sources", icon = icon("table"),
                         miniContentPanel(
                             fillCol(flex = NA,
                                     fillRow(flex = c(NA, .025, NA, NA),
                                             selectizeInput(
                                                 "source",
                                                 label = "Sample source:",
                                                 choices = c("", sources),
                                                 options = list(create = TRUE),
                                                 width = "300px"
                                             ),
                                             tags$style(type='text/css', ".selectize-dropdown-content {max-height: 500px; }")
                                     ),
                                     hr(),
                                     DT::dataTableOutput("srcdf", height = srcht)
                             )
                         ),
                         miniButtonBlock(
                             actionButton("add_source", "Add"),
                             actionButton("del_source", "Delete")
                         )

            ),

            miniTabPanel("Pairs", icon = icon("list-ol"),
                         miniContentPanel(
                             fillCol(flex = NA,
                                     fillRow(flex = c(NA, .025, NA, NA),
                                             selectizeInput(
                                                 "paired",
                                                 label = "Paired sources:",
                                                 multiple = TRUE,
                                                 choices = added_src,
                                                 width = "300px"
                                             )
                                     ),
                                     hr(),
                                     DT::dataTableOutput("prsdf", height = prsht)
                             )
                         ),
                         miniButtonBlock(
                             actionButton("add_pair", "Add"),
                             actionButton("del_pair", "Delete")
                         )
            )
        )
    )



    # ------------------------- server


    server <- function(input, output, session) {


        output$prsdf <- DT::renderDataTable({

            DT::datatable(
                prsdf,
                # rownames = FALSE,
                options = list(
                    scrollY = FALSE,
                    paging = FALSE,
                    searching = FALSE,
                    bInfo = 0
                )
            )
            })




        # show source data
        output$srcdf <- DT::renderDataTable({

            DT::datatable(
                srcdf,
                # rownames = FALSE,
                options = list(
                    scrollY = FALSE,
                    paging = FALSE,
                    bInfo = 0
                ),
                escape = FALSE # need for HTML links
            )
        })


        src_proxy = DT::dataTableProxy('srcdf')
        prs_proxy = DT::dataTableProxy('prsdf')


        # clicked 'Add Source'
        observeEvent(input$add_source, {

            rows  <- input$srcdf_rows_selected
            nrow  <- length(row)
            src   <- input$source


            if (is.null(rows) & src != "") {
                message('Select contrast(s) to add source.')

            } else if (!is.null(rows) & src == "") {
                message('Select a sample source.')

            } else if (!is.null(rows) & src != "") {
                srcdf[rows, 'Source'] <<- input$source
                DT::replaceData(src_proxy, srcdf, resetPaging = FALSE)

                # update added sources
                added_src <<- setdiff(srcdf$Source, NA)
                updateSelectizeInput(session, 'paired', choices = added_src)
                updateSelectizeInput(session, 'source', selected = "")
            }

        })

        # clicked 'Delete Source'
        observeEvent(input$del_source, {

            rows  <- input$srcdf_rows_selected

            if (length(row) == 0) {
                message('Select contrast(s) to delete source.')

            } else {
                srcdf[rows, 'Source'] <<- NA
                DT::replaceData(src_proxy, srcdf, resetPaging = FALSE)

                # update added sources
                added_src <<- setdiff(srcdf$Source, NA)
                updateSelectizeInput(session, 'paired', choices = added_src)
            }

        })

        # clicked 'Add Pair'
        observeEvent(input$add_pair, {

            prd  <- input$paired

            if (length(prd) < 2) {
                message('Select two or more sources to pair.')

            } else {

                # first pairing
                if (length(added_prs) == 0) {
                    added_prs[[length(added_prs)+1]] <<- prd

                } else {

                    # determine if previously paired contain just paired sources
                    have_prd <- sapply(added_prs, function(added_pr) any(prd %in% added_pr))

                    # merge if any
                    if (any(have_prd)) {
                        mrg_prd <- unique(c(prd, unlist(added_prs[have_prd])))
                        added_prs[have_prd] <<- NULL
                        added_prs[[length(added_prs)+1]] <<- mrg_prd

                    } else {
                        # add new otherwise
                        added_prs[[length(added_prs)+1]] <<- prd
                    }
                }

                # update pairs data.frame
                prsvec <- sapply(added_prs, paste, collapse = ', ')
                prsdf <<- data.frame(Pairs = prsvec, stringsAsFactors = FALSE)

                DT::replaceData(prs_proxy, prsdf)

                # unselect paired
                updateSelectizeInput(session, 'paired', selected = "")
            }
        })


        # clicked 'Delete Pair'
        observeEvent(input$del_pair, {


            rows  <- input$prsdf_rows_selected

            if (length(rows) == 0) {
                message('Select row(s) to delete pairs.')

            } else {
                added_prs[[rows]] <<- NULL

                # update pairs data.frame
                if (length(added_prs) == 0) {
                    prsdf <<- data.frame(Pairs = character(0), stringsAsFactors = FALSE)

                } else {
                    prsvec <- sapply(added_prs, paste, collapse = ', ')
                    prsdf <<- data.frame(Pairs = prsvec, stringsAsFactors = FALSE)

                }

                DT::replaceData(prs_proxy, prsdf,  resetPaging = FALSE)
            }
        })


        # clicked 'Done'
        observeEvent(input$done, {
            if (anyNA(srcdf$Source)) {
                message("Contrast(s) without source.")

            }  else {
                srcdf$GSE <- orig_names
                stopApp(list(srcdf = srcdf, added_prs = added_prs))
            }
        })
    }

    runGadget(shinyApp(ui, server), viewer = paneViewer())
}
