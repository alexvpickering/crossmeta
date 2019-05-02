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

select_contrasts <- function(gse_name, eset, data_dir = getwd()) {

    # ------------------- Setup

    # objects we will update
    previous <- list()

    pairs <- 0

    pdata <- data.frame(Accession = sampleNames(eset),
                        Title = pData(eset)$title,
                        Pair = NA,
                        row.names = 1:ncol(eset))

    # remove accession numbers if Illumina
    if ('illum' %in% colnames(pData(eset)))
        pdata$Accession <- NULL


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
                                                 "Previous selections:",
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

                if (group == contrasts[nrow(contrasts), "Control"]) {
                    message("Group name in use for control group.")

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


# ---------------

#' Add sample source information for meta-analysis.
#'
#' User selects a tissue source for each contrast and indicates any sources that
#' should be paired. This step is required if you would like to perform source-specific
#' effect-size/pathway meta-analyses.
#'
#'
#' The \strong{Sources} tab is used to add a source for each contrast. To do so: click the
#' relevant contrast rows, search for a source in the \emph{Sample source} dropdown box,
#' and then click the \emph{Add} button.
#'
#' The \strong{Pairs} tab is used to indicate sources that should be paired
#' (treated as the same source for subsequent effect-size and pathway meta-analyses). To do
#' so: select at least two sources from the \emph{Paired sources} dropdown box,
#' and then click the \emph{Add} button.
#'
#' For each GSE, analysis results with added sources/pairs are saved in the corresponding GSE
#' folder (in \code{data_dir}) that was created by \code{\link{get_raw}}.
#'
#' @import shiny miniUI
#'
#' @param diff_exprs Previous result of \code{\link{diff_expr}}, which can
#'    be reloaded using \code{\link{load_diff}}.
#' @param data_dir String specifying directory of GSE folders.
#'
#' @return Same as \code{\link{diff_expr}} with added slots for each GSE in \code{diff_exprs}:
#'    \item{sources}{Named vector specifying selected sample source for each contrast.
#'       Vector names identify the contrast.}
#'    \item{pairs}{List of character vectors indicating tissue sources that should be
#'       treated as the same source for subsequent effect-size and pathway meta-analyses.}
#'
#' @export
#'
#' @examples
#' library(lydata)
#'
#' # load result of previous call to diff_expr:
#' data_dir  <- system.file("extdata", package = "lydata")
#' gse_names <- c("GSE9601", "GSE34817")
#' anals     <- load_diff(gse_names, data_dir)
#'
#' # run shiny GUI to add tissue sources
#' # anals <- add_sources(anals, data_dir)
#'
add_sources <- function(diff_exprs, data_dir = getwd()) {

    # get source info for each contrast
    srclist <- list('GSE' = character(0),
                    'Contrast' = character(0),
                    'Supplied' = character(0),
                    'Source'   = character(0))

    # pairs info
    added_prs <- lapply(diff_exprs, function(anal) anal$pairs)
    added_prs <- unique(unlist(added_prs, recursive = FALSE, use.names = FALSE))

    # setup inital source list
    for (i in seq_along(diff_exprs)) {

        gse_name  <- names(diff_exprs[i])
        supld_src <- c()
        anal      <- diff_exprs[[i]]
        pdata     <- anal$pdata

        # get contrast names
        contrasts <- colnames(anal$ebayes_sv$contrasts)

        # find tissue/cell type column
        samplecol <- as.character(t(pdata[1, ]))
        is_src    <- grepl("tissue:|cell type:", samplecol)

        # user added sources
        added_src <- unname(anal$source)
        if (is.null(added_src))
            added_src <- rep(NA, length(contrasts))

        # get submitter supplied sources
        # if available and not Illumina (can only guarantee pdata$title)
        if (sum(is_src) == 1 & !'illum' %in% colnames(pdata)) {
            for (con in contrasts) {
                # contrast levels
                groups <- c(gsub('^.+?-', '', con),
                            gsub('-.+?$', '', con))

                # supplied sources
                con_src <- unique(as.character(pdata[pdata$group %in% groups, is_src]))
                con_src <- gsub("tissue: |cell type: ", "", con_src)
                con_src <- paste(con_src, collapse = ', ')

                supld_src <- c(supld_src, con_src)

            }
        } else {
            supld_src <- rep("N/A", length(contrasts))
        }

        # add info to srclist
        srclist$GSE <- c(srclist$GSE, rep(gse_name, length(contrasts)))
        srclist$Contrast <- c(srclist$Contrast, contrasts)
        srclist$Supplied <- c(srclist$Supplied, supld_src)
        srclist$Source   <- c(srclist$Source,   added_src)
    }

    # get sources/pairs info from user
    srcdf  <- as.data.frame(srclist, stringsAsFactors = FALSE)
    selres <- select_sources(srcdf, added_prs)

    srcdf     <- selres$srcdf
    added_prs <- selres$added_prs

    # add to diff_exprs
    for (i in seq_along(diff_exprs)) {

        # get sources
        gse_name <- names(diff_exprs[i])
        gse_rows <- srcdf$GSE == gse_name

        anal_srcs <- srcdf$Source[gse_rows]
        names(anal_srcs) <- paste(gse_name, srcdf$Contrast[gse_rows], sep='_')

        # check if sources have any paired
        have_prd <- sapply(added_prs, function(pairs) any(anal_srcs %in% pairs))

        # add info to diff_exprs
        diff_exprs[[i]]$sources <- anal_srcs
        diff_exprs[[i]]$pairs   <- added_prs[have_prd]

        # save diff_exprs
        gse_folder <- strsplit(gse_name, "\\.")[[1]][1]  # name can be "GSE.GPL"
        gse_dir <- file.path(data_dir, gse_folder)

        save_name <- paste(gse_name, "diff_expr", tolower(diff_exprs[[i]]$annot), sep = "_")
        save_name <- paste0(save_name, ".rds")

        saveRDS(diff_exprs[[i]], file.path(gse_dir, save_name))
    }
    return(diff_exprs)
}



# Shiny GUI used by add_sources
#
# @param srcdf
# @param added_prs
#
# @return
#
# @examples
select_sources <- function(srcdf, added_prs) {

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
                                     DT::dataTableOutput("srcdf", height = "100%")
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
