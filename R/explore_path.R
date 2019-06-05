#' Explore pathway meta analyses.
#'
#' Shiny app for interactively exploring the results of effect-size and pathway meta-analyses.
#' The app also interfaces with the ccmap package in order to explore drugs that are predicted
#' to reverse or mimic your signature.
#'
#' For a given tissue source (top left dropdown box) and KEGG pathway (bottom left dropdown box, ordered
#' by increasing false discovery rate), effect-sizes (y-axis) are plotted for each gene in the pathway
#' (x-axis, ordered by decreasing asbsolute effect size).
#'
#' For each gene, open circles give the effect-sizes for each contrast. The transparency of the open
#' circles is proportional to the standard deviation of the effect-size for each contrast.
#' For each gene, error bars give one standard deviation above and below the the overall meta-analysis
#' effect-size.
#'
#' The top drugs for the full signature in a given tissue (top right dropdown box, red points) and
#' just the pathway genes (bottom right dropdown box, blue points) are orderered by decreasing
#' (if \code{type} is 'both' or 'mimic') or increasing (if \code{type} is 'reverse') similarity.
#' Positive and negative pearson correlations correspond to drugs that, respectively, mimic and
#' reverse the query signature.
#'
#' Drug effect sizes can be made visible by either clicking the legend entries (top left of plot) or
#' selecting a new drug in the dropdown boxes.
#'
#' When a new tissue source or pathway is selected, the top drug and pathway dropdown boxes
#' are approriately updated.
#'
#'
#' @import plotly
#'
#' @param es_res Result of call to \code{\link{es_meta}}.
#' @param path_res Result of call to \code{\link{path_meta}}.
#' @param drug_info Matrix of differential expression values for drugs (rows are genes, columns are drugs).
#'    If NULL (default), \code{\link[ccdata]{cmap_es}} is used.
#' @param type Desired direction of drug action on query signature (see details). One of either \code{'both'} (Default), 
#'   \code{'mimic'}, or \code{'reverse'}.
#'
#' @return None
#' @export
#'
#' @examples
#' library(lydata)
#'
#' data_dir  <- system.file("extdata", package = "lydata")
#' gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#'
#' # load result of previous call to diff_expr:
#' es_anals <- load_diff(gse_names, data_dir)
#'
#' # run shiny GUI to add tissue sources
#' # es_anals <- add_sources(es_anals, data_dir)
#'
#' # perform effect-size meta-analyses for each tissue source
#' es_res <- es_meta(es_anals, by_source = TRUE)
#'
#' # load result of previous call to diff_path:
#' # path_anals <- load_path(gse_names, data_dir)
#'
#' # perform pathway meta-analyses for each tissue source
#' # path_res <- path_meta(path_anals, ncores = 1, nperm = 100, by_source = TRUE)
#'
#' # explore pathway meta-analyses
#' # explore_paths(es_res, path_res)
#'
explore_paths <- function(es_res, path_res, drug_info = NULL, type = c('both', 'mimic', 'reverse')) {
    # global binding to pass CHK
    cmap_es = NULL

    if (is.null(drug_info)) {
        utils::data('cmap_es', package = "ccdata", envir = environment())
        drug_info <- cmap_es
        rm(cmap_es)
    }

    i <- 1

    # visibility
    restart <- FALSE
    vis1 <- FALSE
    vis2 <- FALSE

    # get top drugs for full signature and pathway genes only
    fullpath <- function(query_genes, drug_info, path, type, all_full = NULL) {

        # all drugs (full signature)
        if (is.null(all_full))
            all_full <- ccmap::query_drugs(query_genes, drug_info)

        # all drugs (pathway genes only)
        all_path <- ccmap::query_drugs(query_genes, drug_info, path = path)

        # limit to top 50
        if (type == 'mimic')   {
            top_full <- utils::head(all_full, 50)
            top_path <- utils::head(all_path, 50)
        } else if (type == 'reverse') {
            top_full <- utils::head(rev(all_full), 50)
            top_path <- utils::head(rev(all_path), 50)
        } else {
            top_full <- c(utils::head(all_full, 25), utils::tail(all_full, 25))
            top_path <- c(utils::head(all_path, 25), utils::tail(all_path, 25))
        }

        # construct names as Drug (full cor, path cor)
        full_names <- paste0(names(top_full), ' (',
                             round(top_full, 2), ', ',
                             round(all_path[names(top_full)], 2), ')')

        path_names <- paste0(names(top_path), ' (',
                             round(all_full[names(top_path)], 2), ', ',
                             round(top_path, 2), ')')

        # values as drug names
        top_full <- names(top_full)
        names(top_full) <- full_names

        top_path <- names(top_path)
        names(top_path) <- path_names

        return(list(full=top_full, path=top_path))

    }
    dprimes  <- ccmap::get_dprimes(es_res)

    sources      <- names(es_res)
    paths        <- row.names(path_res[[1]])
    names(paths) <- paste0(paths, ' (', format_fdrs(path_res[[1]][, 'fdr']), ')')
    top_drugs    <- fullpath(dprimes[[1]]$meta, drug_info, paths[1], type[1])


    ui <- shinyUI(fluidPage(
        tags$head(
            shiny::tags$style("#kegg {border-top-left-radius: 0; border-bottom-left-radius: 0; position: relative; z-index: 2; margin-left: -8px; margin-top: -26px;}")
        ),
        title = "Pathway Explorer",
        br(),

        fluidRow(
            column(6, align = 'center',
                   strong("Source:"),
                   selectInput('source', '', choices = sources, selected = sources[1], width = '400px'),
                   tags$style(type='text/css', ".selectize-dropdown-content {max-height: 180px; }"),
                   br(),
                   hr(),
                   br(),
                   strong("Pathway (FDR):"),
                   tags$div(
                       tags$div(style = "display:inline-block;", selectInput('path', '', choices = paths, selected = paths[1], width = '359.5px')),
                       shinyBS::bsButton('kegg', label = '', icon = icon('external-link-alt', 'fa-fw'), style='default', title = 'Go to KEGG')
                   )
            ),
            column(6, align = 'center',
                   HTML("<strong>Top Drugs for Full Signature (&rho;<sub>full</sub>, &rho;<sub>path</sub>):</strong>"),
                   selectInput('drug1', '', top_drugs$full, top_drugs$full[1], width = '400px'),
                   br(),
                   hr(),
                   br(),
                   HTML("<strong>Top Drugs for Pathway Genes Only (&rho;<sub>full</sub>, &rho;<sub>path</sub>):</strong>"),
                   selectInput('drug2', '', top_drugs$path, top_drugs$path[1], width = '400px')
            )
        ),

        hr(),

        plotly::plotlyOutput('trendPlot', height = "650px")
    ))


    server <- function(input, output, session) {


        # update source
        observeEvent(input$source, {

            # cat(i, ': source changing\n')
            # i <<- i + 1

            # reset drug visibility
            restart <<- TRUE

            new_src   <- input$source
            old_path  <- input$path
            old_drug2 <- input$drug2

            # update top drugs for new source
            paths        <<- row.names(path_res[[new_src]])
            names(paths) <<- paste0(paths, ' (', format_fdrs(path_res[[new_src]][, 'fdr']), ')')
            top_drugs    <<- fullpath(dprimes[[new_src]]$meta, drug_info, paths[1], type[1])

            # don't re-update top drugs
            state$source <<- TRUE

            # update path selection (updates drug selections as well)
            updateSelectInput(session, 'path',  choices = paths, selected = paths[1])

            # if path doesn't change with source
            if (old_path == paths[1]) {

                # cat(i, ': no path change (source)\n')
                # i <<- i + 1

                # need to update drugs
                updateSelectInput(session, 'drug1', choices = top_drugs$full, selected = top_drugs$full[1])
                updateSelectInput(session, 'drug2', choices = top_drugs$path, selected = top_drugs$path[1])

                # reset source
                # cat(i, ': resetting source state (source)\n')
                # i <<- i + 1

                state$source <<- FALSE

                # if drug2 doesn't change
                if (old_drug2 == top_drugs$path[1]) {

                    # reset visibility
                    restart <<- FALSE
                    vis1 <<- FALSE
                    vis2 <<- FALSE

                    # trigger plot update
                    state$update <<- !state$update

                    # cat(i, ': no drug2 change (source)\n')
                    # i <<- i + 1

                }
            }
        }, ignoreInit = TRUE)




        # update path
        observeEvent(input$path, {

            # cat(i, ': path change\n')
            # i <<- i + 1

            # reset drug visibility
            restart <<- TRUE

            old_src   <- input$source
            old_drug2 <- input$drug2

            # if source changed, already updated top path drugs
            if (!state$source) {
                # cat(i, ': updating top path drugs\n')
                # i <<- i + 1

                top_drugs <<- fullpath(dprimes[[old_src]]$meta, drug_info, input$path, type[1], top_drugs$all_full)

            } else {
                # cat(i, ': resetting source state (path)\n')
                # i <<- i + 1

                # reset source
                state$source <<- FALSE
            }

            # update drug selections
            updateSelectInput(session, 'drug1', choices = top_drugs$full, selected = top_drugs$full[1])
            updateSelectInput(session, 'drug2', choices = top_drugs$path, selected = top_drugs$path[1])


            # if drug2 doesn't change
            if (old_drug2 == top_drugs$path[1]) {

                # reset visibility
                restart <<- FALSE
                vis1 <<- FALSE
                vis2 <<- FALSE

                # trigger plot update
                state$update <<- !state$update

                # cat(i, ': no drug2 change (path)\n')
                # i <<- i + 1
            }

        }, ignoreInit = TRUE)




        # drug 1 changed
        observeEvent(input$drug1, {

            if (!restart) {
                # cat(i, ': drug1 change (not restart)\n')
                # i <<- i + 1

                vis1  <<- TRUE
                vis2  <<- FALSE

                if (input$drug1 == input$drug2)
                    vis2 <<- TRUE

                # trigger plot update
                state$update <<- !state$update

            } else {

                # cat(i, ': drug1 change (restart)\n')
                # i <<- i + 1
            }
        }, ignoreInit = TRUE)




        # drug 2 changed
        observeEvent(input$drug2, {

            if (restart) {
                # cat(i, ': drug2 change (restart)\n')
                # i <<- i + 1

                # reset visibility
                restart <<- FALSE
                vis1 <<- FALSE
                vis2 <<- FALSE

            } else {
                # cat(i, ': drug2 change (not restart)\n')
                # i <<- i + 1

                vis1 <<- FALSE
                vis2 <<- TRUE

                if (input$drug1 == input$drug2)
                    vis1 <<- TRUE
            }

            # trigger plot update
            state$update <<- !state$update

        }, ignoreInit = TRUE)


        state <- reactiveValues(update = FALSE, source = FALSE)


        dataset <- eventReactive(state$update, {

            # cat(i, ': updating dataset\n\n')
            # i <<- 1

            get_dfs(input$path, es_res[[input$source]], drug_info, c(input$drug1, input$drug2), gslist, gs.names)
        })


        output$trendPlot <- renderPlotly({
            
            dfs <- dataset()
            drugs <- levels(dfs$drug_df$drug)
            drug_ids <- seq_along(unique(drugs))+2
            
            g <-  ggplot2::ggplot(data = dfs$sumry_df,
                                  ggplot2::aes_string(x = 'gene', y = 'mus')) +
                ggplot2::geom_point(ggplot2::aes_string('gene', 'dprime', alpha = 'sdinv'),
                                    dfs$query_df, shape=1, colour = '#666666', na.rm=TRUE) +
                ggplot2::geom_errorbar(ggplot2::aes_string(ymin = 'low', ymax = 'high'),
                                       colour = 'black', width = 1.3) +
                ggplot2::geom_point(ggplot2::aes_string(x = 'gene',  y = 'dprime', colour = 'drug'),
                           dfs$drug_df, na.rm=TRUE) +
                ggplot2::ylab("Dprime") +
                ggplot2::xlab("") +
                ggplot2::geom_hline(yintercept = 0, colour = '#999999') +
                ggplot2::scale_color_manual(values = c("#E41A1C", "#6BAED6")) +
                ggplot2::scale_y_continuous(breaks = function(ylims) floor(ylims[1]):ceiling(ylims[2])) +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust=0.5),
                               legend.title = ggplot2::element_blank())


            pl <- plotly_build(g)


            # adjust margins and ranges
            # ymin <- max(-8, min(dfs$query_df$dprime, na.rm = TRUE))
            # ymax <- min(8, max(dfs$query_df$dprime, na.rm = TRUE))

            pl$x$layout$margin$b <- 90
            nsyms <- nrow(dfs$sumry_df)
            pl$x$layout$xaxis$range <- c(0.5, min(50.5, nsyms+0.5))
            # pl$x$layout$yaxis$range <- c(ymin, ymax)

            # hide modebar
            pl$x$config$displayModeBar <- FALSE

            # query hoverinfo
            query_genes <- paste0('Gene: ', as.character(dfs$query_df$gene))
            query_dps   <- paste0('Dprime: ', round(dfs$query_df$dprime, 1))
            query_sds   <- paste0('SD: ', round(1/dfs$query_df$sdinv, 2))

            pl$x$data[[1]]$text <- paste(query_genes, query_dps, query_sds, sep = '<br>')

            # drug hoverinfo
            nsyms <- nrow(dfs$drug_df)/2
            dr1_ind <- 1:nsyms
            dr2_ind <- (nsyms+1):(nsyms*2)
            drug_genes  <- paste0('Gene: ', as.character(dfs$drug_df$gene[dr1_ind]))
            drug_names1 <- paste0('Drug: ', as.character(dfs$drug_df$drug[dr1_ind]))
            drug_names2 <- paste0('Drug: ', as.character(dfs$drug_df$drug[dr2_ind]))
            drug_dps1   <- paste0('Dprime: ', round(dfs$drug_df$dprime[dr1_ind], 1))
            drug_dps2   <- paste0('Dprime: ', round(dfs$drug_df$dprime[dr2_ind], 1))

            pl$x$data[[min(drug_ids)]]$text <- paste(drug_names1, drug_genes, drug_dps1, sep = '<br>')
            pl$x$data[[max(drug_ids)]]$text <- paste(drug_names2, drug_genes, drug_dps2, sep = '<br>')

            # drug visibility
            vis1_status <- ifelse(vis1, "", 'legendonly')
            vis2_status <- ifelse(vis2, "", 'legendonly')

            pl$x$data[[min(drug_ids)]]$visible <- vis1_status
            pl$x$data[[max(drug_ids)]]$visible <- vis2_status

            # sumry hoverinfo
            sumry_genes <- paste0('Gene: ', as.character(dfs$sumry_df$gene))
            sumry_mus   <- paste0('Mu: ', round(dfs$sumry_df$mus, 1))
            sumry_sds   <- paste0('SD: ', round(dfs$sumry_df$high - dfs$sumry_df$mus, 2))

            pl$x$data[[2]]$text <- paste(sumry_genes, sumry_mus, sumry_sds, sep = '<br>')

            pl$x$layout$legend <- list(orientation = 'h', x=0, y=100)

            pl
        })
        
        # link to KEGG pathway ----
        kegg_link <- shiny::reactive({
            path_num <- names(gs.names[which(gs.names == input$path)])
            kegg_link <- paste0('https://www.genome.jp/kegg-bin/show_pathway?map', path_num)
        })
        
        # Click link out to Wikipedia ----
        shiny::observeEvent(input$kegg, {
            utils::browseURL(kegg_link())
        })
    }
    shinyApp(ui, server)
}


# used by explore_paths
get_dfs <- function(path, es_res, drug_info, drugs, gslist, gs.names, topn = 50) {

    drugs <- unique(drugs)

    # dprime and vardprime columns
    isdp  <- grepl('^dprime', colnames(es_res$filt))
    isvar <- grepl('^vardprime', colnames(es_res$filt))
    ndp   <- sum(isdp)

    # path symbols
    path_num <- names(gs.names)[gs.names == path]
    path_sym <- unique(names(gslist[[path_num]]))
    path_sym <- path_sym[path_sym %in% row.names(es_res$filt)]
    nsym <- length(path_sym)


    # order by decreasing absolute mu
    mus <- es_res$filt[path_sym, 'mu']
    names(mus) <- path_sym
    mus <- mus[order(abs(mus), decreasing = TRUE)]
    path_sym <- names(mus)


    if (!is.null(drugs)) {
        # path symbols also must be in drug_info
        path_sym <- path_sym[path_sym %in% row.names(drug_info)]
        nsym <- length(path_sym)

        drug_df <- reshape::melt(drug_info[path_sym, drugs, drop = FALSE])
        colnames(drug_df) <- c('gene', 'drug', 'dprime')
        drug_df$drug <- factor(drug_df$drug, levels = drugs)
    } else {
        path <- paste0(path, '\n')
    }

    # only keep path symbols
    qes <- es_res$filt[path_sym, ]

    # query data.frame
    query_df <- data.frame('dprime' = c(t(qes[, isdp])),
                           'sdinv'  = 1/sqrt(c(t(qes[, isvar]))),
                           'gene'   = factor(rep(path_sym, each=ndp), levels = path_sym))


    # summary data.frame
    sds <- sqrt(qes$var)
    sumry_df <- data.frame('gene' = factor(path_sym, levels = path_sym),
                           'mus'  = qes$mu,
                           'high' = qes$mu + sds,
                           'low'  = qes$mu - sds)
    
    
    dfs <- list(drug_df  = drug_df[drug_df$gene %in% head(query_df$gene, topn), ],
                query_df = head(query_df, topn),
                sumry_df = head(sumry_df, topn))
    
    return(dfs)
}

#' Make FDR values readable
#' 
#' Used by \code{\link{explore_paths}} to make pathway FDR values readable
#'
#' @param fdrs Numeric vector of false discovery rates
#' @keywords internal
#'
#' @return Character vector of formated false discovery rates.
#' @export
#'
#' @examples
#' 
#' fdrs <- c(0.011, 0.00111, 0.000001111, 0.10000001, 1, 1.00000)
#' format_fdrs(fdrs)
format_fdrs <- function(fdrs) {
    
    # two significant digits
    fdrs <- signif(fdrs, 2)
    fdrs <- format(fdrs, scientific = FALSE)
    
    # dont get too specific when very small
    lt_min <- fdrs < 0.001
    fdrs[lt_min] <- '< 0.001'
    
    # remove trailing zeros
    fdrs <- gsub('\\.?0+$', '', fdrs)
    
    return(fdrs)
}

