#' Explore pathway meta analysis.
#'
#' @import ggplot2 plotly shiny miniUI
#'
#' @param es Result of call to \code{\link{es_meta}}.
#' @param path_res Result of call to \code{\link{path_meta}}.
#' @param drug_info Character vector specifying which dataset to query
#'   (either 'cmap' or 'l1000'). Can also provide a matrix of differential expression
#'   values for drugs or drug combinations (rows are genes, columns are drugs).
#' @param type Desired direction of drug action on query signature. If \code{'both'},
#'   drugs that mimic (positive cosines) and reverse (negative cosines) the query signature
#'   are at the top and bottom of drug selection boxes respectively.
#'
#' @return NULL
#' @export
#'
#' @examples
#' blah <- c()
explore_paths <- function(es, path_res, drug_info, type = c('both', 'mimic', 'reverse')) {

    # bindings to pass check
    gslist = gs.names = NULL
    utils::data("gslist", "gs.names", package = "crossmeta", envir = environment())


    # visibility
    start <- TRUE
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
            top_full <- head(all_full, 50)
            top_path <- head(all_path, 50)
        } else if (type == 'reverse') {
            top_full <- head(rev(all_full), 50)
            top_path <- head(rev(all_path), 50)
        } else {
            top_full <- c(head(all_full, 25), tail(all_full, 25))
            top_path <- c(head(all_path, 25), tail(all_path, 25))
        }

        # construct names as Drug (full cos, path cos)
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
    dprimes  <- ccmap::get_dprimes(es)

    sources      <- names(es)
    paths        <- row.names(path_res[[1]])
    names(paths) <- paste0(paths, ' (', format(signif(path_res[[1]][, 'fdr'], 2), scientific = TRUE), ')')
    top_drugs    <- fullpath(dprimes[[1]]$meta, drug_info, paths[1], type[1])


    ui <- shinyUI(fluidPage(

        title = "Pathway Explorer",

        plotlyOutput('trendPlot', height = "600px"),

        hr(),

        fluidRow(
            column(6, align = 'center',
                   strong("Source:"),
                   selectInput('source', '', choices = sources, selected = sources[1], width = '75%'),
                   tags$style(type='text/css', ".selectize-dropdown-content {max-height: 180px; }"),
                   hr(),
                   strong("Pathway (FDR):"),
                   selectInput('path', '', choices = paths, selected = paths[1], width = '75%')
            ),
            column(6, align = 'center',
                   HTML("<strong>Top Drugs for Full Signature (cos &theta;<sub>full</sub>, cos &theta;<sub>path</sub>):</strong>"),
                   selectInput('drug1', '', top_drugs$full, top_drugs$full[1], width = '75%'),
                   hr(),
                   HTML("<strong>Top Drugs for Pathway Genes Only (cos &theta;<sub>full</sub>, cos &theta;<sub>path</sub>):</strong>"),
                   selectInput('drug2', '', top_drugs$path, top_drugs$path[1], width = '75%')
            )
        )
    ))


    server <- function(input, output, session) {


        observeEvent(input$source, {

            src <- input$source

            # update global top drugs
            paths        <<- row.names(path_res[[src]])
            names(paths) <<- paste0(paths, ' (', format(signif(path_res[[src]][, 'fdr'], 2), scientific = TRUE), ')')
            top_drugs    <<- fullpath(dprimes[[src]]$meta, drug_info, paths[1], type[1])

            # update selections
            updateSelectInput(session, 'path',  choices = paths, selected = paths[1])

        })

        observeEvent(input$path, {

            src <- input$source
            start <<- TRUE

            # update global top drugs
            top_drugs <<- fullpath(dprimes[[src]]$meta, drug_info, input$path, type[1], top_drugs$all_full)

            # update selections
            updateSelectInput(session, 'drug1', choices = top_drugs$full, selected = top_drugs$full[1])
            updateSelectInput(session, 'drug2', choices = top_drugs$path, selected = top_drugs$path[1])

            vis1 <<- FALSE
            vis2 <<- FALSE
        })

        observeEvent(input$drug1, {

            if (!start) {
                vis1 <<- TRUE
                vis2 <<- FALSE

                if (input$drug1 == input$drug2)
                    vis2 <<- TRUE
            }
        })

        observeEvent(input$drug2, {

            if (start) {
                start <<- FALSE

            } else {
                vis1 <<- FALSE
                vis2 <<- TRUE

                if (input$drug1 == input$drug2)
                    vis1 <<- TRUE
            }
        })




        #add reactive data information. Dataset = built in diamonds data
        dataset <- reactive({
            src <- input$source
            get_dfs(input$path, es[[src]], drug_info, c(input$drug1, input$drug2), gslist, gs.names)
        })


        output$trendPlot <- renderPlotly({

            dfs <- dataset()

            drugs <- levels(dfs$drug_df$drug)
            drug_ids <- seq_along(unique(drugs))+1
            sumry_id <- max(drug_ids)+1


            g <-  ggplot(data = dfs$sumry_df, aes_string(x = 'gene', y = 'mus')) +
                geom_point(aes_string('gene', 'dprime', alpha = 'sdinv'), dfs$query_df,
                           shape=1, colour = '#666666', na.rm=TRUE) +
                geom_point(aes_string(x = 'gene',  y = 'dprime', colour = 'drug'), dfs$drug_df,
                           na.rm=TRUE) +
                geom_errorbar(aes_string(ymin = 'low', ymax = 'high'),
                              colour = 'black', width = 0.7) +
                ylab("Dprime") +
                xlab("") +
                geom_hline(yintercept = 0, colour = '#999999') +
                scale_color_manual(values = c("#E41A1C", "#377EB8")) +
                scale_y_continuous(breaks = function(ylims) floor(ylims[1]):ceiling(ylims[2])) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, vjust=0.5),
                      legend.title = element_blank())


            pl <- plotly_build(g)


            # adjust margins and ranges
            ymin <- max(-8, min(dfs$query_df$dprime, na.rm = TRUE))
            ymax <- min(8, max(dfs$query_df$dprime, na.rm = TRUE))

            pl$x$layout$margin$b <- 90
            nsyms <- nrow(dfs$sumry_df)
            pl$x$layout$xaxis$range <- c(0.5, min(30, nsyms+0.5))
            pl$x$layout$yaxis$range <- c(ymin, ymax)

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

            pl$x$data[[sumry_id]]$text <- paste(sumry_genes, sumry_mus, sumry_sds, sep = '<br>')


            pl %>% layout(legend = list(orientation = 'h', x=0, y=100))
        })


    }

    shinyApp(ui, server)

}


# used by explore_paths
get_dfs <- function(path, es, drug_info, drugs, gslist, gs.names) {

    drugs <- unique(drugs)

    # dprime and vardprime columns
    isdp  <- grepl('^dprime', colnames(es$filt))
    isvar <- grepl('^vardprime', colnames(es$filt))
    ndp   <- sum(isdp)

    # path symbols
    path_num <- names(gs.names)[gs.names == path]
    path_sym <- unique(names(gslist[[path_num]]))
    path_sym <- path_sym[path_sym %in% row.names(es$filt)]
    nsym <- length(path_sym)


    # order by decreasing absolute mu
    mus <- es$filt[path_sym, 'mu']
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
    qes <- es$filt[path_sym, ]

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

    return(list(drug_df  = drug_df,
                query_df = query_df,
                sumry_df = sumry_df))
}

