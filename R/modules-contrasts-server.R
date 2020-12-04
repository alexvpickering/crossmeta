#' Logic Bulk Data page
#' @export
#' @keywords internal
bulkPage <- function(input, output, session, eset, gse_name, prev) {
  
  
  
  bulkForm <- callModule(bulkForm, 'form',
                         pdata = bulkTable$pdata,
                         prev = prev)
  
  
  dataset_name <- reactive(gse_name)
  
  up_annot <- callModule(bulkAnnot, 'anal',
                         dataset_name = dataset_name,
                         pdata = bulkTable$pdata)
  
  
  bulkTable <- callModule(bulkTable, 'explore',
                          eset = eset,
                          prev = prev,
                          up_annot = up_annot)
  
  return(list(
    pdata = bulkTable$pdata,
    contrasts = bulkForm$contrasts
  ))
  
}


#' Logic for Bulk Data form
#' @export
#' @keywords internal
bulkForm <- function(input, output, session,  pdata, prev) {
  
  contrasts <- reactiveVal(colnames(prev$ebayes_sv$contrasts))
  
  group_levels <- reactive({
    get_group_levels(pdata())
  })
  
  callModule(addContrast, 'add_contrast', group_levels, contrasts)
  callModule(delContrasts, 'del_contrasts', group_levels, contrasts)
  
  
  
  return(list(
    contrasts = contrasts
  ))
  
}

addContrast <- function(input, output, session, group_levels, contrasts) {
  
  contrast_options <- list(render = I('{option: addContrastOptions, item: addContrastItem}'))
  reset_sel <- reactiveVal(FALSE)
  
  group_colors <- reactive(get_palette(group_levels()))
  
  group_choices <- reactive({
    
    data.frame(
      name = group_levels(),
      value = group_levels(),
      color = group_colors(), stringsAsFactors = FALSE
    )
  })
  
  observe({
    reset_sel()
    updateSelectizeInput(session,
                         'select_groups', 
                         choices = group_choices(),
                         selected = NULL,
                         server = TRUE, 
                         options = contrast_options)
  })
  
  full_contrast <- reactive(length(input$select_groups) == 2)
  
  observeEvent(input$add_contrast, {
    req(full_contrast())
    
    # add contrast to previous contrasts
    prev <- contrasts()
    sel <- input$select_groups
    new <- paste0(sel[1], '-', sel[2])
    
    req(!new %in% prev)
    contrasts(c(prev, new))
    
    reset_sel(!reset_sel())
  })
}

delContrasts <- function(input, output, session, group_levels, contrasts) {
  contrast_options <- list(render = I('{option: delContrastOptions, item: delContrastItem}'))
  
  contrast_choices <- reactive({
    contrasts <- contrasts()
    group_levels <- group_levels()
    req(contrasts)
    get_contrast_choices(contrasts, group_levels)
  })
  
  observe({
    updateSelectizeInput(session,
                         'select_contrasts', 
                         choices = rbind(NA, contrast_choices()),
                         server = TRUE, 
                         options = contrast_options)
  })
  
  observeEvent(input$del_contrasts, {
    prev <- contrasts()
    del <- input$select_contrasts
    contrasts(setdiff(prev, del))
  })
  
  
  
}

get_contrast_choices <- function(contrasts, group_levels) {
  group_colors <- get_palette(group_levels)
  names(group_colors) <- group_levels
  
  cons <- strsplit(contrasts, '-')
  test <- sapply(cons, `[[`, 1)
  ctrl <- sapply(cons, `[[`, 2)
  
  data.frame(test,
             ctrl,
             testColor = group_colors[test],
             ctrlColor = group_colors[ctrl],
             value = contrasts,
             stringsAsFactors = FALSE)
}

#' Get group levels for bulk data plots
#'
#' @param pdata Data.frame of phenotype data
#' @export
#'
#' @keywords internal
get_group_levels <- function(pdata) {
  pdata <- pdata[!is.na(pdata$Group), ]
  group <- pdata$`Group name`
  group_order <- order(unique(pdata$Group))
  unique(group)[group_order]
}


#' Logic for pdata table
#' @export
#' @keywords internal
bulkTable <- function(input, output, session, eset, prev, up_annot) {
  
  pdata <- reactive({
    pdata <- init_pdata(eset, prev)
    annot <- up_annot()
    if (is.null(annot)) return(pdata)
    else return(annot)
  })
  
  html_pdata <- reactive({
    # things that trigger update
    pdata <- pdata()
    
    group <- pdata$Group
    name  <- pdata$`Group name`
    
    # update pdata Group column
    not.na <- !is.na(group)
    group_nums  <- unique(group[not.na])
    group_names <- unique(name[not.na])
    ind <- order(group_nums)
    
    group_nums  <- group_nums[ind]
    group_names <- group_names[ind]
    
    group_colors <- get_palette(group_nums)
    for (i in seq_along(group_nums)) {
      group_num <- group_nums[i]
      group_name <- group_names[i]
      
      color <- group_colors[group_num]
      
      rows <- which(group == group_num)
      pdata[rows, 'Group'] <- paste('<div style="background-color:', color, ';"></div>')
      pdata[rows, 'Group name'] <- group_name
    }
    return(pdata)
  })
  
  # redraw table when eset changes (otherwise update data using proxy)
  output$pdata <- DT::renderDataTable({
    DT::datatable(
      html_pdata(),
      class = 'cell-border dt-fake-height',
      rownames = FALSE,
      escape = FALSE, # to allow HTML in table
      options = list(
        columnDefs = list(list(className = 'dt-nopad', targets = 0)),
        scrollY = FALSE,
        scrollX = TRUE,
        paging = FALSE,
        bInfo = 0
      )
    )
  })
  
  return(list(
    pdata = pdata
  ))
  
}

#' Logic for downloading and uploading bulk annotation
#' @export
#' @keywords internal
bulkAnnot <- function(input, output, session, dataset_name, pdata) {
  
  observeEvent(input$click_up, {
    shinyjs::click('up_annot')
    error_msg(NULL)
  })
  
  observeEvent(input$click_dl, {
    shinyjs::click('dl_annot')
    error_msg(NULL)
  })
  
  fname <- reactive(paste0(dataset_name(), '_annot.csv'))
  
  output$dl_annot <- downloadHandler(
    filename = fname,
    content = function(con) {
      
      write.csv(format_dl_annot(pdata()), con, row.names = FALSE)
    }
  )
  
  # uploaded annotation
  up_annot <- reactiveVal()
  error_msg <- reactiveVal()
  
  
  observe({
    msg <- error_msg()
    html('error_msg', html = msg)
    toggleClass('validate-up', 'has-error', condition = isTruthy(msg))
  })
  
  observeEvent(input$up_annot, {
    
    ref <- pdata()
    
    infile <- input$up_annot
    if (!isTruthy(infile)){
      res <- msg <- NULL
      
    } else {
      res <- read.csv(infile$datapath, check.names = FALSE, stringsAsFactors = FALSE)
      msg <- validate_up_annot(res, ref)
      
      res <- if (is.null(msg)) format_up_annot(res, ref) else NULL
    }
    
    error_msg(msg)
    up_annot(res)
  })
  
  return(up_annot)
  
}

#' Format downloaded annotation
#' @export
#' @keywords internal
format_dl_annot <- function(annot) {
  
  add_pair <- function(df) {
    pair <- df$Pair
    if (is.null(pair)) pair <- df$pair
    if (is.null(pair)) pair <- NA
    
    df$pair <- df$Pair <- NULL
    tibble::add_column(df, Pair = pair, .after = 'Group name')
  }
  
  annot <- add_pair(annot)
  annot <- annot[, colnames(annot) != 'Group']
  return(annot)
  
}


#' Validate uploaded bulk annotation
#' @export
#' @keywords internal
validate_up_annot <- function(up, ref) {
  msg <- NULL
  
  req_cols <- c('Group name', 'Pair', 'Accession')
  miss_cols <- req_cols[!req_cols %in% colnames(up)]
  
  group <- up$`Group name`
  group <- group[!is.na(group)]
  ngroups <- length(unique(group))
  
  
  if (length(miss_cols)) {
    msg <- paste('Missing columns:', paste(miss_cols, collapse = ', '))
    
  } else if (!all(up$Accession %in% ref$Accession)) {
    msg <- 'Do not change Accession column'
    
  } else if (ngroups < 2) {
    msg <- 'Need at least 2 groups'
    
  } else if (length(group) < 3) {
    msg <- 'Need at least 3 grouped samples'
    
  } else if (!is_invertible(up)) {
    msg <- 'Group name and Pair combination not solvable'
  }
  
  return(msg)
}

#' Check uploaded bulk pdata to make sure the study design is invertible
#' @export
#' @keywords internal
is_invertible <- function(pdata) {
  pdata <- pdata[!is.na(pdata$`Group name`), ]
  pdata$group <- pdata$`Group name`
  
  mod <- get_sva_mods(pdata)$mod
  is(try(solve.default(t(mod) %*% mod),silent=T), 'matrix')
}

init_pdata <- function(eset, prev) {
  
  pcols <- colnames(eset@phenoData)
  pdata <- data.frame('Group' = NA,
                      'Group name' = NA,
                      'Pair' = NA,
                      'Accession' = sampleNames(eset), 
                      row.names = sampleNames(eset), check.names = FALSE)
  
  
  # title not helpful if two-channel
  ch2 <- any(grepl('_red|_green', colnames(eset)))
  
  if (ch2) {
    pdata$Source <- eset$source_name
    pdata$Label <- eset$label
    used_cols <- c('source_name', 'label', 'geo_accession')
    
  } else {
    pdata$Title <- pData(eset)$title
    used_cols <- c('title', 'geo_accession')
  }
  
  # remove accession numbers if Illumina (match with raw data not guaranteed)
  if ('illum' %in% pcols) {
    warning(
      "Unmatched Illumina raw samples and GEO annotation. ",
      "Use 'Title' (not 'Accession') to decide groups.", call. = FALSE)
    
  } else {
    add_cols <- setdiff(pcols, used_cols)
    pdata <- cbind(pdata, pData(eset)[, add_cols])
  }
  
  if (!is.null(prev)) {
    pdata <- format_prev_pdata(prev$pdata, pdata)
  }
  
  return(pdata)
  
}

format_prev_pdata <- function(prev, pdata) {
  matches <- match(row.names(prev), row.names(pdata))
  group <- prev$group
  levels <- unique(group[!is.na(group)])
  
  pdata[matches, 'Group name'] <- group
  pdata[matches, 'Group'] <- as.numeric(factor(group, levels = levels))
  
  # workaround as previously saved 'pairs' column
  pair_col <- grep('^pairs?$', colnames(prev))[1]
  if (!is.na(pair_col)) 
    pdata[matches, 'Pair'] <- prev[[pair_col]]
  
  return(pdata)
}

#' Format uploaded annotation
#' @export
#' @keywords internal
format_up_annot <- function(up, ref) {
  row.names(up) <- up$Accession
  up[up == ''] <- NA
  
  # Group in order of Group name
  # allows changing color of groups by changing order of samples
  group <- up$`Group name`
  levels <- unique(group[!is.na(group)])
  group <- as.numeric(factor(group, levels = levels))
  up <- tibble::add_column(up, Group = group, .before = 1)
  
  up$Pair <- factor(up$Pair)
  
  # in case order of sample was changed
  up <- up[row.names(ref), ]
  
  return(up)
  
}


#' Get a pallete for cluster plots
#'
#' @param levs Character vector of levels to get colour pallete for.
#'
#' @return Character vector with colour codes of \code{length(levs)}.
#' @export
#' @keywords internal
get_palette <- function(levs, dark = FALSE) {
  
  # palettes from scater
  tableau20 <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
                 "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
                 "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
                 "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")
  
  tableau20_dark <- c("#004771", "#4075A5", "#984802", "#A06705", "#036003",
                      "#33870E", "#821919", "#CF1701", "#5B3979", "#7D5E91",
                      "#55342D", "#7B574F", "#A22283", "#CE308A", "#4B4B4B",
                      "#727272", "#6D6D01", "#7F7F07", "#056F79", "#028491")
  
  tableau10medium <- c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                       "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                       "#CDCC5D", "#6DCCDA")
  
  tableau10medium_dark <- c("#365D83", "#9C5800", "#33702A", "#A22616",
                            "#6D4B86", "#644741", "#B62B8A", "#5E5E5E",
                            "#777600", "#097984")
  
  
  nlevs <- length(levs)
  if (nlevs == 2) {
    values <- c("#729ECE", "#FF9E4A")
    
  } else if (nlevs <= 10) {
    pal <- if(dark) tableau10medium_dark else tableau10medium
    values <- head(pal, nlevs)
    
  } else if (nlevs <= 20) {
    pal <- if(dark) tableau20_dark else tableau20
    values <- head(pal, nlevs)
    
  } else {
    set.seed(0)
    values <- randomcoloR::distinctColorPalette(nlevs)
    
  }
  return(values)
}

