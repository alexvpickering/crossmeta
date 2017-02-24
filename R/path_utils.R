
# KEGG database
process_kegg <- function() {

    mapkKGML <- list.files('/home/alex/R/library/SPIA/extdata/keggxml/hsa/', full.names = TRUE)
    gslist <- lapply(mapkKGML, KEGGgraph::parseKGML)

    names(gslist) <- lapply(gslist, function(path) path@pathwayInfo@number)
    gs.names <- sapply(gslist, function(path) path@pathwayInfo@title)

    gslist <- lapply(gslist, function(path) {
        # nodes and node types
        nodes <- path@nodes
        types <- sapply(nodes, function(node) node@type)

        # get genes in pathway
        unlist(lapply(nodes[types == 'gene'], function(node) gsub('^hsa:', '', node@name)), use.names = FALSE)
    })
    return(list(gslist = gslist, gs.names=gs.names))
}


#' Differential expression of KEGG pathways.
#'
#' User selects contrasts, then surrogate variable analysis (sva) and
#' differential expression analysis (limma) is performed.
#'
#' @param esets List of annotated esets. Created by \code{load_raw}.
#' @param prev_anals Previous result of \code{diff_expr}.
#' @param data_dir String specifying directory of GSE folders.
#'
#' @return PADOG result for each GSE.
#' @export
#'
#' @examples
#' dssdfsd <- 4
diff_path <- function(esets, prev_anals, data_dir = getwd()) {

    prev_anals <- prev_anals[names(esets)]
    anals <- list()

    # load KEGG data
    utils::data(list = gslist, package = "crossmeta", envir = environment())
    utils::data(list = gs.names, package = "crossmeta", envir = environment())

    for (i in seq_along(esets)) {

        eset <- esets[[i]]
        gse_name <- names(esets)[i]
        prev_anal <- prev_anals[[i]]

        cat('Working on', paste0(gse_name, ':'), '\n')

        gse_folder <- strsplit(gse_name, "\\.")[[1]][1]  # name can be "GSE.GPL"
        gse_dir <- paste(data_dir, gse_folder, sep = "/")

        # select contrasts
        cons <- crossmeta:::add_contrasts(eset, gse_name, prev_anal)
        contrasts <- cons$contrasts

        # remove replicates/duplicates and annotate with human ENTREZID
        dups <- crossmeta:::iqr_replicates(cons$eset, annot = 'ENTREZID_HS', rm.dup = TRUE)
        eset <- dups$eset

        groups <- pData(eset)$group

        conl <- list()
        # padog for each contrast
        for (con in contrasts) {

            cat('    ', paste0(con, '...', '\n'))

            # contrast levels
            ctrl <- gsub('^.+?-', '', con)
            test <- gsub('-.+?$', '', con)

            # group
            incon <- groups %in% c(ctrl, test)
            group <- ifelse(groups[incon] == ctrl, 'c', 'd')

            # other padog inputs
            esetm  <- exprs(eset[, incon])
            pairs  <- pData(eset[, incon])$pairs
            paired <- !anyNA(pairs)
            block  <- NULL
            if (paired)
                block <- pairs

            # run padog
            conl[[con]] <- PADOG::padog(esetm, group, paired, block, gslist, gs.names = gs.names,
                                        parallel = TRUE, ncr = parallel::detectCores(),
                                        verbose = FALSE)
        }
        anals[[gse_name]] <- conl
    }
    return (anals)
}

