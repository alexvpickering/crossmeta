#' Differential expression of list of expression sets.
#'
#' User selects contrasts, then surrogate variable analysis (sva) and
#' differential expression (limma) is run.
#'
#' @param esets list from load_affy, load_illum, or load_agil.
#' @param data_dir Character, data directory containing GSE folders.
#' @param annot Character, one of "PROBE" or "SYMBOL" for probe or gene level
#'        analysis respectively. For duplicated genes symbols,the feature with
#'        the highest interquartile range across selected samples will be kept.
#' @param prev_anals previous output to re-use previous contrast selections.
#' @export
#' @return List of lists (one per GSE), each containing:
#'         \item{eset}{expression set with selected samples and non-duplicate
#'         features named by either gene symbol or probe. pData slot of eset also
#'         contains treatment, tissue, and group columns.}
#'         \item{top_tables}{list with results of \link{topTable} call (one per
#'         contrast). These results model the effect of \link{sva} discovered
#'         surrogate variables.}
#'         \item{ebayes_sv}{results of call to \link{eBayes} with surrogate
#'         variables included in the model matrix.}
#'         \item{ebayes}{results of call to \link{eBayes} with surrogate
#'         variables NOT included in the model matrix.}
#'         \item{mama_data}{list used by \link{make_ma}.}

diff_expr <- function (esets, data_dir, annot="SYMBOL", prev_anals=list(NULL)) {

  prev_anals <- prev_anals[names(esets)]
  anals <- list()
  for (i in seq_along(esets)) {

    eset <- esets[[i]]
    gse_name <- names(esets)[i]
    prev_anal <- prev_anals[[i]]

    gse_folder <- strsplit(gse_name, "\\.")[[1]][1]  #name can be "GSE.GPL"
    gse_dir <- paste(data_dir, gse_folder, sep="/")

    #select contrasts
    cons <- add_contrasts(eset, gse_name, annot, prev_anal)

    #differential expression
    setup <- diff_setup(cons$eset, cons$levels)

    anal <- diff_anal(cons$eset, cons$contrasts, cons$levels, cons$mama_data,
                      setup$mod, setup$modsv, setup$svobj,
                      gse_dir, gse_name, annot)

    anals[[gse_name]] <- anal
  }
  return (anals)
}


#' Reuse contrast selections from previous analysis.
#'
#' Transfers user-supplied selections from previous call of diff_expr.
#' @importFrom Biobase sampleNames pData
#'
#' @param eset expression set passed to diff_expr.
#' @param prev_anal result of previous call to diff_expr.
#'
#' @seealso \code{\link{diff_expr}}
#' @return expression set with samples and pData as in prev_anal.
#'
match_prev_eset <- function(eset, prev_anal) {

    #retain previously selected samples only
    selected_samples <- sampleNames(prev_anal$eset)
    eset <- eset[, selected_samples]

    #transfer previous treatment, tissue & group to eset
    pData(eset)$treatment <- pData(prev_anal$eset)$treatment
    pData(eset)$tissue    <- pData(prev_anal$eset)$tissue
    pData(eset)$group     <- pData(prev_anal$eset)$group

    return (eset)
}


#------------------------


#' Select contrasts for each GSE.
#'
#' Function is used by \code{diff_expr} to get sample selections for each
#' contrast from user. After selecting samples, duplicated features
#' (probes or genes) are removed.
#'
#' @importFrom Biobase sampleNames pData
#' @importFrom tcltk tk_select.list
#' @importFrom stringr str_extract
#'
#' @param eset expression set with "SYMBOL" column in \code{fData}.
#' @param gse_name Character, GSE name of eset.
#' @param annot Character, one of "PROBE" or "SYMBOL" for probe or gene level
#'        analysis respectively. For duplicated gene symbols,the feature with
#'        the highest interquartile range across selected samples will be kept.
#' @param prev_anal result of previous call to diff_expr
#'
#' @seealso \code{\link{diff_expr}}
#' @return list needed for \code{\link{diff_setup}}
#'

add_contrasts <- function (eset, gse_name, annot="SYMBOL", prev_anal) {

    if (!is.null(prev_anal)) {
        #re-use selections/sample labels from previous analysis
        eset <- match_prev_eset(eset, prev_anal)

        #get contrast info from previous analysis
        mama_samples   <- prev_anal$mama_data$sample_names
        mama_clinicals <- prev_anal$mama_data$clinicals

        contrasts    <- colnames(prev_anal$ebayes$contrasts)
        group_levels <- colnames(prev_anal$ebayes$design)
        selected_samples <- sampleNames(prev_anal$eset)

    } else {
        #get contrast info from user input
        choices <- paste(sampleNames(eset), pData(eset)$title)
        contrasts <- c()
        group_levels <- c()
        selected_samples <- c()

        #for MAMA data
        mama_samples <- list()
        mama_clinicals <- list()

        #repeat until all contrasts selected
        while (TRUE) {
            #select ctrl samples
            ctrl <- tk_select.list(choices, multiple=T,
                                   title="Control samples for contrast")
            ctrl <- str_extract(ctrl, "GSM[0-9]+")
            if (length(ctrl) == 0) {break}

            #select test samples
            test <- tk_select.list(choices, multiple=T,
                                   title="Test samples for contrast")
            test <- str_extract(test, "GSM[0-9]+")

            #add treatment to pheno
            pData(eset)[ctrl, "treatment"] <- "ctrl"
            pData(eset)[test, "treatment"] <- "test"

            #add group names & tissue to pheno
            tissue <- inputs("Tissue source", box1="eg. liver", def1="")
            group_names <- paste(inputs("Group names", two=T), tissue, sep=".")
            pData(eset)[c(ctrl, test), "tissue"] <- tissue
            pData(eset)[ctrl, "group"] <- group_names[1]
            pData(eset)[test, "group"] <- group_names[2]

            #add to contrasts
            contrast <- paste(group_names[2], group_names[1], sep="-")
            contrasts <- c(contrasts, contrast)

            #add to group_levels
            group_levels <- unique(c(group_levels,
                                     group_names[1],
                                     group_names[2]))

            #add to selected_samples
            selected_samples <- unique(c(selected_samples, ctrl, test))

            #store sample names for each contrast
            contrast_name <- paste(gse_name, contrast, sep="_")
            mama_samples[[contrast_name]] <- c(ctrl, test)

            #add labels for contrast to mama_clinicals
            labels <- factor(pData(eset)[c(ctrl, test), "treatment"],
                             levels = c("test", "ctrl"))
            clinical <- data.frame(treatment = labels,
                                   row.names = c(ctrl, test))
            mama_clinicals[[contrast_name]] <- clinical
        }
        #retain selected samples only
        eset <- eset[, selected_samples]
    }
    #remove rows with duplicated/NA annot (SYMBOL or PROBE)
    eset <- iqr_duplicates(eset, annot)

    #subset eset for each contrast (for MAMA)
    mama_esets <- get_contrast_esets(mama_samples, eset)

    #put data together
    mama_data <- list(esets=mama_esets, sample_names=mama_samples,
                      clinicals=mama_clinicals)

    contrast_data <- list(eset=eset, contrasts=contrasts, levels=group_levels,
                          samples=selected_samples, mama_data=mama_data)
    return (contrast_data)
}


#------------------------


#' Removes features with duplicated annotation.
#'
#' For rows with duplicated annot, highested IQR retained.
#'
#' @import dplyr
#' @importFrom Biobase fData sampleNames exprs
#' @importFrom matrixStats rowIQRs
#'
#' @param eset expression set with fData column "SYMBOL" or "PROBE".
#' @param annot Character, one of "PROBE" or "SYMBOL" for probe or gene level
#'        analysis respectively. For duplicated gene symbols,the feature with
#'        the highest interquartile range across selected samples will be kept.
#' @return expression set with unique features at probe or gene level.
#'
iqr_duplicates <- function (eset, annot="SYMBOL") {

    #remove rows with NA annot (occurs if annot is SYMBOL)
    eset <- eset[!is.na(fData(eset)[, annot]),]

    data <- as.data.frame(exprs(eset))
    #add inter-quartile ranges to data
    data$IQR <- rowIQRs(exprs(eset))
    #add feature data
    data[, colnames(fData(eset))] <- fData(eset)

    #for rows with same annot, keep highest IQR
    data %>%
        group_by_(annot) %>%
        arrange(desc(IQR)) %>%
        slice(1) %>%
        ungroup ->
        data

    #seperate exprs and fdata columns
    exprs <- as.matrix(data[, sampleNames(eset)])
    fdata <- as.matrix(data[, colnames(fData(eset))])
    row.names(exprs) <- fdata[, annot]
    row.names(fdata) <- fdata[, annot]

    #put exprs and fdata into eset
    exprs(eset) <- exprs
    fData(eset) <- as.data.frame(fdata)

    return (eset)
}


#------------------------


#' Generate model matrix with surrogate variables.
#'
#' Used by \code{\link{diff_expr}} to get model matrix with surrogate variabled
#' in order to run \code{\link{diff_anal}}.
#'
#' @importFrom Biobase pData exprs
#' @importFrom sva sva
#'
#' @param eset eset with "group" column in phenoData.
#' @param group_levels unique group names created by \code{\link{add_contrasts}}.
#'
#' @seealso \code{\link{add_contrasts}}, \code{\link{diff_expr}}.
#' @return list with model matrix(mod), model matrix with surrogate
#'         variables(modsv), and result of \code{\link{sva}} function.
#'
diff_setup <- function(eset, group_levels){

    #make full model matrix
    group <- factor(pData(eset)$group, levels=group_levels)
    mod <- model.matrix(~0+group)
    colnames(mod) <- group_levels

    #make null model matrix (sva)
    mod0 <- model.matrix(~1, data=pData(eset))

    #surrogate variable analysis
    capture.output(svobj <- sva(exprs(eset), mod, mod0))
    modsv <- cbind(mod, svobj$sv)
    colnames(modsv) <- c(colnames(mod), paste("SV", 1:svobj$n.sv, sep=""))

    return (list("mod"=mod, "modsv"=modsv, "svobj"=svobj))
}


#------------------------


#' Run limma analysis.
#'
#' Runs limma differential expression analysis on all contrasts selected by
#' \code{\link{add_contrasts}}. Analysis performed with and without surrogate
#' variables discovered by \code{\link{diff_setup}}. Also print MDS plot and saves
#' results to disc.
#'
#' @importFrom Biobase pData exprs
#' @importFrom limma topTable plotMDS
#' @importFrom RColorBrewer brewer.pal
#'
#' @param eset eset, generated by \code{\link{add_contrasts}}.
#' @param contrasts Character vector, generated by \code{\link{add_contrasts}}.
#' @param group_levels Character vector, generated by \code{\link{add_contrasts}}.
#' @param mama_data List, generated by \code{\link{add_contrasts}}.
#' @param mod model matrix, generated by \code{\link{diff_setup}}.
#' @param modsv model matrix with surrogate variables,
#'        generated by \code{\link{diff_setup}}.
#' @param svobj result from \code{\link{sva}} function called during
#'        \code{\link{diff_setup}}.
#' @param gse_dir Character, path to directory with GSE data.
#' @param gse_name Character, name of GSE.
#' @param annot one of "SYMBOL" or "PROBE". Different save names used for each
#'        so that can run analysis both at gene or probe level.
#'
#' @seealso \code{\link{diff_expr}}.
#' @return List, final result of \code{diff_expr}. Used for subsequent meta-analysis.

diff_anal <- function(eset, contrasts, group_levels, mama_data,
                      mod, modsv, svobj, gse_dir, gse_name, annot="SYMBOL"){

    #differential expression (surrogate variables modeled and not)
    ebayes_sv <- fit_ebayes(eset, contrasts, modsv)
    ebayes <- fit_ebayes(eset, contrasts, mod)

    #annotate/store results
    top_tables <- list()
    contrast_names <- names(mama_data$esets)

    for (i in seq_along(contrast_names)) {
        top_genes <- topTable(ebayes_sv, coef=i, n=Inf)
        num_sig <- dim(filter(top_genes, adj.P.Val <0.05))[1]
        top_tables[[contrast_names[i]]] <- top_genes
        cat (contrast_names[i], "(n significant):", num_sig, "\n")
    }
    cat("\n")

    #plot MDS
    group <- factor(pData(eset)$group, levels=group_levels)
    palette <- brewer.pal(12, "Paired")
    colours <- palette[group]

    sva_exprs <- clean_y(exprs(eset), mod, svobj$sv)

    plotMDS(sva_exprs, pch=19, main = gse_name, col = colours)
    legend("center", legend=group_levels, fill=unique(colours))

    #add sva adjusted exprs to mama_data
    mama_data$esets_sva <- get_contrast_esets(mama_data$sample_names,
                                              eset, sva_exprs)

    #save to disk
    diff_expr <- list(eset=eset, top_tables=top_tables,
                      ebayes_sv=ebayes_sv, ebayes=ebayes, mama_data=mama_data)

    if (annot == "SYMBOL") save_name <- paste (gse_name,
                                               "diff_expr.rds", sep="_")
    if (annot == "PROBE")  save_name <- paste (gse_name,
                                               "diff_expr_probe.rds", sep="_")

    saveRDS(diff_expr, file = paste(gse_dir, save_name, sep="/"))
    return (diff_expr)
}


#------------------------


#' Perform eBayes analysis from limma.
#'
#' Generates contrast matrix then runs eBayes analysis from limma.
#'
#' @importFrom Biobase exprs
#' @importFrom limma makeContrasts lmFit eBayes contrasts.fit
#'
#' @param eset eset, generated by \code{\link{add_contrasts}}.
#' @param contrasts Character vector, generated by \code{\link{add_contrasts}}.
#' @param mod model matrix (with or without surrogate variables modeled),
#'        generated by \code{\link{diff_setup}}.
#'
#' @return result from called to limma \code{\link{eBayes}}.

fit_ebayes <- function(eset, contrasts, mod) {
    contrast_matrix <- makeContrasts(contrasts=contrasts, levels=mod)
    fit <- contrasts.fit (lmFit(exprs(eset),mod), contrast_matrix)
    return (eBayes(fit))
}


#------------------------


#' Generate eset for each contrast.
#'
#' Sets up esets for \link[MAMA]{MetaArray-class}. An eset is generated using
#' only the supplied sample names. Used to generate eset for each contrast
#' within a GSE.
#'
#' @importFrom Biobase exprs
#'
#' @param sample_names Character vector, sample names to include in eset.
#' @param eset eset to subset.
#' @param data expression data for eset (sva-adjusted or not).
#'
#' @seealso  \code{\link{clean_y}}.
#' @return eset including only the supplied sample names.

get_contrast_esets <- function(sample_names, eset, data=exprs(eset)) {
    #used to setup esets for MAMA
    #data = exprs_sva will return sva-adjusted eset for each contrast
    exprs(eset) <- data
    lapply(sample_names, function(x) eset[,x])
}


#------------------------


#' Adjustes expression data for surrogate variables.
#'
#' Factors out effect of surrogate variables discovered during analysis by
#' \link{sva}.
#'
#' @param y expression data of eset.
#' @param mod full model matrix supplied to \link{sva}.
#' @param svs surrogate variables returned by \link{sva} (svobj$sv).
#'
#' @seealso \code{\link{get_contrast_esets}}.
#' @return expression data with effect of svs removed.

clean_y <- function(y, mod, svs) {

    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    rm(Hat)
    gc()
    P = ncol(mod)
    return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}
