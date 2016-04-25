#hack to fool R CMD CHK
#issue caused by dplyr in diff_anal
globalVariables(c("x", "y", "adj.P.Val"))


#' Differential expression of esets.
#'
#' User selects contrasts, then surrogate variable analysis (sva) and
#' differential expression analysis (limma) is performed.
#'
#' For each GSE, analysis results are saved in the corresponding GSE
#' folder (in \code{data_dir}) that was created by \code{get_raw}. If analysis
#' needs to be repeated, previous results can be reloaded with \code{load_diff}
#' and supplied to the \code{prev_anals} parameter. In this case, previous
#' selections/names will be reused.
#'
#' @import magrittr Biobase
#' @importFrom BiocGenerics annotation
#'
#' @param esets List of annotated esets. Created by \code{load_raw}.
#' @param data_dir String specifying directory of GSE folders.
#' @param annot String, either "PROBE" or "SYMBOL" for probe or gene level
#'   analysis respectively. For duplicated genes symbols,the feature with
#'   the highest interquartile range across selected samples will be kept.
#' @param prev_anals Previous result of \code{diff_expr}. If Present, previous
#'   selections and names will be reused.
#'
#' @export
#' @seealso \code{\link{get_raw}}, \code{\link{load_raw}}, and
#'   \code{\link{load_diff}}.
#'
#' @return List of lists (one per GSE), each containing:
#'   \item{eset}{Expression set with selected samples and non-duplicate
#'      features named by either gene symbol or probe. Treatment (ctl or test),
#'      tissue, and group columns have been added to the \code{pData} slot.}
#'   \item{top_tables}{List with results of \code{\link{topTable}} call (one per
#'      contrast). These results account for the effects of nuissance variables
#'      discovered by surrogate variable analysis.}
#'   \item{ebayes_sv, ebayes}{Results of call to \code{\link{eBayes}}, with
#'      surrogate variables included and not included in the model matrix.}
#'   \item{mama_data}{List used by \code{\link{make_ma}} to generate MetaArray
#'      object.}
#'
#' @examples
#' library(lydata)
#'
#' #location of raw data
#' data_dir <- system.file("extdata", package = "lydata")
#'
#' #gather GSE names
#' affy_names  <- c("GSE9601", "GSE15069")
#' illum_names <- c("GSE50841", "GSE34817", "GSE29689")
#' gse_names   <- c(affy_names, illum_names)
#'
#' #load and commonize esets
#' esets <- load_raw(gse_names, data_dir)
#' com_esets <- commonize(esets)
#'
#' #run analysis
#' anals <- diff_expr(com_esets, data_dir)
#'
#' #re-run analysis
#' prev <- load_diff(gse_names, data_dir)
#' anals <- diff_expr(com_esets, data_dir, prev_anals=prev)


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


# Reuse contrast selections from previous analysis.
#
# Transfers user-supplied selections from previous call of diff_expr.
#
# @param eset Annotated eset. Created by \code{load_raw}.
# @param prev_anal Previous result of \code{diff_expr}.
#
# @seealso \code{\link{diff_expr}}
# @return Expression set with samples and pData as in prev_anal.

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


# Select contrasts for each GSE.
#
# Function is used by \code{diff_expr} to get sample selections for each
# contrast from user. After selecting samples, duplicated features
# (probes or genes) are removed.
#
# @param eset Annotated eset. Created by \code{load_raw}.
# @param gse_name String specifying GSE name for eset.
# @param annot String, either "PROBE" or "SYMBOL" for probe or gene level
#   analysis respectively. For duplicated genes symbols,the feature with
#   the highest interquartile range across selected samples will be kept.
# @param prev_anal Previous result of \code{diff_expr}. If present, previous
#   selections and names will be reused.
#
# @seealso \code{\link{diff_expr}}
# @return list needed for \code{diff_setup}.

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
            ctrl <- tcltk::tk_select.list(choices, multiple=TRUE,
                                          title="Control samples for contrast")
            ctrl <- stringr::str_extract(ctrl, "GSM[0-9]+")
            if (length(ctrl) == 0) {break}

            #select test samples
            test <- tcltk::tk_select.list(choices, multiple=TRUE,
                                          title="Test samples for contrast")
            test <- stringr::str_extract(test, "GSM[0-9]+")

            #add treatment to pheno
            pData(eset)[ctrl, "treatment"] <- "ctrl"
            pData(eset)[test, "treatment"] <- "test"

            #add group names & tissue to pheno
            tissue <- inputs("Tissue source", box1="eg. liver", def1="")
            group_names <- paste(inputs("Group names", two=TRUE),
                                 tissue, sep=".")
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


# Removes features with duplicated annotation.
#
# For rows with duplicated annot, highested IQR retained.
#
# @param eset Annotated eset. Created by \code{load_raw}.
# @param annot String, either "PROBE" or "SYMBOL" for probe or gene level
#   analysis respectively. For duplicated genes symbols,the feature with
#   the highest interquartile range across selected samples will be kept.
#
# @return Expression set with unique features at probe or gene level.

iqr_duplicates <- function (eset, annot="SYMBOL") {

    #remove rows with NA annot (occurs if annot is SYMBOL)
    eset <- eset[!is.na(fData(eset)[, annot]),]

    data <- as.data.frame(exprs(eset))
    #add inter-quartile ranges to data
    data$IQR <- matrixStats::rowIQRs(exprs(eset))
    #add feature data
    data[, colnames(fData(eset))] <- fData(eset)

    #for rows with same annot, keep highest IQR
    data %>%
        dplyr::group_by_(annot) %>%
        dplyr::arrange(dplyr::desc(IQR)) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() ->
        data

    #seperate exprs and fData columns
    exprs <- as.matrix(data[, sampleNames(eset)])
    fData <- as.matrix(data[, colnames(fData(eset))])
    row.names(exprs) <- fData[, annot]
    row.names(fData) <- fData[, annot]

    #put exprs and fData into eset
    exprs(eset) <- exprs
    fData(eset) <- as.data.frame(fData)

    return (eset)
}


#------------------------


# Generate model matrix with surrogate variables.
#
# Used by \code{diff_expr} to create model matrix with surrogate variables
# in order to run \code{diff_anal}.
#
# @param eset Annotated eset. Created by \code{load_raw}.
# @param group_levels Unique group names created by \code{add_contrasts}.
#
# @seealso \code{\link{add_contrasts}}, \code{\link{diff_expr}}.
# @return List with model matrix(mod), model matrix with surrogate
#         variables(modsv), and result of \code{sva} function.

diff_setup <- function(eset, group_levels){

    #make full model matrix
    group <- factor(pData(eset)$group, levels=group_levels)
    mod <- model.matrix(~0+group)
    colnames(mod) <- group_levels

    #make null model matrix (sva)
    mod0 <- model.matrix(~1, data=pData(eset))

    #surrogate variable analysis
    capture.output(svobj <- sva::sva(exprs(eset), mod, mod0))
    modsv <- cbind(mod, svobj$sv)
    colnames(modsv) <- c(colnames(mod), paste("SV", 1:svobj$n.sv, sep=""))

    return (list("mod"=mod, "modsv"=modsv, "svobj"=svobj))
}


#------------------------


# Run limma analysis.
#
# Runs limma differential expression analysis on all contrasts selected by
# \code{add_contrasts}. Analysis performed with and without surrogate
# variables discovered by \code{diff_setup}. Also prints MDS plot and saves
# results.
#
# @param eset Annotated eset created by \code{load_raw}. Duplicate features and
#   non-selected samples removed by \code{add_contrasts}.
# @param contrasts Character vector generated by \code{add_contrasts}.
# @param group_levels Character vector generated by \code{add_contrasts}.
# @param mama_data List generated by \code{add_contrasts}.
# @param mod, modsv Model matrix generated by \code{diff_setup}. With
#   and without surrogate variables.
# @param svobj Result from \code{sva} function called during \code{diff_setup}.
# @param gse_dir String, path to directory with GSE folders.
# @param gse_name String, name of GSE.
# @param annot String, either "PROBE" or "SYMBOL" for probe or gene level
#   analysis respectively. If "PROBE", appends "_probe.rds" to save name.
#
# @seealso \code{\link{diff_expr}}.
# @return List, final result of \code{diff_expr}. Used for subsequent
#   meta-analysis.

diff_anal <- function(eset, contrasts, group_levels, mama_data,
                      mod, modsv, svobj, gse_dir, gse_name, annot="SYMBOL"){

    #differential expression (surrogate variables modeled and not)
    ebayes_sv <- fit_ebayes(eset, contrasts, modsv)
    ebayes <- fit_ebayes(eset, contrasts, mod)

    #annotate/store results
    top_tables <- list()
    contrast_names <- names(mama_data$esets)

    for (i in seq_along(contrast_names)) {
        top_genes <- limma::topTable(ebayes_sv, coef=i, n=Inf)
        num_sig <- dim(dplyr::filter(top_genes, adj.P.Val < 0.05))[1]
        top_tables[[contrast_names[i]]] <- top_genes
        cat (contrast_names[i], "(n significant):", num_sig, "\n")
    }
    cat("\n")

    #plot MDS
    group <- factor(pData(eset)$group, levels=group_levels)
    palette <- RColorBrewer::brewer.pal(12, "Paired")
    colours <- palette[group]

    sva_exprs <- clean_y(exprs(eset), mod, svobj$sv)

    limma::plotMDS(sva_exprs, pch=19, main = gse_name, col = colours)
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


# Perform eBayes analysis from limma.
#
# Generates contrast matrix then runs eBayes analysis from limma.
#
# @param eset Annotated eset created by \code{load_raw}. Duplicate features and
#   non-selected samples removed by \code{add_contrasts}.
# @param contrasts Character vector generated by \code{add_contrasts}.
# @param mod Model matrix generated by \code{diff_setup}. With
#   or without surrogate variables.
#
# @return result from call to limma \code{eBayes}.

fit_ebayes <- function(eset, contrasts, mod) {
    contrast_matrix <- limma::makeContrasts(contrasts=contrasts, levels=mod)
    fit <- limma::contrasts.fit(limma::lmFit(exprs(eset), mod), contrast_matrix)
    return (limma::eBayes(fit))
}


#------------------------


# Create eset for each contrast.
#
# Expression sets are generated using sample names for each contrast within a
# GSE. These esets are stored in the \code{mama_data} slot created during
# \code{diff_expr}.
#
# @param sample_names List of character vectors specifying the sample names
#   for each contrast.
# @param eset Annotated eset created by \code{load_raw}. Duplicate features and
#   non-selected samples removed by \code{add_contrasts}.
# @param data Expression data for eset. Used to supply sva adjusted data created
#   by \code{clean_y}.
#
# @seealso  \code{\link{clean_y}}.
# @return List of expression sets (one per contrast).

get_contrast_esets <- function(sample_names, eset, data=exprs(eset)) {
    #used to setup esets for MAMA
    #data = exprs_sva will return sva-adjusted eset for each contrast
    exprs(eset) <- data
    lapply(sample_names, function(x) eset[,x])
}


#------------------------


# Adjustes expression data for surrogate variables.
#
# Factors out effect of surrogate variables discovered during surrogate variable
# analysis.
#
# @param y Expression data of eset.
# @param mod Full model matrix supplied to \code{sva}.
# @param svs Surrogate variables returned by \code{sva} (svobj$sv).
#
# @seealso \code{\link{get_contrast_esets}}.
# @return Expression data with effects of svs removed.

clean_y <- function(y, mod, svs) {

    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    rm(Hat)
    gc()
    P = ncol(mod)
    return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}
