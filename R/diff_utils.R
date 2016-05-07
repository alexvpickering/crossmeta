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
#' @import magrittr Biobase shiny miniUI
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
        cons <- add_contrasts(eset, gse_name, prev_anal)

        #setup for differential expression
        setup <- diff_setup(cons$eset, cons$levels)

        #remove rows with duplicated/NA annot (SYMBOL or PROBE)
        dups <- iqr_duplicates(cons$eset, setup$mod, setup$svobj, annot)

        #add esets (sva adjusted and not) to mama_data
        cons$mama_data <- get_contrast_esets(cons$mama_data, dups)

        #differential expression
        anal <- diff_anal(dups$eset, dups$exprs_sva,
                          cons$contrasts, cons$levels, cons$mama_data,
                          setup$mod, setup$modsv, setup$svobj,
                          gse_dir, gse_name, annot)

        anals[[gse_name]] <- anal
    }
    return (anals)
}



#---------------------


# Reuse contrast selections from previous analysis.
#
# Transfers user-supplied selections from previous call of diff_expr.
#
# @param eset Annotated eset. Created by \code{load_raw}.
# @param prev_anal One item (for eset) from previous result of \code{diff_expr}.
#    If present, previous selections and names will be reused.
#
# @seealso \code{\link{diff_expr}}
# @return Expression set with samples and pData as in prev_anal.

match_prev_eset <- function(eset, prev_anal) {

    #retain previously selected samples only
    selected_samples <- sampleNames(prev_anal$eset)
    eset <- eset[, selected_samples]

    #transfer previous treatment, tissue & group to eset
    pData(eset)$treatment <- pData(prev_anal$eset)$treatment
    pData(eset)$group     <- pData(prev_anal$eset)$group

    return (eset)
}


#------------------------


# Select contrasts for each GSE.
#
# Function is used by \code{diff_expr} to get sample selections for each
# contrast from user.
#
# @inheritParams match_prev_eset
# @param gse_name String specifying GSE name for eset.
#
# @seealso \code{\link{diff_expr}}
# @return List with
#    \item{eset}{Expression set with selected samples only.}
#    \item{contrasts}{Character vector of contrast names.}
#    \item{levels}{Character vector of group variable names.
#    \item{mama_data}{List with :
#       \item{sample_names}{Named list  of character vectors (one per contrast)
#          with sample names.}
#       \item{clinicals}{Named list of dataframes (one per contrast) with rows
#          equal to sample names and 'treatment' column a factor with levels
#          equal to 'test' and 'ctrl'.

add_contrasts <- function (eset, gse_name, prev_anal) {

    if (!is.null(prev_anal)) {
        #re-use selections/sample labels from previous analysis
        eset <- match_prev_eset(eset, prev_anal)

        #get contrast info from previous analysis
        mama_samples   <- prev_anal$mama_data$sample_names
        mama_clinicals <- prev_anal$mama_data$clinicals

        contrasts    <- colnames(prev_anal$ebayes$contrasts)
        group_levels <- colnames(prev_anal$ebayes$design)

    } else {
        #get contrast info from user input
        pdata <- data.frame(Accession=sampleNames(eset),
                            Title=pData(eset)$title)
        sels <- select_contrasts(gse_name, pdata)

        #for MAMA data
        mama_samples <- list()
        mama_clinicals <- list()

        #setup data for each contrast
        for (i in 1:nrow(sels$cons)) {
            #get group names
            cgrp <- sels$cons[i, "Control"]
            tgrp <- sels$cons[i, "Test"]

            #get sample names
            ctrl <- sampleNames(eset)[ sels$rows[[cgrp]] ]
            test <- sampleNames(eset)[ sels$rows[[tgrp]] ]

            #add treatment/group labels to pheno
            pData(eset)[ctrl, "treatment"] <- "ctrl"
            pData(eset)[test, "treatment"] <- "test"

            pData(eset)[ctrl, "group"] <- cgrp
            pData(eset)[test, "group"] <- tgrp


            #add sample names for contrast to mama_samples
            con_name <- paste(tgrp, cgrp, sep="-")
            con_name <- paste(gse_name, con_name, sep="_")
            mama_samples[[con_name]] <- c(ctrl, test)

            #add treatment labels for contrast to mama_clinicals
            labels <- factor(pData(eset)[c(ctrl, test), "treatment"],
                             levels = c("test", "ctrl"))
            clinical <- data.frame(treatment=labels, row.names=c(ctrl, test))
            mama_clinicals[[con_name]] <- clinical
        }
        #create contrast names
        contrasts <- paste(sels$cons$Test, sels$cons$Control, sep="-")

        #store levels for group name variable
        group_levels <- names(sels$rows)

        #retain selected samples only
        eset <- eset[, unique(unlist(sels$rows))]
    }

    #put data together
    mama_data <- list(sample_names=mama_samples, clinicals=mama_clinicals)
    contrast_data <- list(eset=eset, contrasts=contrasts,
                          levels=group_levels, mama_data=mama_data)

    return (contrast_data)
}



#------------------------


# Generate model matrix with surrogate variables.
#
# Used by \code{diff_expr} to create model matrix with surrogate variables
# in order to run \code{diff_anal}.
#
# @param eset Annotated eset with samples selected during \code{add_contrasts}.
# @param group_levels Character vector of unique group names created by
#    \code{add_contrasts}.
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
    svobj <- tryCatch (
        {capture.output(svobj <- sva::sva(exprs(eset), mod, mod0))
            svobj},

        error = function(cond) {
            message("sva failed - continuing without. Could also try to: \n",
                    " - re-select samples (if made mistake) \n",
                    " - select more/fewer contrasts \n",
                    " - remove offending eset \n")
            return(list("sv"=NULL))
        })

    if (is.null(svobj$sv)) {
        modsv <- mod
    } else {
        modsv <- cbind(mod, svobj$sv)
        colnames(modsv) <- c(colnames(mod), paste("SV", 1:svobj$n.sv, sep=""))
    }
    return (list("mod"=mod, "modsv"=modsv, "svobj"=svobj))
}


#------------------------



# Removes features with duplicated annotation.
#
# For rows with duplicated annot, highested IQR retained.
#
# @inheritParams diff_expr
# @inheritParams diff_setup
# @param mod Model matrix without surrogate variables. generated by \code{diff_setup}.
# @param svobj Result from \code{sva} function called during \code{diff_setup}.
#
# @return List with:
#    \item{eset}{Expression set with unique features at probe or gene level.}
#    \item{exprs_sva}{Expression data from eset with effect of surrogate
#       variable removed.}

iqr_duplicates <- function (eset, mod, svobj, annot="SYMBOL") {

    #get eset with surrogate variables modeled out
    exprs_sva <- clean_y(exprs(eset), mod, svobj$sv)

    #add inter-quartile ranges, row, and feature data to exprs data
    data <- as.data.frame(exprs_sva)
    data$IQR <- matrixStats::rowIQRs(exprs_sva)
    data$row <- 1:nrow(data)
    data[, colnames(fData(eset))] <- fData(eset)

    #remove rows with NA annot (occurs if annot is SYMBOL)
    data <- data[!is.na(data[, annot]), ]

    #for rows with same annot, keep highest IQR
    data %>%
        dplyr::group_by_(annot) %>%
        dplyr::arrange(dplyr::desc(IQR)) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() ->
        data


    #use row number to keep selected features
    eset <- eset[data$row, ]
    exprs_sva <- exprs_sva[data$row, ]

    #use annot for feature names
    featureNames(eset) <- fData(eset)[, annot]
    row.names(exprs_sva)   <- fData(eset)[, annot]

    return (list(eset=eset, exprs_sva=exprs_sva))
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
#   non-selected samples removed by \code{iqr_duplicates}.
# @param exprs_sva Expression data with surrogate variables removed. Created by
#    \code{iqr_duplicates}
# @param contrasts Character vector generated by \code{add_contrasts}.
# @param group_levels Character vector of group names generated by
#    \code{add_contrasts}.
# @param mama_data List generated by \code{add_contrasts} then
#    \code{get_contrast_esets}.
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

diff_anal <- function(eset, exprs_sva, contrasts, group_levels, mama_data,
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
        cat (contrast_names[i], "(# p < 0.05):", num_sig, "\n")
    }
    cat("\n")

    #setup plot items
    group <- factor(pData(eset)$group, levels=group_levels)
    palette <- RColorBrewer::brewer.pal(12, "Paired")
    colours <- palette[group]

    #Add extra space to right of plot area
    par(mar = c(5, 4, 2, 6))


    #plot MDS
    limma::plotMDS(exprs_sva, pch=19, main = gse_name, col = colours)
    legend("topright", inset=c(-0.4, 0), legend=group_levels,
           fill=unique(colours), xpd=TRUE, bty="n", cex=0.65)


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
# @param eset Annotated eset created by \code{load_raw}. Non-selected samples
#    and duplicate features removed by \code{add_contrasts} and
#    \code{iqr_duplicates}.
# @param contrasts Character vector of contrast names generated by
#    \code{add_contrasts}.
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
# @param mama_data result of call to \code{add_contrasts}.
# @param dups result of call to \code{iqr_duplicates}.
# @param data Expression data for eset. Used to supply sva adjusted data created
#   by \code{clean_y}.
#
# @seealso  \code{\link{clean_y}}.
# @return A list with:
#    \item{sample_names}{Named list  of character vectors (one per contrast)
#       with sample names.}
#    \item{clinicals}{Named list of dataframes (one per contrast) with rows
#       equal to sample names and 'treatment' column a factor with levels
#       equal to 'test' and 'ctrl'.
#    \item{esets_sva, esets}{List of esets with samples for each contrast only.
#       Expression data has effects of surrogate variables removed or not.}
#

get_contrast_esets <- function(mama_data, dups) {

    sample_names <- mama_data$sample_names
    eset <- dups$eset
    mama_data$esets <- lapply(sample_names, function(x) eset[,x])

    exprs(eset) <- dups$exprs_sva
    mama_data$esets_sva <- lapply(sample_names, function(x) eset[,x])

    return(mama_data)
}


#------------------------


# Adjusts expression data for surrogate variables.
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
