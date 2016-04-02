library(dplyr)
library(matrixStats)
library(RColorBrewer)
library(tcltk)
library(sva)
library(limma)

source("~/Documents/Batcave/GEO/1-meta/load_utils.R")

#------------------------

diff_expr <- function (eset, data_dir, annot="gene") {
    #wrapper for diff_expr_one
    diff_expr <- mapply(diff_expr_one, eset, names(eset), data_dir, annot, SIMPLIFY=F)
    names(diff_expr) <- names(eset)
    return (diff_expr)
}


diff_expr_one <- function (eset, gse_name, data_dir, annot) {
    #INPUT:
    #OUTPUT:
    gse_dir <- paste(data_dir, gse_name, sep="/")

    #select contrasts
    cons <- add_contrasts(eset, gse_name, annot)

   #differential expression
    setup <- diff_setup(cons$eset, cons$levels)

    anal <- diff_anal(cons$eset, cons$contrasts, cons$levels, cons$mama_data,
                      setup$mod, setup$modsv, setup$svobj,
                      gse_dir, gse_name, annot)

    return (anal)
}


#------------------------

add_contrasts <- function (eset, gse_name, annot="gene") {

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
        ctrl <- tk_select.list(choices, multiple=T, title="Control samples for contrast")
        ctrl <- str_extract(ctrl, "GSM[0-9]+")
        if (length(ctrl) == 0) {break}

        #select test samples
        test <- tk_select.list(choices, multiple=T, title="Test samples for contrast")
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
        group_levels <- unique(c(group_levels, group_names[1], group_names[2]))

        #add to selected_samples
        selected_samples <- unique(c(selected_samples, ctrl, test))

        #store sample names for each contrast
        contrast_name <- paste(gse_name, contrast, sep="_")
        mama_samples[[contrast_name]] <- c(ctrl, test)

        #add labels for contrast to mama_clinicals
        labels <- factor(pData(eset)[c(ctrl, test), "treatment"], levels = c("test", "ctrl"))   #  (test, ref)
        clinical <- data.frame(treatment = labels, row.names = c(ctrl, test))
        mama_clinicals[[contrast_name]] <- clinical
    }
    #retain selected samples only
    eset <- eset[, selected_samples]

    if (annot == "gene") {
        #remove rows with duplicated SYMBOL
        eset <- iqr_duplicates(eset)
    } else if (annot == "probe") {
        #remove rows with duplicated probes
        eset <- eset[unique(featureNames(eset)), ]
    }

    #subset eset for each contrast (for MAMA)
    mama_esets <- get_contrast_esets(mama_samples, eset)

    #put data together
    mama_data <- list(esets=mama_esets, sample_names=mama_samples, clinicals=mama_clinicals)

    contrast_data <- list("eset"=eset, "contrasts"=contrasts, "levels"=group_levels,
                       "samples"=selected_samples, "mama_data"=mama_data)

    return (contrast_data)
}

#------------------------

iqr_duplicates <- function (eset) {
    #for rows with duplicated SYMBOL highested IQR retained

    data <- as.data.frame(exprs(eset))
    #add inter-quartile ranges to data
    data$IQR <- rowIQRs(data)
    #add feature data
    data[, colnames(fData(eset))] <- fData(eset)

    #for rows with same SYMBOL, keep highest IQR
    data %>%
        group_by(SYMBOL) %>%
        arrange(desc(IQR)) %>%
        slice(1) %>%
        ungroup ->
        data

    #seperate exprs and fdata columns
    exprs <- as.matrix(data[, sampleNames(eset)])
    fdata <- as.matrix(data[, colnames(fData(eset))])
    row.names(exprs) <- data$SYMBOL
    row.names(fdata) <- data$SYMBOL

    #put exprs and fdata into eset
    exprs(eset) <- exprs
    fData(eset) <- as.data.frame(fdata)

    return (eset)
}

#------------------------

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

diff_anal <- function(eset, contrasts, group_levels, mama_data,
                      mod, modsv, svobj, gse_dir, gse_name, annot="gene"){

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
    mama_data$esets_sva <- get_contrast_esets(mama_data$sample_names, eset, sva_exprs)

    #save to disk
    diff_expr <- list(eset=eset, top_tables=top_tables, 
                      ebayes_sv=ebayes_sv, ebayes=ebayes, mama_data=mama_data)

    if (annot == "gene") {
        save_name <- paste (gse_name, "diff_expr.rds", sep="_")

    } else if (annot == "probe") {
        save_name <- paste (gse_name, "diff_expr_probe.rds", sep="_")
    }

    saveRDS(diff_expr, file = paste(gse_dir, save_name, sep="/"))
    return (diff_expr)
}

#------------------------

fit_ebayes <- function(eset, contrasts, mod) {
    contrast_matrix <- makeContrasts(contrasts=contrasts, levels=mod)
    fit <- contrasts.fit (lmFit(exprs(eset),mod), contrast_matrix)  
    return (eBayes(fit))
}

#------------------------

get_contrast_esets <- function(sample_names, eset, data=exprs(eset)) {
    #used to setup esets for MAMA 
    #data = exprs_sva will return sva-adjusted eset for each contrast
    exprs(eset) <- data
    lapply(sample_names, function(x) eset[,x])
}

#------------------------

clean_y <- function(y, mod, svs) {
    #IN:  y   = exprs(eset), 
    #     mod = full model matrix supplied to sva, 
    #     svs = surrogate variables returned by sva (svobj$sv)
    #OUT: exprs(eset) adjusted for surrogate variables

    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    rm(Hat)
    gc()
    P = ncol(mod)
    return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}