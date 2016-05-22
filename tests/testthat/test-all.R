library(crossmeta)
library(testthat)
library(Biobase)
library(org.Hs.eg.db)



# ------------------- Setup Toy Data


hgene <- data.frame(ENTREZID    = c(469356, 455237, 451528, 470565, 468403, 100608083),
                    ENTREZID_HS = c(34, 37, 38, 90, 6442, 158))

dfdat <- data.frame(row.names = paste(1:6, "at", sep = "_"),
                    ProbeName = paste(1:6, "at", sep = "_"),
                    Symbol = c("ACADM", "ACADVL", "ACAT1", "ACVR1", "SGCA", "ADSL"),
                    Entrez_Gene_ID = hgene$ENTREZID,
                    stringsAsFactors = FALSE)


efdat <- data.frame(row.names = row.names(dfdat)[3:1],
                    ID = row.names(dfdat)[3:1],
                    GeneSymbol = dfdat$Symbol[3:1],
                    Gene_ID = dfdat$Entrez_Gene_ID[3:1],
                    stringsAsFactors = FALSE)

# set last two samples (test group) to have higher expression
set.seed(0)
expr <- matrix(abs(rnorm(12, 0, 0.2)), nrow = 3, ncol = 4)
expr[, 3:4] <- expr[, 3:4] + 1

eset <- ExpressionSet(expr, featureData = AnnotatedDataFrame(efdat))
data <- list(genes = dfdat)


# set annotation to platform without bioconductor annotation data package
annotation(eset) <- "GPL4032"



# ------------------- Fix Agil Features


test_that("fix_agil_features: eset rownames to best match for data ProbeName", {

    #change eset row names to gene symbols
    row.names(eset) <- fData(eset)$GeneSymbol
    eset <- crossmeta:::fix_agil_features(eset, data)

    expect_equal(row.names(eset), fData(eset)$ID)


    #change data ProbeName to Symbol
    data$genes$ProbeName <- data$genes$Symbol
    eset <- crossmeta:::fix_agil_features(eset, data)

    expect_equal(row.names(eset), fData(eset)$GeneSymbol)
})



# ------------------- Fix Illum Features


test_that("fix_illum_features: maps data rownames to eset rownames", {

    #change eset row names to gene
    row.names(eset) <- fData(eset)$GeneSymbol
    dfdat <- crossmeta:::fix_illum_features(eset, dfdat)

    expect_equal(sum(row.names(dfdat) %in% row.names(eset)), 3)

    #restore eset row names to probe id
    row.names(eset) <- fData(eset)$ID
    dfdat <- crossmeta:::fix_illum_features(eset, dfdat)

    expect_equal(sum(row.names(dfdat) %in% row.names(eset)), 3)

})



# ------------------- Merge Feature Data


test_that("merge_fdata merges on rownames, preserving data nrow and order", {

    fdat <- crossmeta:::merge_fdata(fData(eset), dfdat)

    # expected order and dimensions
    expect_equal(row.names(fdat), row.names(dfdat))
    expect_equal(ncol(fdat), ncol(efdat) + ncol(dfdat))

    # NAs only in rows where no eset feature data
    expect_equal(sum(is.na(fdat$GeneSymbol)), nrow(dfdat) - nrow(efdat))

    # matched rows aligned
    expect_equal(fdat[row.names(efdat), "GeneSymbol"],
                 fdat[row.names(efdat), "Symbol"])


})


# ------------------- Symbol Annotation


test_that("symbol_annot finds entrez in eset fData and maps to hgnc symbol", {

    # gets correct hgnc symbol
    eset <- crossmeta:::symbol_annot(eset, hgene, "GSE1")

    expect_equal(fData(eset)$GeneSymbol, fData(eset)$SYMBOL)


    # Gene_ID is column being used
    fData(eset)$Gene_ID <- fData(eset)$Gene_ID[3:1]
    eset <- crossmeta:::symbol_annot(eset, hgene, "GSE1")

    expect_equal(fData(eset)$GeneSymbol, fData(eset)$SYMBOL[3:1])


    # mapping is through homologene
    hgene$ENTREZID[1:3] <- hgene$ENTREZID[3:1]
    eset <- crossmeta:::symbol_annot(eset, hgene)

    expect_equal(fData(eset)$GeneSymbol, fData(eset)$SYMBOL)
    expect_equal(efdat$Gene_ID, fData(eset)$Gene_ID[3:1])
})



# ------------------- Differential Expression

# Setup

dir.create("GSE1")

# annotated eset
pData(eset)$title <- c("VEH", "VEH", "FOO", "FOO")
eset  <- crossmeta:::symbol_annot(eset, hgene)
esets <- list(GSE1 = eset)

# previous analysis mimic
pData(eset)$treatment <- c("ctrl", "ctrl", "test", "test")
pData(eset)$group     <- c("VEH", "VEH", "FOO", "FOO")
pData(eset)$pairs     <- NA

contrasts <- limma::makeContrasts("FOO-VEH", levels = c("VEH", "FOO"))

prev <- list(GSE1 = list(eset = eset,
                         ebayes = list(contrasts = contrasts)))

# run analysis on eset, using previous mimic
anals <- diff_expr(esets, prev_anals = prev)



# Tests

test_that("diff_expr: test group has higher expression", {

    # check t-statistics
    t <- anals$GSE1$top_tables$`GSE1_FOO-VEH`$t
    expect_equal(sign(t), c(1, 1, 1))
})


test_that("sva not affected by replicate rows", {

    # replicate each feature
    esets[[1]] <- esets[[1]][rep(1:3, 100), ]

    # re-run and test
    anals2 <- diff_expr(esets, prev_anals = prev)
    expect_equal(anals, anals2)
})

test_that("diff_expr removes duplicates based on IQR", {

    # duplicate each feature
    esets[[1]] <- esets[[1]][rep(1:3, 2), ]

    # add random noise to duplicates
    set.seed(0)
    expr <- matrix(abs(rnorm(12, 0, 0.2)), nrow = 3, ncol = 4)
    exprs(esets[[1]])[1:3, ] <- expr

    # re-run and test expressions data
    anals2 <- diff_expr(esets, prev_anals = prev)
    expect_equal(exprs(anals$GSE1$eset), exprs(anals2$GSE1$eset))
})



# ----------------- Cleanup

rm(hgene, dfdat, efdat, expr, eset, data, esets, prev, anals, contrasts)

unlink("GSE1", recursive = TRUE)

