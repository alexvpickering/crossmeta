library(crossmeta)
library(testthat)
library(Biobase)
library(org.Hs.eg.db)

# bugs:
# -----

# 1) load_raw(overwrite = "SYMBOL") interpreted as overwrite = TRUE


# Setup Toy Data-------------------


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

eset <- ExpressionSet(expr, featureData = AnnotatedDataFrame(efdat), phenoData = AnnotatedDataFrame(data.frame(organism = rep('Homo sapiens', 4))))
data <- list(genes = dfdat)


# set annotation to platform without bioconductor annotation data package
annotation(eset) <- "GPL4032"



# Merge Feature Data  -------------------


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
