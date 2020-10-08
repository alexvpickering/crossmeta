library(testthat)
library(Biobase)
library(limma)

# RG normalized as in load_agil for two-channel arrays
apoa1_path <- system.file('testdata', 'Apoa1.RData', package = 'crossmeta')
load(apoa1_path)
elist <- RG
elist <- limma::backgroundCorrect(elist, method="normexp", offset=50)
elist <- limma::normalizeWithinArrays(elist, method="loess")
elist <- limma::normalizeBetweenArrays(elist, method="Aquantile")

# setup dummy eset with required pdata columns needed for phenoData.ch2
targets <- elist$targets
colnames(elist) <- targets$geo_accession <- row.names(targets)
elist$genes$SYMBOL <- elist$genes$NAME

targets$label_ch1 <- 'Cy3'
targets$label_ch2 <- 'Cy5'
targets$source_name_ch1 <- targets$Cy3
targets$source_name_ch2 <- targets$Cy5
targets$Cy3 <- targets$Cy5 <- NULL

eset <- ExpressionSet(elist$M,
                      phenoData = as(targets, 'AnnotatedDataFrame'),
                      featureData = as(elist$genes, 'AnnotatedDataFrame'))

pdata <- crossmeta::phenoData.ch2(eset)

test_that("phenoData.ch2 orders all reds first then all greens", {
  
  res <- pdata@data$label
  target <- rep(c('Cy5', 'Cy3'), each = ncol(RG))
 
  expect_equal(res, target)
})

test_that("phenoData.ch2 produces twice as many rows as arrays", {
  
  res <- length(pdata@data$label)
  target <- ncol(RG)*2
  
  expect_equal(res, target)
})

test_that("crossmeta produces similar results to limma", {
  # setups eset for diff_expr
  E <- crossmeta:::exprs.MA(elist)
  eset <- ExpressionSet(E,
                        phenoData = pdata,
                        featureData = as(elist$genes, 'AnnotatedDataFrame'))
  
  # avoid GUI selection
  eset$group <- make.names(eset$source_name)
  eset <- list(Apoa1 = eset)
  prev <- setup_prev(eset, 'ApoAI...-C57BL.6')
  
  data_dir = tempdir()
  dir.create(file.path(data_dir, names(eset)))
  res <- diff_expr(eset, prev_anals = prev, data_dir = data_dir, svanal = FALSE)
  
  # as in limma user guide
})





# setup expression sets
pdata <- data.frame(col1 = letters[1:5], col2 = rep('mouse', 5), row.names = paste0('SAMPLE', 1:5))
mat   <- matrix(rnorm(50), ncol = 5, dimnames = list(1:10, paste0('SAMPLE', 1:5)))
eset <- ExpressionSet(mat, as(pdata, 'AnnotatedDataFrame'))

# make data with sample names equal to col1 of pdata
data <- eset
sampleNames(data) <- pdata$col1


test_that("match_samples correctly orders data when sample names are eset pdata column", {
  
  # trivial order is already correct
  res <- crossmeta:::match_samples(eset, data)
  
  # reorder and rename
  eset <- eset[, res$eset_order]
  data <- data[, res$elist_order]
  colnames(data) <- colnames(eset)
  expect_equal(eset, data)
  
  # shuffle data order
  sampleNames(data) <- pdata$col1
  data <- data[, sample(5)]
  
  # see if recovers
  res <- crossmeta:::match_samples(eset, data)
  eset <- eset[, res$eset_order]
  data <- data[, res$elist_order]
  colnames(data) <- colnames(eset)
  expect_equal(eset, data)
})