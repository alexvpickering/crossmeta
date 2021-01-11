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
elist$genes$rn <- seq_len(nrow(elist))

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
  prev <- crossmeta::setup_prev(eset, 'ApoAI...-C57BL.6')
  
  data_dir = tempdir()
  dir.create(file.path(data_dir, names(eset)))
  res <- crossmeta::diff_expr(eset, prev_anals = prev, data_dir = data_dir, svanal = FALSE, annot='rn')
  
  # cleanup
  unlink('Rplots.pdf')
  res <- res$Apoa1$top_tables$`Apoa1_ApoAI...-C57BL.6`
  
  # as in limma user guide
  MA <- backgroundCorrect(RG, method="normexp", offset=50)
  MA <- normalizeWithinArrays(MA, method = 'loess')
  MA <- normalizeBetweenArrays(MA, method = 'Aquantile')
  design <- cbind("Control-Ref"=1,"KO-Control"=MA$targets$Cy5=="ApoAI-/-")
  fit <- lmFit(MA, design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=2, n = Inf)
  tt <- tt[!is.na(tt$NAME), ]
  
  # annotation by row name so should be unchanged
  expect_equal(nrow(tt), nrow(res))
  
  # correlation between logFCs
  lfc1 <- res$logFC
  names(lfc1) <- row.names(res)
  
  lfc2 <- tt$logFC
  names(lfc2) <- row.names(tt)
  lfc2 <- lfc2[names(lfc1)]
  
  cor12 <- cor(lfc1, lfc2)
  expect_gt(cor12, 0.96)
})

