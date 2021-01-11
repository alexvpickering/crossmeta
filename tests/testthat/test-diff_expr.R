context("sanity check diff_expr")
library(Biobase)
library(limma)
library(tximportData)


# Simulate gene expression data for 100 probes and 6 microarrays
# Microarray are in two groups
# First two probes are differentially expressed in second group
# Std deviations vary between genes with prior df=4
sd <- 0.3*sqrt(4/rchisq(1000,df=4))
y <- matrix(rnorm(1000*6,sd=sd),1000,6)
rownames(y) <- paste("Gene",1:1000)
y[1:2,4:6] <- y[1:2,4:6] + 2

y[100:200, c(1,2,6)] <- y[100:200, c(1,2,6)] + 0.5

pdata <- data.frame(group = c(rep('healthy', 3), rep('disease', 3)))
fdata <- data.frame(SYMBOL = row.names(y),
                    PROBE = row.names(y), 
                    row.names = row.names(y))

eset <- ExpressionSet(y,
                      phenoData = as(pdata, 'AnnotatedDataFrame'),
                      featureData = as(fdata, 'AnnotatedDataFrame'))


eset <- list(GSE1 = eset)
data_dir <- tempdir()
dir.create(file.path(data_dir, 'GSE1'))

test_that("diff_expr runs without sva", {
  
  # run diff_expr without sva
  # mock previous analysis (to skip UI)
  prev_anal <- setup_prev(eset, 'disease-healthy')
  expect_error(diff_expr(eset, data_dir, svanal = FALSE, prev_anal = prev_anal), NA)
  
  # cleanup
  unlink(list.files(file.path(data_dir, 'GSE1'), full.names = TRUE))
  
})

test_that("diff_expr runs with sva", {
  
  # mock previous analysis (to skip UI)
  prev_anal <- setup_prev(eset, 'disease-healthy')
  expect_error(diff_expr(eset, data_dir, svanal = TRUE, prev_anal = prev_anal), NA)
  
  # cleanup
  unlink(list.files(file.path(data_dir, 'GSE1'), full.names = TRUE))
})
