library(testthat)
library(Biobase)

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
    data <- data[, res$data_order]
    colnames(data) <- colnames(eset)
    expect_equal(eset, data)

    # shuffle data order
    sampleNames(data) <- pdata$col1
    data <- data[, sample(5)]

    # see if recovers
    res <- crossmeta:::match_samples(eset, data)
    eset <- eset[, res$eset_order]
    data <- data[, res$data_order]
    colnames(data) <- colnames(eset)
    expect_equal(eset, data)
})


test_that("match_samples correctly order data when fewer data samples than eset samples", {

    # trivial case where sample names are the same
    data <- data[, sample(3)]
    res <- crossmeta:::match_samples(eset, data)
    eset <- eset[, res$eset_order]
    data <- data[, res$data_order]
    colnames(data) <- colnames(eset)
    expect_equal(eset, data)

    #  data sample names are pdata col1 and shuffled order
    eset <- data <- ExpressionSet(mat, as(pdata, 'AnnotatedDataFrame'))
    sampleNames(data) <- pdata$col1
    data <- data[, sample(3)]
    res <- crossmeta:::match_samples(eset, data)
    eset <- eset[, res$eset_order]
    data <- data[, res$data_order]
    colnames(data) <- colnames(eset)
    expect_equal(eset, data)
})


test_that("match_samples correctly order data when fewer eset samples than data samples", {



    # trivial case where sample names are the same
    eset <- data <- ExpressionSet(mat, as(pdata, 'AnnotatedDataFrame'))
    eset <- eset[, sample(3)]
    res <- crossmeta:::match_samples(eset, data)
    eset <- eset[, res$eset_order]
    data <- data[, res$data_order]
    colnames(data) <- colnames(eset)
    expect_equal(eset, data)

    #  data sample names are pdata col1 and shuffled order
    eset <- data <- ExpressionSet(mat, as(pdata, 'AnnotatedDataFrame'))
    sampleNames(data) <- pdata$col1
    eset <- eset[, sample(3)]
    res <- crossmeta:::match_samples(eset, data)
    eset <- eset[, res$eset_order]
    data <- data[, res$data_order]
    colnames(data) <- colnames(eset)
    expect_equal(eset, data)
})

testthat("match_samples correctly order data when sample names are contained in pdata column", {

    # sample sample number
    pdata$col1 <- paste('condition: treatment', 1:5)
    eset <- data <- ExpressionSet(mat, as(pdata, 'AnnotatedDataFrame'))
    colnames(data) <- paste('TREATMENT', 1:5)
    data <- data[, sample(5)]
    res <- crossmeta:::match_samples(eset, data)
    eset <- eset[, res$eset_order]
    data <- data[, res$data_order]
    colnames(data) <- colnames(eset)
    expect_equal(eset, data)

    # fewer data
    eset <- data <- ExpressionSet(mat, as(pdata, 'AnnotatedDataFrame'))
    colnames(data) <- paste('TREATMENT', 1:5)
    data <- data[, sample(3)]
    res <- crossmeta:::match_samples(eset, data)
    eset <- eset[, res$eset_order]
    data <- data[, res$data_order]
    colnames(data) <- colnames(eset)
    expect_equal(eset, data)

    # fewer eset
    eset <- data <- ExpressionSet(mat, as(pdata, 'AnnotatedDataFrame'))
    colnames(data) <- paste('TREATMENT', 1:5)
    eset <- eset[, sample(3)]
    res <- crossmeta:::match_samples(eset, data)
    eset <- eset[, res$eset_order]
    data <- data[, res$data_order]
    colnames(data) <- colnames(eset)
    expect_equal(eset, data)
})
