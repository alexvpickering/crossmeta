library(testthat)
library(Biobase)

# setup expression sets
pdata <- data.frame(col1 = letters[1:5], col2 = rep('mouse', 5), row.names = paste0('SAMPLE', 1:5))
mat   <- matrix(rnorm(50), ncol = 5, dimnames = list(1:10, paste0('SAMPLE', 1:5)))
fdata <- data.frame(SYMBOL = c('a', 'a', 'b', 'b', letters[3:8]), PROBE = 1:10, stringsAsFactors = FALSE)

eset <- ExpressionSet(mat, as(pdata, 'AnnotatedDataFrame'), as(fdata, 'AnnotatedDataFrame'))

test_that("which_max_iqr returns one row per unique value of group_by", {
  
  # same length as unique group_by values
  iqr_rows <- which_max_iqr(eset, 'SYMBOL')
  expect_length(iqr_rows, length(unique(fdata$SYMBOL)))
  
  # same unique group_by values
  expect_equal(unique(fdata$SYMBOL), fdata$SYMBOL[iqr_rows])
  
})