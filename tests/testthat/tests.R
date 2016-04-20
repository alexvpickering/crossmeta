library(crossmeta)


#------
#SETUP
#------

setwd("~/Documents/Batcave/GEO/1-meta")
data_dir <- paste(getwd(), "testdat", sep="/")

#load esets
afy_name<- "GSE9601"
ilm_name <- "GSE34817"

afy <- load_affy(afy_name,  data_dir)
ilm <- load_illum(ilm_name, data_dir)

#commonize
esets <- commonize(c(afy, ilm))

#load and re-run analysis
prev <- load_diff(c(afy_name, ilm_name), data_dir)
re_anal <- diff_expr(esets, data_dir, prev_anals=prev)

#-------
# TESTS
#-------


test_that("get_raw downloads/unpacks affy, illum, and agil", {
  skip_on_bioc()

  expect_null(get_raw("GSE54558", data_dir))  #affy
  expect_null(get_raw("GSE41845", data_dir))  #illum
  expect_null(get_raw("GSE66133", data_dir))  #agil
})

#cleanup
unlink(paste(data_dir, c("GSE54558", "GSE41845", "GSE66133"), sep="/"),
       recursive=T)

#---------


test_that("load_affy/illum loads and annotates raw data", {
  skip_on_bioc()

  #SYMBOL annotation?
  expect_equal(fData(afy[[1]])$SYMBOL[1], "MAPK3")
  expect_equal(fData(ilm[[1]])$SYMBOL[1], "EEF1A1")

  #PROBE annotation?
  expect_equal(fData(afy[[1]])$PROBE[1], "1000_at")
  expect_equal(fData(ilm[[1]])$PROBE[1], "ILMN_1343291")

  #map 1:many?
  expect_gt(length(featureNames(afy[[1]])), length(unique(featureNames(afy[[1]]))))
  expect_gt(length(featureNames(ilm[[1]])), length(unique(featureNames(ilm[[1]]))))

  #different features?
  expect_false(same_features(c(afy, ilm)))
})


#---------


test_that("load_diff and diff_expr work", {
  skip_on_bioc()

  #GSEs all load?
  expect_equal(names(anals), gse_names)

  #same features?
  expect_true(same_features(lapply(prev, function(x) x$eset)))

  #rerun gives same result?
  expect_equal(re_anal, prev)
})
