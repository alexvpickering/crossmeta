library(crossmeta)
library(lydata)
library(testthat)


#------
#SETUP
#------

setwd("~/Documents/Batcave/GEO/1-meta")
data_dir <- system.file("extdata", package = "lydata")

#load esets
affy_names  <- c("GSE9601", "GSE15069")
illum_names <- c("GSE50841", "GSE34817", "GSE29689")
gse_names   <- c(affy_names, illum_names)

esets <- load_raw(gse_names, data_dir)
afy <- esets[1]
ilm <- esets[3]

#commonize
esets <- commonize(esets)

#load and re-run analysis
prev <- load_diff(gse_names, data_dir)
re_anal <- diff_expr(esets, data_dir, prev_anals=prev)


#-------
# TESTS
#-------


test_that("get_raw downloads/unpacks affy, illum, and agil", {
    skip_on_bioc()

    #download a small GSE from each manufacturer
    expect_null(get_raw("GSM1318806", data_dir))  #affy
    expect_null(get_raw("GSE41845", data_dir))  #illum
    expect_null(get_raw("GSE66133", data_dir))  #agil
})

#cleanup
unlink(paste(data_dir, c("GSM1318806", "GSE41845", "GSE66133"), sep="/"),
       recursive=TRUE)

#---------


test_that("load_raw loads and annotates raw data", {
    skip_on_bioc()

    #SYMBOL annotation?
    expect_equal(fData(afy[[1]])$SYMBOL[1], "MAPK3")
    expect_equal(fData(ilm[[1]])$SYMBOL[1], "EEF1A1")

    #PROBE annotation?
    expect_equal(fData(afy[[1]])$PROBE[1], "1000_at")
    expect_equal(fData(ilm[[1]])$PROBE[1], "ILMN_1343291")

    #map 1:many?
    expect_gt(length(featureNames(afy[[1]])),
              length(unique(featureNames(afy[[1]]))))
    expect_gt(length(featureNames(ilm[[1]])),
              length(unique(featureNames(ilm[[1]]))))

    #different features?
    expect_false(same_features(c(afy, ilm)))
})


#---------

# Check for commonize.
#
# @param esets List of expression sets.
#
# @return TRUE if featureNames of all esets are the same. Otherwise, FALSE.
#

same_features <- function(esets) {

    #same number?
    n_features <- sapply(esets, function(eset) nrow(eset))
    if (diff(range(n_features) >= 1)) return(FALSE)

    #same names?
    feature_names <- lapply(esets, function(eset) featureNames(eset))
    com <- Reduce(intersect, feature_names)
    if (length(com) != length(feature_names[[1]])) return(FALSE)

    return(TRUE)
}


test_that("load_diff and diff_expr work", {
    skip_on_bioc()

    #GSEs all load?
    expect_equal(names(prev), gse_names)

    #same features?
    expect_true(same_features(lapply(prev, function(x) x$eset)))

    #rerun gives same result?
    expect_equal(re_anal, prev)
})


#-----------
# CLEANUP
#-----------
rm(afy, ilm, affy_names, illum_names, data_dir,
   esets, gse_names, prev, re_anal, same_features)
