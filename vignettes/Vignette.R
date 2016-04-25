## ------------------------------------------------------------------------
library(crossmeta)

#specify where data will be downloaded
data_dir <- file.path(getwd(), "data", "LY")

#gather GSE names
affy_names  <- c("GSE9601", "GSE15069")
illum_names <- c("GSE50841", "GSE34817", "GSE29689")
gse_names   <- c(affy_names, illum_names)

## ---- eval=FALSE---------------------------------------------------------
#  #download raw data
#  get_raw(gse_names, data_dir)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(lydata)

#location of raw data
data_dir <- system.file("extdata", package = "lydata")

## ---- message=FALSE, warning=FALSE, results='hide'-----------------------
esets <- load_raw(gse_names, data_dir)

## ------------------------------------------------------------------------
com_esets <- commonize(esets)

## ---- eval=FALSE---------------------------------------------------------
#  anals <- diff_expr(com_esets, data_dir)

## ---- fig.height=5, fig.width=5, message=FALSE, warning=FALSE------------
#load auto-saved results of previous call to diff_expr
#(In this case, all of gse_names)
prev <- load_diff(gse_names, data_dir)

#supply prev to diff_expr
#omit "[1]" to analyze all esets
anals <- diff_expr(com_esets[1], data_dir, prev_anals=prev)

## ---- message=FALSE, warning=FALSE---------------------------------------
#re-load previous analyses if need to
anals <- load_diff(gse_names, data_dir)

#MetaArray object (effects of surrogate variables removed)
ma_sva <- make_ma(anals, sva=TRUE)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(MAMA)

es <- ES.GeneMeta(ma_sva, "treatment", nperm=10)

#save signature for further analysis (e.g. by ccmap - see below)
#saveRDS(es, "es.rds")

