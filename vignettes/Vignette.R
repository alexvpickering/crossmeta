## ------------------------------------------------------------------------
library(crossmeta)

# specify where data will be downloaded
data_dir <- file.path(getwd(), "data", "LY")

# gather GSE names; also gather Illumina GSE names in seperate vector 
# (see 'Checking Raw Illumina Data')
gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
illum_names <- c("GSE50841", "GSE34817", "GSE29689")

## ---- eval=FALSE---------------------------------------------------------
#  # download raw data
#  get_raw(gse_names, data_dir)

## ---- eval=FALSE---------------------------------------------------------
#  # this is why we kept track of Illumina names in seperate vector
#  crossmeta:::open_raw_illum(illum_names, data_dir)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(lydata)

# location of raw data
data_dir <- system.file("extdata", package = "lydata")

## ---- message=FALSE, warning=FALSE, results='hide'-----------------------
esets <- load_raw(gse_names, data_dir = data_dir)

## ---- eval=FALSE---------------------------------------------------------
#  anals <- diff_expr(esets, data_dir)

## ---- fig.height=5, fig.width=5, message=FALSE, warning=FALSE------------
# load auto-saved results of previous call to diff_expr
# (In this case, all of gse_names)
prev <- load_diff(gse_names, data_dir)

# supply prev to diff_expr
# omit "[1]" to analyze all esets
anals <- diff_expr(esets[1], data_dir, prev_anals=prev)

## ---- message=FALSE, warning=FALSE---------------------------------------
# re-load previous analyses if need to
anals <- load_diff(gse_names, data_dir)

## ---- warning=FALSE------------------------------------------------------
# perform meta analysis
es <- es_meta(anals)

