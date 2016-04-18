## ---- eval=FALSE---------------------------------------------------------
#  library(crossmeta)
#  
#  #specify where data will be downloaded
#  data_dir <- file.path(getwd(), "data", "LY294002")
#  
#  #gather GSE names into a vector for each manufacturer
#  affy_names  <- c("GSE9601", "GSE15069")
#  illum_names <- c("GSE50841", "GSE34817", "GSE29689")
#  
#  #obtain the raw data(supplementary files).
#  get_raw_affy(affy_names, data_dir)
#  get_raw_illum(illum_names, data_dir)
#  

## ---- eval=FALSE---------------------------------------------------------
#  illum_names <- open_raw_illum(illum_names, data_dir)

## ---- eval=FALSE---------------------------------------------------------
#  affy_esets <- load_affy(affy_names, data_dir)
#  illum_esets <- load_illum(illum_names, data_dir)
#  
#  #put esets together
#  esets <- c(affy_esets, illum_esets)

## ---- eval=FALSE---------------------------------------------------------
#  com_esets <- commonize(esets)

## ---- eval=FALSE---------------------------------------------------------
#  anals <- diff_expr(com_esets, data_dir)

## ---- eval=FALSE---------------------------------------------------------
#  #load auto-saved results of diff_expr that you previously selected samples for
#  prev_names <- c("GSE9601", "GSE15069")
#  prev <- load_diff(prev_names, data_dir)
#  
#  #supply prev to diff_expr
#  anals <- diff_expr(com_esets, data_dir, prev_anals=prev)

## ---- eval=FALSE---------------------------------------------------------
#  #re-load previous analyses if need to
#  gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")
#  anals <- load_diff(gse_names, data_dir)
#  
#  #MetaArray object (effect of SVs removed)
#  ma_sva <- make_ma(anals, sva=T)

## ---- eval=FALSE---------------------------------------------------------
#  library(MAMA)
#  
#  es <- ES.GeneMeta(ma_sva, "treatment", nperm=100)

## ---- eval=FALSE---------------------------------------------------------
#  saveRDS(es, "LY_es.rds")

