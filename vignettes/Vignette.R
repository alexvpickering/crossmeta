## ------------------------------------------------------------------------
library(crossmeta)

# specify where data will be downloaded
data_dir <- file.path(getwd(), "data", "LY")

# gather all GSEs
gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")

# gather Illumina GSEs (see 'Checking Raw Illumina Data')
illum_names <- c("GSE50841", "GSE34817", "GSE29689")

# download raw data
# get_raw(gse_names, data_dir)

## ---- eval=FALSE---------------------------------------------------------
#  # this is why we gathered Illumina GSEs
#  open_raw_illum(illum_names, data_dir)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(lydata)

# location of raw data
data_dir <- system.file("extdata", package = "lydata")

## ---- message=FALSE, warning=FALSE, results='hide'-----------------------
# reloads if previously called
esets <- load_raw(gse_names, data_dir)

## ----eval = FALSE--------------------------------------------------------
#  library(Biobase)
#  library(AnnotationDbi)
#  
#  # check feature data to see what columns are available
#  head(fData(esets$GSE15069))
#  
#  # if using RStudio
#  # View(fData(esets$GSE15069))
#  
#  # annotation package for appropriate species
#  library(org.Mm.eg.db)
#  
#  # map from accession number to entrez gene ids
#  acnums  <- as.character(fData(esets$GSE15069)$GB_ACC)
#  enids   <- mapIds(org.Mm.eg.db, acnums, "ENTREZID", "ACCNUM")
#  
#  # add 'GENE_ID' column with entrez ids
#  fData(esets$GSE15069)$GENE_ID <- enids
#  
#  # use crossmeta to map from entrez gene ids to homologous hgnc symbol
#  esets$GSE15069 <- symbol_annot(esets$GSE15069)
#  
#  # to overwrite saved eset (to avoid repeating above)
#  saveRDS(esets$GSE15069, file.path(data_dir, "GSE15069", "GSE15069_eset.rds"))

## ---- eval=FALSE---------------------------------------------------------
#  anals <- diff_expr(esets, data_dir)

## ------------------------------------------------------------------------
# load auto-saved results of previous call to diff_expr
prev <- load_diff(gse_names, data_dir)

# supply prev to diff_expr
# anals <- diff_expr(esets, data_dir, prev_anals=prev)

## ---- message=FALSE, results='hide'--------------------------------------
# re-load previous analyses if need to
anals <- load_diff(gse_names, data_dir)

# perform meta analysis
es <- es_meta(anals)

# for explanation of values
# ?es_meta

## ----eval=FALSE----------------------------------------------------------
#  # subject is the focus of the meta-analysis (e.g. drug/disease name)
#  contribute(anals, subject = "LY294002")
#  
#  # Thank you!

