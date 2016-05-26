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
#  crossmeta:::open_raw_illum(illum_names, data_dir)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(lydata)

# location of raw data
data_dir <- system.file("extdata", package = "lydata")

## ---- message=FALSE, warning=FALSE, results='hide'-----------------------
# for initial loading:
# homologene_path <- "path/to/homologene.data"
# esets <- load_raw(gse_names, homologene_path, data_dir)

# to reload:
esets <- load_raw(gse_names, data_dir = data_dir)

## ----eval = FALSE--------------------------------------------------------
#  # check feature data to see what columns are available
#  library(Biobase)
#  View(fData(esets$GSE29689))
#  View(fData(esets$GSE15069))
#  
#  # if human platform and hgnc symbol present in fData:
#  # ---------------------------------------------------
#  
#  # use column with hgnc symbol for annotation
#  fData(esets$GSE29689)$SYMBOL <- toupper(fData(esets$GSE29689)$Symbol)
#  
#  # to overwrite saved eset (avoids repeating above)
#  saveRDS(esets$GSE29689, file.path(data_dir, "GSE29689", "GSE29689_eset.rds"))
#  
#  
#  # if non-human platform and accession numbers present in fData:
#  # -------------------------------------------------------------
#  
#  # load annotation package for appropriate species
#  library(org.Mm.eg.db)
#  
#  # map from accession number to entrez gene ids
#  ac_nums <- as.character(fData(esets$GSE15069)$GB_ACC)
#  map <- AnnotationDbi::select(org.Mm.eg.db, ac_nums, "ENTREZID", "ACCNUM")
#  
#  # add entrez gene ids to fData 'GENE_ID' column
#  fData(esets$GSE15069)$GENE_ID <- map$ENTREZID
#  
#  # use crossmeta to map from entrez gene ids to homologous hgnc symbol
#  homologene <- crossmeta:::get_homologene(homologene_path)
#  esets$GSE15069 <- crossmeta:::symbol_annot(esets$GSE15069, homologene)
#  
#  # to overwrite saved eset (avoids repeating above)
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

# to see results of meta-analysis
View(es$filt)

# for explanation of values
# ?es_meta

