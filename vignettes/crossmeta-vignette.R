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

## ---- message=FALSE, warning=FALSE, results='hide', fig.keep='none'------
library(Biobase)

# load eset
gse_name  <- c("GSE34817")
eset <- load_raw(gse_name, data_dir)

# inspect pData of eset
# View(pData(eset$GSE34817))  # if using RStudio
head(pData(eset$GSE34817))    # otherwise

# get group info from pData (differs based on eset)
group <- pData(eset$GSE34817)$characteristics_ch1.1

# make group names concise and valid
group <- gsub("treatment: ", "", group)
group <- make.names(group)

# add group to eset pData
pData(eset$GSE34817)$group <- group

# setup selections
sel <- setup_prev(eset, contrasts = "LY-DMSO")

# run differential expression analysis
anal <- diff_expr(eset, data_dir, prev_anal = sel)

## ---- message=FALSE, results='hide'--------------------------------------
# re-load previous analyses if need to
anals <- load_diff(gse_names, data_dir)

# perform meta analysis
es <- es_meta(anals)

# for explanation of values
# ?es_meta

## ---- message=FALSE, results='hide'--------------------------------------
library(ccmap)
library(ccdata)
data(cmap_es)

# get meta-analysis effect size values
dprimes <- get_dprimes(es)

# query using entire transcriptional profile
top_drugs <- query_drugs(dprimes$meta, cmap_es)

# drug with greatest transcriptional similarity
topd <- names(top_drugs)[1]

## ----eval=FALSE----------------------------------------------------------
#  # subject is the focus of the meta-analysis (e.g. drug/disease name)
#  contribute(anals, subject = "LY294002")
#  
#  # Thank you!

## ---- message=FALSE, results='hide'--------------------------------------
# perform pathway analysis for each contrast
# path_anals <- diff_path(esets, anals, data_dir)

# load previous pathway analyses
path_anals <- load_path(gse_names, data_dir)

## ---- message=FALSE, results='hide'--------------------------------------
# perform pathway meta analysis 
# (default ncores and nperm will be faster and give more reliable fdr/p-values)
path_res <- path_meta(path_anals, ncores = 1, nperm = 100)

# highest ranking pathway
path <- row.names(path_res)[1]

## ---- message=FALSE, results='hide'--------------------------------------
# plot transcript expression values and summary error bars for pathway
# plot_path(path, es, cmap_es)

# drug with most similar transcriptional profile for top pathway

# query using only transcripts for genes in path
top_drugsp <- query_drugs(dprimes$meta, cmap_es, path=path)

# drug with most similar transcriptional profile for top pathway (milrinone)
topdp <- names(top_drugsp)[1] 

# plot meta-analysis of pathway alongside top overall drug and top drug for pathway
plot_path(path, es, cmap_es, c(topd, topdp))

## ------------------------------------------------------------------------
sessionInfo()

