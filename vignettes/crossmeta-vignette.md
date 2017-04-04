---
title: "Cross-Platform Meta Analysis"
author: "Alex Pickering"
date: "2017-04-03"
output:
  html_document:
    theme: flatly
---
<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{crossmeta vignette}
-->


`crossmeta` streamlines the cross-platform effect-size and pathway meta-analysis of 
microarray data. For the analysis, you will need a list of Affymetrix, Illumina,
and/or Agilent GSE numbers from [GEO](http://www.ncbi.nlm.nih.gov/geo/). All 21
species in the current [homologene](http://1.usa.gov/1TGoIy7) build are supported.

-----------------

### Obtaining Raw Data

Search [GEO](http://www.ncbi.nlm.nih.gov/geo/) to find relevant microarray 
data for the meta-analysis. For this example, I searched for the PI3K 
inhibitor LY-294002. The search was filtered as follows:

* **Entry type**: Series
* **Organism**: Homo sapiens and Mus musculus
* **Study type**: Expression profiling by array

This search produced 35 hits from which five were chosen. In practice, it is 
best to select all GSEs matching the necessary criteria (supplementary raw data
files available, multiple samples for each treatment, and from either Affymetrix,
Illumina, or Agilent single channel arrays).

After identifying GSEs for the meta-analysis, download and decompress the raw
data as follows. For raw Illumina data only, you must check that it has the 
correct format (and edit it if needed).
  

```r
library(crossmeta)

# specify where data will be downloaded
data_dir <- file.path(getwd(), "data", "LY")

# gather all GSEs
gse_names  <- c("GSE9601", "GSE15069", "GSE50841", "GSE34817", "GSE29689")

# gather Illumina GSEs (see 'Checking Raw Illumina Data')
illum_names <- c("GSE50841", "GSE34817", "GSE29689")

# download raw data
# get_raw(gse_names, data_dir)
```

-----------------

### Checking Raw Illumina Data

To format raw Illumina data, I recommend that you download and set Sublime Text 
2 as your default text editor. It has very nice regular expression capabilities.
[Here](https://www.cheatography.com/davechild/cheat-sheets/regular-expressions/)
is a good regular expression cheat-sheat.

Raw illumina files will be in `data_dir` in a seperate folder for each GSE. They
are usually _.txt_ files and include _non-normalized_ in their name. Ensure the
following:

* __Detection p-values__: present (usually every second column)
* __File format__: tab seperated _.txt_ file 
* __File name__: includes _non-normalized_

Also ensure that column names have the following format:

* __Probe ID__: _ID_REF_
* __Expression values__: _AVG_Signal-sample_name_
* __Detection p-values__: _Detection-sample_name_


To open these files one at a time with your default text editor:


```r
# this is why we gathered Illumina GSEs
open_raw_illum(illum_names, data_dir)
```


For GSE50841:

* click _Find_ then _Replace_ (or _Ctrl + h_)
* _Find What_: _(SAMPLE \\d+)\\tDetection Pval_
* _Replace With_: _AVG_Signal-\\1\\tDetection-\\1_


For GSE34817:

* click _Find_ then _Replace_ (or _Ctrl + h_)
* _Find What_: _(MDA.+?)\\tDetection Pval_
* _Replace With_: _AVG_Signal-\\1\\tDetection-\\1_

For GSE29689:

* change _PROBE\_ID_ to _ID\_REF_.

A bioconductor data package (`lydata`) is available where all the above was 
performed. This data package will be used for subsequent demonstrations.


```r
library(lydata)

# location of raw data
data_dir <- system.file("extdata", package = "lydata")
```

-----------------

### Loading and Annotating Data

After downloading the raw data, it must be loaded and annotated. The necessary 
bioconductor annotation data packages will be downloaded as needed. 


```r
# reloads if previously called
esets <- load_raw(gse_names, data_dir)
```


If `crossmeta` does not know the bioconductor annotation data package for a given
platform, entry will be requested. Use the given platform (GPL) to search online
for the correct package name. For example, entry was requested for GSE29689 
(platform GPL6883). A search for _GPL6883_ on [GEO](http://www.ncbi.nlm.nih.gov/geo/) identified the title for this platform as _Illumina HumanRef-8 v3.0 expression beadchip_.
A subsequent search on [Bioconductor](https://bioconductor.org/) for 
_illuminahuman_ identified the appropriate package as _illuminaHumanv3_.

If you can't find the bioconductor annotation data package (or none exists), type
enter with a blank entry and `crossmeta` will attempt annotation using entrez 
gene ids in the feature data from the study in question. If entrez ids are absent,
annotation will fail. To proceed, I reccommend that you add entrez ids and then
use `crossmeta` to map these to hgnc symbols. As an example, if annotation had
failed for GSE15069:


```r
library(Biobase)
library(AnnotationDbi)

# check feature data to see what columns are available
head(fData(esets$GSE15069))

# if using RStudio
# View(fData(esets$GSE15069))

# annotation package for appropriate species
library(org.Mm.eg.db)

# map from accession number to entrez gene ids
acnums  <- as.character(fData(esets$GSE15069)$GB_ACC)
enids   <- mapIds(org.Mm.eg.db, acnums, "ENTREZID", "ACCNUM")

# add 'GENE_ID' column with entrez ids
fData(esets$GSE15069)$GENE_ID <- enids 

# use crossmeta to map from entrez gene ids to homologous hgnc symbol
esets$GSE15069 <- symbol_annot(esets$GSE15069)

# to overwrite saved eset (to avoid repeating above)
saveRDS(esets$GSE15069, file.path(data_dir, "GSE15069", "GSE15069_eset.rds"))
```

-----------------

### Differential Expression

After loading and annotating the data, you may proceed with differential
expression analysis:


```r
anals <- diff_expr(esets, data_dir)
```

You will be prompted to type group names and select samples for each control and
test group that you wish to compare. Once done for a study, click `Done` and you
will be prompted to do the same for the next study. You can re-select a previous
control group by selecting the group name in the drop down box. You can also 
view and delete current contrasts in the `Contrasts` tab.

If you need to look up experimental details for the study (e.g. group membership)
Click the `GEO` button in the top left corner.

**For Illumina samples**, Only titles are supplied to select groups as the
accession numbers may be incorrect. The supplied titles are the headers from 
the raw data and may differ from those on GEO. As such, you may need to look in 
the individual sample records on GEO to identify which sample belongs to which
group (not always possible).

If you need to re-run the above analysis, you can avoid having to reselect 
samples and rename groups:


```r
# load auto-saved results of previous call to diff_expr
prev <- load_diff(gse_names, data_dir)

# supply prev to diff_expr
# anals <- diff_expr(esets, data_dir, prev_anals=prev)
```

Multi-dimensional scaling plots are generated for each GSE. These plots are 
generated from expression data with the effects of surrogate variables removed.

### Non-GUI Selection

Using the graphical user interface for sample selection can be error-prone and 
time-consuming for GSEs with a large number of samples. For these cases, you may
prefer to specify group membership of samples using sample information from the
study in question. As an example:


```r
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
# anal <- diff_expr(eset, data_dir, prev_anal = sel)
```

-----------------

### Tissue Sources

Pathway and effect-size meta-analyses can be performed seperately for each tissue
source. To do so, first add tissue sources and specify any sources that should
be paired (treated as the same source for subsequent meta-analyses).


```r
# run GUI to add tissue sources
anals <- add_sources(anals, data_dir)


# for further details
?add_sources
```


### Effect Size Meta-Analysis

After tissue sources have been specified for differential expression analyses, 
overall effect sizes can be determined through meta-analysis. The meta-analysis
method determines an overall standardized effect size and false discovery rate
for each gene that is present in a specified fraction of studies (30% by default). 


```r
# re-load previous analyses if need to
anals <- load_diff(gse_names, data_dir)

# perform meta analyses by tissue source
es_res <- es_meta(anals, by_source = TRUE)

# for explanation of values
# ?es_meta
```

The meta-analysis method was adapted from `GeneMeta`. Differences include the
use of moderated unbiased effect sizes calculated by `metaMA`, the calculation 
of false discovery rates by `fdrtool`, and the allowance for genes measured in 
only a fraction of studies.

### Pathway Meta-Analysis

After tissue sources have been specified for differential expression analyses, 
Pathway analysis for each contrast and subsequent meta-analysis is possible.


```r
# pathway analysis for each contrast
path_anals <- diff_path(esets, anals, data_dir)

# pathway meta analysis by tissue source
path_res <- path_meta(path_anals, by_source = TRUE)
```

Pathway analyses are performed using [PADOG](http://bioconductor.org/packages/release/bioc/html/PADOG.html), which outperforms other methods at prioritizing expected pathways ([ref1](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0079217#pone-0079217-g002), [ref2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707541/figure/f5/)).

Pathway meta-analysis is accomplished using Fisher's p-value combination method.

### Exploring Results and Finding Drug Candidates

One application of the signatures generated by meta-analysis is to 
search for a drug, or combination of drugs that is predicted to reverse or mimic
your signature. These drugs might reverse diseases or mimic healthy lifestyles. 

The `crossmeta` function `explore_paths` interfaces with [ccmap](https://bioconductor.org/packages/release/bioc/html/ccmap.html)
to find drug candidates and graphically explore the results of the pathway meta-analysis.



```r
explore_paths(es_res, path_res)

# for details about interface
?explore_paths
```

If `crossmeta` is usefull to you, please contribute your signature! Your
contribution will be used to build a public database of microarray meta-analyses.
To contribute:


```r
# subject is the focus of the meta-analysis (e.g. drug/disease name)
contribute(anals, subject = "LY294002")

# Thank you!
```





```r
sessionInfo()
```

```
## R version 3.3.3 (2017-03-06)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 14393)
## 
## locale:
## [1] LC_COLLATE=C                          
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] Biobase_2.34.0      BiocGenerics_0.20.0 lydata_1.1.1       
## [4] crossmeta_1.1.2    
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-131               bitops_1.0-6              
##  [3] matrixStats_0.51.0         devtools_1.12.0           
##  [5] doParallel_1.0.10          RColorBrewer_1.1-2        
##  [7] httr_1.2.1                 rprojroot_1.2             
##  [9] GenomeInfoDb_1.10.3        backports_1.0.5           
## [11] tools_3.3.3                doRNG_1.6                 
## [13] R6_2.2.0                   DT_0.2                    
## [15] affyio_1.44.0              DBI_0.6                   
## [17] lazyeval_0.2.0             mgcv_1.8-17               
## [19] colorspace_1.3-2           withr_1.0.2               
## [21] bit_1.1-12                 preprocessCore_1.36.0     
## [23] fdrtool_1.2.15             xml2_1.1.1                
## [25] desc_1.1.0                 plotly_4.5.6              
## [27] pkgmaker_0.22              scales_0.4.1              
## [29] genefilter_1.56.0          affy_1.52.0               
## [31] commonmark_1.2             stringr_1.2.0             
## [33] digest_0.6.12              GEOquery_2.40.0           
## [35] XVector_0.14.1             base64enc_0.1-3           
## [37] htmltools_0.3.5            limma_3.30.13             
## [39] htmlwidgets_0.8            RSQLite_1.1-2             
## [41] BiocInstaller_1.24.0       shiny_1.0.1               
## [43] jsonlite_1.3               dplyr_0.5.0               
## [45] RCurl_1.95-4.8             magrittr_1.5              
## [47] Matrix_1.2-8               oligoClasses_1.36.0       
## [49] Rcpp_0.12.10               munsell_0.4.3             
## [51] S4Vectors_0.12.2           stringi_1.1.3             
## [53] oligo_1.38.0               SummarizedExperiment_1.4.0
## [55] zlibbioc_1.20.0            plyr_1.8.4                
## [57] grid_3.3.3                 affxparser_1.46.0         
## [59] crayon_1.3.2               miniUI_0.1.1              
## [61] lattice_0.20-34            Biostrings_2.42.1         
## [63] splines_3.3.3              annotate_1.52.1           
## [65] pander_0.6.0               knitr_1.15.1              
## [67] GenomicRanges_1.26.4       rngtools_1.2.4            
## [69] codetools_0.2-15           stats4_3.3.3              
## [71] rdrop2_0.7.0               XML_3.98-1.6              
## [73] evaluate_0.10              metap_0.8                 
## [75] data.table_1.10.4          httpuv_1.3.3              
## [77] foreach_1.4.3              testthat_1.0.2            
## [79] SMVar_1.3.3                gtable_0.2.0              
## [81] purrr_0.2.2                tidyr_0.6.1               
## [83] reshape_0.8.6              assertthat_0.1            
## [85] ggplot2_2.2.1              mime_0.5                  
## [87] xtable_1.8-2               ff_2.2-13                 
## [89] roxygen2_6.0.1             survival_2.40-1           
## [91] viridisLite_0.2.0          tibble_1.2                
## [93] metaMA_3.1.2               iterators_1.0.8           
## [95] AnnotationDbi_1.36.2       registry_0.3              
## [97] memoise_1.0.0              IRanges_2.8.2             
## [99] sva_3.22.0
```
