##Cross-Platform Meta-Analyis

`crossmeta` streamlines the cross-platform meta-analysis of 
microarray data. For the analysis, you will need a list of Affymetrix, Illumina,
and/or Agilent GSE numbers from [GEO](http://www.ncbi.nlm.nih.gov/geo/). All 21
species in the current [homologene](http://1.usa.gov/1TGoIy7) build are supported. 
See [vignette](http://bit.ly/1P199F9) for usage.


-----------------


A high quality meta-analysis is achieved by addressing the key issues in 
conducting a meta-analysis of microarray data (1):
  
  
  
**Consistently normalize raw data**   

  * Raw data is downloaded from GEO supplementary files.
  * Raw data is norm-exp background corrected, quantile normalized, and log2
    transformed.
  
**Annotate probes** 

  * Probes are mapped to human gene symbol (or homolog).
  * One-to-many (probe-to-symbol): keeps all.
  * Many-to-one (probe-to-symbol): keeps highest interquartile range among
    selected samples.
  
  
**Differential expression analysis**  

  * Controls for unknown batch effects (`sva`).
  * User selects contrasts and any paired samples (`shiny` GUI).
  * Plots group clustering (multidimensional scaling plots).

**Meta-analysis of results**  

  * Extends method in `GeneMeta` to allow for genes that were not measured in 
    all studies.
  * Analysis uses moderated unbiased effect sizes calculated by `metaMA` and
    determines false discovery rates using `fdrtool`.

-------------------------------

In general, the analysis workflow is as follows (see 
[vignette](http://bit.ly/1P199F9) for details):

```R
# studies from GEO
gse_names  <- c("GSE9601", "GSE15069")

# get raw data for specified studies
get_raw(gse_names)

# load and annotate raw data
esets <- load_raw(gse_names)

# perform differential expression analysis
anals <- diff_expr(esets)

# perform meta-analysis
es <- es_meta(anals)

# contribute your results
contribute(anals)
```
  
-----------------


*(1) Ramasamy A, Mondry A, Holmes CC, Altman DG (2008) Key Issues in Conducting a*
*Meta-Analysis of Gene Expression Microarray Datasets. PLoS Med 5(9): e184.* doi:10.1371/journal.pmed.0050184
