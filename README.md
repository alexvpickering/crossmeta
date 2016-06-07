Cross-Platform Meta-Analyis
===========================

`crossmeta` streamlines the cross-platform meta-analysis of 
microarray data. For the analysis, you will need a list of Affymetrix, Illumina,
and/or Agilent GSE numbers from [GEO](http://www.ncbi.nlm.nih.gov/geo/). All 21
species in the current [homologene](http://1.usa.gov/1TGoIy7) build are supported. 
See [vignette](http://bit.ly/1P199F9) for detailed usage.

Basic Workflow
--------------

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

Approach
--------

A high quality meta-analysis is achieved by addressing the key issues in 
conducting a meta-analysis of microarray data (1):
  
  
  
##### Uses raw data

  * Different labs process their raw data differently. These differences in
    data processing may lead to flawed conclusions upon meta-analysis. 
    
  * `crossmeta` starts with raw data, and uses a consistent processing pipeline 
    for all studies.

  
##### Maps probes to human genes

  * Different species use different symbols to reference related genes. These
    differences can make it challenging to compare similar microarray experiments
    in diverse species.
    
  * `crossmeta` maps probes to human gene symbols using homology relationships
    established by [HomoloGene](http://www.ncbi.nlm.nih.gov/homologene).
    
##### Resolves many-to-many mappings

  * Single probes can measure multiple genes. To incorporate these measurements
    appropriately, they are replaced with a new record for each gene.
    
  * Multiple probes can measure the same gene. Averaging these measurements is
    inappropriate because measurement scales will vary with probe affinity.
    `crossmeta` selects the measurement with the highest inter-quartile range as 
    it is the least likely to occur by chance.
    
    
##### Models nuisance variables

  * In addition to variables of interest, there are sources of signal due to 
    factors that are unknown, unmeasured, or too complicated to capture through
    simple models. These factors can either hide true effects or introduce 
    spurious ones.
      
  * `crossmeta` discovers and accounts for these nuissance variables using 
    surrogate variable analysis.
  
  
##### Simplifies model specification

  * Correctly specifying a model and contrast matrix can be challenging.
  
  * `crossmeta` uses an attractive user interface that allows you to simply select
    the samples you want to compare. 
    
  * Paired samples (eg. the same subject before and after treatment) can also be
    selected using the same interface.
  

##### Meta-analyzes genes with missing data

  * Differences in microarray platforms often lead to thousands of genes that are
    not measured in all studies. Most existing meta-analysis software requires 
    that these genes with missing data are discarded.

  * `crossmeta` extends the effect size meta-analysis method in `GeneMeta` to 
    allow for genes that were not measured in all studies. Keep your data, find
    more insights.


-----------------

*(1) Ramasamy A, Mondry A, Holmes CC, Altman DG (2008) Key Issues in Conducting a*
*Meta-Analysis of Gene Expression Microarray Datasets. PLoS Med 5(9): e184.* doi:10.1371/journal.pmed.0050184
