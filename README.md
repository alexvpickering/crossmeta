#crossmeta

`crossmeta` streamlines the cross-platform meta-analysis of 
microarray data. For the analysis, you will need a list of Affymetrix, Illumina,
and/or Agilent GSE numbers from [GEO](http://www.ncbi.nlm.nih.gov/geo/). Mouse 
and human inter-species analyses are supported.


-----------------


A high quality meta-analysis is achieved by addressing the key issues in 
conducting a meta-analysis of microarray data (1):
  
  
  
**Consistently normalize raw data**   

  * Raw data is downloaded from GEO supplementary files.
  * Raw data is norm-exp background corrected, quantile normalized, and log2
    transformed.
  
**Annotate probes** 

  * Uses gene symbols from most recent bioconductor annotation data packages.
  * One-to-many (probe-to-symbol): keeps all.
  * Many-to-one (probe-to-symbol): keeps highest interquartile range among
    selected samples.
  
  
**Differential expression analysis**  

  * Controls for unknown batch effects (`sva`).
  * User selects desired contrasts (interactive).
  * Multidimensional scaling plots.


**Meta-analysis of results**  

  * Results formated to work with `MAMA` package.
  
  
-----------------


*(1) Ramasamy A, Mondry A, Holmes CC, Altman DG (2008) Key Issues in Conducting a*
*Meta-Analysis of Gene Expression Microarray Datasets. PLoS Med 5(9): e184.* doi:10.1371/journal.pmed.0050184
