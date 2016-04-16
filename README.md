#cross-meta

###Description:
These scripts are used to streamline the cross-platform meta-analysis of 
microarray data. The analysis requires a list of Affymetrix, Illumina, and/or
Agilent GSE numbers (from GEO). Only supports mouse and human.

####Features:
1. Download/save raw data from GEO ('Supplementary file')

2. Consistently normalize raw data:
  * __Affymetrix__: RMA (includes norm-exp bg correct, inter array quantile 
  normalization, and log2 transformation)
  * __Illumina/Agilent__: neqc (includes norm-exp bg correct, inter array 
  quantile normalization, and log2 transformation)
  
3. Annotate using Bioconductor annotation data:
  * uses SYMBOL (NAs removed)
  * one-to-many: expands to keep all
  * many-to-one: keeps highest interquartile range of samples selected in (4)
  
4. Differential expression analysis:
  * control for unknown surogate variables (SVA)
  * select desired contrasts (interactive)
  * multidimensional scaling plot (SVA-adjusted)

5. Meta-Analysis of results:
  * format results to work with MAMA package
    * pvalcombination & EScombination edited:
        * use results from (4) 
          * p/t-values reflect SVA/multi-contrast analysis
    * Incorporate RankMerging from GeneExpressionSignature
        
