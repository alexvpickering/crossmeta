# caloric-restriction (CR)

###Description:
These scripts are used to streamline the meta-analysis of microarray data. The analysis requires a list of Affymetrix and/or Illumina GSE numbers (from GEO). Although developed for meta-analysis of CR microarray data, should be appropriate for other microarray data.

####Features:
1. Download/save raw data from GEO ('Supplementary file')

2. Consistently normalize raw data:
  * __Affymetrix__: RMA (includes norm-exp bg correct, inter array quantile normalization, and log2 transformation)
  * __Illumina__: neqc (includes norm-exp bg correct, inter array quantile normalization, and log2 transformation)
  
3. Annotate using Bioconductor annotation data:
  * uses SYMBOL (NAs removed)
  * one-to-many: expands to keep all
  * many-to-one: keeps highest interquartile range samples selected in (4)
  
4. Differential expression analysis:
  * control for unknown surogate variables (SVA)
  * select desired contrasts (CR=test, AL=control)
  * multidimensional scaling plot (SVA-adjusted)

5. Meta-Analysis of results:
  * format results to work with MAMA package
    * pvalcombination & EScombination edited:
        * use results from (4) 
          * p/t-values reflect SVA/multi-contrast analysis
    * VennDiagram for 2, 3, or 4 methods
    * Incorporate RankMerging from GeneExpressionSignature
        