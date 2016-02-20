# caloric-restriction (CR)

###Description:
These scripts are used to streamline the meta-analysis of caloric-restriction (CR) microarray data. The analysis requires are a list of Affymetrix or Illumina GSE numbers (from GEO). Although developed for meta-analysis of CR microarray data, can be used for non-CR microarray data.

#### Features:
1. Download/save raw data from GEO ('Supplementary file')

2. Consistently normalize raw data:
  * __Affymetrix__: RMA (includes norm-exp bg correct, inter array quantile normalization, and log2 transformation)
  * __Illumina__: neqc (includes norm-exp bg correct, inter array quantile normalization, and log2 transformation)
  
3. Differential expression analysis:
  * control for known nuissance variables (e.g. batch, scan date)
    * scan date detected automatically for Affymetrix arrays
  * control for unknown surogate variabled (sva)
  * select desired contrasts (CR=test, AL=control)
  * returns list of probes with adjusted P-values > 0.05 and fold-change > 2.
  * analysis saved (with expression data from contrasted samples, model matrix, contrast matrix, and differentially expressed probes)
