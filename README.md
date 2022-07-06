# LTR7Y Project Analysis

Scripts for bioinfomatics analysis of LTR7Y project.

Summary: 
##### 01preprocess 

| File                   | Description                                                                     |
|------------------------|---------------------------------------------------------------------------------|
| Preprocess_ATACseq.sh  | Preprocess of ATAC-seq data, including QC, trimming, alignment, peak calling    |
| Preprocess_CutTag.sh   | Preprocess of Cut&Tag data, including QC, trimming, alignment, peak calling     |
| Preprocess_ChIPseq.sh  | Preprocess of ChIP-seq data, including QC, trimming, alignment, peak calling    |
| Preprocess_RNAseq.sh   | Preprocess of RNA-seq data, including QC, trimming, alignment, quantification   |
| Preprocess_scRNAseq.sh | Preprocess of scRNA-seq data, including QC, trimming, alignment, quantification |

##### 02downstreamanalysis

| File                        | Destription                                                  |
|-----------------------------|--------------------------------------------------------------|
| DEG.R                       | Differential expressed genes                                 |
| DETE.R                      | Differential expressed TEs                                   |
| Barplot_DETEpct.R           | Barplot of top DETE percentage                               |
| Heatmap_expression.R        | Heatmap of gene/TE expression pattern                        |
| PCA.R                       | PCA plot of KLF5 KOs                                         |
| scRNAseq_analysis.R         | naive hESC scRNA-seq downstream analysis                     |
| Heatmap_signalenrichment.sh | Heatmap of TF/histone/motif signals over TEs or peak regions |

