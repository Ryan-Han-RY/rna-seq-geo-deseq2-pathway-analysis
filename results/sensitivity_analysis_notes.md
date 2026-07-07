# limma-voom sensitivity analysis notes

An independent limma-voom sensitivity analysis was performed to assess whether the major DESeq2 findings were robust to an alternative differential expression framework.

The same filtered count matrix and metadata were used. edgeR was used to create the DGEList object and calculate TMM normalisation factors, followed by limma-voom transformation, linear modelling and empirical Bayes moderation.

## Summary

- Genes compared: 15510
- DESeq2 significant genes: 835
- limma-voom significant genes: 898
- Shared significant genes: 720
- Significant-gene Jaccard index: 0.7108
- Pearson logFC correlation: 0.9658
- Spearman logFC correlation: 0.9762

## Interpretation

Strong agreement in log fold-change direction and positive logFC correlation would support the robustness of the main disease-control transcriptional signal. Differences in exact significant-gene calls are expected because DESeq2 and limma-voom use different mean-variance modelling assumptions and test statistics.
