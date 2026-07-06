# QC and DESeq2 differential expression analysis notes

## Filtering

Low-expression genes were removed by retaining genes with counts >= 10 in at least 17 samples, corresponding to the smallest analysis group size.

The count matrix contained 60683 genes before filtering and 15510 genes after filtering.

## PCA and sample-level QC

PCA was performed on variance-stabilised expression values. PC1 explained 49.6% of variance and PC2 explained 10.9% of variance.

The PCA plot and sample distance heatmap should be inspected for group separation, outlying samples, or unexpected clustering.

## Differential expression

Using DESeq2 with design ~ group and the contrast disease vs control, 835 genes passed the threshold padj < 0.05 and absolute shrunken log2 fold change >= 1.

736 genes were upregulated and 99 genes were downregulated in disease relative to control.

## Top-ranked genes by adjusted p-value

- ENSG00000167900: log2FC = 3.112, padj = 1.77e-49
- ENSG00000171848: log2FC = 3.713, padj = 2.15e-46
- ENSG00000166851: log2FC = 3.773, padj = 2.61e-46
- ENSG00000175063: log2FC = 3.459, padj = 1.02e-41
- ENSG00000145386: log2FC = 3.641, padj = 1.21e-40
- ENSG00000111206: log2FC = 2.982, padj = 2e-40
- ENSG00000165949: log2FC = 8.551, padj = 2.19e-40
- ENSG00000088325: log2FC = 3.062, padj = 2.19e-40
- ENSG00000178999: log2FC = 3.109, padj = 2.71e-39
- ENSG00000117399: log2FC = 3.547, padj = 1.37e-38

## Next interpretation step

The differential expression table should be connected to gene annotation and pathway-level analysis, including GO, KEGG, Reactome and ranked-list enrichment.
