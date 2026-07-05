# Reproducible RNA-seq Differential Expression and Pathway Analysis of COVID-19-associated PBMC Transcriptomic Signatures

This project implements a reproducible bulk RNA-seq differential expression and pathway analysis workflow using public GEO count data.

## Dataset

- GEO accession: GSE152418
- Organism: Homo sapiens
- Sample type: PBMC
- Data type: RNA-seq raw gene count matrix
- Comparison: COVID-19-related PBMC samples versus healthy control PBMC samples
- Samples retained after metadata curation: 17 control samples and 17 disease samples

## Research question

Which genes and biological pathways are transcriptionally altered between COVID-19-related PBMC samples and healthy control PBMC samples in public GEO RNA-seq data?

## Main methods

- GEO metadata curation
- Raw count matrix inspection
- Count-metadata sample alignment
- DESeq2 differential expression analysis
- edgeR/limma-voom sensitivity analysis
- PCA and sample distance analysis
- Volcano plot and clustered heatmap
- GO, KEGG, Reactome and GSEA-style enrichment analysis
- Quarto-based reproducible report

## Current project status

Preprocessing completed:

- Project structure created
- GEO metadata downloaded
- Raw count matrix downloaded
- Metadata cleaned
- Count matrix aligned to metadata
- Dataset decision table created

## Reproducibility

The project uses `renv` to record the R package environment.
