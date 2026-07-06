# Reproducible RNA-seq Differential Expression and Pathway Analysis of COVID-19-associated PBMC Transcriptomic Signatures

## Project overview

This repository contains a reproducible bulk RNA-seq analysis workflow using public GEO count data from **GSE152418**. The project compares COVID-19-related peripheral blood mononuclear cell (PBMC) samples against healthy control PBMC samples to identify disease-associated gene-level and pathway-level transcriptomic signatures.

The workflow starts from a public raw count matrix and cleaned sample metadata, rather than FASTQ files. It focuses on metadata curation, count-matrix quality control, DESeq2 differential expression modelling, gene annotation, pathway enrichment, biological interpretation, and reproducible reporting.

## Research question

Which genes and biological pathways are transcriptionally altered between COVID-19-related PBMC samples and healthy control PBMC samples in public bulk RNA-seq data?

## Dataset

- **GEO accession:** GSE152418
- **Organism:** Homo sapiens
- **Sample type:** PBMC
- **Data type:** Bulk RNA-seq raw gene count matrix
- **Comparison:** COVID-19-related PBMC samples versus healthy control PBMC samples
- **Samples retained after metadata curation:** 17 control samples and 17 disease samples
- **Initial count matrix:** 60,683 genes × 34 samples
- **Filtered count matrix:** 15,510 genes × 34 samples

The dataset was selected because it provides a manageable public RNA-seq count matrix with a clear disease-control design, making it suitable for a focused, reproducible differential expression and pathway interpretation workflow.

## Repository structure

```text
.
├── data/
│   ├── raw/
│   │   ├── counts_raw.tsv
│   │   ├── counts_raw_matrix.rds
│   │   └── gse152418_eset.rds
│   └── processed/
│       ├── metadata_clean.tsv
│       ├── counts_clean.tsv
│       ├── counts_clean_matrix.rds
│       ├── counts_filtered.tsv
│       ├── counts_filtered_matrix.rds
│       └── vst_transformed_counts.tsv
├── metadata/
│   ├── metadata_raw_from_geo.tsv
│   └── dataset_decision.tsv
├── results/
│   ├── figures/
│   ├── objects/
│   ├── tables/
│   ├── analysis_notes.md
│   └── pathway_interpretation_notes.md
├── scripts/
│   ├── 00_setup.R
│   ├── 01_download_data.R
│   ├── 02_metadata_cleaning.R
│   ├── 03_qc_deseq2_analysis.R
│   └── 04_gene_annotation_pathway_analysis.R
├── report/
├── renv.lock
└── README.md
