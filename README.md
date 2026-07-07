# Reproducible RNA-seq Differential Expression and Pathway Analysis of COVID-19-associated PBMC Transcriptomic Signatures

## Overview

This repository contains a reproducible bulk RNA-seq analysis workflow using public GEO count data from **GSE152418**. The project investigates transcriptional differences between COVID-19-related peripheral blood mononuclear cell (PBMC) samples and healthy control PBMC samples.

The analysis starts from a public gene-level count matrix rather than FASTQ files. The focus is on practical bioinformatics workflow execution: metadata curation, count-matrix alignment, quality control, differential expression modelling, gene annotation, pathway-level interpretation, sensitivity analysis, and reproducible reporting.

## Research question

Which genes and biological pathways are transcriptionally altered between COVID-19-related PBMC samples and healthy control PBMC samples?

## Dataset

| Item | Description |
|---|---|
| GEO accession | GSE152418 |
| Organism | Homo sapiens |
| Sample type | Peripheral blood mononuclear cells, PBMCs |
| Data type | Bulk RNA-seq gene count matrix |
| Comparison | COVID-19-related PBMC samples vs healthy control PBMC samples |
| Samples retained | 17 disease samples and 17 healthy controls |
| Initial count matrix after metadata alignment | 60,683 genes × 34 samples |
| Filtered count matrix | 15,510 genes × 34 samples |
| Primary model | DESeq2, design `~ group` |
| Primary contrast | Disease vs control |
| Sensitivity analysis | edgeR + limma-voom |

## Repository structure

```text
.
├── data/
│   ├── raw/
│   └── processed/
│       ├── metadata_clean.tsv
│       ├── counts_clean.tsv
│       ├── counts_clean_matrix.rds
│       ├── counts_filtered.tsv
│       ├── counts_filtered_matrix.rds
│       └── vst_transformed_counts.tsv
├── metadata/
│   ├── dataset_decision.tsv
│   └── metadata_raw_from_geo.tsv
├── report/
│   ├── rna_seq_geo_analysis.qmd
│   └── rna_seq_geo_analysis.html
├── results/
│   ├── figures/
│   ├── objects/
│   ├── tables/
│   ├── analysis_notes.md
│   ├── pathway_interpretation_notes.md
│   ├── sensitivity_analysis_notes.md
│   └── session_info.txt
├── scripts/
│   ├── 00_setup.R
│   ├── 01_download_data.R
│   ├── 02_metadata_cleaning.R
│   ├── 03_qc_deseq2_analysis.R
│   ├── 04_gene_annotation_pathway_analysis.R
│   ├── 05_sensitivity_limma_voom_analysis.R
│   └── 06_render_report.R
├── renv.lock
├── README.md
└── project1.Rproj
```

## Workflow

The analysis was performed as a scripted, reproducible workflow:

1. Set up project structure and package environment using `renv`
2. Download GEO metadata and public count data
3. Curate sample metadata into a clean disease-control comparison
4. Align count matrix columns with cleaned metadata sample IDs
5. Filter low-expression genes
6. Run DESeq2 differential expression analysis
7. Generate PCA, sample distance heatmap, MA plot, volcano plot and top-gene heatmap
8. Map Ensembl gene IDs to gene symbols, Entrez IDs and gene names
9. Perform GO Biological Process, Reactome and KEGG enrichment analysis
10. Perform ranked-list GO enrichment analysis
11. Perform independent limma-voom sensitivity analysis
12. Generate a Quarto HTML report

## Main results

### Metadata and filtering

The cleaned dataset retained 34 samples: 17 healthy controls and 17 COVID-19-related PBMC samples.

Low-expression genes were removed by retaining genes with counts `>=10` in at least 17 samples, corresponding to the smallest analysis group size.

| Metric | Value |
|---|---:|
| Genes before filtering | 60,683 |
| Genes after filtering | 15,510 |
| Genes removed | 45,173 |
| Smallest group size | 17 |

### DESeq2 differential expression

Differential expression was modelled using DESeq2 with design `~ group`, comparing disease samples against healthy controls. Log2 fold changes were shrunk using `apeglm`.

| Metric | Value |
|---|---:|
| Tested genes | 15,510 |
| Genes with adjusted p-value | 15,510 |
| Significant genes, padj < 0.05 and abs(shrunken log2FC) >= 1 | 835 |
| Upregulated genes in disease | 736 |
| Downregulated genes in disease | 99 |

The disease-control contrast showed a strong transcriptional signal, with most significant genes upregulated in disease samples.

### Gene-level findings

After annotation from Ensembl IDs to gene symbols, the strongest upregulated genes included:

| Gene | Description | Shrunken log2FC | Adjusted p-value |
|---|---|---:|---:|
| TK1 | thymidine kinase 1 | 3.11 | 1.77e-49 |
| RRM2 | ribonucleotide reductase regulatory subunit M2 | 3.71 | 2.15e-46 |
| PLK1 | polo like kinase 1 | 3.77 | 2.61e-46 |
| UBE2C | ubiquitin conjugating enzyme E2 C | 3.46 | 1.02e-41 |
| CCNA2 | cyclin A2 | 3.64 | 1.21e-40 |
| FOXM1 | forkhead box M1 | 2.98 | 2.00e-40 |
| IFI27 | interferon alpha inducible protein 27 | 8.55 | 2.19e-40 |
| TPX2 | TPX2 microtubule nucleation factor | 3.06 | 2.19e-40 |
| AURKB | aurora kinase B | 3.11 | 2.71e-39 |
| CDC20 | cell division cycle 20 | 3.55 | 1.37e-38 |

The leading gene-level signal was dominated by cell-cycle and proliferative markers, together with immune-associated genes such as `IFI27`, `IGHV1-24` and `MZB1`.

### Pathway-level findings

Upregulated genes were strongly enriched for two major biological programmes:

1. **Humoral/B-cell immune response**
2. **Cell-cycle and proliferation**

Top GO Biological Process enrichment results among upregulated genes included:

| Pathway | Adjusted p-value | Gene ratio | Gene count |
|---|---:|---:|---:|
| Immunoglobulin mediated immune response | 3.79e-31 | 59/635 | 59 |
| B cell mediated immunity | 6.50e-31 | 59/635 | 59 |
| Immunoglobulin production | 3.01e-26 | 52/635 | 52 |
| Production of molecular mediator of immune response | 4.43e-20 | 58/635 | 58 |
| Adaptive immune response based on somatic recombination of immune receptors | 2.34e-19 | 65/635 | 65 |
| Lymphocyte mediated immunity | 6.67e-19 | 63/635 | 63 |
| Leukocyte mediated immunity | 1.62e-15 | 66/635 | 66 |
| Nuclear division | 1.01e-14 | 61/635 | 61 |

Reactome and KEGG independently supported the cell-cycle signal:

| Source | Pathway | Adjusted p-value | Gene count |
|---|---|---:|---:|
| Reactome | Cell Cycle Checkpoints | 5.88e-17 | 56 |
| KEGG | Cell cycle | 1.26e-13 | 38 |
| Reactome | DNA Replication | 7.24e-10 | 32 |
| Reactome | M Phase | 6.58e-10 | 53 |

Downregulated genes showed enrichment for chemotaxis and leukocyte migration-related processes, including G protein-coupled receptor signalling pathway, cell chemotaxis, leukocyte chemotaxis and regulation of leukocyte migration.

### Sensitivity analysis

An independent edgeR + limma-voom analysis was performed using the same filtered count matrix and metadata. The sensitivity analysis strongly supported the DESeq2 findings.

| Metric | Value |
|---|---:|
| Genes compared | 15,510 |
| DESeq2 significant genes | 835 |
| limma-voom significant genes | 898 |
| Shared significant genes | 720 |
| DESeq2-only significant genes | 115 |
| limma-voom-only significant genes | 178 |
| Significant-gene Jaccard index | 0.711 |
| Same logFC direction across all tested genes | 0.959 |
| Same logFC direction among shared significant genes | 1.000 |
| Pearson logFC correlation | 0.966 |
| Spearman logFC correlation | 0.976 |

The high log fold-change correlation and complete directional agreement among shared significant genes indicate that the main disease-control transcriptional signal is robust to an alternative differential expression framework.

## Biological interpretation

COVID-19-related PBMC samples showed a robust disease-associated transcriptional signature. The main signal was characterised by three components:

1. **Humoral/B-cell immune response**

   GO enrichment was dominated by immunoglobulin-mediated immune response, B-cell-mediated immunity and immunoglobulin production. This suggests antibody-producing or plasmablast-associated immune activity in disease PBMC samples.

2. **Cell-cycle and proliferation programme**

   Top upregulated genes included `TK1`, `RRM2`, `PLK1`, `UBE2C`, `CCNA2`, `FOXM1`, `AURKB`, `CDC20`, `MKI67` and `BIRC5`. Reactome and KEGG independently highlighted cell-cycle-related pathways, supporting a proliferative or activated immune-cell programme in disease samples.

3. **Interferon-associated antiviral marker signal**

   `IFI27` was one of the strongest upregulated genes, supporting an interferon-associated antiviral component. However, pathway-level results were dominated more strongly by humoral immunity and cell-cycle activation than by interferon response alone.

Because this is bulk PBMC RNA-seq, these findings may reflect both transcriptional activation within immune cells and differences in immune-cell composition. Cell-type-specific conclusions would require single-cell RNA-seq, cytometry, or computational cell-type deconvolution.

## Main figures

| Figure | Description |
|---|---|
| `results/figures/pca_vst_group.png` | PCA of variance-stabilised expression profiles |
| `results/figures/sample_distance_heatmap.png` | Sample-to-sample distance heatmap |
| `results/figures/volcano_plot_deseq2_gene_symbols.png` | DESeq2 volcano plot with gene symbols |
| `results/figures/top50_gene_heatmap_gene_symbols.png` | Top 50 differentially expressed genes |
| `results/figures/go_bp_top_terms_upregulated_barplot.png` | Top GO Biological Process terms among upregulated genes |
| `results/figures/reactome_dotplot_upregulated.png` | Reactome enrichment among upregulated genes |
| `results/figures/deseq2_limma_voom_logfc_comparison.png` | DESeq2 versus limma-voom logFC comparison |
| `results/figures/significant_gene_overlap_deseq2_limma_voom.png` | Significant-gene overlap between DESeq2 and limma-voom |

## Main result tables

| Table | Description |
|---|---|
| `results/tables/deseq2_all_genes_annotated.tsv` | Full annotated DESeq2 result table |
| `results/tables/deseq2_significant_genes_annotated.tsv` | Significant annotated DESeq2 genes |
| `results/tables/top50_annotated_genes.tsv` | Top 50 annotated genes by adjusted p-value |
| `results/tables/pathway_enrichment_summary.tsv` | Combined pathway enrichment summary |
| `results/tables/biological_interpretation_table.tsv` | Gene/pathway interpretation table |
| `results/tables/limma_voom_all_genes.tsv` | Full limma-voom result table |
| `results/tables/sensitivity_overlap_summary.tsv` | DESeq2 versus limma-voom robustness summary |
| `results/tables/shared_significant_genes_deseq2_limma_voom.tsv` | Shared significant genes between methods |

## Reproducibility

This project uses `renv` to record the R package environment.

To restore the environment:

```r
renv::restore()
```

To reproduce the analysis from the project root:

```r
source("scripts/00_setup.R")
source("scripts/01_download_data.R")
source("scripts/02_metadata_cleaning.R")
source("scripts/03_qc_deseq2_analysis.R")
source("scripts/04_gene_annotation_pathway_analysis.R")
source("scripts/05_sensitivity_limma_voom_analysis.R")
source("scripts/06_render_report.R")
```

The rendered report is available at:

```text
report/rna_seq_geo_analysis.html
```

## Limitations

- The analysis uses bulk PBMC RNA-seq data, so observed gene expression differences may reflect both transcriptional regulation and immune-cell composition differences.
- The primary model uses a simple disease-control design, prioritising a clean and reproducible comparison.
- Public GEO metadata can be incomplete or inconsistently annotated, so metadata curation decisions were made explicitly and recorded.
- Pathway enrichment results are hypothesis-generating and should not be interpreted as causal mechanisms without additional validation.
- Downregulated genes were fewer than upregulated genes, so downregulated pathway results should be interpreted more cautiously.

## Software

Main R packages used:

- `GEOquery`
- `tidyverse`
- `DESeq2`
- `apeglm`
- `edgeR`
- `limma`
- `AnnotationDbi`
- `org.Hs.eg.db`
- `clusterProfiler`
- `ReactomePA`
- `enrichplot`
- `pheatmap`
- `ggrepel`
- `quarto`
- `renv`

Package versions are recorded in `renv.lock` and `results/session_info.txt`.

## Project summary

This project demonstrates a reproducible public RNA-seq workflow from GEO count data to biological interpretation. The main COVID-19-associated PBMC signature was dominated by humoral/B-cell immune response and cell-cycle/proliferation programmes, with selected antiviral marker genes such as `IFI27`. An independent limma-voom sensitivity analysis supported the robustness of the DESeq2 findings.