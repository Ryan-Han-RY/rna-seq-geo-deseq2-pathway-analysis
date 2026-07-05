# scripts/00_setup.R
# Purpose: create project folders, install required packages, initialise renv,
# and record the computational environment.

rm(list = ls())

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  stringsAsFactors = FALSE
)

message("Starting project setup...")

# -----------------------------
# 1. Install renv first
# -----------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Give renv explicit permission to create project files.
# This avoids the first-time interactive y/n prompt.
options(renv.consent = TRUE)
renv::consent(provided = TRUE)

# Initialise renv for this project.
# bare = TRUE means: do not try to capture random packages from your old R setup.
if (!file.exists("renv.lock") && !dir.exists("renv")) {
  renv::init(bare = TRUE, restart = FALSE)
}
# -----------------------------
# 2. Helper: install missing packages
# -----------------------------
install_if_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, function(pkg) requireNamespace(pkg, quietly = TRUE), logical(1))]
  if (length(missing) > 0) {
    install.packages(missing)
  }
}

cran_pkgs <- c(
  "tidyverse",
  "here",
  "janitor",
  "fs",
  "glue",
  "knitr",
  "rmarkdown",
  "quarto",
  "pheatmap",
  "ggrepel",
  "patchwork"
)

install_if_missing_cran(cran_pkgs)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "GEOquery",
  "Biobase",
  "DESeq2",
  "edgeR",
  "limma",
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "EnhancedVolcano",
  "fgsea"
)

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, function(pkg) requireNamespace(pkg, quietly = TRUE), logical(1))]

if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# -----------------------------
# 3. Load packages needed today
# -----------------------------
library(fs)
library(glue)

# -----------------------------
# 4. Create project folders
# -----------------------------
dirs <- c(
  "data/raw",
  "data/processed",
  "metadata",
  "results/figures",
  "results/tables",
  "scripts",
  "report"
)

fs::dir_create(dirs)

# -----------------------------
# 5. Create .gitignore
# -----------------------------
gitignore_text <- c(
  ".Rhistory",
  ".RData",
  ".Rproj.user/",
  ".DS_Store",
  "renv/library/",
  "renv/staging/",
  "renv/python/",
  "results/*.html",
  "results/*.pdf"
)

writeLines(gitignore_text, ".gitignore")

# -----------------------------
# 6. Create first README
# -----------------------------
readme_text <- c(
  "# Reproducible RNA-seq Differential Expression and Pathway Analysis of COVID-19-associated PBMC Transcriptomic Signatures",
  "",
  "This project implements a reproducible bulk RNA-seq differential expression and pathway analysis workflow using public GEO count data.",
  "",
  "## Dataset",
  "",
  "- GEO accession: GSE152418",
  "- Organism: Homo sapiens",
  "- Sample type: PBMC",
  "- Data type: RNA-seq raw gene count matrix",
  "- Comparison: COVID-19-related samples versus healthy controls",
  "",
  "## Research question",
  "",
  "Which genes and biological pathways are transcriptionally altered between COVID-19-related PBMC samples and healthy control PBMC samples in public GEO RNA-seq data?",
  "",
  "## Main methods",
  "",
  "- GEO metadata curation",
  "- Raw count matrix inspection",
  "- DESeq2 differential expression analysis",
  "- edgeR/limma-voom sensitivity analysis",
  "- PCA, sample distance analysis, volcano plot and heatmap",
  "- GO, KEGG, Reactome and GSEA-style enrichment analysis",
  "- Quarto-based reproducible report",
  "",
  "## Reproducibility",
  "",
  "The project uses `renv` to record the R package environment."
)

writeLines(readme_text, "README.md")

# -----------------------------
# 7. Record session info
# -----------------------------
sink("results/session_info.txt")
sessionInfo()
sink()

# -----------------------------
# 8. Snapshot package environment
# -----------------------------
renv::snapshot(prompt = FALSE)

message("Project setup complete.")
