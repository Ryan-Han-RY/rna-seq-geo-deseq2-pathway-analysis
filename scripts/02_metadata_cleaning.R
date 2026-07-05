# scripts/02_metadata_cleaning.R
# Purpose: clean GSE152418 metadata, define COVID/control groups,
# align count matrix with metadata, and save clean Initial preprocessing outputs.

rm(list = ls())

library(tidyverse)
library(here)
library(janitor)
library(fs)
library(glue)

here::i_am("scripts/02_metadata_cleaning.R")

metadata_raw_path <- here("metadata", "metadata_raw_from_geo.tsv")
counts_raw_path <- here("data", "raw", "counts_raw_matrix.rds")

if (!file.exists(metadata_raw_path)) {
  stop("Missing metadata_raw_from_geo.tsv. Run scripts/01_download_data.R first.")
}

if (!file.exists(counts_raw_path)) {
  stop("Missing counts_raw_matrix.rds. Run scripts/01_download_data.R first.")
}

metadata_raw <- readr::read_tsv(
  metadata_raw_path,
  show_col_types = FALSE
)

counts_raw <- readRDS(counts_raw_path)

message("Raw count matrix dimension:")
print(dim(counts_raw))

if (ncol(counts_raw) < 30) {
  stop("Raw count matrix has fewer than 30 samples. Your 01_download_data.R output is wrong. Re-run/fix download script first.")
}

# -----------------------------
# 1. Basic metadata inspection
# -----------------------------
message("Metadata columns:")
print(names(metadata_raw))

message("First sample titles:")
print(head(metadata_raw$title, 20))

message("First count matrix sample names:")
print(head(colnames(counts_raw), 20))

# -----------------------------
# 2. Parse characteristics columns where available
# -----------------------------
char_cols <- names(metadata_raw)[stringr::str_detect(names(metadata_raw), "^characteristics_ch1")]

if (length(char_cols) > 0) {
  char_long <- metadata_raw %>%
    dplyr::select(geo_accession, all_of(char_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(char_cols),
      names_to = "characteristics_column",
      values_to = "characteristics_value"
    ) %>%
    dplyr::filter(!is.na(characteristics_value), characteristics_value != "") %>%
    dplyr::mutate(
      key = stringr::str_match(characteristics_value, "^([^:]+):\\s*(.*)$")[, 2],
      value = stringr::str_match(characteristics_value, "^([^:]+):\\s*(.*)$")[, 3],
      key = janitor::make_clean_names(key)
    ) %>%
    dplyr::filter(!is.na(key), !is.na(value))
  
  char_wide <- char_long %>%
    dplyr::select(geo_accession, key, value) %>%
    tidyr::pivot_wider(
      names_from = key,
      values_from = value,
      values_fn = function(x) paste(unique(x), collapse = "; ")
    )
} else {
  char_wide <- tibble::tibble(geo_accession = metadata_raw$geo_accession)
}

metadata_parsed <- metadata_raw %>%
  dplyr::transmute(
    geo_accession = geo_accession,
    sample_title = title,
    source_name = source_name_ch1,
    organism = organism_ch1
  ) %>%
  dplyr::left_join(char_wide, by = "geo_accession")

take_col <- function(df, candidates, default = "not_available") {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    return(rep(default, nrow(df)))
  }
  x <- as.character(df[[hit[1]]])
  x[is.na(x) | x == ""] <- default
  x
}

# -----------------------------
# 3. Define sample_id
# -----------------------------
count_cols <- colnames(counts_raw)

metadata_clean <- metadata_parsed %>%
  dplyr::mutate(
    sample_id = dplyr::case_when(
      sample_title %in% count_cols ~ sample_title,
      geo_accession %in% count_cols ~ geo_accession,
      TRUE ~ NA_character_
    )
  )

if (sum(!is.na(metadata_clean$sample_id)) < 30) {
  message("Matched sample IDs:")
  print(metadata_clean %>% dplyr::select(geo_accession, sample_title, sample_id))
  
  stop("Fewer than 30 metadata samples matched count matrix columns. Need to inspect count column names.")
}

# -----------------------------
# 4. Define group using GSE152418 sample title pattern
# -----------------------------
# In GSE152418, COVID samples have nCoV/nCOV in sample titles.
# Healthy controls are the remaining PBMC samples.
metadata_clean <- metadata_clean %>%
  dplyr::mutate(
    group = dplyr::case_when(
      stringr::str_detect(sample_title, regex("ncov", ignore_case = TRUE)) ~ "disease",
      TRUE ~ "control"
    ),
    group = factor(group, levels = c("control", "disease")),
    
    disease_state = dplyr::case_when(
      group == "disease" ~ "COVID-19-related",
      group == "control" ~ "healthy_control",
      TRUE ~ "not_available"
    ),
    
    severity = take_col(metadata_parsed, c("severity", "disease_severity", "clinical_status")),
    sex = take_col(metadata_parsed, c("gender", "sex")),
    age = take_col(metadata_parsed, c("age", "age_years")),
    batch = take_col(metadata_parsed, c("geographical_location", "batch", "sequencing_batch")),
    cell_type = take_col(metadata_parsed, c("cell_type")),
    days_post_symptom_onset = take_col(metadata_parsed, c("days_post_symptom_onset"))
  ) %>%
  dplyr::filter(!is.na(sample_id), !is.na(group)) %>%
  dplyr::arrange(group, sample_id) %>%
  dplyr::select(
    sample_id,
    geo_accession,
    sample_title,
    group,
    disease_state,
    severity,
    sex,
    age,
    batch,
    cell_type,
    source_name,
    organism,
    days_post_symptom_onset
  )

# -----------------------------
# 5. Hard checks
# -----------------------------
group_tab <- table(metadata_clean$group)
print(group_tab)

if (nrow(metadata_clean) != 34) {
  stop(glue("Expected 34 samples, but retained {nrow(metadata_clean)} samples. Do not continue."))
}

if (!all(c("control", "disease") %in% names(group_tab))) {
  stop("Both control and disease groups must exist.")
}

if (as.integer(group_tab["control"]) != 17 || as.integer(group_tab["disease"]) != 17) {
  stop("Expected 17 control and 17 disease samples. Metadata grouping is wrong.")
}

if (!all(metadata_clean$sample_id %in% colnames(counts_raw))) {
  missing_samples <- setdiff(metadata_clean$sample_id, colnames(counts_raw))
  print(missing_samples)
  stop("Some metadata sample_id values are not present in count matrix columns.")
}

counts_clean <- counts_raw[, metadata_clean$sample_id, drop = FALSE]

stopifnot(identical(colnames(counts_clean), metadata_clean$sample_id))
stopifnot(!anyNA(counts_clean))
stopifnot(all(counts_clean >= 0))

# -----------------------------
# 6. Save clean files
# -----------------------------
fs::dir_create(c(
  here("data", "processed"),
  here("metadata"),
  here("results", "tables"),
  here("results", "figures")
))

readr::write_tsv(
  metadata_clean,
  here("data", "processed", "metadata_clean.tsv")
)

saveRDS(
  counts_clean,
  here("data", "processed", "counts_clean_matrix.rds")
)

counts_clean_out <- as.data.frame(counts_clean, check.names = FALSE) %>%
  tibble::rownames_to_column("gene_id")

readr::write_tsv(
  counts_clean_out,
  here("data", "processed", "counts_clean.tsv")
)

# -----------------------------
# 7. Save Initial preprocessing summary tables
# -----------------------------
group_summary <- metadata_clean %>%
  dplyr::count(group, name = "n_samples")

metadata_summary <- metadata_clean %>%
  dplyr::count(group, disease_state, severity, sex, batch, name = "n_samples") %>%
  dplyr::arrange(group, disease_state, severity, sex)

library_size_table <- tibble::tibble(
  sample_id = colnames(counts_clean),
  library_size = colSums(counts_clean)
) %>%
  dplyr::left_join(
    metadata_clean %>% dplyr::select(sample_id, group, disease_state, severity, sex),
    by = "sample_id"
  ) %>%
  dplyr::arrange(group, sample_id)

matrix_checks <- tibble::tibble(
  metric = c(
    "n_genes",
    "n_samples",
    "n_control_samples",
    "n_disease_samples",
    "counts_are_numeric",
    "no_missing_values",
    "metadata_count_order_identical",
    "min_library_size",
    "median_library_size",
    "max_library_size"
  ),
  value = c(
    as.character(nrow(counts_clean)),
    as.character(ncol(counts_clean)),
    as.character(sum(metadata_clean$group == "control")),
    as.character(sum(metadata_clean$group == "disease")),
    as.character(is.numeric(counts_clean)),
    as.character(!anyNA(counts_clean)),
    as.character(identical(colnames(counts_clean), metadata_clean$sample_id)),
    as.character(min(colSums(counts_clean))),
    as.character(median(colSums(counts_clean))),
    as.character(max(colSums(counts_clean)))
  )
)

dataset_decision <- tibble::tribble(
  ~criterion, ~decision,
  "GEO accession", "GSE152418",
  "Platform", "RNA-seq / expression profiling by high throughput sequencing",
  "Species", "Homo sapiens",
  "Sample type", "PBMC",
  "Data type", "Raw gene count matrix",
  "Comparison", "COVID-19-related PBMC samples vs healthy control PBMC samples",
  "Sample size after metadata cleaning", paste0(nrow(metadata_clean), " samples"),
  "Main DESeq2 design planned", "~ group",
  "Primary grouping rule", "Samples with nCoV/nCOV in the GEO sample title were assigned to disease; all remaining PBMC samples were assigned to control.",
  "Age variable", ifelse(all(metadata_clean$age == "not_available"), "not_available", "available"),
  "Sex variable", ifelse(all(metadata_clean$sex == "not_available"), "not_available", "available"),
  "Batch variable", ifelse(length(unique(metadata_clean$batch)) <= 1, "single/constant or not informative", "available")
)

readr::write_tsv(group_summary, here("results", "tables", "group_summary.tsv"))
readr::write_tsv(metadata_summary, here("results", "tables", "metadata_summary.tsv"))
readr::write_tsv(library_size_table, here("results", "tables", "library_sizes.tsv"))
readr::write_tsv(matrix_checks, here("results", "tables", "count_matrix_checks.tsv"))
readr::write_tsv(dataset_decision, here("metadata", "dataset_decision.tsv"))

# -----------------------------
# 8. Simple Initial preprocessing figure
# -----------------------------
p_group <- ggplot(metadata_clean, aes(x = group)) +
  geom_bar() +
  labs(
    title = "Samples retained after metadata curation",
    x = "Analysis group",
    y = "Number of samples"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = here("results", "figures", "group_counts_by_condition.png"),
  plot = p_group,
  width = 5,
  height = 4,
  dpi = 300
)

# -----------------------------
# 9. Update README
# -----------------------------
n_control <- sum(metadata_clean$group == "control")
n_disease <- sum(metadata_clean$group == "disease")

readme_text <- glue::glue(
  "# Reproducible RNA-seq Differential Expression and Pathway Analysis of COVID-19-associated PBMC Transcriptomic Signatures

This project implements a reproducible bulk RNA-seq differential expression and pathway analysis workflow using public GEO count data.

## Dataset

- GEO accession: GSE152418
- Organism: Homo sapiens
- Sample type: PBMC
- Data type: RNA-seq raw gene count matrix
- Comparison: COVID-19-related PBMC samples versus healthy control PBMC samples
- Samples retained after metadata curation: {n_control} control samples and {n_disease} disease samples

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
"
)

writeLines(readme_text, here("README.md"))

message("Metadata cleaning complete.")
message("Final group summary:")
print(group_summary)

message("Final count matrix checks:")
print(matrix_checks)

message("Clean files saved.")