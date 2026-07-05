# scripts/01_download_data.R
# Purpose: download GEO metadata and raw count matrix for GSE152418.

rm(list = ls())

library(tidyverse)
library(GEOquery)
library(Biobase)
library(here)
library(janitor)
library(fs)

here::i_am("scripts/01_download_data.R")

gse_id <- "GSE152418"

raw_dir <- here("data", "raw")
metadata_dir <- here("metadata")
processed_dir <- here("data", "processed")
table_dir <- here("results", "tables")

fs::dir_create(c(raw_dir, metadata_dir, processed_dir, table_dir))

message("Downloading GEO metadata for ", gse_id, "...")

gse_list <- GEOquery::getGEO(
  GEO = gse_id,
  GSEMatrix = TRUE,
  getGPL = FALSE
)

if (length(gse_list) > 1) {
  message("More than one ExpressionSet found. Using the first one.")
}

eset <- gse_list[[1]]

saveRDS(eset, here("data", "raw", "gse152418_eset.rds"))

metadata_raw_df <- as.data.frame(Biobase::pData(eset), stringsAsFactors = FALSE)

metadata_raw <- metadata_raw_df %>%
  tibble::rownames_to_column("geo_accession_row") %>%
  as_tibble()

if (!"geo_accession" %in% names(metadata_raw)) {
  metadata_raw <- metadata_raw %>%
    rename(geo_accession = geo_accession_row)
}

readr::write_tsv(
  metadata_raw,
  here("metadata", "metadata_raw_from_geo.tsv")
)

message("Raw metadata saved: metadata/metadata_raw_from_geo.tsv")
message("Downloading supplementary raw count file...")

GEOquery::getGEOSuppFiles(
  GEO = gse_id,
  baseDir = raw_dir,
  makeDirectory = TRUE
)

count_file <- list.files(
  path = raw_dir,
  pattern = "RawCounts.*\\.txt\\.gz$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(count_file) == 0) {
  stop("No RawCounts file found. Check GEO supplementary file download.")
}

if (length(count_file) > 1) {
  message("More than one RawCounts file found. Using the first one:")
  print(count_file)
}

count_file <- count_file[1]
message("Using count file: ", count_file)

counts_df_raw <- readr::read_tsv(
  file = count_file,
  col_types = readr::cols(.default = "c"),
  name_repair = "minimal",
  show_col_types = FALSE
)

message("Raw count file dimensions:")
print(dim(counts_df_raw))

sample_titles <- metadata_raw$title
geo_ids <- metadata_raw$geo_accession

sample_cols <- intersect(names(counts_df_raw), c(sample_titles, geo_ids))

if (length(sample_cols) < 10) {
  message("Column names in count file:")
  print(names(counts_df_raw))
  message("Sample titles from GEO metadata:")
  print(sample_titles)
  stop("Could not match enough count columns to GEO metadata. Need manual inspection.")
}

id_cols <- setdiff(names(counts_df_raw), sample_cols)

if (length(id_cols) < 1) {
  stop("No gene ID column detected.")
}

gene_id_col <- id_cols[1]

message("Detected gene ID column: ", gene_id_col)
message("Matched sample columns: ", length(sample_cols))

counts_mat <- counts_df_raw %>%
  dplyr::select(all_of(gene_id_col), all_of(sample_cols)) %>%
  dplyr::rename(gene_id = all_of(gene_id_col)) %>%
  dplyr::filter(!is.na(gene_id), gene_id != "") %>%
  dplyr::distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    dplyr::across(
      all_of(sample_cols),
      ~ as.integer(round(as.numeric(.x)))
    )
  ) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

if (anyNA(counts_mat)) {
  stop("Count matrix contains NA values after numeric conversion.")
}

if (!all(counts_mat >= 0)) {
  stop("Count matrix contains negative values, which should not happen for raw counts.")
}

saveRDS(
  counts_mat,
  here("data", "raw", "counts_raw_matrix.rds")
)

counts_out <- as.data.frame(counts_mat, check.names = FALSE) %>%
  tibble::rownames_to_column("gene_id")

readr::write_tsv(
  counts_out,
  here("data", "raw", "counts_raw.tsv")
)

count_file_check <- tibble::tibble(
  item = c(
    "geo_accession",
    "raw_count_file",
    "n_genes",
    "n_samples",
    "matched_sample_columns",
    "first_gene_id",
    "first_sample_id"
  ),
  value = c(
    gse_id,
    basename(count_file),
    as.character(nrow(counts_mat)),
    as.character(ncol(counts_mat)),
    as.character(length(sample_cols)),
    rownames(counts_mat)[1],
    colnames(counts_mat)[1]
  )
)

readr::write_tsv(
  count_file_check,
  here("results", "tables", "download_check.tsv")
)

message("Raw count matrix saved:")
message(" - data/raw/counts_raw_matrix.rds")
message(" - data/raw/counts_raw.tsv")
message("Download script complete.")