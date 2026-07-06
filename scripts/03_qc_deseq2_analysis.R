# scripts/03_qc_deseq2_analysis.R
# Purpose: low-count filtering, DESeq2 modelling, QC plots, and primary differential expression outputs.

rm(list = ls())

library(tidyverse)
library(here)
library(DESeq2)
library(apeglm)
library(pheatmap)
library(ggrepel)

here::i_am("scripts/03_qc_deseq2_analysis.R")

# -----------------------------
# 1. Input files
# -----------------------------
metadata_path <- here("data", "processed", "metadata_clean.tsv")
counts_path <- here("data", "processed", "counts_clean_matrix.rds")

if (!file.exists(metadata_path)) {
  stop("Missing data/processed/metadata_clean.tsv. Finish metadata curation first.")
}

if (!file.exists(counts_path)) {
  stop("Missing data/processed/counts_clean_matrix.rds. Finish count matrix cleaning first.")
}

metadata_clean <- readr::read_tsv(metadata_path, show_col_types = FALSE)
counts_clean <- readRDS(counts_path)

# -----------------------------
# 2. Hard checks before modelling
# -----------------------------
metadata_clean <- metadata_clean %>%
  mutate(
    group = factor(group, levels = c("control", "disease"))
  ) %>%
  arrange(group, sample_id)

counts_clean <- counts_clean[, metadata_clean$sample_id, drop = FALSE]

stopifnot(identical(colnames(counts_clean), metadata_clean$sample_id))
stopifnot(all(c("control", "disease") %in% levels(metadata_clean$group)))
stopifnot(!anyNA(counts_clean))
stopifnot(all(counts_clean >= 0))

# DESeq2 expects integer count data.
counts_clean <- round(counts_clean)
storage.mode(counts_clean) <- "integer"

message("Input count matrix:")
print(dim(counts_clean))

message("Group table:")
print(table(metadata_clean$group))

# -----------------------------
# 3. Create output folders
# -----------------------------
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "processed"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 4. Low-count filtering
# -----------------------------
smallest_group_size <- min(table(metadata_clean$group))

keep_genes <- rowSums(counts_clean >= 10) >= smallest_group_size
counts_filtered <- counts_clean[keep_genes, , drop = FALSE]

filtering_summary <- tibble(
  metric = c(
    "genes_before_filtering",
    "genes_after_filtering",
    "genes_removed",
    "smallest_group_size",
    "filtering_rule"
  ),
  value = c(
    as.character(nrow(counts_clean)),
    as.character(nrow(counts_filtered)),
    as.character(nrow(counts_clean) - nrow(counts_filtered)),
    as.character(smallest_group_size),
    "Retained genes with counts >= 10 in at least the smallest group size"
  )
)

readr::write_tsv(
  filtering_summary,
  here("results", "tables", "filtering_summary.tsv")
)

counts_filtered_out <- as.data.frame(counts_filtered, check.names = FALSE) %>%
  rownames_to_column("gene_id")

readr::write_tsv(
  counts_filtered_out,
  here("data", "processed", "counts_filtered.tsv")
)

saveRDS(
  counts_filtered,
  here("data", "processed", "counts_filtered_matrix.rds")
)

message("Filtering summary:")
print(filtering_summary)

# -----------------------------
# 5. Library size QC
# -----------------------------
library_size_tbl <- tibble(
  sample_id = colnames(counts_filtered),
  library_size = colSums(counts_filtered)
) %>%
  left_join(
    metadata_clean %>% select(sample_id, group, disease_state, severity),
    by = "sample_id"
  )

readr::write_tsv(
  library_size_tbl,
  here("results", "tables", "library_sizes_after_filtering.tsv")
)

p_library <- ggplot(library_size_tbl, aes(x = reorder(sample_id, library_size), y = library_size, fill = group)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Library sizes after low-count filtering",
    x = "Sample",
    y = "Total filtered counts"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top"
  )

ggsave(
  filename = here("results", "figures", "library_sizes_after_filtering.png"),
  plot = p_library,
  width = 8,
  height = 7,
  dpi = 300
)

# -----------------------------
# 6. Build DESeq2 object
# -----------------------------
metadata_df <- as.data.frame(metadata_clean)
rownames(metadata_df) <- metadata_df$sample_id

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata_df,
  design = ~ group
)

dds$group <- relevel(dds$group, ref = "control")

# Remove any all-zero rows after construction, just as a final safety check.
dds <- dds[rowSums(counts(dds)) > 0, ]

saveRDS(
  dds,
  here("results", "objects", "deseq2_dataset_prefit.rds")
)

# -----------------------------
# 7. Run DESeq2 model
# -----------------------------
dds <- DESeq(dds)

saveRDS(
  dds,
  here("results", "objects", "deseq2_dataset_fitted.rds")
)

message("DESeq2 result names:")
print(resultsNames(dds))

# Save dispersion plot
png(
  filename = here("results", "figures", "dispersion_estimates.png"),
  width = 1800,
  height = 1400,
  res = 220
)
plotDispEsts(dds)
dev.off()

# -----------------------------
# 8. Variance stabilising transformation
# -----------------------------
vsd <- vst(dds, blind = FALSE)

saveRDS(
  vsd,
  here("results", "objects", "vst_transformed_counts.rds")
)

vst_mat <- assay(vsd)

vst_out <- as.data.frame(vst_mat, check.names = FALSE) %>%
  rownames_to_column("gene_id")

readr::write_tsv(
  vst_out,
  here("data", "processed", "vst_transformed_counts.tsv")
)

# -----------------------------
# 9. PCA
# -----------------------------
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

pca_out <- pca_data %>%
  as_tibble() %>%
  rename(sample_id = name) %>%
  select(sample_id, group, PC1, PC2)

readr::write_tsv(
  pca_out,
  here("results", "tables", "pca_sample_coordinates.tsv")
)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 2.6, max.overlaps = 30) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  labs(
    title = "PCA of variance-stabilised expression profiles"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top"
  )

ggsave(
  filename = here("results", "figures", "pca_vst_group.png"),
  plot = p_pca,
  width = 7,
  height = 6,
  dpi = 300
)

# -----------------------------
# 10. Sample distance heatmap
# -----------------------------
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

rownames(sample_dist_matrix) <- metadata_clean$sample_id
colnames(sample_dist_matrix) <- metadata_clean$sample_id

annotation_col <- metadata_clean %>%
  select(sample_id, group) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$sample_id
annotation_col$sample_id <- NULL

png(
  filename = here("results", "figures", "sample_distance_heatmap.png"),
  width = 2200,
  height = 2000,
  res = 230
)

pheatmap(
  sample_dist_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  main = "Sample-to-sample distance heatmap",
  fontsize = 8
)

dev.off()

sample_distance_out <- as.data.frame(sample_dist_matrix, check.names = FALSE) %>%
  rownames_to_column("sample_id")

readr::write_tsv(
  sample_distance_out,
  here("results", "tables", "sample_distance_matrix.tsv")
)

# -----------------------------
# 11. Extract DESeq2 results
# -----------------------------
res_raw <- results(
  dds,
  contrast = c("group", "disease", "control"),
  alpha = 0.05
)

coef_name <- "group_disease_vs_control"

if (!coef_name %in% resultsNames(dds)) {
  stop(
    paste0(
      "Expected coefficient '", coef_name, "' was not found. Available coefficients: ",
      paste(resultsNames(dds), collapse = ", ")
    )
  )
}

res_shrunk <- lfcShrink(
  dds,
  coef = coef_name,
  type = "apeglm"
)

saveRDS(
  res_raw,
  here("results", "objects", "deseq2_results_raw.rds")
)

saveRDS(
  res_shrunk,
  here("results", "objects", "deseq2_results_shrunken.rds")
)

res_tbl <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(
    significant = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >= 1 ~ "upregulated",
      !is.na(padj) & padj < 0.05 & log2FoldChange <= -1 ~ "downregulated",
      TRUE ~ "not_significant"
    )
  )

readr::write_tsv(
  res_tbl,
  here("results", "tables", "deseq2_all_genes.tsv")
)

sig_genes <- res_tbl %>%
  filter(significant != "not_significant")

readr::write_tsv(
  sig_genes,
  here("results", "tables", "deseq2_significant_genes.tsv")
)

result_summary <- tibble(
  metric = c(
    "tested_genes",
    "genes_with_adjusted_pvalue",
    "significant_genes_fdr_0_05_abs_lfc_1",
    "upregulated_genes",
    "downregulated_genes",
    "contrast",
    "lfc_shrinkage"
  ),
  value = c(
    as.character(nrow(res_tbl)),
    as.character(sum(!is.na(res_tbl$padj))),
    as.character(nrow(sig_genes)),
    as.character(sum(sig_genes$significant == "upregulated")),
    as.character(sum(sig_genes$significant == "downregulated")),
    "disease vs control",
    "apeglm"
  )
)

readr::write_tsv(
  result_summary,
  here("results", "tables", "deseq2_result_summary.tsv")
)

message("DESeq2 result summary:")
print(result_summary)

# -----------------------------
# 12. MA plot
# -----------------------------
png(
  filename = here("results", "figures", "ma_plot_shrunken_lfc.png"),
  width = 1800,
  height = 1400,
  res = 220
)

plotMA(
  res_shrunk,
  ylim = c(-5, 5),
  main = "MA plot with shrunken log2 fold changes"
)

dev.off()

# -----------------------------
# 13. Volcano plot
# -----------------------------
volcano_tbl <- res_tbl %>%
  mutate(
    padj_for_plot = if_else(
      is.na(padj),
      NA_real_,
      pmax(padj, .Machine$double.xmin)
    ),
    neg_log10_padj = -log10(padj_for_plot)
  )

top_labels <- volcano_tbl %>%
  filter(!is.na(padj), significant != "not_significant") %>%
  arrange(padj) %>%
  slice_head(n = 15)

p_volcano <- ggplot(
  volcano_tbl,
  aes(x = log2FoldChange, y = neg_log10_padj, color = significant)
) +
  geom_point(alpha = 0.65, size = 1.2, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = gene_id),
    size = 2.6,
    max.overlaps = 30
  ) +
  labs(
    title = "Differential expression volcano plot",
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top"
  )

ggsave(
  filename = here("results", "figures", "volcano_plot_deseq2.png"),
  plot = p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

# -----------------------------
# 14. Top gene heatmap
# -----------------------------
top_genes <- res_tbl %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(gene_id)

top_gene_tbl <- res_tbl %>%
  filter(gene_id %in% top_genes) %>%
  arrange(padj)

readr::write_tsv(
  top_gene_tbl,
  here("results", "tables", "top50_genes_for_heatmap.tsv")
)

heatmap_mat <- assay(vsd)[top_genes, , drop = FALSE]
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
heatmap_mat_scaled[is.na(heatmap_mat_scaled)] <- 0

png(
  filename = here("results", "figures", "top50_gene_heatmap.png"),
  width = 2200,
  height = 2400,
  res = 240
)

pheatmap(
  heatmap_mat_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top 50 differentially expressed genes",
  fontsize_row = 6,
  fontsize_col = 6
)

dev.off()

# -----------------------------
# 15. Initial interpretation notes
# -----------------------------
n_sig <- nrow(sig_genes)
n_up <- sum(sig_genes$significant == "upregulated")
n_down <- sum(sig_genes$significant == "downregulated")

top_gene_text <- res_tbl %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  mutate(
    line = paste0(
      "- ", gene_id,
      ": log2FC = ", round(log2FoldChange, 3),
      ", padj = ", signif(padj, 3)
    )
  ) %>%
  pull(line)

analysis_notes <- c(
  "# QC and DESeq2 differential expression analysis notes",
  "",
  "## Filtering",
  "",
  paste0(
    "Low-expression genes were removed by retaining genes with counts >= 10 in at least ",
    smallest_group_size,
    " samples, corresponding to the smallest analysis group size."
  ),
  "",
  paste0(
    "The count matrix contained ",
    nrow(counts_clean),
    " genes before filtering and ",
    nrow(counts_filtered),
    " genes after filtering."
  ),
  "",
  "## PCA and sample-level QC",
  "",
  paste0(
    "PCA was performed on variance-stabilised expression values. PC1 explained ",
    percent_var[1],
    "% of variance and PC2 explained ",
    percent_var[2],
    "% of variance."
  ),
  "",
  "The PCA plot and sample distance heatmap should be inspected for group separation, outlying samples, or unexpected clustering.",
  "",
  "## Differential expression",
  "",
  paste0(
    "Using DESeq2 with design ~ group and the contrast disease vs control, ",
    n_sig,
    " genes passed the threshold padj < 0.05 and absolute shrunken log2 fold change >= 1."
  ),
  "",
  paste0(n_up, " genes were upregulated and ", n_down, " genes were downregulated in disease relative to control."),
  "",
  "## Top-ranked genes by adjusted p-value",
  "",
  top_gene_text,
  "",
  "## Next interpretation step",
  "",
  "The differential expression table should be connected to gene annotation and pathway-level analysis, including GO, KEGG, Reactome and ranked-list enrichment."
)

writeLines(
  analysis_notes,
  here("results", "analysis_notes.md")
)

message("QC and DESeq2 analysis complete.")
message("Main outputs:")
message(" - data/processed/counts_filtered.tsv")
message(" - results/tables/deseq2_all_genes.tsv")
message(" - results/tables/deseq2_significant_genes.tsv")
message(" - results/figures/pca_vst_group.png")
message(" - results/figures/sample_distance_heatmap.png")
message(" - results/figures/ma_plot_shrunken_lfc.png")
message(" - results/figures/volcano_plot_deseq2.png")
message(" - results/figures/top50_gene_heatmap.png")
message(" - results/analysis_notes.md")
