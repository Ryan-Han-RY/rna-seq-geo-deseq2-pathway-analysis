# scripts/05_sensitivity_limma_voom_analysis.R
# Purpose: assess robustness of DESeq2 findings using an independent limma-voom sensitivity analysis.

rm(list = ls())

library(tidyverse)
library(here)
library(edgeR)
library(limma)
library(ggrepel)

here::i_am("scripts/05_sensitivity_limma_voom_analysis.R")

# -----------------------------
# 1. Input files
# -----------------------------
counts_path <- here("data", "processed", "counts_filtered_matrix.rds")
metadata_path <- here("data", "processed", "metadata_clean.tsv")
deseq2_path <- here("results", "tables", "deseq2_all_genes_annotated.tsv")

required_files <- c(counts_path, metadata_path, deseq2_path)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

counts_filtered <- readRDS(counts_path)
metadata_clean <- readr::read_tsv(metadata_path, show_col_types = FALSE)
deseq2_res <- readr::read_tsv(deseq2_path, show_col_types = FALSE)

dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "objects"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Hard checks
# -----------------------------
metadata_clean <- metadata_clean %>%
  mutate(group = factor(group, levels = c("control", "disease"))) %>%
  arrange(group, sample_id)

counts_filtered <- counts_filtered[, metadata_clean$sample_id, drop = FALSE]

stopifnot(identical(colnames(counts_filtered), metadata_clean$sample_id))
stopifnot(all(counts_filtered >= 0))
stopifnot(!anyNA(counts_filtered))

storage.mode(counts_filtered) <- "integer"

message("Input filtered count matrix:")
print(dim(counts_filtered))

message("Group table:")
print(table(metadata_clean$group))

# -----------------------------
# 3. limma-voom model
# -----------------------------
dge <- edgeR::DGEList(
  counts = counts_filtered,
  samples = as.data.frame(metadata_clean)
)

dge <- edgeR::calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ 0 + group, data = metadata_clean)
colnames(design) <- gsub("^group", "", colnames(design))

contrast_matrix <- limma::makeContrasts(
  disease_vs_control = disease - control,
  levels = design
)

# Save voom mean-variance plot
png(
  filename = here("results", "figures", "limma_voom_mean_variance_trend.png"),
  width = 1800,
  height = 1400,
  res = 220
)

v <- limma::voom(
  dge,
  design = design,
  plot = TRUE
)

dev.off()

fit <- limma::lmFit(v, design)
fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2, robust = TRUE)

saveRDS(dge, here("results", "objects", "edgeR_dge_list.rds"))
saveRDS(v, here("results", "objects", "limma_voom_object.rds"))
saveRDS(fit2, here("results", "objects", "limma_voom_fitted_model.rds"))

limma_tbl <- limma::topTable(
  fit2,
  coef = "disease_vs_control",
  number = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  rename(
    limma_logFC = logFC,
    limma_AveExpr = AveExpr,
    limma_t = t,
    limma_pvalue = P.Value,
    limma_padj = adj.P.Val,
    limma_B = B
  ) %>%
  mutate(
    ensembl_gene_id = stringr::str_remove(gene_id, "\\..*$"),
    limma_significant = case_when(
      limma_padj < 0.05 & limma_logFC >= 1 ~ "upregulated",
      limma_padj < 0.05 & limma_logFC <= -1 ~ "downregulated",
      TRUE ~ "not_significant"
    )
  )

readr::write_tsv(
  limma_tbl,
  here("results", "tables", "limma_voom_all_genes.tsv")
)

limma_sig <- limma_tbl %>%
  filter(limma_significant != "not_significant")

readr::write_tsv(
  limma_sig,
  here("results", "tables", "limma_voom_significant_genes.tsv")
)

# -----------------------------
# 4. Compare DESeq2 and limma-voom
# -----------------------------
deseq2_compare <- deseq2_res %>%
  mutate(
    ensembl_gene_id = stringr::str_remove(ensembl_gene_id, "\\..*$"),
    deseq2_significant = significant,
    deseq2_log2FC = log2FoldChange,
    deseq2_padj = padj
  ) %>%
  select(
    gene_id,
    ensembl_gene_id,
    gene_symbol,
    entrez_id,
    gene_name,
    deseq2_log2FC,
    deseq2_padj,
    deseq2_significant
  )

comparison_tbl <- deseq2_compare %>%
  inner_join(
    limma_tbl %>%
      select(
        gene_id,
        limma_logFC,
        limma_AveExpr,
        limma_pvalue,
        limma_padj,
        limma_significant
      ),
    by = "gene_id"
  ) %>%
  mutate(
    both_tested = TRUE,
    same_logfc_direction = sign(deseq2_log2FC) == sign(limma_logFC),
    both_significant = deseq2_significant != "not_significant" & limma_significant != "not_significant",
    both_upregulated = deseq2_significant == "upregulated" & limma_significant == "upregulated",
    both_downregulated = deseq2_significant == "downregulated" & limma_significant == "downregulated",
    method_agreement = case_when(
      both_upregulated ~ "both_upregulated",
      both_downregulated ~ "both_downregulated",
      both_significant ~ "both_significant_opposite_or_mixed",
      deseq2_significant != "not_significant" & limma_significant == "not_significant" ~ "deseq2_only",
      deseq2_significant == "not_significant" & limma_significant != "not_significant" ~ "limma_voom_only",
      TRUE ~ "not_significant_in_both"
    )
  )

readr::write_tsv(
  comparison_tbl,
  here("results", "tables", "deseq2_limma_voom_gene_level_comparison.tsv")
)

# Spearman/Pearson correlation across all tested genes
cor_pearson <- cor(
  comparison_tbl$deseq2_log2FC,
  comparison_tbl$limma_logFC,
  use = "complete.obs",
  method = "pearson"
)

cor_spearman <- cor(
  comparison_tbl$deseq2_log2FC,
  comparison_tbl$limma_logFC,
  use = "complete.obs",
  method = "spearman"
)

deseq2_sig_ids <- comparison_tbl %>%
  filter(deseq2_significant != "not_significant") %>%
  pull(gene_id)

limma_sig_ids <- comparison_tbl %>%
  filter(limma_significant != "not_significant") %>%
  pull(gene_id)

shared_sig_ids <- intersect(deseq2_sig_ids, limma_sig_ids)
union_sig_ids <- union(deseq2_sig_ids, limma_sig_ids)

overlap_summary <- tibble(
  metric = c(
    "genes_compared",
    "deseq2_significant_genes",
    "limma_voom_significant_genes",
    "shared_significant_genes",
    "deseq2_only_significant_genes",
    "limma_voom_only_significant_genes",
    "jaccard_index_significant_genes",
    "same_logfc_direction_all_tested",
    "same_logfc_direction_shared_significant",
    "pearson_logfc_correlation",
    "spearman_logfc_correlation"
  ),
  value = c(
    as.character(nrow(comparison_tbl)),
    as.character(length(deseq2_sig_ids)),
    as.character(length(limma_sig_ids)),
    as.character(length(shared_sig_ids)),
    as.character(length(setdiff(deseq2_sig_ids, limma_sig_ids))),
    as.character(length(setdiff(limma_sig_ids, deseq2_sig_ids))),
    round(length(shared_sig_ids) / length(union_sig_ids), 4),
    round(mean(comparison_tbl$same_logfc_direction, na.rm = TRUE), 4),
    round(
      comparison_tbl %>%
        filter(gene_id %in% shared_sig_ids) %>%
        summarise(x = mean(same_logfc_direction, na.rm = TRUE)) %>%
        pull(x),
      4
    ),
    round(cor_pearson, 4),
    round(cor_spearman, 4)
  )
)

readr::write_tsv(
  overlap_summary,
  here("results", "tables", "sensitivity_overlap_summary.tsv")
)

shared_significant_tbl <- comparison_tbl %>%
  filter(gene_id %in% shared_sig_ids) %>%
  arrange(deseq2_padj)

readr::write_tsv(
  shared_significant_tbl,
  here("results", "tables", "shared_significant_genes_deseq2_limma_voom.tsv")
)

top100_deseq2 <- comparison_tbl %>%
  filter(!is.na(deseq2_padj)) %>%
  arrange(deseq2_padj) %>%
  slice_head(n = 100) %>%
  mutate(deseq2_top100 = TRUE) %>%
  select(gene_id, deseq2_top100)

top100_limma <- comparison_tbl %>%
  filter(!is.na(limma_padj)) %>%
  arrange(limma_padj) %>%
  slice_head(n = 100) %>%
  mutate(limma_top100 = TRUE) %>%
  select(gene_id, limma_top100)

top100_comparison <- comparison_tbl %>%
  left_join(top100_deseq2, by = "gene_id") %>%
  left_join(top100_limma, by = "gene_id") %>%
  mutate(
    deseq2_top100 = replace_na(deseq2_top100, FALSE),
    limma_top100 = replace_na(limma_top100, FALSE)
  ) %>%
  filter(deseq2_top100 | limma_top100) %>%
  arrange(deseq2_padj)

readr::write_tsv(
  top100_comparison,
  here("results", "tables", "top100_method_comparison.tsv")
)

message("Sensitivity overlap summary:")
print(overlap_summary)

# -----------------------------
# 5. LogFC comparison figure
# -----------------------------
label_tbl <- comparison_tbl %>%
  filter(
    gene_id %in% shared_sig_ids,
    !is.na(deseq2_padj),
    !is.na(limma_padj)
  ) %>%
  arrange(deseq2_padj) %>%
  slice_head(n = 15)

p_logfc <- ggplot(
  comparison_tbl,
  aes(x = deseq2_log2FC, y = limma_logFC)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(color = method_agreement), alpha = 0.55, size = 1.2) +
  ggrepel::geom_text_repel(
    data = label_tbl,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = 30
  ) +
  labs(
    title = "DESeq2 versus limma-voom log fold-change comparison",
    subtitle = paste0(
      "Pearson r = ", round(cor_pearson, 3),
      "; Spearman rho = ", round(cor_spearman, 3)
    ),
    x = "DESeq2 shrunken log2 fold change",
    y = "limma-voom log2 fold change",
    color = "Method agreement"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right"
  )

ggsave(
  filename = here("results", "figures", "deseq2_limma_voom_logfc_comparison.png"),
  plot = p_logfc,
  width = 8,
  height = 6.5,
  dpi = 300
)

# -----------------------------
# 6. Agreement category bar plot
# -----------------------------
agreement_tbl <- comparison_tbl %>%
  dplyr::group_by(method_agreement) %>%
  dplyr::summarise(n_genes = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_genes))

readr::write_tsv(
  agreement_tbl,
  here("results", "tables", "method_agreement_categories.tsv")
)

p_agreement <- ggplot(
  agreement_tbl,
  aes(x = reorder(method_agreement, n_genes), y = n_genes)
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "DESeq2 and limma-voom agreement categories",
    x = "Agreement category",
    y = "Number of genes"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = here("results", "figures", "method_agreement_categories.png"),
  plot = p_agreement,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

# -----------------------------
# 7. Notes for report
# -----------------------------
notes <- c(
  "# limma-voom sensitivity analysis notes",
  "",
  "An independent limma-voom sensitivity analysis was performed to assess whether the major DESeq2 findings were robust to an alternative differential expression framework.",
  "",
  "The same filtered count matrix and metadata were used. edgeR was used to create the DGEList object and calculate TMM normalisation factors, followed by limma-voom transformation, linear modelling and empirical Bayes moderation.",
  "",
  "## Summary",
  "",
  paste0("- Genes compared: ", overlap_summary$value[overlap_summary$metric == "genes_compared"]),
  paste0("- DESeq2 significant genes: ", overlap_summary$value[overlap_summary$metric == "deseq2_significant_genes"]),
  paste0("- limma-voom significant genes: ", overlap_summary$value[overlap_summary$metric == "limma_voom_significant_genes"]),
  paste0("- Shared significant genes: ", overlap_summary$value[overlap_summary$metric == "shared_significant_genes"]),
  paste0("- Significant-gene Jaccard index: ", overlap_summary$value[overlap_summary$metric == "jaccard_index_significant_genes"]),
  paste0("- Pearson logFC correlation: ", overlap_summary$value[overlap_summary$metric == "pearson_logfc_correlation"]),
  paste0("- Spearman logFC correlation: ", overlap_summary$value[overlap_summary$metric == "spearman_logfc_correlation"]),
  "",
  "## Interpretation",
  "",
  "Strong agreement in log fold-change direction and positive logFC correlation would support the robustness of the main disease-control transcriptional signal. Differences in exact significant-gene calls are expected because DESeq2 and limma-voom use different mean-variance modelling assumptions and test statistics."
)

writeLines(
  notes,
  here("results", "sensitivity_analysis_notes.md")
)

message("limma-voom sensitivity analysis complete.")
message("Main outputs:")
message(" - results/tables/limma_voom_all_genes.tsv")
message(" - results/tables/sensitivity_overlap_summary.tsv")
message(" - results/tables/shared_significant_genes_deseq2_limma_voom.tsv")
message(" - results/figures/deseq2_limma_voom_logfc_comparison.png")
message(" - results/figures/method_agreement_categories.png")
message(" - results/sensitivity_analysis_notes.md")

# -----------------------------
# 8. Significant-gene agreement plot for final report
# -----------------------------
significant_agreement_tbl <- tibble::tibble(
  category = c(
    "Shared significant",
    "DESeq2 only",
    "limma-voom only"
  ),
  n_genes = c(
    length(shared_sig_ids),
    length(setdiff(deseq2_sig_ids, limma_sig_ids)),
    length(setdiff(limma_sig_ids, deseq2_sig_ids))
  )
)

readr::write_tsv(
  significant_agreement_tbl,
  here("results", "tables", "significant_gene_overlap_categories.tsv")
)

p_sig_agreement <- ggplot(
  significant_agreement_tbl,
  aes(x = reorder(category, n_genes), y = n_genes)
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Significant gene overlap between DESeq2 and limma-voom",
    x = NULL,
    y = "Number of significant genes"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = here("results", "figures", "significant_gene_overlap_deseq2_limma_voom.png"),
  plot = p_sig_agreement,
  width = 7,
  height = 4.5,
  dpi = 300
)