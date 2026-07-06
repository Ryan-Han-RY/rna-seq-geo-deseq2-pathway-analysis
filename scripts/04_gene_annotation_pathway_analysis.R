# scripts/04_gene_annotation_pathway_analysis.R
# Purpose: annotate DESeq2 results, improve gene-level visualisations, and perform GO/KEGG/Reactome/GSEA-style pathway analysis.

rm(list = ls())

library(tidyverse)
library(here)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(ggrepel)
library(pheatmap)
library(SummarizedExperiment)

here::i_am("scripts/04_gene_annotation_pathway_analysis.R")

# -----------------------------
# 1. Input files
# -----------------------------
res_path <- here("results", "tables", "deseq2_all_genes.tsv")
sig_path <- here("results", "tables", "deseq2_significant_genes.tsv")
vsd_path <- here("results", "objects", "vst_transformed_counts.rds")
metadata_path <- here("data", "processed", "metadata_clean.tsv")

required_files <- c(res_path, sig_path, vsd_path, metadata_path)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

res_tbl <- readr::read_tsv(res_path, show_col_types = FALSE)
sig_tbl <- readr::read_tsv(sig_path, show_col_types = FALSE)
vsd <- readRDS(vsd_path)
metadata_clean <- readr::read_tsv(metadata_path, show_col_types = FALSE)

dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "objects"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Clean Ensembl IDs
# -----------------------------
res_tbl <- res_tbl %>%
  mutate(
    ensembl_gene_id = str_remove(gene_id, "\\..*$")
  )

sig_tbl <- sig_tbl %>%
  mutate(
    ensembl_gene_id = str_remove(gene_id, "\\..*$")
  )

all_ensembl_ids <- unique(res_tbl$ensembl_gene_id)

# -----------------------------
# 3. Gene annotation: ENSEMBL -> SYMBOL / ENTREZID / GENENAME
# -----------------------------
annotation_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_ensembl_ids,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL", "ENTREZID", "GENENAME")
)

annotation_table <- annotation_raw %>%
  as_tibble() %>%
  rename(
    ensembl_gene_id = ENSEMBL,
    gene_symbol = SYMBOL,
    entrez_id = ENTREZID,
    gene_name = GENENAME
  ) %>%
  arrange(ensembl_gene_id, gene_symbol, entrez_id) %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    gene_symbol = first(na.omit(gene_symbol)),
    entrez_id = first(na.omit(entrez_id)),
    gene_name = first(na.omit(gene_name)),
    .groups = "drop"
  ) %>%
  mutate(
    gene_symbol = if_else(is.na(gene_symbol), ensembl_gene_id, gene_symbol),
    gene_label = gene_symbol
  )

readr::write_tsv(
  annotation_table,
  here("results", "tables", "gene_annotation_table.tsv")
)

annotation_summary <- tibble(
  metric = c(
    "genes_in_deseq2_table",
    "genes_with_symbol",
    "genes_with_entrez_id",
    "symbol_mapping_rate",
    "entrez_mapping_rate"
  ),
  value = c(
    as.character(n_distinct(res_tbl$ensembl_gene_id)),
    as.character(sum(!is.na(annotation_table$gene_symbol) & annotation_table$gene_symbol != annotation_table$ensembl_gene_id)),
    as.character(sum(!is.na(annotation_table$entrez_id))),
    round(mean(!is.na(annotation_table$gene_symbol) & annotation_table$gene_symbol != annotation_table$ensembl_gene_id), 4),
    round(mean(!is.na(annotation_table$entrez_id)), 4)
  )
)

readr::write_tsv(
  annotation_summary,
  here("results", "tables", "gene_annotation_summary.tsv")
)

message("Gene annotation summary:")
print(annotation_summary)

# -----------------------------
# 4. Annotated DESeq2 result tables
# -----------------------------
res_annotated <- res_tbl %>%
  left_join(annotation_table, by = "ensembl_gene_id") %>%
  mutate(
    gene_symbol = if_else(is.na(gene_symbol), ensembl_gene_id, gene_symbol),
    gene_label = if_else(is.na(gene_label), gene_symbol, gene_label)
  ) %>%
  select(
    gene_id,
    ensembl_gene_id,
    gene_symbol,
    entrez_id,
    gene_name,
    baseMean,
    log2FoldChange,
    lfcSE,
    pvalue,
    padj,
    significant
  ) %>%
  arrange(padj)

sig_annotated <- res_annotated %>%
  filter(significant != "not_significant")

readr::write_tsv(
  res_annotated,
  here("results", "tables", "deseq2_all_genes_annotated.tsv")
)

readr::write_tsv(
  sig_annotated,
  here("results", "tables", "deseq2_significant_genes_annotated.tsv")
)

top_annotated_genes <- res_annotated %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 50)

readr::write_tsv(
  top_annotated_genes,
  here("results", "tables", "top50_annotated_genes.tsv")
)

# -----------------------------
# 5. Improved volcano plot with gene symbols
# -----------------------------
volcano_tbl <- res_annotated %>%
  mutate(
    padj_for_plot = if_else(is.na(padj), NA_real_, pmax(padj, .Machine$double.xmin)),
    neg_log10_padj = -log10(padj_for_plot),
    direction_for_plot = factor(
      significant,
      levels = c("downregulated", "not_significant", "upregulated")
    )
  )

top_labels <- volcano_tbl %>%
  filter(!is.na(padj), significant != "not_significant") %>%
  arrange(padj) %>%
  slice_head(n = 15)

p_volcano_symbol <- ggplot(
  volcano_tbl,
  aes(x = log2FoldChange, y = neg_log10_padj, color = direction_for_plot)
) +
  geom_point(alpha = 0.6, size = 1.15, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = 50
  ) +
  labs(
    title = "Differential expression volcano plot",
    subtitle = "COVID-19-related PBMC samples versus healthy controls",
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p-value",
    color = "DESeq2 category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top"
  )

ggsave(
  filename = here("results", "figures", "volcano_plot_deseq2_gene_symbols.png"),
  plot = p_volcano_symbol,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

# -----------------------------
# 6. Improved top gene heatmap with symbols
# -----------------------------
metadata_clean <- metadata_clean %>%
  mutate(group = factor(group, levels = c("control", "disease"))) %>%
  arrange(group, sample_id)

annotation_col <- metadata_clean %>%
  select(sample_id, group) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$sample_id
annotation_col$sample_id <- NULL

top_gene_ids <- top_annotated_genes$gene_id

heatmap_mat <- SummarizedExperiment::assay(vsd)[top_gene_ids, metadata_clean$sample_id, drop = FALSE]

row_labels <- top_annotated_genes$gene_symbol
row_labels[is.na(row_labels) | row_labels == ""] <- top_annotated_genes$ensembl_gene_id
rownames(heatmap_mat) <- make.unique(row_labels)

heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
heatmap_mat_scaled[is.na(heatmap_mat_scaled)] <- 0

png(
  filename = here("results", "figures", "top50_gene_heatmap_gene_symbols.png"),
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
  fontsize_row = 7,
  fontsize_col = 6
)

dev.off()

# -----------------------------
# 7. Prepare gene sets for enrichment analysis
# -----------------------------
universe_entrez <- res_annotated %>%
  filter(!is.na(entrez_id)) %>%
  pull(entrez_id) %>%
  unique()

up_entrez <- sig_annotated %>%
  filter(significant == "upregulated", !is.na(entrez_id)) %>%
  pull(entrez_id) %>%
  unique()

down_entrez <- sig_annotated %>%
  filter(significant == "downregulated", !is.na(entrez_id)) %>%
  pull(entrez_id) %>%
  unique()

enrichment_input_summary <- tibble(
  gene_set = c("universe", "upregulated", "downregulated"),
  n_entrez_genes = c(length(universe_entrez), length(up_entrez), length(down_entrez))
)

readr::write_tsv(
  enrichment_input_summary,
  here("results", "tables", "enrichment_input_summary.tsv")
)

message("Enrichment input summary:")
print(enrichment_input_summary)

# Helper to save enrichment result tables safely
save_enrichment_table <- function(enrich_obj, path) {
  if (is.null(enrich_obj)) {
    readr::write_tsv(tibble(), path)
    return(invisible(NULL))
  }
  
  out <- as.data.frame(enrich_obj) %>% as_tibble()
  
  readr::write_tsv(out, path)
  invisible(out)
}

save_dotplot <- function(enrich_obj, path, title, n_show = 15) {
  if (is.null(enrich_obj)) {
    return(invisible(NULL))
  }
  
  enrich_df <- as.data.frame(enrich_obj)
  
  if (nrow(enrich_df) == 0) {
    return(invisible(NULL))
  }
  
  p <- enrichplot::dotplot(enrich_obj, showCategory = min(n_show, nrow(enrich_df))) +
    ggtitle(title) +
    theme_bw(base_size = 11)
  
  ggsave(
    filename = path,
    plot = p,
    width = 8.5,
    height = 6.5,
    dpi = 300
  )
  
  invisible(p)
}

# -----------------------------
# 8. GO Biological Process enrichment
# -----------------------------
ego_up <- enrichGO(
  gene = up_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_down <- enrichGO(
  gene = down_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_up_tbl <- save_enrichment_table(
  ego_up,
  here("results", "tables", "go_bp_enrichment_upregulated.tsv")
)

go_down_tbl <- save_enrichment_table(
  ego_down,
  here("results", "tables", "go_bp_enrichment_downregulated.tsv")
)

saveRDS(ego_up, here("results", "objects", "go_bp_enrichment_upregulated.rds"))
saveRDS(ego_down, here("results", "objects", "go_bp_enrichment_downregulated.rds"))

save_dotplot(
  ego_up,
  here("results", "figures", "go_bp_dotplot_upregulated.png"),
  "GO Biological Process enrichment: upregulated genes"
)

save_dotplot(
  ego_down,
  here("results", "figures", "go_bp_dotplot_downregulated.png"),
  "GO Biological Process enrichment: downregulated genes"
)

# -----------------------------
# 9. Reactome pathway enrichment
# -----------------------------
reactome_up <- enrichPathway(
  gene = up_entrez,
  universe = universe_entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  readable = TRUE
)

reactome_down <- enrichPathway(
  gene = down_entrez,
  universe = universe_entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  readable = TRUE
)

reactome_up_tbl <- save_enrichment_table(
  reactome_up,
  here("results", "tables", "reactome_enrichment_upregulated.tsv")
)

reactome_down_tbl <- save_enrichment_table(
  reactome_down,
  here("results", "tables", "reactome_enrichment_downregulated.tsv")
)

saveRDS(reactome_up, here("results", "objects", "reactome_enrichment_upregulated.rds"))
saveRDS(reactome_down, here("results", "objects", "reactome_enrichment_downregulated.rds"))

save_dotplot(
  reactome_up,
  here("results", "figures", "reactome_dotplot_upregulated.png"),
  "Reactome pathway enrichment: upregulated genes"
)

save_dotplot(
  reactome_down,
  here("results", "figures", "reactome_dotplot_downregulated.png"),
  "Reactome pathway enrichment: downregulated genes"
)

# -----------------------------
# 10. KEGG enrichment
# -----------------------------
# KEGG queries can occasionally fail if the online KEGG service is temporarily unavailable.
# The script keeps running and records empty outputs if that happens.
kegg_up <- tryCatch(
  enrichKEGG(
    gene = up_entrez,
    universe = universe_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2
  ),
  error = function(e) {
    message("KEGG upregulated enrichment failed: ", e$message)
    NULL
  }
)

kegg_down <- tryCatch(
  enrichKEGG(
    gene = down_entrez,
    universe = universe_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2
  ),
  error = function(e) {
    message("KEGG downregulated enrichment failed: ", e$message)
    NULL
  }
)

kegg_up_tbl <- save_enrichment_table(
  kegg_up,
  here("results", "tables", "kegg_enrichment_upregulated.tsv")
)

kegg_down_tbl <- save_enrichment_table(
  kegg_down,
  here("results", "tables", "kegg_enrichment_downregulated.tsv")
)

saveRDS(kegg_up, here("results", "objects", "kegg_enrichment_upregulated.rds"))
saveRDS(kegg_down, here("results", "objects", "kegg_enrichment_downregulated.rds"))

save_dotplot(
  kegg_up,
  here("results", "figures", "kegg_dotplot_upregulated.png"),
  "KEGG enrichment: upregulated genes"
)

save_dotplot(
  kegg_down,
  here("results", "figures", "kegg_dotplot_downregulated.png"),
  "KEGG enrichment: downregulated genes"
)

# -----------------------------
# 11. Ranked-list GSEA-style GO analysis
# -----------------------------
rank_tbl <- res_annotated %>%
  filter(
    !is.na(entrez_id),
    !is.na(log2FoldChange),
    !is.na(pvalue)
  ) %>%
  mutate(
    # Signed ranking metric: direction from log2FC, strength from p-value.
    rank_metric = sign(log2FoldChange) * -log10(pmax(pvalue, .Machine$double.xmin))
  ) %>%
  group_by(entrez_id) %>%
  slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(rank_metric))

gene_rank <- rank_tbl$rank_metric
names(gene_rank) <- rank_tbl$entrez_id
gene_rank <- sort(gene_rank, decreasing = TRUE)

readr::write_tsv(
  rank_tbl %>% select(entrez_id, gene_symbol, gene_name, log2FoldChange, pvalue, padj, rank_metric),
  here("results", "tables", "ranked_gene_list_for_gsea.tsv")
)

gsea_go_bp <- gseGO(
  geneList = gene_rank,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

gsea_go_tbl <- save_enrichment_table(
  gsea_go_bp,
  here("results", "tables", "gsea_go_bp_ranked_results.tsv")
)

saveRDS(
  gsea_go_bp,
  here("results", "objects", "gsea_go_bp_ranked_results.rds")
)

save_dotplot(
  gsea_go_bp,
  here("results", "figures", "gsea_go_bp_dotplot.png"),
  "Ranked-list GSEA-style GO Biological Process analysis"
)

# -----------------------------
# 12. Create pathway summary table
# -----------------------------
make_directional_enrichment_summary <- function(tbl, source_name, direction_name, n = 10) {
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(tibble())
  }
  
  tbl %>%
    slice_head(n = n) %>%
    transmute(
      source = source_name,
      direction = direction_name,
      pathway_id = ID,
      pathway = Description,
      adjusted_pvalue = p.adjust,
      gene_ratio = GeneRatio,
      gene_count = Count,
      leading_genes = geneID
    )
}

pathway_summary <- bind_rows(
  make_directional_enrichment_summary(go_up_tbl, "GO Biological Process", "upregulated", 10),
  make_directional_enrichment_summary(go_down_tbl, "GO Biological Process", "downregulated", 10),
  make_directional_enrichment_summary(reactome_up_tbl, "Reactome", "upregulated", 10),
  make_directional_enrichment_summary(reactome_down_tbl, "Reactome", "downregulated", 10),
  make_directional_enrichment_summary(kegg_up_tbl, "KEGG", "upregulated", 10),
  make_directional_enrichment_summary(kegg_down_tbl, "KEGG", "downregulated", 10)
) %>%
  arrange(adjusted_pvalue)

readr::write_tsv(
  pathway_summary,
  here("results", "tables", "pathway_enrichment_summary.tsv")
)

# -----------------------------
# 13. Biological interpretation table
# -----------------------------
interpret_theme <- function(x) {
  x_lower <- str_to_lower(x)
  
  case_when(
    str_detect(x_lower, "interferon|viral|virus|defense response|innate immune|cytokine|chemokine") ~
      "Immune and antiviral response",
    str_detect(x_lower, "neutrophil|granulocyte|myeloid|leukocyte|phagocyt") ~
      "Innate immune cell activation",
    str_detect(x_lower, "antigen|mhc|t cell|lymphocyte|adaptive immune") ~
      "Antigen presentation or adaptive immune response",
    str_detect(x_lower, "ribosome|translation|protein synthesis|rna processing") ~
      "Protein synthesis or RNA-processing programme",
    str_detect(x_lower, "mitochond|oxidative|respiratory|metabolic|metabolism") ~
      "Metabolic or mitochondrial regulation",
    str_detect(x_lower, "cell cycle|mitotic|dna replication|chromosome") ~
      "Cell-cycle or proliferative activity",
    TRUE ~ "Biological pathway-level signal"
  )
}

why_it_matters <- function(theme, direction) {
  case_when(
    theme == "Immune and antiviral response" & direction == "upregulated" ~
      "Upregulation of antiviral and inflammatory programmes is consistent with systemic immune activation in disease samples.",
    theme == "Innate immune cell activation" & direction == "upregulated" ~
      "This suggests activation or altered abundance of innate immune cell populations in PBMC disease samples.",
    theme == "Antigen presentation or adaptive immune response" & direction == "upregulated" ~
      "This may reflect immune recognition and coordination of adaptive immune responses in disease.",
    theme == "Metabolic or mitochondrial regulation" ~
      "Changes in metabolic pathways can reflect altered immune-cell activation state and cellular stress.",
    theme == "Cell-cycle or proliferative activity" ~
      "Cell-cycle signals may indicate changes in proliferating immune-cell subsets or sample composition.",
    direction == "downregulated" ~
      "Downregulation suggests that this programme is relatively suppressed in disease compared with healthy controls.",
    TRUE ~
      "This pathway provides a candidate biological axis for interpreting the disease-control transcriptomic contrast."
  )
}

top_gene_interpretation <- sig_annotated %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%
  transmute(
    item_type = "gene",
    gene_or_pathway = gene_symbol,
    direction = significant,
    statistical_evidence = paste0(
      "padj=", signif(padj, 3),
      "; shrunken log2FC=", round(log2FoldChange, 3)
    ),
    biological_meaning = if_else(
      log2FoldChange > 0,
      "Among the strongest upregulated gene-level signals in disease samples.",
      "Among the strongest downregulated gene-level signals in disease samples."
    ),
    why_it_matters_in_context = "This gene contributes to the gene-level signature separating disease and healthy control PBMC samples."
  )

top_pathway_interpretation <- pathway_summary %>%
  slice_head(n = 20) %>%
  mutate(
    biological_theme = interpret_theme(pathway),
    why_context = why_it_matters(biological_theme, direction)
  ) %>%
  transmute(
    item_type = "pathway",
    gene_or_pathway = pathway,
    direction = direction,
    statistical_evidence = paste0(
      source,
      "; adjusted p-value=", signif(adjusted_pvalue, 3),
      "; gene ratio=", gene_ratio,
      "; count=", gene_count
    ),
    biological_meaning = biological_theme,
    why_it_matters_in_context = why_context
  )

biological_interpretation_table <- bind_rows(
  top_gene_interpretation,
  top_pathway_interpretation
)

readr::write_tsv(
  biological_interpretation_table,
  here("results", "tables", "biological_interpretation_table.tsv")
)

# -----------------------------
# 14. Interpretation notes
# -----------------------------
n_sig <- nrow(sig_annotated)
n_up <- sum(sig_annotated$significant == "upregulated")
n_down <- sum(sig_annotated$significant == "downregulated")

top_symbols <- sig_annotated %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  pull(gene_symbol) %>%
  paste(collapse = ", ")

top_pathways <- pathway_summary %>%
  slice_head(n = 10) %>%
  mutate(line = paste0("- ", source, " / ", direction, ": ", pathway, " (adjusted p-value=", signif(adjusted_pvalue, 3), ")")) %>%
  pull(line)

interpretation_notes <- c(
  "# Gene annotation and pathway interpretation notes",
  "",
  "## Gene annotation",
  "",
  paste0(
    "DESeq2 results were annotated by mapping Ensembl gene identifiers to human gene symbols, Entrez IDs and gene names using org.Hs.eg.db."
  ),
  "",
  paste0(
    "A total of ", n_sig, " genes passed the threshold padj < 0.05 and absolute shrunken log2 fold change >= 1, including ",
    n_up, " upregulated and ", n_down, " downregulated genes in disease relative to control."
  ),
  "",
  paste0("The top-ranked annotated genes included: ", top_symbols, "."),
  "",
  "## Pathway-level interpretation",
  "",
  "GO Biological Process, Reactome and KEGG over-representation analyses were performed separately for upregulated and downregulated genes using the tested genes with Entrez IDs as the background universe.",
  "",
  "A ranked-list GO Biological Process analysis was also performed using a signed statistic based on log2 fold change direction and nominal p-value strength.",
  "",
  "## Top pathway findings",
  "",
  if (length(top_pathways) > 0) top_pathways else "- No enriched pathways passed the selected reporting threshold.",
  "",
  "## Interpretation caution",
  "",
  "These pathway results should be interpreted as hypothesis-generating biological summaries of the disease-control contrast. They do not establish cell-type causality or clinical mechanism without additional validation."
)

writeLines(
  interpretation_notes,
  here("results", "pathway_interpretation_notes.md")
)

message("Gene annotation and pathway analysis complete.")
message("Main outputs:")
message(" - results/tables/gene_annotation_table.tsv")
message(" - results/tables/deseq2_all_genes_annotated.tsv")
message(" - results/tables/deseq2_significant_genes_annotated.tsv")
message(" - results/tables/biological_interpretation_table.tsv")
message(" - results/tables/go_bp_enrichment_upregulated.tsv")
message(" - results/tables/reactome_enrichment_upregulated.tsv")
message(" - results/tables/kegg_enrichment_upregulated.tsv")
message(" - results/tables/gsea_go_bp_ranked_results.tsv")
message(" - results/figures/volcano_plot_deseq2_gene_symbols.png")
message(" - results/figures/top50_gene_heatmap_gene_symbols.png")
message(" - results/pathway_interpretation_notes.md")