# scripts/06_render_report.R
# Purpose: render the Quarto report and record final session information.

rm(list = ls())

library(here)

here::i_am("scripts/06_render_report.R")

if (!requireNamespace("quarto", quietly = TRUE)) {
  install.packages("quarto")
}

dir.create(here("report"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)

report_qmd <- here("report", "rna_seq_geo_analysis.qmd")

if (!file.exists(report_qmd)) {
  stop("Missing report/rna_seq_geo_analysis.qmd")
}

quarto::quarto_render(
  input = report_qmd,
  output_format = "html"
)

sink(here("results", "session_info.txt"))
sessionInfo()
sink()

if (requireNamespace("renv", quietly = TRUE)) {
  renv::snapshot(prompt = FALSE)
}

message("Report rendering complete.")
message("Main report:")
message(" - report/rna_seq_geo_analysis.html")
message("Session info:")
message(" - results/session_info.txt")