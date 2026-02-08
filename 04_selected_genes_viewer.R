## ============================================================
## 04_selected_genes_viewer.R
## Quick viewer for specific genes in GA-adjusted continuous benzene model
## - Barplot of log2FC per 1 SD benzene (Up with benzene / Down)
## - Expression vs benzene scatter for each gene (VST)
##
## Inputs (produced by 01_... script):
## - results/DESeq2_results_continuous_benzene_GAadjusted.csv
## - results/vsd_continuous_benzene_GAadjusted.rds
##
## Outputs (written to results/):
## - barplot_selectedGenes_log2FC.png
## - scatter_<GENE>_expr_vs_adductz.png (one per gene found)
## - selectedGenes_DESeq2_subset.csv
## - sessionInfo_selected_genes.txt
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(SummarizedExperiment)
})

## -------------------------------
## 0) User-editable parameters
## -------------------------------

in_results_csv <- "results/DESeq2_results_continuous_benzene_GAadjusted.csv"
in_vsd_rds     <- "results/vsd_continuous_benzene_GAadjusted.rds"

out_dir <- "results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

genes_of_interest <- c("TNFAIP3", "KCNIP3")  # edit this

## -------------------------------
## 1) Helpers
## -------------------------------

assert_file_exists <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
}

safe_neglog10 <- function(p) {
  p2 <- ifelse(is.na(p) | p <= 0, 1, p)
  -log10(p2)
}

## -------------------------------
## 2) Load inputs
## -------------------------------

assert_file_exists(in_results_csv)
assert_file_exists(in_vsd_rds)

res_adduct_tbl <- read.csv(in_results_csv, stringsAsFactors = FALSE)
vsd <- readRDS(in_vsd_rds)

needed_cols <- c("gene", "stat", "log2FoldChange", "pvalue", "padj")
if (!all(needed_cols %in% colnames(res_adduct_tbl))) {
  stop(
    "DESeq2 results CSV must contain columns: ",
    paste(needed_cols, collapse = ", "),
    call. = FALSE
  )
}

vsd_mat <- SummarizedExperiment::assay(vsd)
meta_df <- as.data.frame(SummarizedExperiment::colData(vsd))

req_meta <- c("sex", "adduct_z", "ga_weeks_num")
if (!all(req_meta %in% colnames(meta_df))) {
  stop(
    "VSD colData is missing required columns: ",
    paste(setdiff(req_meta, colnames(meta_df)), collapse = ", "),
    call. = FALSE
  )
}

# Ensure meta_df has sample_id and is aligned to vsd columns
if (!("sample_id" %in% colnames(meta_df))) {
  meta_df <- meta_df %>% tibble::rownames_to_column("sample_id")
}

# Reorder meta_df to match vsd columns explicitly
meta_df <- meta_df %>%
  filter(sample_id %in% colnames(vsd_mat)) %>%
  mutate(sample_id = as.character(sample_id))

meta_df <- meta_df[match(colnames(vsd_mat), meta_df$sample_id), , drop = FALSE]

if (!all(meta_df$sample_id == colnames(vsd_mat))) {
  stop("Metadata sample_id does not align to VSD column names after matching.", call. = FALSE)
}

message("Loaded DESeq2 results: ", nrow(res_adduct_tbl), " genes")
message("Loaded VSD matrix: ", nrow(vsd_mat), " genes x ", ncol(vsd_mat), " samples")
message("Genes requested: ", paste(genes_of_interest, collapse = ", "))

## -------------------------------
## 3) Benzene effect summary for selected genes
## -------------------------------

gene_res <- res_adduct_tbl %>%
  filter(gene %in% genes_of_interest) %>%
  mutate(
    direction   = ifelse(log2FoldChange > 0, "Up with benzene", "Down with benzene"),
    pvalue_raw  = as.numeric(pvalue),
    neg_log10_p = safe_neglog10(pvalue_raw),
    p_label     = ifelse(is.na(pvalue_raw), "p = NA", paste0("p = ", signif(pvalue_raw, 2))),
    # hjust for coord_flip label placement
    hjust_lbl   = ifelse(log2FoldChange >= 0, -0.05, 1.05)
  ) %>%
  arrange(match(gene, genes_of_interest))

if (nrow(gene_res) == 0) {
  stop("None of the requested genes were found in the DESeq2 results table.", call. = FALSE)
}

message("Benzene coefficients for requested genes (GA-adjusted):")
print(gene_res %>% select(gene, log2FoldChange, pvalue_raw, padj, stat))

out_subset <- file.path(out_dir, "selectedGenes_DESeq2_subset.csv")
write.csv(gene_res, file = out_subset, row.names = FALSE)
message("Wrote: ", out_subset)

## Barplot of log2FC per 1 SD benzene
gg_gene_fc <- ggplot(
  gene_res,
  aes(x = reorder(gene, log2FoldChange), y = log2FoldChange, fill = direction)
) +
  geom_col() +
  geom_text(aes(label = p_label, hjust = hjust_lbl), size = 3) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Log2 fold change per 1 SD benzene (adjusted for sex + GA)",
    title = "Selected genes: benzene association (GA-adjusted)"
  ) +
  theme_classic(base_size = 12)

print(gg_gene_fc)

out_fc_png <- file.path(out_dir, "barplot_selectedGenes_log2FC.png")
ggsave(out_fc_png, gg_gene_fc, width = 8, height = 4.8, dpi = 300)
message("Wrote: ", out_fc_png)

## -------------------------------
## 4) Expression vs benzene scatterplots (VST)
## -------------------------------

genes_in_vsd <- intersect(gene_res$gene, rownames(vsd_mat))

if (length(genes_in_vsd) == 0) {
  message("Selected genes not found in VST matrix, skipping scatterplots.")
} else {

  plot_gene_vs_benzene <- function(gene_symbol) {
    expr <- as.numeric(vsd_mat[gene_symbol, ])

    df <- data.frame(
      sample_id    = meta_df$sample_id,
      expr         = expr,
      adduct_z     = meta_df$adduct_z,
      ga_weeks_num = meta_df$ga_weeks_num,
      sex          = meta_df$sex
    ) %>%
      filter(!is.na(expr), !is.na(adduct_z), !is.na(sex))

    ggplot(df, aes(x = adduct_z, y = expr, color = sex)) +
      geom_point(size = 3) +
      geom_smooth(method = "lm", se = FALSE) +
      labs(
        x = "Benzene adduct_z (per SD)",
        y = paste0("VST expression: ", gene_symbol),
        title = paste0(gene_symbol, " expression vs benzene")
      ) +
      theme_classic(base_size = 14)
  }

  for (g in genes_in_vsd) {
    p <- plot_gene_vs_benzene(g)
    print(p)

    out_scatter <- file.path(out_dir, paste0("scatter_", g, "_expr_vs_adductz.png"))
    ggsave(out_scatter, p, width = 7.5, height = 6, dpi = 300)
    message("Wrote: ", out_scatter)
  }
}

## -------------------------------
## 5) Session info
## -------------------------------

out_session <- file.path(out_dir, "sessionInfo_selected_genes.txt")
writeLines(capture.output(sessionInfo()), con = out_session)
message("Wrote: ", out_session)

message("Done.")
