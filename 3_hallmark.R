## ============================================================
## 03_hallmark_ifn_alpha_gene_view.R
## Hallmark interferon gene-level view (GA-adjusted DESeq2 model)
## Focus: HALLMARK_INTERFERON_ALPHA_RESPONSE (plus optional IFN gamma)
##
## Inputs (produced by 01_... script):
## - results/DESeq2_results_continuous_benzene_GAadjusted.csv
## - results/vsd_continuous_benzene_GAadjusted.rds
##
## Outputs (written to results/):
## - Hallmark_IFN_alpha_DESeq2_GAadjusted.csv
## - barplot_IFNalpha_log2FC_topN.png
## - heatmap_IFNalpha_topN.png
## - scatter_IFNalpha_signature_vs_adductz.png
## - lm_IFNalpha_signature_summary.txt
## - sessionInfo_ifnalpha.txt
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(msigdbr)
  library(ggplot2)
  library(pheatmap)
  library(stringr)
  library(SummarizedExperiment)  # for assay() and colData()
})

## -------------------------------
## 0) User-editable parameters
## -------------------------------

in_results_csv <- "results/DESeq2_results_continuous_benzene_GAadjusted.csv"
in_vsd_rds     <- "results/vsd_continuous_benzene_GAadjusted.rds"

out_dir <- "results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# How many IFN alpha genes to show in barplot and heatmap (ranked by raw p-value)
top_n_ifna_genes <- 30

# Optional: also compute IFN gamma list/table (no extra plots by default)
also_extract_ifng <- TRUE

## -------------------------------
## 1) Helper functions
## -------------------------------

assert_file_exists <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
}

safe_neglog10 <- function(p) {
  # Replace NA or non-positive p with 1 for plotting metrics
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

# Pull expression and metadata from the VST object
vsd_mat <- SummarizedExperiment::assay(vsd)
meta_df <- as.data.frame(SummarizedExperiment::colData(vsd))

# Ensure required covariates exist (these should come from DESeq2 colData)
req_meta <- c("adduct_z", "ga_weeks_num", "sex")
if (!all(req_meta %in% colnames(meta_df))) {
  stop(
    "VSD colData is missing required columns: ",
    paste(setdiff(req_meta, colnames(meta_df)), collapse = ", "),
    call. = FALSE
  )
}

# Add sample_id column if it does not exist
if (!("sample_id" %in% colnames(meta_df))) {
  meta_df <- meta_df %>% tibble::rownames_to_column("sample_id")
}

message("Loaded DESeq2 results: ", nrow(res_adduct_tbl), " genes")
message("Loaded VSD matrix: ", nrow(vsd_mat), " genes x ", ncol(vsd_mat), " samples")

## -------------------------------
## 3) Pull Hallmark IFN gene sets
## -------------------------------

msig_h <- msigdbr(species = "Homo sapiens", category = "H")

ifna_set <- msig_h %>%
  filter(gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE")

ifng_set <- msig_h %>%
  filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE")

ifna_genes <- unique(ifna_set$gene_symbol)
ifng_genes <- unique(ifng_set$gene_symbol)

message("MSigDB IFN alpha genes: ", length(ifna_genes))
message("MSigDB IFN gamma genes: ", length(ifng_genes))

## -------------------------------
## 4) Intersect with DESeq2 results (GA-adjusted model)
## -------------------------------

ifna_res <- res_adduct_tbl %>%
  filter(gene %in% ifna_genes) %>%
  mutate(
    pvalue_raw  = as.numeric(pvalue),
    neg_log10_p = safe_neglog10(pvalue_raw),
    direction   = ifelse(log2FoldChange > 0, "Up with benzene", "Down with benzene")
  ) %>%
  arrange(pvalue_raw)

out_ifna_csv <- file.path(out_dir, "Hallmark_IFN_alpha_DESeq2_GAadjusted.csv")
write.csv(ifna_res, file = out_ifna_csv, row.names = FALSE)
message("Wrote: ", out_ifna_csv)

message("Top IFN alpha genes by raw p-value (GA-adjusted benzene effect):")
print(ifna_res %>% select(gene, log2FoldChange, pvalue_raw, padj, stat) %>% head(15))

# Optional: IFN gamma table
if (isTRUE(also_extract_ifng)) {
  ifng_res <- res_adduct_tbl %>%
    filter(gene %in% ifng_genes) %>%
    mutate(
      pvalue_raw  = as.numeric(pvalue),
      neg_log10_p = safe_neglog10(pvalue_raw),
      direction   = ifelse(log2FoldChange > 0, "Up with benzene", "Down with benzene")
    ) %>%
    arrange(pvalue_raw)

  out_ifng_csv <- file.path(out_dir, "Hallmark_IFN_gamma_DESeq2_GAadjusted.csv")
  write.csv(ifng_res, file = out_ifng_csv, row.names = FALSE)
  message("Wrote: ", out_ifng_csv)
}

## Pick top N IFN alpha genes for plotting (by raw p-value)
ifna_res_plot <- ifna_res %>% slice_head(n = top_n_ifna_genes)

## -------------------------------
## 5) Barplot: IFN alpha gene log2FC (GA-adjusted)
## -------------------------------

if (nrow(ifna_res_plot) > 0) {

  gg_ifna_fc <- ggplot(
    ifna_res_plot,
    aes(x = reorder(gene, log2FoldChange), y = log2FoldChange, fill = direction)
  ) +
    geom_col() +
    coord_flip() +
    labs(
      x = NULL,
      y = "Log2 fold change per 1 SD benzene (adjusted for sex + GA)",
      title = paste0("Hallmark IFN alpha response genes (top ", nrow(ifna_res_plot), " by raw p)")
    ) +
    theme_classic(base_size = 12)

  out_fc_png <- file.path(out_dir, paste0("barplot_IFNalpha_log2FC_top", top_n_ifna_genes, ".png"))
  ggsave(out_fc_png, gg_ifna_fc, width = 9, height = 7, dpi = 300)
  message("Wrote: ", out_fc_png)

  print(gg_ifna_fc)

} else {
  message("No IFN alpha genes found in DESeq2 results for barplot.")
}

## -------------------------------
## 6) Heatmap: IFN alpha genes across samples (VST)
## -------------------------------

top_ifna_genes <- ifna_res_plot$gene
ifna_genes_in_vsd <- intersect(top_ifna_genes, rownames(vsd_mat))

# Annotation columns (only include what exists)
ann_cols_wanted <- c("sex", "adduct_group", "adduct_z", "ga_weeks_num")
ann_cols_present <- intersect(ann_cols_wanted, colnames(meta_df))

if (length(ifna_genes_in_vsd) > 1) {

  ifna_mat <- vsd_mat[ifna_genes_in_vsd, , drop = FALSE]

  # Row z-score (gene-wise standardization for heatmap visibility)
  ifna_mat_z <- t(scale(t(ifna_mat)))

  ann_col <- meta_df %>%
    dplyr::select(all_of(c("sample_id", ann_cols_present))) %>%
    tibble::column_to_rownames("sample_id")

  out_hm_png <- file.path(out_dir, paste0("heatmap_IFNalpha_top", length(ifna_genes_in_vsd), ".png"))
  png(out_hm_png, width = 1100, height = 900, res = 150)

  pheatmap(
    ifna_mat_z,
    annotation_col = ann_col,
    show_rownames  = TRUE,
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    fontsize_row   = 7,
    main = paste0("Hallmark IFN alpha genes (top ", length(ifna_genes_in_vsd), " by raw p)")
  )

  dev.off()
  message("Wrote: ", out_hm_png)

} else {
  message("Not enough IFN alpha genes in VST matrix to make a heatmap.")
}

## -------------------------------
## 7) Sample-level IFN alpha signature vs benzene
## -------------------------------

if (length(ifna_genes_in_vsd) > 1) {

  # IFN alpha genes x samples (VST)
  ifna_mat <- vsd_mat[ifna_genes_in_vsd, , drop = FALSE]

  # Row z-score then average per sample
  ifna_mat_z <- t(scale(t(ifna_mat)))
  ifna_score <- colMeans(ifna_mat_z, na.rm = TRUE)

  # Build a tidy data frame for plotting and modeling
  ifna_score_df <- meta_df %>%
    filter(sample_id %in% colnames(ifna_mat)) %>%
    mutate(
      IFNA_score = ifna_score[match(sample_id, colnames(ifna_mat))]
    ) %>%
    select(sample_id, sex, adduct_z, ga_weeks_num, IFNA_score) %>%
    filter(!is.na(IFNA_score), !is.na(adduct_z), !is.na(ga_weeks_num), !is.na(sex))

  message("IFN alpha score rows: ", nrow(ifna_score_df))

  if (nrow(ifna_score_df) > 0) {

    # Simple regression to show association adjusted for sex + GA
    fit <- lm(IFNA_score ~ adduct_z + ga_weeks_num + sex, data = ifna_score_df)

    out_lm <- file.path(out_dir, "lm_IFNalpha_signature_summary.txt")
    writeLines(capture.output(summary(fit)), con = out_lm)
    message("Wrote: ", out_lm)

    gg_ifna_score <- ggplot(
      ifna_score_df,
      aes(x = adduct_z, y = IFNA_score, color = sex)
    ) +
      geom_point(size = 3) +
      geom_smooth(method = "lm", se = FALSE) +
      labs(
        x = "Benzene adduct_z (per SD)",
        y = "Mean z-scored VST expression (IFN alpha genes)",
        title = "Hallmark IFN alpha signature vs benzene"
      ) +
      theme_classic(base_size = 14)

    out_scatter <- file.path(out_dir, "scatter_IFNalpha_signature_vs_adductz.png")
    ggsave(out_scatter, gg_ifna_score, width = 7.5, height = 6, dpi = 300)
    message("Wrote: ", out_scatter)

    print(gg_ifna_score)

  } else {
    message("No complete cases for IFN alpha signature vs benzene.")
  }

} else {
  message("Not enough IFN alpha genes to compute a signature score.")
}

## -------------------------------
## 8) Save session info
## -------------------------------

out_session <- file.path(out_dir, "sessionInfo_ifnalpha.txt")
writeLines(capture.output(sessionInfo()), con = out_session)
message("Wrote: ", out_session)

message("Done.")
