## ============================================================
## 02_fgsea_pathways_continuous_benzene.R
## Pathway-level analysis for continuous benzene model
## - Uses DESeq2 Wald statistic ("stat") as ranking metric for GSEA
## - Global GO BP, immune-focused GO BP subset, and Hallmark immune sets
##
## Input:  results/DESeq2_results_continuous_benzene_GAadjusted.csv
## Output: fgsea CSVs + barplot PNGs in results/
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
  library(stringr)
})

## -------------------------------
## 0) User-editable parameters
## -------------------------------

in_results_csv <- "results/DESeq2_results_continuous_benzene_GAadjusted.csv"
out_dir <- "results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

set.seed(123)

# fgsea parameters
min_size <- 10
max_size <- 500
nperm <- 10000

# plotting
top_n <- 15

# Hallmark immune sets to test
hallmark_immune_names <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

## -------------------------------
## 1) Helper functions
## -------------------------------

assert_file_exists <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
}

# Flatten fgsea leadingEdge list column for CSV export
flatten_fgsea <- function(fgsea_res) {
  fgsea_res %>%
    mutate(
      leadingEdge = vapply(
        leadingEdge,
        function(x) paste(x, collapse = ";"),
        FUN.VALUE = character(1)
      )
    ) %>%
    as.data.frame()
}

# Build ranked stats vector; resolves duplicate gene symbols
# Strategy: keep the entry with largest absolute stat per gene
build_ranked_stats <- function(res_tbl, gene_col = "gene", stat_col = "stat") {
  req <- c(gene_col, stat_col)
  if (!all(req %in% colnames(res_tbl))) {
    stop("Missing required columns: ", paste(setdiff(req, colnames(res_tbl)), collapse = ", "), call. = FALSE)
  }

  rank_df <- res_tbl %>%
    dplyr::select(gene = all_of(gene_col), stat = all_of(stat_col)) %>%
    filter(!is.na(gene), !is.na(stat)) %>%
    mutate(gene = as.character(gene), stat = as.numeric(stat)) %>%
    filter(is.finite(stat)) %>%
    group_by(gene) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(stat))

  ranked <- rank_df$stat
  names(ranked) <- rank_df$gene
  ranked
}

# Save fgsea outputs (CSV + barplot PNG)
save_fgsea_outputs <- function(fgsea_res, prefix, top_n = 15, title = NULL, short_fun = NULL) {
  fgsea_res <- fgsea_res %>% arrange(padj)

  out_csv <- file.path(out_dir, paste0(prefix, ".csv"))
  write.csv(flatten_fgsea(fgsea_res), file = out_csv, row.names = FALSE)

  message("Wrote: ", out_csv)

  # Plot top_n by padj (if available)
  plot_df <- fgsea_res %>%
    filter(!is.na(padj)) %>%
    slice_head(n = top_n)

  if (nrow(plot_df) == 0) {
    message("No non-NA padj results for plotting: ", prefix)
    return(invisible(NULL))
  }

  plot_df <- plot_df %>%
    mutate(
      pathway_short = if (!is.null(short_fun)) short_fun(pathway) else pathway
    )

  p <- ggplot(plot_df, aes(x = NES, y = reorder(pathway_short, NES), fill = NES)) +
    geom_col() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      title = ifelse(is.null(title), prefix, title)
    ) +
    theme_classic(base_size = 12)

  out_png <- file.path(out_dir, paste0(prefix, "_top", top_n, "_bar.png"))
  ggsave(out_png, p, width = 10, height = 7, dpi = 300)
  message("Wrote: ", out_png)

  print(p)
  invisible(list(csv = out_csv, png = out_png))
}

## -------------------------------
## 2) Load DESeq2 results table
## -------------------------------

assert_file_exists(in_results_csv)

res_adduct_tbl <- read.csv(in_results_csv, stringsAsFactors = FALSE)

needed_cols <- c("gene", "stat", "log2FoldChange", "pvalue")
if (!all(needed_cols %in% colnames(res_adduct_tbl))) {
  stop(
    "Input results CSV must contain columns: ",
    paste(needed_cols, collapse = ", "),
    call. = FALSE
  )
}

message("Loaded DESeq2 results: ", nrow(res_adduct_tbl), " genes")
message("Input: ", in_results_csv)

## -------------------------------
## 3) Build ranked statistics
## -------------------------------

ranked_stats <- build_ranked_stats(res_adduct_tbl, gene_col = "gene", stat_col = "stat")

message("Ranked stats length (unique genes with finite stat): ", length(ranked_stats))

if (length(ranked_stats) < 50) {
  stop("Too few ranked genes after filtering. Check input results table.", call. = FALSE)
}

## ============================================================
## PART 1 - GLOBAL GO BP GSEA
## ============================================================

message("\nRunning global GO BP fgsea...")

msig_go_bp <- msigdbr(
  species     = "Homo sapiens",
  category    = "C5",
  subcategory = "GO:BP"
)

go_bp_sets <- split(msig_go_bp$gene_symbol, msig_go_bp$gs_name)

common_genes_go <- intersect(names(ranked_stats), unique(msig_go_bp$gene_symbol))
stats_go <- ranked_stats[common_genes_go]

message("Overlap genes (GO BP): ", length(stats_go))

if (length(stats_go) < 50) {
  stop("Too few overlapping genes for GO BP fgsea. Check gene symbols.", call. = FALSE)
}

fgsea_go <- fgsea(
  pathways = go_bp_sets,
  stats    = stats_go,
  minSize  = min_size,
  maxSize  = max_size,
  nperm    = nperm
)

save_fgsea_outputs(
  fgsea_res = fgsea_go,
  prefix    = "fgsea_GO_BP_continuous_benzene",
  top_n     = top_n,
  title     = "Top GO BP pathways enriched with continuous benzene (global)",
  short_fun = function(x) {
    x %>%
      gsub("^GOBP_", "", .) %>%
      gsub("_", " ", .)
  }
)

## ============================================================
## PART 2 - IMMUNE FOCUSED GSEA
##   2A) Immune-related GO BP subset
##   2B) Hallmark immune sets
## ============================================================

## 2A) Immune-related GO BP terms (simple keyword filter)
message("\nRunning immune-filtered GO BP fgsea...")

go_immune_bp <- msig_go_bp %>%
  filter(
    grepl("IMMUNE", gs_name, ignore.case = TRUE) |
      grepl("INFLAMMATORY", gs_name, ignore.case = TRUE) |
      grepl("CYTOKINE", gs_name, ignore.case = TRUE)
  )

immune_bp_sets <- split(go_immune_bp$gene_symbol, go_immune_bp$gs_name)

common_genes_immune <- intersect(names(ranked_stats), unique(go_immune_bp$gene_symbol))
stats_immune_bp <- ranked_stats[common_genes_immune]

message("Immune GO BP sets: ", length(immune_bp_sets))
message("Overlap genes (immune GO BP): ", length(stats_immune_bp))

if (length(stats_immune_bp) >= 50 && length(immune_bp_sets) >= 5) {
  fgsea_immune_bp <- fgsea(
    pathways = immune_bp_sets,
    stats    = stats_immune_bp,
    minSize  = min_size,
    maxSize  = max_size,
    nperm    = nperm
  )

  save_fgsea_outputs(
    fgsea_res = fgsea_immune_bp,
    prefix    = "fgsea_GO_BP_immune_continuous_benzene",
    top_n     = top_n,
    title     = "Top immune GO BP pathways enriched with continuous benzene",
    short_fun = function(x) {
      x %>%
        gsub("^GOBP_", "", .) %>%
        gsub("_", " ", .)
    }
  )
} else {
  message("Skipping immune GO BP fgsea: too few overlapping genes or sets.")
}

## 2B) Hallmark immune sets
message("\nRunning Hallmark immune fgsea...")

msig_h <- msigdbr(species = "Homo sapiens", category = "H")

hallmark_immune <- msig_h %>%
  filter(gs_name %in% hallmark_immune_names)

hallmark_immune_sets <- split(hallmark_immune$gene_symbol, hallmark_immune$gs_name)

common_genes_hallmark <- intersect(names(ranked_stats), unique(hallmark_immune$gene_symbol))
stats_hallmark <- ranked_stats[common_genes_hallmark]

message("Hallmark immune sets: ", length(hallmark_immune_sets))
message("Overlap genes (Hallmark immune): ", length(stats_hallmark))

if (length(stats_hallmark) >= 50 && length(hallmark_immune_sets) >= 2) {
  fgsea_hallmark <- fgsea(
    pathways = hallmark_immune_sets,
    stats    = stats_hallmark,
    minSize  = min_size,
    maxSize  = max_size,
    nperm    = nperm
  )

  # Plot all hallmark immune sets (small number)
  save_fgsea_outputs(
    fgsea_res = fgsea_hallmark,
    prefix    = "fgsea_Hallmark_immune_continuous_benzene",
    top_n     = length(hallmark_immune_sets),
    title     = "Hallmark immune pathways enriched with continuous benzene",
    short_fun = function(x) {
      x %>%
        gsub("^HALLMARK_", "", .) %>%
        gsub("_", " ", .)
    }
  )
} else {
  message("Skipping Hallmark immune fgsea: too few overlapping genes or sets.")
}

## -------------------------------
## 4) Save session info
## -------------------------------

out_session <- file.path(out_dir, "sessionInfo_fgsea.txt")
writeLines(capture.output(sessionInfo()), con = out_session)
message("\nWrote: ", out_session)

message("Done.")
