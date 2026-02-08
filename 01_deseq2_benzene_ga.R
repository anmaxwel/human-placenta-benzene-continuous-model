## ============================================================
## 01_deseq2_benzene_ga.R
## Human placenta bulk RNA-seq + benzene adducts + EMR covariates
## DESeq2: continuous exposure (z-scored adducts), adjusted for sex + GA
## Optional exploratory High vs Low split (median)
## Outputs: CSV tables, PCA plot, session info
## ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(stringr)
})

## -------------------------------
## 0) User-editable parameters
## -------------------------------

# Input files (recommend using relative paths in a repo)
rna_file  <- "data/raw/BenezeneAll.xlsx"
benz_file <- "data/raw/PPIH_plasma_Benzene_adduct_countsKG.xlsx"
emr_file  <- "data/raw/PPIH_EMR_data_deidentified.xlsx"

# Excel sheet names
rna_sheet  <- "BenzeneAllw1179"
benz_sheet <- "CLEAR_adducts"
emr_sheet  <- "Sheet1"

# Output directory
out_dir <- "results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Prefilter parameters (DESeq2 best practice for sparse genes)
min_count_per_sample <- 10
min_samples_with_min <- 3

# Exploratory median split flag
run_high_low <- TRUE

## -------------------------------
## 1) Helper functions
## -------------------------------

# Robust gestational age parser to numeric weeks.
# Handles common variants: "40w3d", "40 w 3 d", "40w", "40+3", "40.4"
parse_ga_to_weeks <- function(x) {
  x0 <- as.character(x)
  x1 <- str_trim(tolower(x0))

  # remove spaces
  x1 <- str_replace_all(x1, "\\s+", "")

  # Case A: 40w3d
  m <- str_match(x1, "^(\\d+)w(\\d+)d$")
  out <- suppressWarnings(as.numeric(m[, 2]) + as.numeric(m[, 3]) / 7)

  # Case B: 40w
  m2 <- str_match(x1, "^(\\d+)w$")
  out2 <- suppressWarnings(as.numeric(m2[, 2]))

  # Case C: 40+3
  m3 <- str_match(x1, "^(\\d+)\\+(\\d+)$")
  out3 <- suppressWarnings(as.numeric(m3[, 2]) + as.numeric(m3[, 3]) / 7)

  # Case D: 40.4 (assume already weeks as decimal)
  out4 <- suppressWarnings(as.numeric(x1))
  # But only trust decimal if it contains a dot and is plausible
  out4[!str_detect(x1, "\\.")] <- NA_real_

  # Combine with priority: A > B > C > D
  out_final <- out
  out_final[is.na(out_final)] <- out2[is.na(out_final)]
  out_final[is.na(out_final)] <- out3[is.na(out_final)]
  out_final[is.na(out_final)] <- out4[is.na(out_final)]

  return(out_final)
}

# Safe file existence checks
assert_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
}

# If duplicate gene symbols exist, collapse by summing counts
collapse_duplicate_genes <- function(count_mat) {
  dup <- duplicated(rownames(count_mat))
  if (!any(dup)) return(count_mat)

  message("Warning: duplicate gene symbols detected. Collapsing by sum.")
  df <- as.data.frame(count_mat) %>%
    rownames_to_column("gene") %>%
    group_by(gene) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  out <- as.matrix(df[, -1, drop = FALSE])
  rownames(out) <- df$gene
  storage.mode(out) <- "integer"
  return(out)
}

## -------------------------------
## 2) Validate inputs
## -------------------------------

assert_file_exists(rna_file)
assert_file_exists(benz_file)
assert_file_exists(emr_file)

message("Inputs:")
message("  RNA counts:  ", rna_file, " (sheet: ", rna_sheet, ")")
message("  Adducts:     ", benz_file, " (sheet: ", benz_sheet, ")")
message("  EMR:         ", emr_file, " (sheet: ", emr_sheet, ")")
message("  Output dir:  ", out_dir)

## -------------------------------
## 3) Read RNA-seq count matrix
## -------------------------------

rna_df <- readxl::read_excel(rna_file, sheet = rna_sheet) %>%
  as.data.frame()

if (ncol(rna_df) < 2) stop("RNA sheet has <2 columns. Expect gene + samples.", call. = FALSE)

# First column must be gene symbols
if (colnames(rna_df)[1] != "gene") {
  stop("First column of RNA sheet must be named 'gene'. Found: ", colnames(rna_df)[1], call. = FALSE)
}

# Set rownames and convert to integer matrix
rownames(rna_df) <- rna_df$gene
rna_df$gene <- NULL

count_mat <- as.matrix(rna_df)
mode(count_mat) <- "numeric"

# Coerce NA to 0 and round to integer counts
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)

# Sanity checks
if (any(count_mat < 0, na.rm = TRUE)) stop("Count matrix contains negative values.", call. = FALSE)

# Ensure integer storage for DESeq2
storage.mode(count_mat) <- "integer"

# Collapse duplicate gene symbols if present
count_mat <- collapse_duplicate_genes(count_mat)

samples_in_counts <- colnames(count_mat)
message("Counts matrix loaded: ", nrow(count_mat), " genes x ", ncol(count_mat), " samples")

## -------------------------------
## 4) Read benzene adduct data
## -------------------------------

benz_df <- readxl::read_excel(benz_file, sheet = benz_sheet) %>%
  dplyr::rename(
    sample_id    = sample,
    adduct_count = `adduct count`
  ) %>%
  dplyr::mutate(
    sample_id    = as.character(sample_id),
    adduct_count = as.numeric(adduct_count)
  )

if (!all(c("sample_id", "adduct_count") %in% colnames(benz_df))) {
  stop("Adduct sheet must contain columns: 'sample' and 'adduct count'.", call. = FALSE)
}

message("Adduct rows: ", nrow(benz_df), " | Unique samples: ", length(unique(benz_df$sample_id)))

## -------------------------------
## 5) Read EMR data
## -------------------------------

emr_df <- readxl::read_excel(emr_file, sheet = emr_sheet) %>%
  dplyr::rename(
    redcap_id       = `REDCap ID`,
    sex             = `Infant Sex`,
    birth_weight_g  = `Birth Weight (g)`,
    birth_length_cm = `Birth Length (cm)`,
    ga              = `Gestational Age`,
    delivery        = `Mode of Delivery`
  ) %>%
  dplyr::mutate(
    redcap_id = as.character(redcap_id),
    sample_id = paste0("s", redcap_id)
  )

message("EMR rows: ", nrow(emr_df))

## -------------------------------
## 6) Build sample metadata
##    Intersection: counts + adducts + EMR
## -------------------------------

meta <- benz_df %>%
  # Only keep adduct samples that exist in RNA count matrix
  dplyr::filter(sample_id %in% samples_in_counts) %>%
  dplyr::left_join(emr_df, by = "sample_id") %>%
  # Require sex + gestational age for the GA-adjusted model
  dplyr::filter(!is.na(sex), !is.na(ga))

if (nrow(meta) == 0) {
  stop("No overlapping samples between counts, adducts, and EMR (with sex + GA present).", call. = FALSE)
}

# Remove duplicated sample_id rows if present (keep first, but warn)
if (any(duplicated(meta$sample_id))) {
  message("Warning: duplicated sample_id in metadata. Keeping first instance per sample.")
  meta <- meta %>% dplyr::distinct(sample_id, .keep_all = TRUE)
}

## -------------------------------
## 7) Parse gestational age to numeric weeks
## -------------------------------

meta <- meta %>%
  dplyr::mutate(
    ga_weeks_num = parse_ga_to_weeks(ga)
  )

# Fail loudly if parsing fails
if (any(is.na(meta$ga_weeks_num))) {
  bad <- meta %>% dplyr::filter(is.na(ga_weeks_num)) %>% dplyr::select(sample_id, ga)
  message("Failed to parse GA for these rows:")
  print(bad)
  stop("Gestational age parsing produced NA values. Fix GA formatting or extend parser.", call. = FALSE)
}

message("Gestational age (weeks) summary:")
print(summary(meta$ga_weeks_num))

## -------------------------------
## 8) Define covariates and align counts to metadata
## -------------------------------

meta <- meta %>%
  dplyr::mutate(
    # Normalize sex values and factorize
    sex = factor(str_trim(as.character(sex))),
    adduct_count = as.numeric(adduct_count),
    adduct_z     = as.numeric(scale(adduct_count)),
    # Exploratory median split (information-losing; treat as exploratory)
    adduct_group = ifelse(
      adduct_count <= median(adduct_count, na.rm = TRUE),
      "Low", "High"
    ),
    adduct_group = factor(adduct_group, levels = c("Low", "High"))
  )

# Set rownames for DESeq2
rownames(meta) <- meta$sample_id

# Align count matrix columns to metadata sample order
meta$sample_id <- as.character(meta$sample_id)

# Ensure all meta sample_ids exist in counts
missing_in_counts <- setdiff(meta$sample_id, colnames(count_mat))
if (length(missing_in_counts) > 0) {
  stop("Metadata contains sample_ids not in count matrix: ",
       paste(missing_in_counts, collapse = ", "), call. = FALSE)
}

count_mat <- count_mat[, meta$sample_id, drop = FALSE]

# Final alignment check: column order equals rownames(meta)
stopifnot(all(colnames(count_mat) == rownames(meta)))

message("Final analysis samples: ", ncol(count_mat))
message("Final metadata columns: ", paste(colnames(meta), collapse = ", "))

## -------------------------------
## 9) DESeq2: continuous exposure model
##    design: ~ sex + ga_weeks_num + adduct_z
## -------------------------------

dds_cont <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ sex + ga_weeks_num + adduct_z
)

# Prefilter low-count genes
keep <- rowSums(counts(dds_cont) >= min_count_per_sample) >= min_samples_with_min
dds_cont <- dds_cont[keep, ]

message("Genes kept after prefiltering: ", nrow(dds_cont))

# Run DESeq2
dds_cont <- DESeq(dds_cont)

message("Available coefficients (continuous model):")
print(resultsNames(dds_cont))

## -------------------------------
## 10) Results: benzene association (continuous)
## -------------------------------

# Coefficient name should match "adduct_z" exactly
# If it does not, inspect resultsNames(dds_cont) and update below
if (!("adduct_z" %in% resultsNames(dds_cont))) {
  stop("Coefficient 'adduct_z' not found in resultsNames(dds_cont).", call. = FALSE)
}

res_adduct <- results(dds_cont, name = "adduct_z")
message("Continuous model summary:")
print(summary(res_adduct))

res_adduct_tbl <- as.data.frame(res_adduct) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

out_cont_csv <- file.path(out_dir, "DESeq2_results_continuous_benzene_GAadjusted.csv")
write.csv(res_adduct_tbl, file = out_cont_csv, row.names = FALSE)

message("Wrote: ", out_cont_csv)
message("Top benzene-associated genes (continuous model, GA-adjusted):")
print(head(res_adduct_tbl, 20))

## -------------------------------
## 11) Optional: exploratory High vs Low model (median split)
##    design: ~ sex + ga_weeks_num + adduct_group
## -------------------------------

if (isTRUE(run_high_low)) {

  dds_group <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = meta,
    design    = ~ sex + ga_weeks_num + adduct_group
  )

  # Use the same keep filter for comparability
  dds_group <- dds_group[keep, ]
  dds_group <- DESeq(dds_group)

  message("Available coefficients (High vs Low model):")
  print(resultsNames(dds_group))

  res_high_vs_low <- results(
    dds_group,
    contrast = c("adduct_group", "High", "Low")
  )

  message("High vs Low model summary:")
  print(summary(res_high_vs_low))

  res_high_vs_low_tbl <- as.data.frame(res_high_vs_low) %>%
    rownames_to_column("gene") %>%
    arrange(padj)

  out_hl_csv <- file.path(out_dir, "DESeq2_results_high_vs_low_benzene_GAadjusted.csv")
  write.csv(res_high_vs_low_tbl, file = out_hl_csv, row.names = FALSE)

  message("Wrote: ", out_hl_csv)
  message("Top genes High vs Low (exploratory, GA-adjusted):")
  print(head(res_high_vs_low_tbl, 20))
}

## -------------------------------
## 12) Quick PCA for QC
## -------------------------------

vsd <- vst(dds_cont, blind = TRUE)

p_pca <- plotPCA(vsd, intgroup = c("sex", "adduct_group")) +
  ggtitle("PCA: VST counts (colored by sex and adduct group)")

print(p_pca)

out_pca <- file.path(out_dir, "PCA_vst_sex_adductGroup.png")
ggsave(out_pca, p_pca, width = 9, height = 6, dpi = 300)
message("Wrote: ", out_pca)

## -------------------------------
## 13) Save session info for reproducibility
## -------------------------------

out_session <- file.path(out_dir, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), con = out_session)
message("Wrote: ", out_session)

message("Done.")
