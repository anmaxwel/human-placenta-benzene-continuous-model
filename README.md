# Human placenta benzene RNA-seq (continuous exposure, GA-adjusted)

This repository demonstrates a reproducible R workflow for differential expression and pathway analysis of human placental bulk RNA-seq data in relation to **continuous benzene exposure**, quantified by plasma adduct levels and modeled as a z-scored covariate. The primary DESeq2 model adjusts for **infant sex** and **gestational age (numeric weeks)**.

This public repository contains **code and synthetic/example data only**. Real human RNA-seq counts and clinical/EMR data are **not** included due to IRB and data-use restrictions.

---

## What this repo demonstrates

- Reproducible analysis structure with numbered scripts
- Continuous exposure modeling in DESeq2 using adduct z-scores
- Covariate adjustment for gestational age and sex
- Exportable results tables and QC plots
- Pathway-level analysis using fgsea (GO BP global, immune-filtered GO BP, Hallmark immune)
- Gene-level view for Hallmark IFN alpha response
- Quick gene viewer for targeted genes (example: TNFAIP3, KCNIP3)

---

## Human subjects and data policy

- Do not commit any real human data to this repository.
- Keep real inputs in `data/raw/` (ignored by git).
- If you want the repo runnable publicly, add synthetic example files in `data/example/` that mimic the same column names and formats.

Suggested `.gitignore` rules (safe default):

~~~gitignore
data/raw/
results/
*.xlsx
*.csv
*.tsv
*.rds
~~~

If you want to commit example outputs from synthetic data:
- do not ignore `results/` globally
- ignore `results/real/`
- track `results/example/`

---

## Repository layout

Recommended structure:

~~~text
.
├── analysis/
│   ├── 01_deseq2_benzene_ga.R
│   ├── 02_fgsea_pathways_continuous_benzene.R
│   ├── 03_hallmark_ifn_alpha_gene_view.R
│   └── 04_selected_genes_viewer.R
├── data/
│   ├── raw/          # NOT tracked (real data lives here locally)
│   └── example/      # optional: synthetic example inputs
├── results/          # outputs written here (not tracked by default)
├── .gitignore
└── README.md
~~~

---

## Inputs expected (real data not included)

### 1) RNA-seq counts (Excel)
- Path: `data/raw/BenezeneAll.xlsx`
- Sheet: `BenzeneAllw1179`
- Format:
  - First column named `gene` (gene symbols)
  - Remaining columns are sample IDs (example: `s1043`, `s1081`, ...)

### 2) Benzene adducts (Excel)
- Path: `data/raw/PPIH_plasma_Benzene_adduct_countsKG.xlsx`
- Sheet: `CLEAR_adducts`
- Required columns:
  - `sample`
  - `adduct count`

### 3) EMR metadata (Excel)
- Path: `data/raw/PPIH_EMR_data_deidentified.xlsx`
- Sheet: `Sheet1`
- Expected columns (as used in the scripts):
  - `REDCap ID`
  - `Infant Sex`
  - `Gestational Age`
  - optional: `Birth Weight (g)`, `Birth Length (cm)`, `Mode of Delivery`

The workflow constructs `sample_id` as `paste0("s", redcap_id)` to match RNA-seq sample IDs.

---

## Modeling overview

### Primary model: continuous benzene, GA-adjusted
DESeq2 design:

`~ sex + ga_weeks_num + adduct_z`

- `adduct_z` is the exposure variable (1 SD increase in adduct counts).
- `ga_weeks_num` is numeric gestational age in weeks (parsed from GA strings).
- `sex` is included as a covariate.

### Exploratory model: High vs Low split
A secondary model is included as an exploratory comparison using a median split:

`~ sex + ga_weeks_num + adduct_group`

This is included as a convenience check and should be treated as exploratory due to information loss from dichotomization.

---

## Requirements

- R (>= 4.1 recommended)
- CRAN packages: `readxl`, `dplyr`, `tibble`, `ggplot2`, `stringr`, `msigdbr`, `fgsea`, `pheatmap`
- Bioconductor packages: `DESeq2`, `SummarizedExperiment`

Install packages:

~~~r
install.packages(c("readxl","dplyr","tibble","ggplot2","stringr","msigdbr","fgsea","pheatmap"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","SummarizedExperiment"), ask = FALSE, update = FALSE)
~~~

---

## How to run

### 1) Put real data in `data/raw/` (locally only)
Place your input Excel files in:

- `data/raw/BenezeneAll.xlsx`
- `data/raw/PPIH_plasma_Benzene_adduct_countsKG.xlsx`
- `data/raw/PPIH_EMR_data_deidentified.xlsx`

If your filenames or sheet names differ, edit the parameters at the top of `analysis/01_deseq2_benzene_ga.R`.

### 2) Run scripts in order
From the repo root:

~~~bash
Rscript analysis/01_deseq2_benzene_ga.R
Rscript analysis/02_fgsea_pathways_continuous_benzene.R
Rscript analysis/03_hallmark_ifn_alpha_gene_view.R
Rscript analysis/04_selected_genes_viewer.R
~~~

---

## Outputs

All outputs are written to `results/`.

### `01_deseq2_benzene_ga.R`
- `results/DESeq2_results_continuous_benzene_GAadjusted.csv`
- `results/DESeq2_results_high_vs_low_benzene_GAadjusted.csv` (optional exploratory)
- `results/PCA_vst_sex_adductGroup.png`
- `results/vsd_continuous_benzene_GAadjusted.rds`
- `results/sessionInfo.txt`

### `02_fgsea_pathways_continuous_benzene.R`
- `results/fgsea_GO_BP_continuous_benzene.csv`
- `results/fgsea_GO_BP_continuous_benzene_top15_bar.png`
- `results/fgsea_GO_BP_immune_continuous_benzene.csv`
- `results/fgsea_GO_BP_immune_continuous_benzene_top15_bar.png`
- `results/fgsea_Hallmark_immune_continuous_benzene.csv`
- `results/fgsea_Hallmark_immune_continuous_benzene_top<N>_bar.png`
- `results/sessionInfo_fgsea.txt`

Note: the immune GO BP analysis uses a keyword filter on GO BP gene set names (IMMUNE, INFLAMMATORY, CYTOKINE). Treat this as a practical subset, not a formal immune ontology.

### `03_hallmark_ifn_alpha_gene_view.R`
- `results/Hallmark_IFN_alpha_DESeq2_GAadjusted.csv`
- `results/Hallmark_IFN_gamma_DESeq2_GAadjusted.csv` (optional)
- `results/barplot_IFNalpha_log2FC_top30.png` (N is user-set)
- `results/heatmap_IFNalpha_top<N>.png`
- `results/scatter_IFNalpha_signature_vs_adductz.png`
- `results/lm_IFNalpha_signature_summary.txt`
- `results/sessionInfo_ifnalpha.txt`

### `04_selected_genes_viewer.R`
- `results/selectedGenes_DESeq2_subset.csv`
- `results/barplot_selectedGenes_log2FC.png`
- `results/scatter_<GENE>_expr_vs_adductz.png` (one per gene found)
- `results/sessionInfo_selected_genes.txt`

---

## Notes on interpretation

- The DESeq2 coefficient for `adduct_z` represents the estimated log2 fold change in expression per 1 SD increase in benzene adduct counts, adjusted for sex and gestational age.
- The PCA plot is based on VST-transformed counts for QC and visualization, not for inference.
- In the selected gene viewer, scatterplots are descriptive. The DESeq2 inference comes from the model coefficient table, not the scatter trend.

---

## Troubleshooting

- "No overlapping samples": confirm sample IDs match across counts, adducts, and EMR. The workflow expects RNA-seq sample IDs like `s####` and constructs `sample_id` from `REDCap ID` as `s<REDCap ID>`.
- GA parsing failures: confirm gestational age formatting. Common formats like `40w3d`, `40w`, `40+3`, or numeric weeks should work. If not, extend the GA parser in `01_deseq2_benzene_ga.R`.
- Missing genes in VST matrix: some genes may be filtered out or not present after prefiltering. Use the DESeq2 results table as the source of truth for modeled effects.

---

## Citation

If you use or adapt this workflow, please cite the repository in your methods or acknowledgments and cite the primary tools used (DESeq2, fgsea, MSigDB).

---


