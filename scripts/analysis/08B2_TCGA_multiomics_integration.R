#!/usr/bin/env Rscript
# ==============================================================================
# Step 08B2: TCGA Multi-Omics Integration (Mutation, CNV, Methylation, RPPA)
# ==============================================================================
# Integrates genomic alterations with validated program scores from Step 08B1-QC2
# using curated TCGA multi-omics data (TCGAbiolinks v2.1.1)
#
# Inputs:
#   - Approved programs from 08B1-QC2 (50 programs, evidence levels B/C)
#   - TCGA LUAD/LUSC MultiAssayExperiment objects (mutation, CNV, RPPA)
#   - TCGA methyl450 beta-value matrices
#   - ssGSEA program scores from 08B1
#
# Outputs:
#   - Mutation burden & gene-level mutation associations per program
#   - CNV (GISTIC) associations per program
#   - Methylation (beta) associations per program
#   - RPPA protein associations per program
#   - Cross-omics concordance summary
#   - Figures and tables in 03_results/step08_TCGA/B2/
#   - 04_figures/step08_TCGA/B2/
#
# IMPORTANT RULES:
#   - Never download new data
#   - Never re-run ssGSEA or GSVA
#   - Never modify Step 07/08A/08B1/QC results
#   - Use /usr/local/bin/Rscript
# ==============================================================================

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Step 08B2: TCGA Multi-Omics Integration\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ==============================================================================
# Section I - Setup and Package Loading
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION I: Setup and Package Loading\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

.libPaths(c(path.expand("~/Library/R/arm64/4.6/library"), .libPaths()))
options(stringsAsFactors = FALSE, warn = 1)

required_pkgs <- c("survival", "data.table", "dplyr", "tidyr", "tibble",
                    "stringr", "ggplot2", "pheatmap", "RColorBrewer",
                    "MultiAssayExperiment", "SummarizedExperiment",
                    "S4Vectors", "RaggedExperiment", "metafor")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("RaggedExperiment")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      stop(paste("Package", pkg, "not available."))
    }
  }
}

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(RaggedExperiment)
})

cat("All packages loaded successfully.\n\n")

# ==============================================================================
# Section II - Create Directories
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION II: Create Directories\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

dirs <- c(
  "03_results/step08_TCGA/B2",
  "03_results/step08_TCGA/B2/mutation",
  "03_results/step08_TCGA/B2/cnv",
  "03_results/step08_TCGA/B2/methylation",
  "03_results/step08_TCGA/B2/rppa",
  "03_results/step08_TCGA/B2/cross_omics",
  "03_results/step08_TCGA/B2/final",
  "04_figures/step08_TCGA/B2",
  "04_figures/step08_TCGA/B2/mutation",
  "04_figures/step08_TCGA/B2/cnv",
  "04_figures/step08_TCGA/B2/methylation",
  "04_figures/step08_TCGA/B2/rppa",
  "04_figures/step08_TCGA/B2/cross_omics"
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
cat("Output directories created.\n\n")

# ==============================================================================
# Section III - Freeze Inputs
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION III: Freeze Inputs\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

input_files <- c(
  "03_results/GSE243013_step08B1_QC2_COMPLETE.txt",
  "03_results/GSE243013_step08B1_VALIDATED_FOR_B2.txt",
  "03_results/step08_TCGA/B1_QC2/final/GSE243013_programs_approved_for_step08B2.csv",
  "03_results/step08_TCGA/B1_QC2/cox/GSE243013_canonical_Cox_results.csv.gz",
  "03_results/step08_TCGA/B1_QC2/meta/GSE243013_canonical_fixed_effect_meta.csv",
  "02_data/tcga/clinical/TCGA_LUAD_program_scores_ssGSEA.rds",
  "02_data/tcga/clinical/TCGA_LUSC_program_scores_ssGSEA.rds",
  "02_data/tcga/clinical/GSE243013_TCGA_LUAD_LUSC_patient_manifest.csv",
  "02_data/tcga/curated/LUAD/TCGA_LUAD_core_primary_tumors_v2.1.1.rds",
  "02_data/tcga/curated/LUSC/TCGA_LUSC_core_primary_tumors_v2.1.1.rds",
  "02_data/tcga/curated/LUAD/TCGA_LUAD_methyl450_primary_tumors_v2.1.1.rds",
  "02_data/tcga/curated/LUSC/TCGA_LUSC_methyl450_primary_tumors_v2.1.1.rds"
)

freeze_manifest <- data.frame(
  file_id = seq_along(input_files),
  file_path = input_files,
  md5 = NA_character_,
  size_bytes = NA_real_,
  exists = FALSE,
  stringsAsFactors = FALSE
)

for (i in seq_along(input_files)) {
  fp <- input_files[i]
  freeze_manifest$exists[i] <- file.exists(fp)
  if (freeze_manifest$exists[i]) {
    freeze_manifest$size_bytes[i] <- file.size(fp)
    md5_out <- tryCatch(
      system(paste("md5 -q", shQuote(fp)), intern = TRUE),
      error = function(e) NA_character_
    )
    freeze_manifest$md5[i] <- md5_out
  }
  cat(sprintf("  [%d/%d] %s: exists=%s\n",
              i, length(input_files), basename(fp), freeze_manifest$exists[i]))
}

write.csv(freeze_manifest,
          "03_results/step08_TCGA/B2/GSE243013_B2_frozen_input_manifest.csv",
          row.names = FALSE)
cat("\nFrozen input manifest saved.\n\n")

if (!all(freeze_manifest$exists)) {
  stop("Not all required input files exist. Cannot proceed.")
}
cat("All required input files verified.\n\n")

# ==============================================================================
# Section IV - Load Approved Programs & Clinical Data
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION IV: Load Approved Programs & Clinical Data\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

approved <- read.csv(
  "03_results/step08_TCGA/B1_QC2/final/GSE243013_programs_approved_for_step08B2.csv",
  stringsAsFactors = FALSE
)
cat(sprintf("Approved programs: %d rows\n", nrow(approved)))

# Deduplicate by program_id (take unique program IDs)
approved_unique <- approved[!duplicated(approved$program_id), ]
cat(sprintf("Unique program IDs: %d\n", nrow(approved_unique)))
cat("Program IDs:\n")
for (pid in approved_unique$program_id) {
  cat("  -", pid, "\n")
}

patient_manifest <- read.csv(
  "02_data/tcga/clinical/GSE243013_TCGA_LUAD_LUSC_patient_manifest.csv",
  stringsAsFactors = FALSE
)
cat(sprintf("\nPatient manifest: %d patients\n", nrow(patient_manifest)))
cat(sprintf("Cohorts: %s\n", paste(sort(unique(patient_manifest$cohort)), collapse = ", ")))

cox_results <- data.table::fread(
  "03_results/step08_TCGA/B1_QC2/cox/GSE243013_canonical_Cox_results.csv.gz",
  stringsAsFactors = FALSE
)
cat(sprintf("Canonical Cox results: %d rows\n", nrow(cox_results)))

meta_results <- read.csv(
  "03_results/step08_TCGA/B1_QC2/meta/GSE243013_canonical_fixed_effect_meta.csv",
  stringsAsFactors = FALSE
)
cat(sprintf("Meta-analysis results: %d rows\n", nrow(meta_results)))

# ==============================================================================
# Section V - Load Program Scores
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION V: Load Program Scores\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

score_list <- list()
for (cohort in c("LUAD", "LUSC")) {
  path <- sprintf("02_data/tcga/clinical/TCGA_%s_program_scores_ssGSEA.rds", cohort)
  score_list[[cohort]] <- readRDS(path)
  cat(sprintf("[OK] %s ssGSEA scores: %d programs x %d patients\n",
              cohort, nrow(score_list[[cohort]]), ncol(score_list[[cohort]])))
}

# Filter to approved program IDs
approved_ids <- approved_unique$program_id
score_filtered <- list()
for (cohort in c("LUAD", "LUSC")) {
  scores <- score_list[[cohort]]
  common <- intersect(rownames(scores), approved_ids)
  score_filtered[[cohort]] <- scores[common, , drop = FALSE]
  cat(sprintf("[OK] %s filtered to %d programs\n", cohort, nrow(score_filtered[[cohort]])))
}

# ==============================================================================
# Section VI - Map Sample IDs to Patient IDs
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION VI: Map Sample IDs to Patient IDs\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Load sample-patient mapping
sample_map <- data.table::fread(
  "02_data/tcga/manifests/GSE243013_TCGA_sample_patient_map.csv.gz",
  stringsAsFactors = FALSE
)
cat(sprintf("Sample-patient map: %d rows\n", nrow(sample_map)))
cat("Columns:", paste(colnames(sample_map), collapse = ", "), "\n")

# Create patient ID extraction function (first 12 chars of TCGA barcode)
extract_patient_id <- function(barcode) {
  substr(barcode, 1, 12)
}

# ==============================================================================
# Section VII - Load TCGA Multi-Assay Data
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION VII: Load TCGA Multi-Assay Data\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

mae_data <- list()
for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s ---\n", cohort))
  
  mae_path <- sprintf(
    "02_data/tcga/curated/%s/TCGA_%s_core_primary_tumors_v2.1.1.rds",
    cohort, cohort
  )
  mae <- readRDS(mae_path)
  mae_data[[cohort]] <- mae
  
  exp_names <- names(experiments(mae))
  cat(sprintf("[INFO] Assays: %s\n", paste(exp_names, collapse = ", ")))
  cat(sprintf("[INFO] Patients: %d, Samples: %d\n", nrow(colData(mae)), ncol(mae)))
}

# ==============================================================================
# Section VIII - Extract Mutation Data
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION VIII: Extract Mutation Data\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Define known driver genes in NSCLC
driver_genes <- c(
  "TP53", "KRAS", "EGFR", "ALK", "BRAF", "PIK3CA", "PTEN", "RB1",
  "STK11", "KEAP1", "NF1", "MAP2K1", "ARID1A", "RBM10", "SMARCA4",
  "U2AF1", "SRSF2", "ATM", "CDKN2A", "ERBB2", "RET", "ROS1", "MET",
  "HER2", "FGFR1", "FGFR2", "FGFR3", "AKT1", "MAP3K1", "STAT3"
)

mutation_data <- list()

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s Mutations ---\n", cohort))
  
  mae <- mae_data[[cohort]]
  exp_names <- names(experiments(mae))
  mut_name <- grep("Mutation", exp_names, value = TRUE)
  
  if (length(mut_name) == 0) {
    cat(sprintf("[WARN] No Mutation assay found in %s\n", cohort))
    next
  }
  
  mut_se <- experiments(mae)[[mut_name[1]]]
  cat(sprintf("[INFO] Mutation class: %s\n", class(mut_se)))
  
  # Get sample barcodes and patient IDs
  mut_samples <- colnames(mut_se)
  mut_patients <- extract_patient_id(mut_samples)
  unique_patients <- sort(unique(mut_patients))
  
  # Get Hugo_Symbol from rowData (gene names)
  rd <- rowData(mut_se)
  hugo_symbols <- rd$Hugo_Symbol
  
  # Get assay matrix (values are Hugo_Symbol or NA)
  assay_mat <- assay(mut_se)
  
  cat(sprintf("[INFO] Mutation data: %d genomic positions x %d samples\n",
              nrow(assay_mat), ncol(assay_mat)))
  cat(sprintf("[INFO] Unique Hugo Symbols: %d\n", length(unique(hugo_symbols[!is.na(hugo_symbols)]))))
  
  # Create mutation matrix (gene x patient): 1 = mutated, 0 = wild-type
  # First, get unique gene names
  all_genes <- sort(unique(hugo_symbols[!is.na(hugo_symbols) & hugo_symbols != "Unknown"]))
  
  patient_mut <- matrix(0, nrow = length(all_genes), ncol = length(unique_patients),
                        dimnames = list(all_genes, unique_patients))
  
  # For each sample, check which genes have mutations
  for (j in seq_len(length(mut_samples))) {
    sample <- mut_samples[j]
    patient <- mut_patients[j]
    
    # Get gene names with mutations in this sample
    sample_genes <- assay_mat[, j]
    mutated_genes <- unique(sample_genes[!is.na(sample_genes) & sample_genes != "Unknown"])
    
    # Mark mutated genes
    common_genes <- intersect(mutated_genes, all_genes)
    if (length(common_genes) > 0) {
      patient_mut[common_genes, patient] <- 1
    }
  }
  
  # Remove genes with no mutations across all patients
  gene_mut_count <- rowSums(patient_mut)
  patient_mut <- patient_mut[gene_mut_count > 0, , drop = FALSE]
  
  cat(sprintf("[OK] Mutation matrix: %d genes mutated across %d patients\n",
              nrow(patient_mut), ncol(patient_mut)))
  
  # Calculate mutation burden per patient (total mutated genes)
  mut_burden <- colSums(patient_mut)
  
  # Calculate mutation burden per patient (total mutated genes)
  mut_burden <- colSums(patient_mut)
  
  mutation_data[[cohort]] <- list(
    matrix = patient_mut,
    burden = mut_burden,
    n_patients = ncol(patient_mut),
    n_genes = nrow(patient_mut)
  )
  
  cat(sprintf("[OK] %s: %d patients x %d genes mutated\n",
              cohort, ncol(patient_mut), nrow(patient_mut)))
  cat(sprintf("[OK] Mutation burden: mean=%.1f, median=%.1f, range=[%d, %d]\n",
              mean(mut_burden), median(mut_burden), min(mut_burden), max(mut_burden)))
}

# ==============================================================================
# Section IX - Extract CNV Data (GISTIC)
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION IX: Extract CNV Data (GISTIC)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cnv_data <- list()

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s CNV ---\n", cohort))
  
  mae <- mae_data[[cohort]]
  exp_names <- names(experiments(mae))
  
  # Use GISTIC_ThresholdedByGene for categorical CNV (-2, -1, 0, 1, 2)
  cnv_name <- grep("GISTIC_ThresholdedByGene", exp_names, value = TRUE)
  
  if (length(cnv_name) == 0) {
    cat(sprintf("[WARN] No GISTIC CNV assay found in %s\n", cohort))
    next
  }
  
  cnv_se <- experiments(mae)[[cnv_name[1]]]
  cnv_mat <- assay(cnv_se)
  
  # Map sample to patient
  cnv_samples <- colnames(cnv_mat)
  cnv_patients <- extract_patient_id(cnv_samples)
  
  # Aggregate to patient level (take mean of multiple samples per patient)
  patient_ids <- unique(cnv_patients)
  patient_cnv <- matrix(NA, nrow = nrow(cnv_mat), ncol = length(patient_ids),
                        dimnames = list(rownames(cnv_mat), patient_ids))
  
  for (pat in patient_ids) {
    idx <- which(cnv_patients == pat)
    if (length(idx) == 1) {
      patient_cnv[, pat] <- as.numeric(cnv_mat[, idx])
    } else {
      patient_cnv[, pat] <- rowMeans(cnv_mat[, idx, drop = FALSE])
    }
  }
  
  # Convert to numeric matrix
  storage.mode(patient_cnv) <- "double"
  
  # Remove genes with zero variance
  gene_var <- apply(patient_cnv, 1, var, na.rm = TRUE)
  patient_cnv <- patient_cnv[gene_var > 0, , drop = FALSE]
  
  # Ensure numeric matrix
  storage.mode(patient_cnv) <- "double"
  
  cnv_data[[cohort]] <- list(
    matrix = patient_cnv,
    n_patients = ncol(patient_cnv),
    n_genes = nrow(patient_cnv),
    cnv_levels = sort(unique(as.vector(patient_cnv)))
  )
  
  cat(sprintf("[OK] %s: %d patients x %d genes CNV\n",
              cohort, ncol(patient_cnv), nrow(patient_cnv)))
  cat(sprintf("[OK] CNV levels: %s\n",
              paste(cnv_data[[cohort]]$cnv_levels, collapse = ", ")))
}

# ==============================================================================
# Section X - Extract Methylation Data (450K)
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION X: Extract Methylation Data (450K)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

methylation_data <- list()

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s Methylation ---\n", cohort))
  
  methyl_path <- sprintf(
    "02_data/tcga/curated/%s/TCGA_%s_methyl450_primary_tumors_v2.1.1.rds",
    cohort, cohort
  )
  
  methyl_obj <- readRDS(methyl_path)
  cat(sprintf("[INFO] Methylation class: %s\n", class(methyl_obj)))
  
  # Check if it's a MultiAssayExperiment or a matrix
  if (is(methyl_obj, "MultiAssayExperiment")) {
    cat(sprintf("[INFO] Methylation is MultiAssayExperiment with %d experiments\n",
                length(experiments(methyl_obj))))
    exp_names <- names(experiments(methyl_obj))
    methyl_name <- grep("Methylation|methyl450", exp_names, value = TRUE)
    if (length(methyl_name) == 0) methyl_name <- exp_names[1]
    methyl_se <- experiments(methyl_obj)[[methyl_name[1]]]
    methyl_mat <- assay(methyl_se)
  } else if (is.matrix(methyl_obj) || is(methyl_obj, "SummarizedExperiment")) {
    methyl_mat <- if (is.matrix(methyl_obj)) methyl_obj else assay(methyl_obj)
  } else {
    cat(sprintf("[WARN] Unknown methylation class: %s\n", class(methyl_obj)))
    next
  }
  
  cat(sprintf("[INFO] Methylation dim: %d CpGs x %d samples\n",
              nrow(methyl_mat), ncol(methyl_mat)))
  
  # Map sample to patient
  methyl_samples <- colnames(methyl_mat)
  methyl_patients <- extract_patient_id(methyl_samples)
  
  # Aggregate to patient level (mean beta values)
  patient_ids <- unique(methyl_patients)
  
  # For large matrices, compute in chunks
  n_cpgs <- nrow(methyl_mat)
  chunk_size <- 10000
  patient_methyl <- matrix(NA, nrow = n_cpgs, ncol = length(patient_ids),
                           dimnames = list(rownames(methyl_mat), patient_ids))
  
  for (start in seq(1, n_cpgs, by = chunk_size)) {
    end <- min(start + chunk_size - 1, n_cpgs)
    chunk <- methyl_mat[start:end, , drop = FALSE]
    
    for (pat in patient_ids) {
      idx <- which(methyl_patients == pat)
      if (length(idx) == 1) {
        patient_methyl[start:end, pat] <- chunk[, idx]
      } else {
        patient_methyl[start:end, pat] <- rowMeans(chunk[, idx, drop = FALSE])
      }
    }
    cat(sprintf("  Processed CpGs %d-%d\n", start, end))
  }
  
  # Remove CpGs with zero variance
  cpg_var <- apply(patient_methyl, 1, var, na.rm = TRUE)
  patient_methyl <- patient_methyl[cpg_var > 0, , drop = FALSE]
  
  methylation_data[[cohort]] <- list(
    matrix = patient_methyl,
    n_patients = ncol(patient_methyl),
    n_cpgs = nrow(patient_methyl)
  )
  
  cat(sprintf("[OK] %s: %d patients x %d CpGs\n",
              cohort, ncol(patient_methyl), nrow(patient_methyl)))
}

# ==============================================================================
# Section XI - Extract RPPA Data
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION XI: Extract RPPA Data\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

rppa_data <- list()

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s RPPA ---\n", cohort))
  
  mae <- mae_data[[cohort]]
  exp_names <- names(experiments(mae))
  
  rppa_name <- grep("RPPA", exp_names, value = TRUE)
  
  if (length(rppa_name) == 0) {
    cat(sprintf("[WARN] No RPPA assay found in %s\n", cohort))
    next
  }
  
  rppa_se <- experiments(mae)[[rppa_name[1]]]
  rppa_mat <- assay(rppa_se)
  
  # Map sample to patient
  rppa_samples <- colnames(rppa_mat)
  rppa_patients <- extract_patient_id(rppa_samples)
  
  # Aggregate to patient level
  patient_ids <- unique(rppa_patients)
  patient_rppa <- matrix(NA, nrow = nrow(rppa_mat), ncol = length(patient_ids),
                         dimnames = list(rownames(rppa_mat), patient_ids))
  
  for (pat in patient_ids) {
    idx <- which(rppa_patients == pat)
    if (length(idx) == 1) {
      patient_rppa[, pat] <- rppa_mat[, idx]
    } else {
      patient_rppa[, pat] <- rowMeans(rppa_mat[, idx, drop = FALSE])
    }
  }
  
  # Remove proteins with zero variance
  prot_var <- apply(patient_rppa, 1, var, na.rm = TRUE)
  patient_rppa <- patient_rppa[prot_var > 0, , drop = FALSE]
  
  rppa_data[[cohort]] <- list(
    matrix = patient_rppa,
    n_patients = ncol(patient_rppa),
    n_proteins = nrow(patient_rppa),
    proteins = rownames(patient_rppa)
  )
  
  cat(sprintf("[OK] %s: %d patients x %d proteins\n",
              cohort, ncol(patient_rppa), nrow(patient_rppa)))
  cat(sprintf("[OK] Proteins: %s\n",
              paste(head(rownames(patient_rppa), 10), collapse = ", ")))
}

# ==============================================================================
# Section XII - Build Shared Patient Sets
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION XII: Build Shared Patient Sets\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

shared_patients <- list()

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s ---\n", cohort))
  
  # Get patient IDs with scores
  score_pats <- colnames(score_filtered[[cohort]])
  
  # Get patient IDs with each omics
  mut_pats <- if (cohort %in% names(mutation_data)) {
    colnames(mutation_data[[cohort]]$matrix)
  } else character(0)
  
  cnv_pats <- if (cohort %in% names(cnv_data)) {
    colnames(cnv_data[[cohort]]$matrix)
  } else character(0)
  
  methyl_pats <- if (cohort %in% names(methylation_data)) {
    colnames(methylation_data[[cohort]]$matrix)
  } else character(0)
  
  rppa_pats <- if (cohort %in% names(rppa_data)) {
    colnames(rppa_data[[cohort]]$matrix)
  } else character(0)
  
  # Find intersection of available omics (at minimum: scores + CNV + RPPA)
  available_omics <- list(scores = score_pats)
  if (length(mut_pats) > 0) available_omics$mutations <- mut_pats
  if (length(cnv_pats) > 0) available_omics$cnv <- cnv_pats
  if (length(methyl_pats) > 0) available_omics$methylation <- methyl_pats
  if (length(rppa_pats) > 0) available_omics$rppa <- rppa_pats
  
  all_pats <- Reduce(intersect, available_omics)
  
  if (length(all_pats) == 0) {
    # Try pairwise intersections
    cat("[WARN] No patients shared across ALL omics. Using pairwise intersections.\n")
    
    # For now, use intersection of scores, CNV, RPPA (most common)
    base_omics <- list(scores = score_pats)
    if (length(cnv_pats) > 0) base_omics$cnv <- cnv_pats
    if (length(rppa_pats) > 0) base_omics$rppa <- rppa_pats
    all_pats <- Reduce(intersect, base_omics)
    cat(sprintf("[INFO] Using %d patients with scores + CNV + RPPA\n", length(all_pats)))
  }
  
  shared_patients[[cohort]] <- list(
    all_omics = all_pats,
    scores = score_pats,
    mutations = mut_pats,
    cnv = cnv_pats,
    methylation = methyl_pats,
    rppa = rppa_pats
  )
  
  cat(sprintf("[OK] %s shared patients (all omics): %d\n", cohort, length(all_pats)))
  cat(sprintf("  Scores: %d, Mutations: %d, CNV: %d, Methylation: %d, RPPA: %d\n",
              length(score_pats), length(mut_pats), length(cnv_pats),
              length(methyl_pats), length(rppa_pats)))
}

# ==============================================================================
# Section XIII - Association Testing: Mutation
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION XIII: Association Testing - Mutation\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

mutation_results <- list()

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s ---\n", cohort))
  
  if (!(cohort %in% names(mutation_data))) {
    cat(sprintf("[SKIP] No mutation data for %s\n", cohort))
    next
  }
  
  shared_pats <- shared_patients[[cohort]]$all_omics
  if (length(shared_pats) < 30) {
    cat(sprintf("[SKIP] Insufficient shared patients (%d) for %s mutation analysis\n",
                length(shared_pats), cohort))
    next
  }
  
  scores <- score_filtered[[cohort]][, shared_pats, drop = FALSE]
  mut_mat <- mutation_data[[cohort]]$matrix[, shared_pats, drop = FALSE]
  mut_burden <- mutation_data[[cohort]]$burden[shared_pats]
  
  # Test driver genes with sufficient mutation frequency
  driver_in_data <- intersect(driver_genes, rownames(mut_mat))
  mut_freq <- rowSums(mut_mat[driver_in_data, , drop = FALSE]) / length(shared_pats)
  drivers_test <- names(mut_freq[mut_freq >= 0.05 & mut_freq <= 0.95])
  
  cat(sprintf("[INFO] Testing %d driver genes (MAF 5-95%%)\n", length(drivers_test)))
  
  results_list <- list()
  
  # Test mutation burden vs program scores
  for (prog in rownames(scores)) {
    score_vals <- scores[prog, shared_pats]
    burden_vals <- mut_burden[shared_pats]
    
    # Correlation test
    cor_test <- tryCatch(
      cor.test(score_vals, burden_vals, method = "spearman"),
      error = function(e) NULL
    )
    
    if (!is.null(cor_test)) {
      results_list[[length(results_list) + 1]] <- data.frame(
        program_id = prog,
        analysis = "mutation_burden",
        feature = "mutation_burden",
        statistic = "rho",
        estimate = cor_test$estimate,
        p_value = cor_test$p.value,
        n = length(shared_pats),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Test each driver gene vs each program
  for (gene in drivers_test) {
    mut_status <- as.numeric(mut_mat[gene, shared_pats])
    
    for (prog in rownames(scores)) {
      score_vals <- scores[prog, shared_pats]
      
      # t-test: mutated vs wild-type
      mut_group <- score_vals[mut_status == 1]
      wt_group <- score_vals[mut_status == 0]
      
      if (length(mut_group) >= 5 && length(wt_group) >= 5) {
        tt <- tryCatch(
          t.test(mut_group, wt_group),
          error = function(e) NULL
        )
        
        if (!is.null(tt)) {
          # Cohen's d
          pooled_sd <- sqrt(((length(mut_group) - 1) * var(mut_group) +
                             (length(wt_group) - 1) * var(wt_group)) /
                            (length(mut_group) + length(wt_group) - 2))
          cohen_d <- if (!is.na(pooled_sd) && pooled_sd > 0) {
            (mean(mut_group) - mean(wt_group)) / pooled_sd
          } else 0
          
          results_list[[length(results_list) + 1]] <- data.frame(
            program_id = prog,
            analysis = "driver_mutation",
            feature = gene,
            statistic = "cohen_d",
            estimate = cohen_d,
            p_value = tt$p.value,
            mean_mutated = mean(mut_group),
            mean_wildtype = mean(wt_group),
            n_mutated = length(mut_group),
            n_wildtype = length(wt_group),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (length(results_list) > 0) {
    mutation_results[[cohort]] <- rbindlist(results_list, fill = TRUE)
    
    # FDR correction per program
    mutation_results[[cohort]][, FDR := p.adjust(p_value, method = "BH"), by = program_id]
    
    # Save results
    fwrite(mutation_results[[cohort]],
           sprintf("03_results/step08_TCGA/B2/mutation/GSE243013_%s_mutation_associations.csv",
                   cohort))
    
    n_sig <- sum(mutation_results[[cohort]]$FDR < 0.05, na.rm = TRUE)
    cat(sprintf("[OK] %s: %d tests, %d significant (FDR<0.05)\n",
                cohort, nrow(mutation_results[[cohort]]), n_sig))
  }
}

# ==============================================================================
# Section XIV - Association Testing: CNV
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION XIV: Association Testing - CNV\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cnv_results <- list()

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s ---\n", cohort))
  
  if (!(cohort %in% names(cnv_data))) {
    cat(sprintf("[SKIP] No CNV data for %s\n", cohort))
    next
  }
  
  shared_pats <- shared_patients[[cohort]]$all_omics
  if (length(shared_pats) < 30) {
    cat(sprintf("[SKIP] Insufficient shared patients for %s CNV analysis\n", cohort))
    next
  }
  
  scores <- score_filtered[[cohort]][, shared_pats, drop = FALSE]
  cnv_mat <- cnv_data[[cohort]]$matrix[, shared_pats, drop = FALSE]
  
  # Test driver genes and frequently altered genes
  # GISTIC thresholded values are -2,-1,0,1,2; count non-zero as "altered"
  cnv_freq <- rowSums(cnv_mat != 0, na.rm = TRUE) / length(shared_pats)
  genes_test <- names(cnv_freq[cnv_freq >= 0.05 & cnv_freq <= 0.95])
  
  # Also include known driver genes
  genes_test <- union(genes_test, intersect(driver_genes, rownames(cnv_mat)))
  
  cat(sprintf("[INFO] Testing %d genes with CNV (freq 5-95%%)\n", length(genes_test)))
  
  results_list <- list()
  
  for (gene in genes_test) {
    cnv_vals <- as.numeric(cnv_mat[gene, shared_pats])
    
    for (prog in rownames(scores)) {
      score_vals <- scores[prog, shared_pats]
      
      # Correlation test (Spearman)
      cor_test <- tryCatch(
        cor.test(score_vals, cnv_vals, method = "spearman"),
        error = function(e) NULL
      )
      
      if (!is.null(cor_test)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          program_id = prog,
          analysis = "cnv_correlation",
          feature = gene,
          statistic = "rho",
          estimate = cor_test$estimate,
          p_value = cor_test$p.value,
          n = length(shared_pats),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(results_list) > 0) {
    cnv_results[[cohort]] <- rbindlist(results_list)
    cnv_results[[cohort]][, FDR := p.adjust(p_value, method = "BH"), by = program_id]
    
    fwrite(cnv_results[[cohort]],
           sprintf("03_results/step08_TCGA/B2/cnv/GSE243013_%s_cnv_associations.csv", cohort))
    
    n_sig <- sum(cnv_results[[cohort]]$FDR < 0.05, na.rm = TRUE)
    cat(sprintf("[OK] %s: %d tests, %d significant (FDR<0.05)\n",
                cohort, nrow(cnv_results[[cohort]]), n_sig))
  }
}

# ==============================================================================
# Section XV - Association Testing: Methylation
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION XV: Association Testing - Methylation\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

methylation_results <- list()

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s ---\n", cohort))
  
  if (!(cohort %in% names(methylation_data))) {
    cat(sprintf("[SKIP] No methylation data for %s\n", cohort))
    next
  }
  
  shared_pats <- shared_patients[[cohort]]$all_omics
  if (length(shared_pats) < 30) {
    cat(sprintf("[SKIP] Insufficient shared patients for %s methylation analysis\n", cohort))
    next
  }
  
  scores <- score_filtered[[cohort]][, shared_pats, drop = FALSE]
  methyl_mat <- methylation_data[[cohort]]$matrix[, shared_pats, drop = FALSE]
  
  cat(sprintf("[INFO] Testing %d CpGs x %d programs\n",
              nrow(methyl_mat), nrow(scores)))
  
  # For computational efficiency, test top 1000 most variable CpGs
  cpg_var <- apply(methyl_mat, 1, var, na.rm = TRUE)
  cpg_var[is.na(cpg_var)] <- 0
  n_variable <- sum(cpg_var > 0)
  top_cpgs <- if (n_variable > 0) {
    names(sort(cpg_var, decreasing = TRUE))[1:min(1000, n_variable)]
  } else {
    character(0)
  }
  
  results_list <- list()
  
  for (cpg in top_cpgs) {
    methyl_vals <- as.numeric(methyl_mat[cpg, shared_pats])
    
    for (prog in rownames(scores)) {
      score_vals <- scores[prog, shared_pats]
      
      cor_test <- tryCatch(
        cor.test(score_vals, methyl_vals, method = "spearman"),
        error = function(e) NULL
      )
      
      if (!is.null(cor_test)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          program_id = prog,
          analysis = "methylation_correlation",
          feature = cpg,
          statistic = "rho",
          estimate = cor_test$estimate,
          p_value = cor_test$p.value,
          n = length(shared_pats),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(results_list) > 0) {
    methylation_results[[cohort]] <- rbindlist(results_list)
    methylation_results[[cohort]][, FDR := p.adjust(p_value, method = "BH"), by = program_id]
    
    fwrite(methylation_results[[cohort]],
           sprintf("03_results/step08_TCGA/B2/methylation/GSE243013_%s_methylation_associations.csv",
                   cohort))
    
    n_sig <- sum(methylation_results[[cohort]]$FDR < 0.05, na.rm = TRUE)
    cat(sprintf("[OK] %s: %d tests, %d significant (FDR<0.05)\n",
                cohort, nrow(methylation_results[[cohort]]), n_sig))
  }
}

# ==============================================================================
# Section XVI - Association Testing: RPPA
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION XVI: Association Testing - RPPA\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

rppa_results <- list()

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s ---\n", cohort))
  
  if (!(cohort %in% names(rppa_data))) {
    cat(sprintf("[SKIP] No RPPA data for %s\n", cohort))
    next
  }
  
  shared_pats <- shared_patients[[cohort]]$all_omics
  if (length(shared_pats) < 30) {
    cat(sprintf("[SKIP] Insufficient shared patients for %s RPPA analysis\n", cohort))
    next
  }
  
  scores <- score_filtered[[cohort]][, shared_pats, drop = FALSE]
  rppa_mat <- rppa_data[[cohort]]$matrix[, shared_pats, drop = FALSE]
  
  cat(sprintf("[INFO] Testing %d proteins x %d programs\n",
              nrow(rppa_mat), nrow(scores)))
  
  results_list <- list()
  
  for (prot in rownames(rppa_mat)) {
    prot_vals <- as.numeric(rppa_mat[prot, shared_pats])
    
    for (prog in rownames(scores)) {
      score_vals <- scores[prog, shared_pats]
      
      cor_test <- tryCatch(
        cor.test(score_vals, prot_vals, method = "spearman"),
        error = function(e) NULL
      )
      
      if (!is.null(cor_test)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          program_id = prog,
          analysis = "rppa_correlation",
          feature = prot,
          statistic = "rho",
          estimate = cor_test$estimate,
          p_value = cor_test$p.value,
          n = length(shared_pats),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(results_list) > 0) {
    rppa_results[[cohort]] <- rbindlist(results_list)
    rppa_results[[cohort]][, FDR := p.adjust(p_value, method = "BH"), by = program_id]
    
    fwrite(rppa_results[[cohort]],
           sprintf("03_results/step08_TCGA/B2/rppa/GSE243013_%s_rppa_associations.csv",
                   cohort))
    
    n_sig <- sum(rppa_results[[cohort]]$FDR < 0.05, na.rm = TRUE)
    cat(sprintf("[OK] %s: %d tests, %d significant (FDR<0.05)\n",
                cohort, nrow(rppa_results[[cohort]]), n_sig))
  }
}

# ==============================================================================
# Section XVII - Cross-Omics Concordance Summary
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION XVII: Cross-Omics Concordance Summary\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cross_omics_summary <- list()

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s ---\n", cohort))
  
  summary_list <- list()
  
  # Mutation summary
  if (cohort %in% names(mutation_results) && nrow(mutation_results[[cohort]]) > 0) {
    mut_summ <- mutation_results[[cohort]][, .(
      n_tests = .N,
      n_sig_FDR05 = sum(FDR < 0.05, na.rm = TRUE),
      n_sig_FDR01 = sum(FDR < 0.01, na.rm = TRUE),
      min_FDR = min(FDR, na.rm = TRUE)
    ), by = program_id]
    mut_summ$omics_type <- "mutation"
    summary_list[["mutation"]] <- mut_summ
  }
  
  # CNV summary
  if (cohort %in% names(cnv_results) && nrow(cnv_results[[cohort]]) > 0) {
    cnv_summ <- cnv_results[[cohort]][, .(
      n_tests = .N,
      n_sig_FDR05 = sum(FDR < 0.05, na.rm = TRUE),
      n_sig_FDR01 = sum(FDR < 0.01, na.rm = TRUE),
      min_FDR = min(FDR, na.rm = TRUE)
    ), by = program_id]
    cnv_summ$omics_type <- "cnv"
    summary_list[["cnv"]] <- cnv_summ
  }
  
  # Methylation summary
  if (cohort %in% names(methylation_results) && nrow(methylation_results[[cohort]]) > 0) {
    methyl_summ <- methylation_results[[cohort]][, .(
      n_tests = .N,
      n_sig_FDR05 = sum(FDR < 0.05, na.rm = TRUE),
      n_sig_FDR01 = sum(FDR < 0.01, na.rm = TRUE),
      min_FDR = min(FDR, na.rm = TRUE)
    ), by = program_id]
    methyl_summ$omics_type <- "methylation"
    summary_list[["methylation"]] <- methyl_summ
  }
  
  # RPPA summary
  if (cohort %in% names(rppa_results) && nrow(rppa_results[[cohort]]) > 0) {
    rppa_summ <- rppa_results[[cohort]][, .(
      n_tests = .N,
      n_sig_FDR05 = sum(FDR < 0.05, na.rm = TRUE),
      n_sig_FDR01 = sum(FDR < 0.01, na.rm = TRUE),
      min_FDR = min(FDR, na.rm = TRUE)
    ), by = program_id]
    rppa_summ$omics_type <- "rppa"
    summary_list[["rppa"]] <- rppa_summ
  }
  
  if (length(summary_list) > 0) {
    cross_omics_summary[[cohort]] <- rbindlist(summary_list)
    
    fwrite(cross_omics_summary[[cohort]],
           sprintf("03_results/step08_TCGA/B2/cross_omics/GSE243013_%s_cross_omics_summary.csv",
                   cohort))
    
    cat(sprintf("[OK] %s cross-omics summary:\n", cohort))
    print(cross_omics_summary[[cohort]])
  }
}

# ==============================================================================
# Section XVIII - Top Associations Heatmap
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION XVIII: Top Associations Heatmap\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("\n--- %s ---\n", cohort))
  
  all_results <- rbindlist(list(
    if (cohort %in% names(mutation_results)) mutation_results[[cohort]],
    if (cohort %in% names(cnv_results)) cnv_results[[cohort]],
    if (cohort %in% names(methylation_results)) methylation_results[[cohort]],
    if (cohort %in% names(rppa_results)) rppa_results[[cohort]]
  ), fill = TRUE)
  
  if (nrow(all_results) == 0) {
    cat("[SKIP] No results to plot\n")
    next
  }
  
  # Get top associations per program (lowest FDR)
  top_hits <- all_results[!is.na(FDR), .SD[which.min(FDR)], by = program_id]
  top_hits <- top_hits[order(FDR)]
  
  # Save top hits
  fwrite(top_hits,
         sprintf("03_results/step08_TCGA/B2/cross_omics/GSE243013_%s_top_associations_per_program.csv",
                 cohort))
  
  cat(sprintf("[OK] %s: Top associations saved\n", cohort))
  cat(sprintf("  Programs with FDR<0.05: %d / %d\n",
              sum(top_hits$FDR < 0.05, na.rm = TRUE), nrow(top_hits)))
  cat(sprintf("  Programs with FDR<0.10: %d / %d\n",
              sum(top_hits$FDR < 0.10, na.rm = TRUE), nrow(top_hits)))
  
  # Create heatmap of top associations
  tryCatch({
    # Create matrix for heatmap
    sig_programs <- top_hits$program_id[top_hits$FDR < 0.10]
    
    if (length(sig_programs) > 0) {
      # Build heatmap matrix: programs x omics types
      omics_types <- c("mutation", "cnv", "methylation", "rppa")
      heatmap_data <- matrix(NA, nrow = length(sig_programs), ncol = length(omics_types),
                             dimnames = list(sig_programs, omics_types))
      
      for (ot in omics_types) {
        ot_results <- all_results[analysis == ot & program_id %in% sig_programs]
        if (nrow(ot_results) > 0) {
          # Get min FDR per program
          ot_fdr <- ot_results[, .(min_FDR = min(FDR, na.rm = TRUE)), by = program_id]
          for (i in seq_len(nrow(ot_fdr))) {
            if (ot_fdr$program_id[i] %in% rownames(heatmap_data)) {
              heatmap_data[ot_fdr$program_id[i], ot] <- -log10(ot_fdr$min_FDR[i])
            }
          }
        }
      }
      
      # Remove programs with no data
      heatmap_data <- heatmap_data[rowSums(!is.na(heatmap_data)) > 0, , drop = FALSE]
      
      if (nrow(heatmap_data) > 0) {
        pdf(sprintf("04_figures/step08_TCGA/B2/%s_top_associations_heatmap.pdf", cohort),
            width = 8, height = max(6, nrow(heatmap_data) * 0.4))
        pheatmap(heatmap_data,
                 main = sprintf("%s: Top Multi-Omics Associations", cohort),
                 xlab = "Omics Type",
                 ylab = "Program",
         cluster_cols = FALSE,
                 cluster_rows = TRUE,
                 color = colorRampPalette(c("white", "yellow", "red"))(100),
                 breaks = seq(0, 3, length.out = 101),
                 fontsize = 8)
        dev.off()
        cat(sprintf("[OK] %s heatmap saved\n", cohort))
      }
    }
  }, error = function(e) {
    cat(sprintf("[WARN] Could not create heatmap: %s\n", conditionMessage(e)))
    tryCatch(dev.off(), error = function(e2) NULL)
  })
}

# ==============================================================================
# Section XIX - Completion Marker
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION XIX: Completion Marker\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Save session info
session_info <- data.frame(
  key = c("R_version", "os", "timestamp", "script", "status"),
  value = c(R.version.string, Sys.info()["sysname"], as.character(Sys.time()),
            "08B2_TCGA_multiomics_integration.R", "COMPLETE"),
  stringsAsFactors = FALSE
)
write.csv(session_info,
          "03_results/step08_TCGA/B2/final/GSE243013_B2_session_info.csv",
          row.names = FALSE)

# Create completion marker
completion_text <- paste(
  "Step 08B2: TCGA Multi-Omics Integration\n",
  "Status: COMPLETE\n",
  "Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n",
  "Input files:", sum(freeze_manifest$exists), "/", nrow(freeze_manifest), "\n",
  "Programs analyzed:", nrow(approved_unique), "\n",
  "Cohorts:", paste(c("LUAD", "LUSC"), collapse = ", "), "\n",
  "Omics types: mutation, CNV, methylation, RPPA\n",
  "Output in: 03_results/step08_TCGA/B2/\n"
)
writeLines(completion_text, "03_results/GSE243013_step08B2_COMPLETE.txt")
cat("Completion marker created.\n\n")

# ==============================================================================
# Section XX - Final Report
# ==============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("SECTION XX: Final Report\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("\n--- Summary ---\n")
cat(sprintf("Programs analyzed: %d\n", nrow(approved_unique)))
cat(sprintf("Cohorts: LUAD, LUSC\n"))
cat("\n")

for (cohort in c("LUAD", "LUSC")) {
  cat(sprintf("--- %s Results ---\n", cohort))
  
  # Shared patients
  cat(sprintf("  Shared patients (all omics): %d\n",
              length(shared_patients[[cohort]]$all_omics)))
  
  # Mutation results
  if (cohort %in% names(mutation_results) && nrow(mutation_results[[cohort]]) > 0) {
    cat(sprintf("  Mutation tests: %d, Significant (FDR<0.05): %d\n",
                nrow(mutation_results[[cohort]]),
                sum(mutation_results[[cohort]]$FDR < 0.05, na.rm = TRUE)))
  } else {
    cat("  Mutation: No results\n")
  }
  
  # CNV results
  if (cohort %in% names(cnv_results) && nrow(cnv_results[[cohort]]) > 0) {
    cat(sprintf("  CNV tests: %d, Significant (FDR<0.05): %d\n",
                nrow(cnv_results[[cohort]]),
                sum(cnv_results[[cohort]]$FDR < 0.05, na.rm = TRUE)))
  } else {
    cat("  CNV: No results\n")
  }
  
  # Methylation results
  if (cohort %in% names(methylation_results) && nrow(methylation_results[[cohort]]) > 0) {
    cat(sprintf("  Methylation tests: %d, Significant (FDR<0.05): %d\n",
                nrow(methylation_results[[cohort]]),
                sum(methylation_results[[cohort]]$FDR < 0.05, na.rm = TRUE)))
  } else {
    cat("  Methylation: No results\n")
  }
  
  # RPPA results
  if (cohort %in% names(rppa_results) && nrow(rppa_results[[cohort]]) > 0) {
    cat(sprintf("  RPPA tests: %d, Significant (FDR<0.05): %d\n",
                nrow(rppa_results[[cohort]]),
                sum(rppa_results[[cohort]]$FDR < 0.05, na.rm = TRUE)))
  } else {
    cat("  RPPA: No results\n")
  }
  
  cat("\n")
}

# Save final results summary
final_summary <- data.frame(
  cohort = character(0),
  omics_type = character(0),
  n_tests = integer(0),
  n_sig_FDR05 = integer(0),
  stringsAsFactors = FALSE
)

for (cohort in c("LUAD", "LUSC")) {
  for (omics in c("mutation", "cnv", "methylation", "rppa")) {
    res_list <- switch(omics,
      "mutation" = mutation_results,
      "cnv" = cnv_results,
      "methylation" = methylation_results,
      "rppa" = rppa_results
    )
    
    if (cohort %in% names(res_list) && nrow(res_list[[cohort]]) > 0) {
      final_summary <- rbind(final_summary, data.frame(
        cohort = cohort,
        omics_type = omics,
        n_tests = nrow(res_list[[cohort]]),
        n_sig_FDR05 = sum(res_list[[cohort]]$FDR < 0.05, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }
}

write.csv(final_summary,
          "03_results/step08_TCGA/B2/final/GSE243013_B2_final_summary.csv",
          row.names = FALSE)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Step 08B2: COMPLETED\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
