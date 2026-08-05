# SUPERSEDED SCRIPT: FAIL_FDR_CALCULATION - 0% pass rate
# Original: 01_scripts/08B1_QC_audit_TCGA_survival_meta_results.R
# This file is kept for historical reference only.
# DO NOT USE for reproduction. 
# --- 
#############################################################################
# 08B1_QC_audit_TCGA_survival_meta_results.R
# Comprehensive QC audit for Step 08B1 NSCLC multi-omics survival analysis
#############################################################################

.libPaths(c(path.expand("~/Library/R/arm64/4.6/library"), .libPaths()))

cat("================================================================\n")
cat("08B1 QC Audit: TCGA Survival Meta-Results\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# =========================================================================
# Directory Setup
# =========================================================================
dir.create("03_results/step08_TCGA/B1_QC", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/step08_TCGA/B1_QC/cox", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/step08_TCGA/B1_QC/meta", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/step08_TCGA/B1_QC/scoring", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/step08_TCGA/B1_QC/redundancy", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/step08_TCGA/B1_QC/negative_control", recursive = TRUE, showWarnings = FALSE)
dir.create("04_figures/step08_TCGA/B1_QC", recursive = TRUE, showWarnings = FALSE)

cat("[DIR] Created QC output directories\n")

# =========================================================================
# Verify prerequisite completion file
# =========================================================================
if (!file.exists("03_results/GSE243013_step08B1_COMPLETE.txt")) {
  stop("FATAL: 03_results/GSE243013_step08B1_COMPLETE.txt not found. Run step08B1 first.")
}
cat("[CHECK] Prerequisite completion file exists\n")

# =========================================================================
# Package Loading
# =========================================================================
required_pkgs <- c("survival", "GSVA", "BiocParallel", "MultiAssayExperiment",
                    "SummarizedExperiment", "S4Vectors", "data.table", "dplyr",
                    "tidyr", "tibble", "stringr", "ggplot2", "pheatmap", "metafor")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("[INSTALL] Installing missing package:", pkg, "\n")
    if (pkg %in% c("GSVA", "BiocParallel", "MultiAssayExperiment",
                    "SummarizedExperiment", "S4Vectors")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

suppressPackageStartupMessages({
  library(survival)
  library(GSVA)
  library(BiocParallel)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(metafor)
})

cat("[LOAD] All packages loaded successfully\n\n")

# =========================================================================
# File Freeze (md5 checksums + dimensions)
# =========================================================================
cat("================================================================\n")
cat("SECTION I - File Freeze (md5 + dimensions)\n")
cat("================================================================\n")

freeze_files <- c(
  "02_data/tcga/clinical/TCGA_LUAD_program_scores_ssGSEA.rds",
  "02_data/tcga/clinical/TCGA_LUSC_program_scores_ssGSEA.rds",
  "02_data/tcga/clinical/TCGA_LUAD_program_scores_GSVA.rds",
  "02_data/tcga/clinical/TCGA_LUSC_program_scores_GSVA.rds",
  "02_data/tcga/clinical/GSE243013_TCGA_LUAD_LUSC_patient_manifest.csv",
  "03_results/step08_TCGA/clinical_models/LUAD/TCGA_LUAD_program_OS_Cox_results.csv.gz",
  "03_results/step08_TCGA/clinical_models/LUSC/TCGA_LUSC_program_OS_Cox_results.csv.gz",
  "03_results/step08_TCGA/clinical_models/GSE243013_TCGA_Cox_PH_assumption_results.csv.gz",
  "03_results/step08_TCGA/meta_analysis/GSE243013_TCGA_program_OS_fixed_effect_meta.csv",
  "03_results/step08_TCGA/combined/GSE243013_TCGA_clinically_validated_programs.csv",
  "03_results/step08_TCGA/scoring/GSE243013_TCGA_ssGSEA_vs_GSVA_concordance.csv",
  "03_results/step08_TCGA/programs/GSE243013_TCGA_program_manifest.csv",
  "03_results/step08_TCGA/programs/GSE243013_TCGA_program_gene_membership.csv.gz"
)

freeze_list <- list()
for (f in freeze_files) {
  if (file.exists(f)) {
    md5 <- tools::md5sum(f)
    info <- file.info(f)
    freeze_list[[f]] <- list(md5 = md5, size = info$size, modified = info$mtime)
    cat("  [FREEZE]", basename(f), "- md5:", md5, "\n")
  } else {
    cat("  [FREEZE]", basename(f), "- FILE NOT FOUND\n")
    freeze_list[[f]] <- list(md5 = NA, size = NA, modified = NA)
  }
}

freeze_df <- data.frame(
  file = names(freeze_list),
  md5 = sapply(freeze_list, `[[`, "md5"),
  size_bytes = sapply(freeze_list, `[[`, "size"),
  modified = sapply(freeze_list, `[[`, "modified"),
  stringsAsFactors = FALSE
)
write.csv(freeze_df, "03_results/step08_TCGA/B1_QC/file_freeze_manifest.csv", row.names = FALSE)
cat("[FREEZE] Saved file_freeze_manifest.csv\n\n")

# =========================================================================
# Load All Input Data
# =========================================================================
cat("================================================================\n")
cat("SECTION II - Load Input Data\n")
cat("================================================================\n")

# ssGSEA scores
ssgsea_luad <- readRDS("02_data/tcga/clinical/TCGA_LUAD_program_scores_ssGSEA.rds")
ssgsea_lusc <- readRDS("02_data/tcga/clinical/TCGA_LUSC_program_scores_ssGSEA.rds")
cat("[LOAD] ssGSEA LUAD:", nrow(ssgsea_luad), "programs x", ncol(ssgsea_luad), "patients\n")
cat("[LOAD] ssGSEA LUSC:", nrow(ssgsea_lusc), "programs x", ncol(ssgsea_lusc), "patients\n")

# GSVA scores
gsva_luad <- readRDS("02_data/tcga/clinical/TCGA_LUAD_program_scores_GSVA.rds")
gsva_lusc <- readRDS("02_data/tcga/clinical/TCGA_LUSC_program_scores_GSVA.rds")
cat("[LOAD] GSVA LUAD:", nrow(gsva_luad), "programs x", ncol(gsva_luad), "patients\n")
cat("[LOAD] GSVA LUSC:", nrow(gsva_lusc), "programs x", ncol(gsva_lusc), "patients\n")

# Patient manifest
patient_manifest <- read.csv("02_data/tcga/clinical/GSE243013_TCGA_LUAD_LUSC_patient_manifest.csv",
                              stringsAsFactors = FALSE)
cat("[LOAD] Patient manifest:", nrow(patient_manifest), "rows\n")

# Cox results (original)
cox_luad <- read.csv("03_results/step08_TCGA/clinical_models/LUAD/TCGA_LUAD_program_OS_Cox_results.csv.gz",
                      stringsAsFactors = FALSE)
cox_lusc <- read.csv("03_results/step08_TCGA/clinical_models/LUSC/TCGA_LUSC_program_OS_Cox_results.csv.gz",
                      stringsAsFactors = FALSE)
cat("[LOAD] Cox LUAD:", nrow(cox_luad), "rows\n")
cat("[LOAD] Cox LUSC:", nrow(cox_lusc), "rows\n")

# PH assumption results
ph_results <- read.csv("03_results/step08_TCGA/clinical_models/GSE243013_TCGA_Cox_PH_assumption_results.csv.gz",
                        stringsAsFactors = FALSE)
cat("[LOAD] PH results:", nrow(ph_results), "rows\n")

# Meta-analysis results (original)
meta_orig <- read.csv("03_results/step08_TCGA/meta_analysis/GSE243013_TCGA_program_OS_fixed_effect_meta.csv",
                       stringsAsFactors = FALSE)
cat("[LOAD] Meta-analysis:", nrow(meta_orig), "rows\n")

# Clinically validated programs
validated_df <- read.csv("03_results/step08_TCGA/combined/GSE243013_TCGA_clinically_validated_programs.csv",
                          stringsAsFactors = FALSE)
cat("[LOAD] Validated programs:", nrow(validated_df), "rows\n")

# Concordance
concordance_df <- read.csv("03_results/step08_TCGA/scoring/GSE243013_TCGA_ssGSEA_vs_GSVA_concordance.csv",
                            stringsAsFactors = FALSE)
cat("[LOAD] Concordance:", nrow(concordance_df), "rows\n")

# Program manifest
program_manifest <- read.csv("03_results/step08_TCGA/programs/GSE243013_TCGA_program_manifest.csv",
                              stringsAsFactors = FALSE)
cat("[LOAD] Program manifest:", nrow(program_manifest), "rows\n")

# Gene membership
gene_membership <- read.csv("03_results/step08_TCGA/programs/GSE243013_TCGA_program_gene_membership.csv.gz",
                             stringsAsFactors = FALSE)
cat("[LOAD] Gene membership:", nrow(gene_membership), "rows\n\n")


# =========================================================================
# SECTION III - Program ID Integrity
# =========================================================================
cat("================================================================\n")
cat("SECTION III - Program ID Integrity\n")
cat("================================================================\n")

program_ids_manifest <- unique(program_manifest$program_id)
n_unique <- length(program_ids_manifest)
cat("[III] Unique program_ids in manifest:", n_unique, "\n")

# Check no duplicates in manifest
dup_mask <- duplicated(program_manifest$program_id)
if (any(dup_mask)) {
  dup_ids <- program_manifest$program_id[dup_mask]
  cat("[III] WARNING: Duplicate program_ids in manifest:", length(dup_ids), "\n")
  cat("      Duplicates:", paste(head(dup_ids, 10), collapse = ", "), "\n")
} else {
  cat("[III] No duplicate program_ids in manifest: PASS\n")
}

# Check no NAs
na_mask <- is.na(program_manifest$program_id)
if (any(na_mask)) {
  cat("[III] WARNING: NA program_ids found:", sum(na_mask), "\n")
} else {
  cat("[III] No NA program_ids: PASS\n")
}

# Check validated_df has exactly 145 rows
validated_ids <- unique(validated_df$program_id)
n_validated <- length(validated_ids)
cat("[III] Unique program_ids in validated_df:", n_validated, "\n")
if (n_validated == 145) {
  cat("[III] Validated_df has exactly 145 programs: PASS\n")
} else {
  cat("[III] WARNING: Validated_df has", n_validated, "programs (expected 145)\n")
}

# Check 1:1 match between manifest and validated
if (nrow(validated_df) != 145) {
  cat("[III] FAIL: validated_df row count != 145 (got", nrow(validated_df), ")\n")
} else {
  manifest_set <- sort(program_ids_manifest)
  validated_set <- sort(validated_ids)
  if (identical(manifest_set, validated_set)) {
    cat("[III] Manifest and validated_df match 1:1: PASS\n")
  } else {
    only_manifest <- setdiff(manifest_set, validated_set)
    only_validated <- setdiff(validated_set, manifest_set)
    cat("[III] FAIL: Mismatch between manifest and validated_df\n")
    cat("      In manifest only:", length(only_manifest), "\n")
    if (length(only_manifest) > 0) cat("        ", paste(head(only_manifest, 10), collapse = ", "), "\n")
    cat("      In validated only:", length(only_validated), "\n")
    if (length(only_validated) > 0) cat("        ", paste(head(only_validated, 10), collapse = ", "), "\n")
  }
}

# Check for identical gene sets with different IDs
cat("[III] Checking for duplicate gene sets...\n")
gm_split <- split(gene_membership$gene, gene_membership$program_id)
gm_frozen <- lapply(gm_split, function(x) paste(sort(unique(x)), collapse = "|"))
gm_df <- data.frame(
  program_id = names(gm_frozen),
  gene_key = unlist(gm_frozen, use.names = FALSE),
  stringsAsFactors = FALSE
)
gene_dup <- gm_df[duplicated(gm_df$gene_key) | duplicated(gm_df$gene_key, fromLast = TRUE), ]
if (nrow(gene_dup) > 0) {
  cat("[III] WARNING: Identical gene sets found for different program_ids:\n")
  gene_dup_groups <- split(gene_dup$program_id, gene_dup$gene_key)
  for (i in seq_along(gene_dup_groups)) {
    cat("      Set", i, ":", paste(gene_dup_groups[[i]], collapse = " = "), "\n")
  }
} else {
  cat("[III] No identical gene sets with different IDs: PASS\n")
}

# Save Section III results
s3_result <- data.frame(
  check = c("n_unique_manifest", "n_validated", "duplicates_in_manifest",
            "na_in_manifest", "manifest_validated_match", "identical_gene_sets"),
  value = c(n_unique, n_validated, any(dup_mask), any(na_mask),
            identical(sort(manifest_set), sort(validated_set)), nrow(gene_dup) > 0),
  status = c(ifelse(n_unique == 145, "PASS", "FAIL"),
             ifelse(n_validated == 145, "PASS", "FAIL"),
             ifelse(any(dup_mask), "FAIL", "PASS"),
             ifelse(any(na_mask), "FAIL", "PASS"),
             ifelse(identical(sort(manifest_set), sort(validated_set)), "PASS", "FAIL"),
             ifelse(nrow(gene_dup) > 0, "WARN", "PASS")),
  stringsAsFactors = FALSE
)
write.csv(s3_result, "03_results/step08_TCGA/B1_QC/section_III_program_id_integrity.csv",
          row.names = FALSE)
cat("[III] Saved section_III_program_id_integrity.csv\n\n")


# =========================================================================
# SECTION IV - Concordance Denominator
# =========================================================================
cat("================================================================\n")
cat("SECTION IV - Concordance Denominator Audit\n")
cat("================================================================\n")

# Check matrix dimensions
cat("[IV] ssGSEA LUAD matrix:", nrow(ssgsea_luad), "x", ncol(ssgsea_luad), "\n")
cat("[IV] ssGSEA LUSC matrix:", nrow(ssgsea_lusc), "x", ncol(ssgsea_lusc), "\n")
cat("[IV] GSVA LUAD matrix:", nrow(gsva_luad), "x", ncol(gsva_luad), "\n")
cat("[IV] GSVA LUSC matrix:", nrow(gsva_lusc), "x", ncol(gsva_lusc), "\n")

# Check program/patient set overlap
luad_prog_set <- sort(rownames(ssgsea_luad))
lusc_prog_set <- sort(rownames(ssgsea_lusc))
luad_pat_set <- sort(colnames(ssgsea_luad))
lusc_pat_set <- sort(colnames(ssgsea_lusc))

cat("[IV] LUAD programs in ssGSEA:", length(luad_prog_set), "\n")
cat("[IV] LUSC programs in ssGSEA:", length(lusc_prog_set), "\n")
cat("[IV] LUAD patients in ssGSEA:", length(luad_pat_set), "\n")
cat("[IV] LUSC patients in ssGSEA:", length(lusc_pat_set), "\n")

luad_prog_in_manifest <- sum(luad_prog_set %in% program_ids_manifest)
lusc_prog_in_manifest <- sum(lusc_prog_set %in% program_ids_manifest)
cat("[IV] LUAD ssGSEA programs in manifest:", luad_prog_in_manifest, "/", length(luad_prog_set), "\n")
cat("[IV] LUSC ssGSEA programs in manifest:", lusc_prog_in_manifest, "/", length(lusc_prog_set), "\n")

# Verify GSVA matrices have same dimensions as ssGSEA
gsva_dims_match <- all(
  dim(ssgsea_luad) == dim(gsva_luad),
  dim(ssgsea_lusc) == dim(gsva_lusc)
)
cat("[IV] GSVA/ssGSEA dimension match:", ifelse(gsva_dims_match, "PASS", "FAIL"), "\n")

# Recalculate Spearman correlation independently
cat("[IV] Recalculating ssGSEA vs GSVA Spearman correlations...\n")
recalc_concordance <- data.frame()

for (cohort in c("LUAD", "LUSC")) {
  if (cohort == "LUAD") {
    ssgsea_mat <- ssgsea_luad
    gsva_mat <- gsva_luad
  } else {
    ssgsea_mat <- ssgsea_lusc
    gsva_mat <- gsva_lusc
  }
  
  common_progs <- intersect(rownames(ssgsea_mat), rownames(gsva_mat))
  cat("[IV]   ", cohort, "- Common programs:", length(common_progs), "\n")
  
  for (pg in common_progs) {
    x <- as.numeric(ssgsea_mat[pg, ])
    y <- as.numeric(gsva_mat[pg, ])
    complete <- !is.na(x) & !is.na(y)
    if (sum(complete) < 3) {
      rho <- NA_real_
    } else {
      rho <- cor(x[complete], y[complete], method = "spearman")
    }
    recalc_concordance <- rbind(recalc_concordance, data.frame(
      cohort = cohort,
      program_id = pg,
      rho_recalc = rho,
      stringsAsFactors = FALSE
    ))
  }
}

# Merge with original concordance
if ("rho" %in% names(concordance_df)) {
  rho_col <- "rho"
} else if ("spearman_rho" %in% names(concordance_df)) {
  rho_col <- "spearman_rho"
} else {
  rho_col <- names(concordance_df)[grepl("rho|cor|spearman", names(concordance_df), ignore.case = TRUE)]
  rho_col <- rho_col[1]
}

merged_conc <- merge(
  concordance_df, recalc_concordance,
  by = c("cohort", "program_id"), all = TRUE, suffixes = c("_orig", "_recalc")
)

if (!is.na(rho_col) && nrow(merged_conc) > 0) {
  conc_diff <- abs(merged_conc[[paste0(rho_col, "_orig")]] - merged_conc$rho_recalc)
  max_diff <- max(conc_diff, na.rm = TRUE)
  cat("[IV] Max Spearman rho difference (orig vs recalc):", max_diff, "\n")
  if (max_diff < 1e-10) {
    cat("[IV] Concordance recalculation matches: PASS\n")
  } else {
    cat("[IV] WARNING: Concordance recalculation differs\n")
  }
}

# Explained 224/225 denominator
luad_count <- sum(luad_prog_set %in% program_ids_manifest)
lusc_count <- sum(lusc_prog_set %in% program_ids_manifest)
theoretical_max <- luad_count + lusc_count
cat("[IV] Theoretical max (denominator):", theoretical_max, "\n")
cat("[IV]  LUAD programs evaluated:", luad_count, "\n")
cat("[IV]  LUSC programs evaluated:", lusc_count, "\n")

if (theoretical_max != 225) {
  missing_from_luad <- setdiff(program_ids_manifest, luad_prog_set)
  missing_from_lusc <- setdiff(program_ids_manifest, lusc_prog_set)
  cat("[IV] Programs missing from LUAD ssGSEA:", length(missing_from_luad), "\n")
  cat("[IV] Programs missing from LUSC ssGSEA:", length(missing_from_lusc), "\n")
  if (length(missing_from_luad) > 0)
    cat("      LUAD missing:", paste(head(missing_from_luad, 10), collapse = ", "), "\n")
  if (length(missing_from_lusc) > 0)
    cat("      LUSC missing:", paste(head(missing_from_lusc, 10), collapse = ", "), "\n")
}

n_rho_ge_06 <- sum(recalc_concordance$rho_recalc >= 0.60, na.rm = TRUE)
cat("[IV] Programs with rho >= 0.60:", n_rho_ge_06, "/", nrow(recalc_concordance), "\n")

s4_result <- data.frame(
  check = c("luad_ssGSEA_dims", "lusc_ssGSEA_dims", "gsva_dims_match",
            "theoretical_denominator", "actual_denominator_225", "n_rho_ge_060"),
  value = c(paste(nrow(ssgsea_luad), "x", ncol(ssgsea_luad)),
            paste(nrow(ssgsea_lusc), "x", ncol(ssgsea_lusc)),
            gsva_dims_match, theoretical_max,
            nrow(recalc_concordance), n_rho_ge_06),
  status = c("INFO", "INFO",
             ifelse(gsva_dims_match, "PASS", "FAIL"),
             ifelse(theoretical_max == 225, "PASS", "INFO"),
             "INFO", "INFO"),
  stringsAsFactors = FALSE
)
write.csv(s4_result, "03_results/step08_TCGA/B1_QC/section_IV_concordance_denominator.csv",
          row.names = FALSE)
cat("[IV] Saved section_IV_concordance_denominator.csv\n\n")


# =========================================================================
# SECTION V - Cox Field Audit
# =========================================================================
cat("================================================================\n")
cat("SECTION V - Cox Field Audit\n")
cat("================================================================\n")

required_cox_fields <- c("program_id", "model_level", "n_complete", "n_events",
                          "logHR", "standard_error", "HR_per_1SD", "lower_95CI", "upper_95CI",
                          "Wald_PValue", "FDR_within_cohort")

for (cohort in c("LUAD", "LUSC")) {
  if (cohort == "LUAD") {
    cx <- cox_luad
  } else {
    cx <- cox_lusc
  }
  cat("[V] Auditing Cox fields for", cohort, "(", nrow(cx), "rows )\n")
  
  # Check required fields
  missing_fields <- setdiff(required_cox_fields, names(cx))
  if (length(missing_fields) > 0) {
    cat("[V]   FAIL: Missing fields:", paste(missing_fields, collapse = ", "), "\n")
  } else {
    cat("[V]   All required fields present: PASS\n")
  }
  
  # Check program_id uniqueness
  if (any(duplicated(cx$program_id))) {
    cat("[V]   FAIL: Duplicate program_ids found\n")
  } else {
    cat("[V]   No duplicate program_ids: PASS\n")
  }
  
  # Check FULL_MODEL count
  full_models <- cx[cx$model_level == "FULL_MODEL", ]
  if (nrow(full_models) == 145) {
    cat("[V]   FULL_MODEL count = 145: PASS\n")
  } else {
    cat("[V]   FAIL: FULL_MODEL count =", nrow(full_models), "(expected 145)\n")
  }
  
  # Check numeric fields are finite
  numeric_fields <- c("logHR", "standard_error", "HR_per_1SD", "lower_95CI", "upper_95CI", "Wald_PValue", "FDR_within_cohort")
  for (nf in numeric_fields) {
    if (nf %in% names(cx)) {
      vals <- cx[[nf]]
      n_na <- sum(is.na(vals))
      n_inf <- sum(is.infinite(vals))
      if (n_na > 0 || n_inf > 0) {
        cat("[V]   WARNING:", nf, "- NAs:", n_na, "Infs:", n_inf, "\n")
      }
    }
  }
  
  # Check SE > 0
  if ("standard_error" %in% names(cx)) {
    n_se_neg <- sum(cx$standard_error <= 0, na.rm = TRUE)
    if (n_se_neg > 0) {
      cat("[V]   FAIL:", n_se_neg, "rows with SE <= 0\n")
    } else {
      cat("[V]   All SE > 0: PASS\n")
    }
  }
  
  # Check P in [0,1] and FDR in [0,1]
  if ("Wald_PValue" %in% names(cx)) {
    p_out <- sum(cx$Wald_PValue < 0 | cx$Wald_PValue > 1, na.rm = TRUE)
    if (p_out > 0) cat("[V]   FAIL:", p_out, "P-values outside [0,1]\n")
    else cat("[V]   All P-values in [0,1]: PASS\n")
  }
  if ("FDR" %in% names(cx)) {
    fdr_out <- sum(cx$FDR_within_cohort < 0 | cx$FDR_within_cohort > 1, na.rm = TRUE)
    if (fdr_out > 0) cat("[V]   FAIL:", fdr_out, "FDR values outside [0,1]\n")
    else cat("[V]   All FDR in [0,1]: PASS\n")
  }
  
  # Verify HR = exp(logHR)
  if (all(c("HR_per_1SD", "logHR") %in% names(cx))) {
    hr_diff <- abs(cx$HR_per_1SD - exp(cx$logHR))
    hr_diff <- hr_diff[!is.na(hr_diff)]
    max_hr_diff <- max(hr_diff)
    if (max_hr_diff < 1e-10) {
      cat("[V]   HR = exp(logHR): PASS\n")
    } else {
      cat("[V]   FAIL: max |HR - exp(logHR)| =", max_hr_diff, "\n")
    }
  }
  
  # Verify lower_CI = exp(logHR - 1.96*SE)
  if (all(c("lower_95CI", "logHR", "standard_error") %in% names(cx))) {
    lower_expected <- exp(cx$logHR - 1.96 * cx$standard_error)
    lower_diff <- abs(cx$lower_95CI - lower_expected)
    lower_diff <- lower_diff[!is.na(lower_diff)]
    max_lower_diff <- max(lower_diff)
    if (max_lower_diff < 1e-8) {
      cat("[V]   lower_CI = exp(logHR - 1.96*SE): PASS\n")
    } else {
      cat("[V]   FAIL: max |lower_CI - expected| =", max_lower_diff, "\n")
    }
  }
  
  # Verify upper_CI = exp(logHR + 1.96*SE)
  if (all(c("upper_95CI", "logHR", "standard_error") %in% names(cx))) {
    upper_expected <- exp(cx$logHR + 1.96 * cx$standard_error)
    upper_diff <- abs(cx$upper_95CI - upper_expected)
    upper_diff <- upper_diff[!is.na(upper_diff)]
    max_upper_diff <- max(upper_diff)
    if (max_upper_diff < 1e-8) {
      cat("[V]   upper_CI = exp(logHR + 1.96*SE): PASS\n")
    } else {
      cat("[V]   FAIL: max |upper_CI - expected| =", max_upper_diff, "\n")
    }
  }
  
  # Check lower <= HR <= upper
  if (all(c("lower_95CI", "HR_per_1SD", "upper_95CI") %in% names(cx))) {
    ci_order <- cx$lower_95CI <= cx$HR_per_1SD & cx$HR_per_1SD <= cx$upper_95CI
    ci_order <- ci_order[!is.na(ci_order)]
    if (all(ci_order)) {
      cat("[V]   lower_CI <= HR <= upper_CI: PASS\n")
    } else {
      cat("[V]   FAIL:", sum(!ci_order), "rows violate CI ordering\n")
    }
  }
}

cat("[V] Cox field audit complete\n\n")


# =========================================================================
# SECTION VI - Independent Cox Refitting
# =========================================================================
cat("================================================================\n")
cat("SECTION VI - Independent Cox Refitting\n")
cat("================================================================\n")

# Prepare clinical covariates from manifest
clin <- patient_manifest
# Map field names
if (!"OS_days" %in% names(clin) && "OS.time" %in% names(clin)) {
  clin$OS_days <- clin$OS.time
}
if (!"OS_event" %in% names(clin) && "OS.status" %in% names(clin)) {
  clin$OS_event <- as.numeric(grepl("[1eE]", clin$OS.status))
}
if (!"age" %in% names(clin)) {
  if ("age_at_diagnosis" %in% names(clin)) {
    clin$age <- clin$age_at_diagnosis / 365.25
  } else if ("age_at_initial_pathologic_diagnosis" %in% names(clin)) {
    clin$age <- clin$age_at_initial_pathologic_diagnosis
  }
}
if (!"sex" %in% names(clin) && "gender" %in% names(clin)) {
  clin$sex <- clin$gender
}
if (!"stage" %in% names(clin)) {
  if ("pathologic_stage_clean" %in% names(clin)) {
    clin$stage <- clin$pathologic_stage_clean
  } else if ("pathologic_stage_original" %in% names(clin)) {
    clin$stage <- clin$pathologic_stage_original
  }
}

# Standardize continuous covariates within cohort
clin$age_z <- as.numeric(scale(clin$age))
clin$sex_f <- as.factor(clin$sex)
clin$stage_f <- as.factor(clin$stage)
clin$stage_f <- droplevels(clin$stage_f)

# Compute standardized score_z within cohort
std_scores_luad <- t(scale(t(ssgsea_luad)))
std_scores_lusc <- t(scale(t(ssgsea_lusc)))

cat("[VI] Refitting Cox models for each program in each cohort...\n")
cat("[VI] Model: Surv(OS_days/365.25, OS_event) ~ score_z + age_z + sex_f + stage_f\n")
cat("[VI] score_z standardized within cohort\n")

cox_recalc_results <- list()

for (cohort in c("LUAD", "LUSC")) {
  if (cohort == "LUAD") {
    std_mat <- std_scores_luad
    orig_cox <- cox_luad
  } else {
    std_mat <- std_scores_lusc
    orig_cox <- cox_lusc
  }
  
  # Align patients
  common_pats <- intersect(colnames(std_mat), clin$patient_id)
  if (length(common_pats) == 0) {
    # Try barcode format
    common_pats <- intersect(colnames(std_mat), clin$barcode)
  }
  
  clin_sub <- clin[match(common_pats, clin$patient_id), ]
  if (all(is.na(clin_sub$patient_id))) {
    clin_sub <- clin[match(common_pats, clin$barcode), ]
  }
  
  cat("[VI]   ", cohort, "- Common patients:", length(common_pats), "\n")
  
  cohort_results <- data.frame()
  
  for (i in seq_len(nrow(std_mat))) {
    pg <- rownames(std_mat)[i]
    score_vec <- as.numeric(std_mat[pg, common_pats])
    
    # Build data frame
    df_i <- data.frame(
      OS_time = clin_sub$OS_days / 365.25,
      OS_event = clin_sub$OS_event,
      score_z = score_vec,
      age_z = as.numeric(clin_sub$age_z),
      sex_f = clin_sub$sex_f,
      stage_f = clin_sub$stage_f,
      stringsAsFactors = FALSE
    )
    
    # Remove rows with any NA
    complete_mask <- complete.cases(df_i)
    df_complete <- df_i[complete_mask, ]
    n_complete <- nrow(df_complete)
    n_events <- sum(df_complete$OS_event == 1, na.rm = TRUE)
    
    if (n_complete < 20 || n_events < 5) {
      cohort_results <- rbind(cohort_results, data.frame(
        program_id = pg,
        model_level = "INSUFFICIENT",
        n_complete = n_complete,
        n_events = n_events,
        logHR = NA, SE = NA, HR = NA, lower_CI = NA, upper_CI = NA,
        Wald_PValue = NA, FDR = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    stage_complete <- droplevels(df_complete$stage_f)
    
    fit <- tryCatch({
      coxph(Surv(OS_time, OS_event) ~ score_z + age_z + sex_f + stage_complete,
            data = df_complete)
    }, error = function(e) NULL)
    
    if (is.null(fit)) {
      cohort_results <- rbind(cohort_results, data.frame(
        program_id = pg,
        model_level = "FIT_ERROR",
        n_complete = n_complete,
        n_events = n_events,
        logHR = NA, SE = NA, HR = NA, lower_CI = NA, upper_CI = NA,
        Wald_PValue = NA, FDR = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    coef_table <- summary(fit)$coefficients
    # Extract score_z coefficient EXACTLY by name, not by row number
    if ("score_z" %in% rownames(coef_table)) {
      score_row <- coef_table["score_z", ]
    } else {
      # Try partial match
      score_idx <- grep("score_z", rownames(coef_table))
      if (length(score_idx) == 1) {
        score_row <- coef_table[score_idx, ]
      } else {
        cat("[VI]   WARNING: score_z not found in coefficients for", pg, "\n")
        next
      }
    }
    
    logHR_i <- score_row["coef"]
    se_i <- score_row["se(coef)"]
    z_val <- score_row["z"]
    p_val <- score_row["Pr(>|z|)"]
    hr_i <- exp(logHR_i)
    lower_i <- exp(logHR_i - 1.96 * se_i)
    upper_i <- exp(logHR_i + 1.96 * se_i)
    
    # Determine ties method
    ties_method <- fit$method
    if (is.null(ties_method)) ties_method <- "efron"
    
    cohort_results <- rbind(cohort_results, data.frame(
      program_id = pg,
      model_level = "FULL_MODEL",
      n_complete = n_complete,
      n_events = n_events,
      logHR = logHR_i,
      SE = se_i,
      HR = hr_i,
      lower_CI = lower_i,
      upper_CI = upper_i,
      Wald_PValue = p_val,
      FDR = NA,
      ties_method = ties_method,
      stringsAsFactors = FALSE
    ))
  }
  
  # Compute FDR
  valid_mask <- !is.na(cohort_results$Wald_PValue) & is.finite(cohort_results$Wald_PValue)
  cohort_results$FDR <- NA
  cohort_results$FDR[valid_mask] <- p.adjust(cohort_results$Wald_PValue[valid_mask], method = "BH")
  
  cox_recalc_results[[cohort]] <- cohort_results
  
  n_full <- sum(cohort_results$model_level == "FULL_MODEL", na.rm = TRUE)
  cat("[VI]   ", cohort, "- FULL_MODEL fitted:", n_full, "/ 145\n")
}

# Save recalculated Cox results
for (cohort in c("LUAD", "LUSC")) {
  outfile <- paste0("03_results/step08_TCGA/B1_QC/cox/",
                     cohort, "_Cox_recalculated.csv.gz")
  write.csv(cox_recalc_results[[cohort]], file = gzfile(outfile), row.names = FALSE)
  cat("[VI] Saved", outfile, "\n")
}

cat("[VI] Independent Cox refitting complete\n\n")


# =========================================================================
# SECTION VII - Compare Cox (Original vs Recalculated)
# =========================================================================
cat("================================================================\n")
cat("SECTION VII - Compare Original vs Recalculated Cox\n")
cat("================================================================\n")

compare_cox_all <- data.frame()

for (cohort in c("LUAD", "LUSC")) {
  if (cohort == "LUAD") {
    orig <- cox_luad
  } else {
    orig <- cox_lusc
  }
  recalc <- cox_recalc_results[[cohort]]
  
  # Rename recalc columns to match original names for clean merge
  names(recalc)[names(recalc) == "SE"] <- "standard_error"
  names(recalc)[names(recalc) == "HR"] <- "HR_per_1SD"
  names(recalc)[names(recalc) == "lower_CI"] <- "lower_95CI"
  names(recalc)[names(recalc) == "upper_CI"] <- "upper_95CI"
  names(recalc)[names(recalc) == "FDR"] <- "FDR_within_cohort"
  
  # Merge by cohort + program_id
  merged <- merge(
    orig, recalc,
    by = "program_id", all = TRUE, suffixes = c("_orig", "_recalc")
  )
  merged$cohort <- cohort
  
  # Calculate differences
  merged$delta_logHR <- merged$logHR_orig - merged$logHR_recalc
  merged$relative_SE_diff <- abs(merged$standard_error_orig - merged$standard_error_recalc) /
    (abs(merged$standard_error_orig) + 1e-15)
  merged$delta_PValue <- merged$Wald_PValue_orig - merged$Wald_PValue_recalc
  
  # Direction match
  merged$direction_match <- sign(merged$logHR_orig) == sign(merged$logHR_recalc)
  
  # Model level match
  merged$model_level_match <- merged$model_level_orig == merged$model_level_recalc
  
  # n_complete match
  merged$n_complete_same <- merged$n_complete_orig == merged$n_complete_recalc
  
  # n_events match
  merged$n_events_same <- merged$n_events_orig == merged$n_events_recalc
  
  # PASS criteria
  merged$PASS <- with(merged,
    abs(delta_logHR) < 1e-8 &
    relative_SE_diff < 1e-6 &
    n_complete_same &
    n_events_same
  )
  
  compare_cox_all <- rbind(compare_cox_all, merged)
  
  n_pass <- sum(merged$PASS, na.rm = TRUE)
  n_total <- sum(!is.na(merged$PASS))
  cat("[VII]   ", cohort, "- PASS:", n_pass, "/", n_total, "\n")
  
  if (n_pass < n_total) {
    failing <- merged[!merged$PASS, ]
    cat("[VII]   FAILURES:", nrow(failing), "programs\n")
    # Show top failures
    failing_sorted <- failing[order(abs(failing$delta_logHR), decreasing = TRUE), ]
    for (j in seq_len(min(5, nrow(failing_sorted)))) {
      cat("         ", failing_sorted$program_id[j],
          " delta_logHR=", failing_sorted$delta_logHR[j],
          " rel_SE_diff=", failing_sorted$relative_SE_diff[j], "\n")
    }
  }
}

write.csv(compare_cox_all,
          "03_results/step08_TCGA/B1_QC/cox/original_vs_recalculated_Cox_comparison.csv",
          row.names = FALSE)
cat("[VII] Saved Cox comparison CSV\n")

# Summary
overall_pass_rate <- mean(compare_cox_all$PASS, na.rm = TRUE)
cat("[VII] Overall PASS rate:", round(overall_pass_rate * 100, 2), "%\n")
cat("[VII] Cox comparison complete\n\n")


# =========================================================================
# SECTION VIII - FDR Recalculation
# =========================================================================
cat("================================================================\n")
cat("SECTION VIII - FDR Recalculation\n")
cat("================================================================\n")

fdr_comparison_all <- data.frame()

for (cohort in c("LUAD", "LUSC")) {
  recalc <- cox_recalc_results[[cohort]]
  
  # Filter: FULL_MODEL + complete + finite P
  filter_mask <- recalc$model_level == "FULL_MODEL" &
    !is.na(recalc$Wald_PValue) &
    is.finite(recalc$Wald_PValue)
  
  fdr_recalc <- rep(NA_real_, nrow(recalc))
  fdr_recalc[filter_mask] <- p.adjust(recalc$Wald_PValue[filter_mask], method = "BH")
  recalc$FDR_recalculated <- fdr_recalc
  
  # Compare with existing FDR from the same refitting
  compare_df <- data.frame(
    program_id = recalc$program_id,
    cohort = cohort,
    FDR_recalculated = recalc$FDR_recalculated,
    PValue = recalc$Wald_PValue,
    model_level = recalc$model_level,
    stringsAsFactors = FALSE
  )
  
  # Check if FDR was re-p.adjusted incorrectly
  if ("FDR_orig" %in% names(compare_df)) {
    fdr_diff <- abs(compare_df$FDR_recalculated - compare_df$FDR_orig)
    n_diff <- sum(fdr_diff > 1e-10, na.rm = TRUE)
    cat("[VIII]   ", cohort, "- FDR differences from original:", n_diff, "\n")
  } else {
    cat("[VIII]   ", cohort, "- Comparing with original FDR from loaded data\n")
    if (cohort == "LUAD") {
      orig_fdr <- cox_luad$FDR[match(recalc$program_id, cox_luad$program_id)]
    } else {
      orig_fdr <- cox_lusc$FDR[match(recalc$program_id, cox_lusc$program_id)]
    }
    compare_df$FDR_orig <- orig_fdr
    fdr_diff <- abs(compare_df$FDR_recalculated - orig_fdr)
    n_diff <- sum(fdr_diff > 1e-10, na.rm = TRUE)
    cat("[VIII]   ", cohort, "- FDR differences:", n_diff, "/", length(fdr_diff), "\n")
  }
  
  # P-value distribution
  p_vals <- recalc$Wald_PValue[filter_mask]
  cat("[VIII]   ", cohort, "- P-value distribution:\n")
  cat("         Min:", min(p_vals, na.rm = TRUE),
      " Median:", median(p_vals, na.rm = TRUE),
      " Max:", max(p_vals, na.rm = TRUE), "\n")
  cat("         P < 0.05:", sum(p_vals < 0.05, na.rm = TRUE), "\n")
  cat("         FDR < 0.05:", sum(fdr_recalc[filter_mask] < 0.05, na.rm = TRUE), "\n")
  cat("         FDR < 0.10:", sum(fdr_recalc[filter_mask] < 0.10, na.rm = TRUE), "\n")
  
  fdr_comparison_all <- rbind(fdr_comparison_all, compare_df)
}

write.csv(fdr_comparison_all,
          "03_results/step08_TCGA/B1_QC/cox/FDR_recalculation_comparison.csv",
          row.names = FALSE)
cat("[VIII] Saved FDR comparison CSV\n\n")


# =========================================================================
# SECTION IX - Meta Recalculation
# =========================================================================
cat("================================================================\n")
cat("SECTION IX - Meta-Analysis Recalculation\n")
cat("================================================================\n")

# Identify programs with FULL_MODEL in both cohorts
luad_full <- cox_recalc_results$LUAD[cox_recalc_results$LUAD$model_level == "FULL_MODEL", ]
lusc_full <- cox_recalc_results$LUSC[cox_recalc_results$LUSC$model_level == "FULL_MODEL", ]
both_full <- intersect(luad_full$program_id, lusc_full$program_id)
cat("[IX] Programs with FULL_MODEL in both cohorts:", length(both_full), "\n")

meta_recalc <- data.frame()

for (pg in both_full) {
  lu_row <- luad_full[luad_full$program_id == pg, ]
  ls_row <- lusc_full[lusc_full$program_id == pg, ]
  
  lu_logHR <- lu_row$logHR
  lu_se <- lu_row$SE
  ls_logHR <- ls_row$logHR
  ls_se <- ls_row$SE
  
  if (any(is.na(c(lu_logHR, lu_se, ls_logHR, ls_se)))) next
  
  # Hand-calculate inverse-variance weighted meta
  w_lu <- 1 / lu_se^2
  w_ls <- 1 / ls_se^2
  meta_logHR <- (w_lu * lu_logHR + w_ls * ls_logHR) / (w_lu + w_ls)
  meta_SE <- sqrt(1 / (w_lu + w_ls))
  meta_Z <- meta_logHR / meta_SE
  meta_P <- 2 * pnorm(-abs(meta_Z))
  
  # metafor::rma.uni comparison
  rma_fit <- tryCatch({
    rma.uni(yi = c(lu_logHR, ls_logHR), sei = c(lu_se, ls_se), method = "EE")
  }, error = function(e) NULL)
  
  if (!is.null(rma_fit)) {
    rma_beta <- rma_fit$beta[[1]]
    rma_se <- rma_fit$se
    rma_pval <- rma_fit$pval
  } else {
    rma_beta <- NA
    rma_se <- NA
    rma_pval <- NA
  }
  
  meta_recalc <- rbind(meta_recalc, data.frame(
    program_id = pg,
    meta_logHR_hand = meta_logHR,
    meta_SE_hand = meta_SE,
    meta_Z_hand = meta_Z,
    meta_P_hand = meta_P,
    meta_logHR_rma = rma_beta,
    meta_SE_rma = rma_se,
    meta_P_rma = rma_pval,
    stringsAsFactors = FALSE
  ))
}

# Compare with original meta results
if ("program_id" %in% names(meta_orig) && "meta_logHR" %in% names(meta_orig)) {
  meta_compare <- merge(meta_orig, meta_recalc, by = "program_id", all = TRUE)
  meta_compare$logHR_diff <- abs(meta_compare$meta_logHR - meta_compare$meta_logHR_hand)
  meta_compare$P_diff <- abs(meta_compare$meta_PValue - meta_compare$meta_P_hand)
  
  cat("[IX] Meta logHR max difference:", max(meta_compare$logHR_diff, na.rm = TRUE), "\n")
  cat("[IX] Meta P-value max difference:", max(meta_compare$P_diff, na.rm = TRUE), "\n")
}

# Compare hand vs rma
rma_diff <- abs(meta_recalc$meta_logHR_hand - meta_recalc$meta_logHR_rma)
cat("[IX] Hand vs rma.uni logHR max diff:", max(rma_diff, na.rm = TRUE), "\n")

rma_p_diff <- abs(meta_recalc$meta_P_hand - meta_recalc$meta_P_rma)
cat("[IX] Hand vs rma.uni P-value max diff:", max(rma_p_diff, na.rm = TRUE), "\n")

# Add heterogeneity (Section X integrated)
cat("[IX] Calculating heterogeneity...\n")
meta_recalc$I2 <- NA
meta_recalc$Cochran_Q <- NA
meta_recalc$heterogeneity_P <- NA

for (i in seq_len(nrow(meta_recalc))) {
  pg <- meta_recalc$program_id[i]
  lu_row <- luad_full[luad_full$program_id == pg, ]
  ls_row <- lusc_full[lusc_full$program_id == pg, ]
  
  lu_logHR <- lu_row$logHR
  lu_se <- lu_row$SE
  ls_logHR <- ls_row$logHR
  ls_se <- ls_row$SE
  
  if (any(is.na(c(lu_logHR, lu_se, ls_logHR, ls_se)))) next
  
  w_lu <- 1 / lu_se^2
  w_ls <- 1 / ls_se^2
  meta_logHR <- (w_lu * lu_logHR + w_ls * ls_logHR) / (w_lu + w_ls)
  
  Q <- w_lu * (lu_logHR - meta_logHR)^2 + w_ls * (ls_logHR - meta_logHR)^2
  meta_recalc$Cochran_Q[i] <- Q
  
  I2 <- max(0, min(100, (Q - 1) / Q * 100))
  meta_recalc$I2[i] <- I2
  
  hetero_p <- pchisq(Q, df = 1, lower.tail = FALSE)
  meta_recalc$heterogeneity_P[i] <- hetero_p
}

write.csv(meta_recalc,
          "03_results/step08_TCGA/B1_QC/meta/GSE243013_TCGA_meta_recalculated.csv",
          row.names = FALSE)
cat("[IX] Saved meta_recalculated.csv\n")
cat("[IX] Meta-analysis recalculation complete\n\n")


# =========================================================================
# SECTION X - Heterogeneity (Cochran Q, I2)
# =========================================================================
cat("================================================================\n")
cat("SECTION X - Heterogeneity Assessment\n")
cat("================================================================\n")

if (nrow(meta_recalc) > 0) {
  cat("[X] Heterogeneity summary:\n")
  cat("    Mean I2:", round(mean(meta_recalc$I2, na.rm = TRUE), 2), "%\n")
  cat("    Median I2:", round(median(meta_recalc$I2, na.rm = TRUE), 2), "%\n")
  cat("    I2 < 25%:", sum(meta_recalc$I2 < 25, na.rm = TRUE), "\n")
  cat("    I2 25-50%:", sum(meta_recalc$I2 >= 25 & meta_recalc$I2 < 50, na.rm = TRUE), "\n")
  cat("    I2 50-75%:", sum(meta_recalc$I2 >= 50 & meta_recalc$I2 < 75, na.rm = TRUE), "\n")
  cat("    I2 >= 75%:", sum(meta_recalc$I2 >= 75, na.rm = TRUE), "\n")
  cat("    Heterogeneity P < 0.05:", sum(meta_recalc$heterogeneity_P < 0.05, na.rm = TRUE), "\n")
}

# Save heterogeneity plot
if (nrow(meta_recalc) > 0 && any(!is.na(meta_recalc$I2))) {
  p_i2 <- ggplot(meta_recalc[!is.na(meta_recalc$I2), ],
                  aes(x = reorder(program_id, I2), y = I2)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 75, linetype = "dashed", color = "darkred") +
    coord_flip() +
    labs(title = "I2 Heterogeneity by Program",
         x = "Program", y = "I2 (%)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6))
  
  ggsave("04_figures/step08_TCGA/B1_QC/X_I2_heterogeneity.pdf",
         p_i2, width = 10, height = max(6, nrow(meta_recalc) * 0.3))
  cat("[X] Saved I2 heterogeneity plot\n")
}

cat("[X] Heterogeneity assessment complete\n\n")


# =========================================================================
# SECTION XI - Meta FDR
# =========================================================================
cat("================================================================\n")
cat("SECTION XI - Meta FDR Recalculation\n")
cat("================================================================\n")

# Apply BH FDR to meta P-values
if ("meta_P_hand" %in% names(meta_recalc)) {
  valid_meta_p <- !is.na(meta_recalc$meta_P_hand) & is.finite(meta_recalc$meta_P_hand)
  meta_recalc$meta_FDR_recalculated <- NA
  meta_recalc$meta_FDR_recalculated[valid_meta_p] <- p.adjust(
    meta_recalc$meta_P_hand[valid_meta_p], method = "BH"
  )
  
  cat("[XI] Meta P-value distribution:\n")
  cat("    Min:", min(meta_recalc$meta_P_hand, na.rm = TRUE), "\n")
  cat("    Median:", median(meta_recalc$meta_P_hand, na.rm = TRUE), "\n")
  cat("    Max:", max(meta_recalc$meta_P_hand, na.rm = TRUE), "\n")
  
  cat("[XI] Counts:\n")
  cat("    P < 0.05:", sum(meta_recalc$meta_P_hand < 0.05, na.rm = TRUE), "\n")
  cat("    FDR < 0.05:", sum(meta_recalc$meta_FDR_recalculated < 0.05, na.rm = TRUE), "\n")
  cat("    FDR < 0.10:", sum(meta_recalc$meta_FDR_recalculated < 0.10, na.rm = TRUE), "\n")
  
  # Compare with original FDR if available
  if ("meta_FDR" %in% names(meta_orig)) {
    fdr_compare <- merge(
      data.frame(program_id = meta_orig$program_id, FDR_orig = meta_orig$meta_FDR),
      data.frame(program_id = meta_recalc$program_id, FDR_recalc = meta_recalc$meta_FDR_recalculated),
      by = "program_id"
    )
    fdr_diff <- abs(fdr_compare$FDR_orig - fdr_compare$FDR_recalc)
    n_diff <- sum(fdr_diff > 1e-10, na.rm = TRUE)
    cat("[XI] FDR differences from original:", n_diff, "/", nrow(fdr_compare), "\n")
  }
}

write.csv(meta_recalc,
          "03_results/step08_TCGA/B1_QC/meta/GSE243013_TCGA_meta_recalculated.csv",
          row.names = FALSE)
cat("[XI] Saved updated meta results with FDR\n\n")


# =========================================================================
# SECTION XII - Column Alignment
# =========================================================================
cat("================================================================\n")
cat("SECTION XII - Column Alignment Check\n")
cat("================================================================\n")

# Check meta_PValue varies across rows (not scalar recycled)
if ("meta_PValue" %in% names(meta_orig)) {
  n_unique_p <- length(unique(meta_orig$meta_PValue))
  if (n_unique_p > 1) {
    cat("[XII] meta_PValue varies across rows:", n_unique_p, "unique values: PASS\n")
  } else {
    cat("[XII] FAIL: meta_PValue appears to be scalar recycled (", n_unique_p, "unique value)\n")
  }
}

# Check meta_FDR varies
if ("meta_FDR" %in% names(meta_orig)) {
  n_unique_fdr <- length(unique(meta_orig$meta_FDR))
  if (n_unique_fdr > 1) {
    cat("[XII] meta_FDR varies across rows:", n_unique_fdr, "unique values: PASS\n")
  } else {
    cat("[XII] FAIL: meta_FDR appears to be scalar recycled\n")
  }
}

# Check values in [0,1]
if ("meta_PValue" %in% names(meta_orig)) {
  p_range <- range(meta_orig$meta_PValue, na.rm = TRUE)
  cat("[XII] meta_PValue range:", p_range[1], "to", p_range[2], "\n")
  if (p_range[1] >= 0 && p_range[2] <= 1) cat("[XII] meta_PValue in [0,1]: PASS\n")
  else cat("[XII] FAIL: meta_PValue outside [0,1]\n")
}

if ("meta_FDR" %in% names(meta_orig)) {
  fdr_range <- range(meta_orig$meta_FDR, na.rm = TRUE)
  cat("[XII] meta_FDR range:", fdr_range[1], "to", fdr_range[2], "\n")
  if (fdr_range[1] >= 0 && fdr_range[2] <= 1) cat("[XII] meta_FDR in [0,1]: PASS\n")
  else cat("[XII] FAIL: meta_FDR outside [0,1]\n")
}

# Align with recalculated values after merge
if (nrow(meta_recalc) > 0) {
  align_merge <- merge(
    meta_orig[, c("program_id", "meta_PValue", "meta_FDR")],
    meta_recalc[, c("program_id", "meta_P_hand", "meta_FDR_recalculated")],
    by = "program_id", all = TRUE
  )
  p_align_diff <- abs(align_merge$meta_PValue - align_merge$meta_P_hand)
  max_p_diff <- max(p_align_diff, na.rm = TRUE)
  cat("[XII] Max meta P-value alignment difference:", max_p_diff, "\n")
  if (max_p_diff < 1e-10) cat("[XII] Meta P-value alignment: PASS\n")
  else cat("[XII] WARNING: Meta P-value alignment differs\n")
  
  fdr_align_diff <- abs(align_merge$meta_FDR - align_merge$meta_FDR_recalculated)
  max_fdr_diff <- max(fdr_align_diff, na.rm = TRUE)
  cat("[XII] Max meta FDR alignment difference:", max_fdr_diff, "\n")
}

cat("[XII] Column alignment check complete\n\n")


# =========================================================================
# SECTION XIII - PH Recalculation
# =========================================================================
cat("================================================================\n")
cat("SECTION XIII - Proportional Hazards Recalculation\n")
cat("================================================================\n")

ph_recalc_results <- list()

for (cohort in c("LUAD", "LUSC")) {
  if (cohort == "LUAD") {
    std_mat <- std_scores_luad
    orig_ph <- ph_results[ph_results$cohort == "LUAD", ]
  } else {
    std_mat <- std_scores_lusc
    orig_ph <- ph_results[ph_results$cohort == "LUSC", ]
  }
  
  common_pats <- intersect(colnames(std_mat), clin$patient_id)
  if (length(common_pats) == 0) {
    common_pats <- intersect(colnames(std_mat), clin$barcode)
  }
  
  clin_sub <- clin[match(common_pats, clin$patient_id), ]
  if (all(is.na(clin_sub$patient_id))) {
    clin_sub <- clin[match(common_pats, clin$barcode), ]
  }
  
  cat("[XIII]   ", cohort, "- Refitting PH for", nrow(std_mat), "programs\n")
  
  cohort_ph <- data.frame()
  
  for (i in seq_len(nrow(std_mat))) {
    pg <- rownames(std_mat)[i]
    score_vec <- as.numeric(std_mat[pg, common_pats])
    
    df_i <- data.frame(
      OS_time = clin_sub$OS_days / 365.25,
      OS_event = clin_sub$OS_event,
      score_z = score_vec,
      age_z = as.numeric(clin_sub$age_z),
      sex_f = clin_sub$sex_f,
      stage_f = clin_sub$stage_f
    )
    
    complete_mask <- complete.cases(df_i)
    df_complete <- df_i[complete_mask, ]
    
    if (nrow(df_complete) < 20 || sum(df_complete$OS_event) < 5) next
    
    stage_complete <- droplevels(df_complete$stage_f)
    
    fit <- tryCatch({
      coxph(Surv(OS_time, OS_event) ~ score_z + age_z + sex_f + stage_complete,
            data = df_complete)
    }, error = function(e) NULL)
    
    if (is.null(fit)) next
    
    zph <- tryCatch({
      cox.zph(fit)
    }, error = function(e) NULL)
    
    if (is.null(zph)) next
    
    zph_table <- zph$table
    
    # Extract score_z PH p-value by name, NOT by row index
    score_ph_p <- NA
    global_ph_p <- NA
    
    if ("score_z" %in% rownames(zph_table)) {
      score_ph_p <- zph_table["score_z", "p"]
    }
    
    # Global test
    if ("GLOBAL" %in% rownames(zph_table)) {
      global_ph_p <- zph_table["GLOBAL", "p"]
    }
    
    cohort_ph <- rbind(cohort_ph, data.frame(
      program_id = pg,
      score_z_PH_p = score_ph_p,
      GLOBAL_PH_p = global_ph_p,
      stringsAsFactors = FALSE
    ))
  }
  
  ph_recalc_results[[cohort]] <- cohort_ph
  
  # Compare with original
  if (nrow(orig_ph) > 0 && nrow(cohort_ph) > 0) {
    ph_merge <- merge(orig_ph, cohort_ph, by = "program_id", all = TRUE)
    if ("score_z_p" %in% names(ph_merge)) {
      ph_diff <- abs(ph_merge$score_z_p - ph_merge$score_z_PH_p)
    } else if ("p" %in% names(ph_merge)) {
      ph_diff <- abs(ph_merge$p - ph_merge$score_z_PH_p)
    } else {
      ph_diff <- NA
    }
    max_ph_diff <- max(ph_diff, na.rm = TRUE)
    cat("[XIII]   ", cohort, "- Max PH p-value difference:", max_ph_diff, "\n")
  }
  
  cat("[XIII]   ", cohort, "- PH recalculation complete:", nrow(cohort_ph), "programs\n")
}

write.csv(do.call(rbind, ph_recalc_results),
          "03_results/step08_TCGA/B1_QC/cox/PH_assumption_recalculated.csv",
          row.names = FALSE)
cat("[XIII] Saved PH recalculation CSV\n\n")


# =========================================================================
# SECTION XIV - Redundancy Analysis
# =========================================================================
cat("================================================================\n")
cat("SECTION XIV - Redundancy Analysis (145x145 Correlation)\n")
cat("================================================================\n")

redundancy_results <- list()

for (cohort in c("LUAD", "LUSC")) {
  if (cohort == "LUAD") {
    std_mat <- std_scores_luad
  } else {
    std_mat <- std_scores_lusc
  }
  
  n_prog <- nrow(std_mat)
  cat("[XIV]   ", cohort, "- Computing", n_prog, "x", n_prog, "Spearman correlation matrix\n")
  
  # 145x145 Spearman correlation (standardized within cohort - already done)
  rho_mat <- cor(t(std_mat), method = "spearman", use = "pairwise.complete.obs")
  
  # Count pairs at various thresholds
  upper_tri <- rho_mat[upper.tri(rho_mat)]
  n_ge_070 <- sum(abs(upper_tri) >= 0.70)
  n_ge_080 <- sum(abs(upper_tri) >= 0.80)
  n_ge_090 <- sum(abs(upper_tri) >= 0.90)
  total_pairs <- length(upper_tri)
  
  cat("[XIV]   ", cohort, " pairs |rho|>=0.70:", n_ge_070, "/", total_pairs,
      " (", round(n_ge_070/total_pairs*100, 2), "%)\n")
  cat("[XIV]   ", cohort, " pairs |rho|>=0.80:", n_ge_080, "/", total_pairs,
      " (", round(n_ge_080/total_pairs*100, 2), "%)\n")
  cat("[XIV]   ", cohort, " pairs |rho|>=0.90:", n_ge_090, "/", total_pairs,
      " (", round(n_ge_090/total_pairs*100, 2), "%)\n")
  
  # Cluster with hclust
  dist_mat <- as.dist(1 - abs(rho_mat))
  hc <- hclust(dist_mat, method = "complete")
  clusters <- cutree(hc, h = 0.20)
  n_clusters <- max(clusters)
  cat("[XIV]   ", cohort, "- Number of clusters (h=0.20):", n_clusters, "\n")
  
  # Identify representatives per cluster
  cluster_representatives <- data.frame()
  for (cl in seq_len(n_clusters)) {
    members <- names(clusters[clusters == cl])
    
    # Priority: tier, meta FDR, gene count
    rep_info <- data.frame(
      program_id = members,
      cluster = cl,
      stringsAsFactors = FALSE
    )
    
    # Add tier info if available
    if ("tier" %in% names(program_manifest)) {
      rep_info <- merge(rep_info,
                        program_manifest[, c("program_id", "tier")],
                        by = "program_id", all.x = TRUE)
    } else {
      rep_info$tier <- NA
    }
    
    # Add meta FDR
    if (nrow(meta_recalc) > 0) {
      rep_info <- merge(rep_info,
                        meta_recalc[, c("program_id", "meta_FDR_recalculated")],
                        by = "program_id", all.x = TRUE)
    } else {
      rep_info$meta_FDR_recalculated <- NA
    }
    
    # Add gene count
    gene_counts <- table(gene_membership$program_id)
    rep_info$n_genes <- gene_counts[rep_info$program_id]
    
    # Select representative: lowest tier number, then lowest meta FDR
    rep_info <- rep_info[order(rep_info$tier, rep_info$meta_FDR_recalculated,
                                -rep_info$n_genes), ]
    best <- rep_info[1, ]
    
    cluster_representatives <- rbind(cluster_representatives, data.frame(
      cluster = cl,
      representative = best$program_id,
      n_members = length(members),
      members = paste(members, collapse = ";"),
      stringsAsFactors = FALSE
    ))
  }
  
  # Save correlation matrix
  write.csv(as.data.frame(rho_mat),
            paste0("03_results/step08_TCGA/B1_QC/redundancy/", cohort,
                   "_145x145_Spearman_rho_matrix.csv"),
            row.names = TRUE)
  
  # Save cluster representatives
  write.csv(cluster_representatives,
            paste0("03_results/step08_TCGA/B1_QC/redundancy/", cohort,
                   "_cluster_representatives.csv"),
            row.names = FALSE)
  
  redundancy_results[[cohort]] <- list(
    rho_mat = rho_mat,
    clusters = clusters,
    n_clusters = n_clusters,
    representatives = cluster_representatives,
    n_ge_070 = n_ge_070, n_ge_080 = n_ge_080, n_ge_090 = n_ge_090
  )
  
  # Save heatmap
  pdf(paste0("04_figures/step08_TCGA/B1_QC/redundancy/", cohort,
             "_145x145_rho_heatmap.pdf"), width = 14, height = 12)
  pheatmap(rho_mat,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-1, 1, length.out = 101),
           cluster_rows = TRUE, cluster_cols = TRUE,
           main = paste(cohort, "145x145 Spearman rho Heatmap"),
           show_rownames = FALSE, show_colnames = FALSE)
  dev.off()
}

cat("[XIV] Redundancy analysis complete\n\n")


# =========================================================================
# SECTION XV - Jaccard Index
# =========================================================================
cat("================================================================\n")
cat("SECTION XV - Jaccard Index (Leading-Edge Genes)\n")
cat("================================================================\n")

# Build gene sets
gene_sets <- split(gene_membership$gene, gene_membership$program_id)
gene_sets <- lapply(gene_sets, unique)
all_programs <- names(gene_sets)

n_comparisons <- length(all_programs) * (length(all_programs) - 1) / 2
cat("[XV] Computing Jaccard for", n_comparisons, "program pairs\n")

jaccard_pairs <- data.frame()
thresholds <- c(0.25, 0.50, 0.75, 1.0)
count_at_threshold <- setNames(rep(0, length(thresholds)), paste0(">=", thresholds))

for (i in seq_len(length(all_programs) - 1)) {
  for (j in seq(i + 1, length(all_programs))) {
    pg1 <- all_programs[i]
    pg2 <- all_programs[j]
    
    genes1 <- gene_sets[[pg1]]
    genes2 <- gene_sets[[pg2]]
    
    intersection <- length(intersect(genes1, genes2))
    union_size <- length(union(genes1, genes2))
    
    if (union_size == 0) next
    jaccard <- intersection / union_size
    
    # Count thresholds
    for (th in thresholds) {
      if (jaccard >= th) count_at_threshold[paste0(">=", th)] <- count_at_threshold[paste0(">=", th)] + 1
    }
    
    # Save pairs at >= 0.25
    if (jaccard >= 0.25) {
      jaccard_pairs <- rbind(jaccard_pairs, data.frame(
        program_1 = pg1,
        program_2 = pg2,
        jaccard = jaccard,
        n_intersection = intersection,
        n_union = union_size,
        genes_intersection = paste(intersect(genes1, genes2), collapse = ";"),
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("[XV] Jaccard threshold counts:\n")
for (th in thresholds) {
  cat("    Jaccard >= ", th, ": ", count_at_threshold[paste0(">=", th)], " pairs\n")
}

write.csv(jaccard_pairs,
          "03_results/step08_TCGA/B1_QC/redundancy/jaccard_pairs_gt_0.25.csv",
          row.names = FALSE)
cat("[XV] Saved jaccard_pairs_gt_0.25.csv\n")

# Save histogram
p_jaccard <- ggplot(jaccard_pairs, aes(x = jaccard)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = c(0.25, 0.50, 0.75), linetype = "dashed", color = c("orange", "red", "darkred")) +
  labs(title = "Jaccard Index Distribution (Program Pairs)",
       x = "Jaccard Index", y = "Count") +
  theme_minimal()
ggsave("04_figures/step08_TCGA/B1_QC/redundancy/XV_Jaccard_distribution.pdf",
       p_jaccard, width = 8, height = 5)

cat("[XV] Jaccard analysis complete\n\n")


# =========================================================================
# SECTION XVI - PCA Analysis
# =========================================================================
cat("================================================================\n")
cat("SECTION XVI - PCA on 145-Program Score Matrix\n")
cat("================================================================\n")

pca_results <- list()

for (cohort in c("LUAD", "LUSC")) {
  if (cohort == "LUAD") {
    std_mat <- std_scores_luad
  } else {
    std_mat <- std_scores_lusc
  }
  
  # Ensure 145 programs
  common_progs <- intersect(rownames(std_mat), program_ids_manifest)
  mat_145 <- std_mat[common_progs, ]
  
  # Handle NAs: remove programs with any NA values
  na_per_program <- rowSums(is.na(mat_145))
  valid_progs <- na_per_program == 0
  mat_valid <- mat_145[valid_progs, ]
  
  cat("[XVI]   ", cohort, "- PCA on", nrow(mat_valid), "programs (after removing", sum(!valid_progs), "with NAs) x", ncol(mat_valid), "patients\n")
  
  # PCA
  pca_fit <- prcomp(t(mat_valid), center = TRUE, scale. = FALSE)
  
  # PC1-3 variance
  pc_var <- summary(pca_fit)$importance[2, 1:3]
  cat("[XVI]   ", cohort, "- PC1 variance:", round(pc_var[1] * 100, 2), "%\n")
  cat("[XVI]   ", cohort, "- PC2 variance:", round(pc_var[2] * 100, 2), "%\n")
  cat("[XVI]   ", cohort, "- PC3 variance:", round(pc_var[3] * 100, 2), "%\n")
  
  # Correlate each program with PC1
  pc1_scores <- pca_fit$x[, 1]
  pc1_correlations <- sapply(seq_len(nrow(mat_145)), function(i) {
    cor(as.numeric(mat_145[i, ]), pc1_scores, method = "spearman")
  })
  names(pc1_correlations) <- rownames(mat_145)
  
  n_high_cor <- sum(abs(pc1_correlations) > 0.5, na.rm = TRUE)
  cat("[XVI]   ", cohort, "- Programs with |rho|>0.5 to PC1:", n_high_cor, "\n")
  
  # Fit Surv ~ PC1 + covariates
  common_pats <- intersect(colnames(mat_145), clin$patient_id)
  if (length(common_pats) == 0) {
    common_pats <- intersect(colnames(mat_145), clin$barcode)
  }
  clin_sub <- clin[match(common_pats, clin$patient_id), ]
  if (all(is.na(clin_sub$patient_id))) {
    clin_sub <- clin[match(common_pats, clin$barcode), ]
  }
  
  pc1_df <- data.frame(
    OS_time = clin_sub$OS_days / 365.25,
    OS_event = clin_sub$OS_event,
    PC1 = pca_fit$x[common_pats, 1],
    age_z = as.numeric(clin_sub$age_z),
    sex = clin_sub$sex_f,
    stage_factor = droplevels(clin_sub$stage_f)
  )
  pc1_df <- pc1_df[complete.cases(pc1_df), ]
  
  pc1_cox <- tryCatch({
    coxph(Surv(OS_time, OS_event) ~ PC1 + age_z + sex + stage_factor, data = pc1_df)
  }, error = function(e) NULL)
  
  if (!is.null(pc1_cox)) {
    pc1_summary <- summary(pc1_cox)$coefficients
    if ("PC1" %in% rownames(pc1_summary)) {
      cat("[XVI]   ", cohort, "- PC1 Cox HR:", round(exp(pc1_summary["PC1", "coef"]), 3),
          " P:", pc1_summary["PC1", "Pr(>|z|)"], "\n")
    }
  }
  
  # Save PC scores
  write.csv(data.frame(program_id = names(pc1_correlations),
                        rho_with_PC1 = pc1_correlations),
            paste0("03_results/step08_TCGA/B1_QC/redundancy/", cohort,
                   "_PC1_program_correlations.csv"),
            row.names = FALSE)
  
  pca_results[[cohort]] <- list(
    fit = pca_fit,
    pc_var = pc_var,
    pc1_correlations = pc1_correlations
  )
  
  # Scree plot
  pdf(paste0("04_figures/step08_TCGA/B1_QC/redundancy/", cohort, "_PCA_scree.pdf"),
      width = 8, height = 5)
  plot(summary(pca_fit)$importance[2, ], type = "b", pch = 19,
       xlab = "Principal Component", ylab = "Proportion of Variance",
       main = paste(cohort, "PCA Scree Plot"))
  dev.off()
}

cat("[XVI] PCA analysis complete\n\n")


# =========================================================================
# SECTION XVII - Negative Controls (Permutation Test)
# =========================================================================
cat("================================================================\n")
cat("SECTION XVII - Negative Controls (Permutation Test)\n")
cat("================================================================\n")

set.seed(2026804)

# Select 20 programs: 10 lowest FDR, 5 median, 5 highest FDR
# Use LUAD as reference cohort
luad_cox <- cox_recalc_results$LUAD
luad_valid <- luad_cox[luad_cox$model_level == "FULL_MODEL" &
                        !is.na(luad_cox$FDR) & is.finite(luad_cox$FDR), ]

luad_valid <- luad_valid[order(luad_valid$FDR), ]
n_valid <- nrow(luad_valid)

# 10 lowest FDR
low_fdr <- head(luad_valid, 10)
# 5 median
mid_start <- floor(n_valid / 2) - 2
mid_fdr <- luad_valid[seq(mid_start, mid_start + 4), ]
# 5 highest FDR
high_fdr <- tail(luad_valid, 5)

selected_programs <- rbind(low_fdr, mid_fdr, high_fdr)
cat("[XVII] Selected", nrow(selected_programs), "programs for permutation testing\n")
cat("[XVII]   Low FDR programs:", paste(low_fdr$program_id, collapse = ", "), "\n")
cat("[XVII]   Mid FDR programs:", paste(mid_fdr$program_id, collapse = ", "), "\n")
cat("[XVII]   High FDR programs:", paste(high_fdr$program_id, collapse = ", "), "\n")

# Prepare data
std_mat <- std_scores_luad
common_pats <- intersect(colnames(std_mat), clin$patient_id)
if (length(common_pats) == 0) {
  common_pats <- intersect(colnames(std_mat), clin$barcode)
}
clin_sub <- clin[match(common_pats, clin$patient_id), ]
if (all(is.na(clin_sub$patient_id))) {
  clin_sub <- clin[match(common_pats, clin$barcode), ]
}

n_perm <- 50
perm_results <- data.frame()

for (pg in selected_programs$program_id) {
  cat("[XVII] Permuting program:", pg, "\n")
  
  score_vec <- as.numeric(std_mat[pg, common_pats])
  orig_logHR <- selected_programs$logHR[selected_programs$program_id == pg]
  orig_P <- selected_programs$Wald_PValue[selected_programs$program_id == pg]
  
  perm_pvals <- numeric(n_perm)
  perm_loghrs <- numeric(n_perm)
  
  for (perm in seq_len(n_perm)) {
    # Shuffle score_z among patients, keep OS and covariates
    shuffled_score <- sample(score_vec)
    
    df_perm <- data.frame(
      OS_time = clin_sub$OS_days / 365.25,
      OS_event = clin_sub$OS_event,
      score_z = shuffled_score,
      age_z = as.numeric(clin_sub$age_z),
      sex_f = clin_sub$sex_f,
      stage_f = clin_sub$stage_f
    )
    
    complete_mask <- complete.cases(df_perm)
    df_perm <- df_perm[complete_mask, ]
    
    if (nrow(df_perm) < 20 || sum(df_perm$OS_event) < 5) {
      perm_pvals[perm] <- NA
      perm_loghrs[perm] <- NA
      next
    }
    
    stage_perm <- droplevels(df_perm$stage_f)
    
    fit_perm <- tryCatch({
      coxph(Surv(OS_time, OS_event) ~ score_z + age_z + sex_f + stage_perm, data = df_perm)
    }, error = function(e) NULL)
    
    if (is.null(fit_perm)) {
      perm_pvals[perm] <- NA
      perm_loghrs[perm] <- NA
      next
    }
    
    coef_perm <- summary(fit_perm)$coefficients
    if ("score_z" %in% rownames(coef_perm)) {
      perm_pvals[perm] <- coef_perm["score_z", "Pr(>|z|)"]
      perm_loghrs[perm] <- coef_perm["score_z", "coef"]
    } else {
      perm_pvals[perm] <- NA
      perm_loghrs[perm] <- NA
    }
  }
  
  perm_results <- rbind(perm_results, data.frame(
    program_id = pg,
    orig_logHR = orig_logHR,
    orig_P = orig_P,
    perm_mean_logHR = mean(perm_loghrs, na.rm = TRUE),
    perm_sd_logHR = sd(perm_loghrs, na.rm = TRUE),
    perm_mean_P = mean(perm_pvals, na.rm = TRUE),
    n_valid_perms = sum(!is.na(perm_pvals)),
    empirical_P = mean(abs(perm_loghrs) >= abs(orig_logHR), na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
}

write.csv(perm_results,
          "03_results/step08_TCGA/B1_QC/negative_control/permutation_test_results.csv",
          row.names = FALSE)

# Save permutation plot
p_perm <- ggplot(perm_results, aes(x = program_id, y = empirical_P)) +
  geom_point(size = 3, color = "steelblue") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "orange") +
  coord_flip() +
  labs(title = "Empirical P-values from Permutation Test",
       x = "Program", y = "Empirical P-value") +
  theme_minimal()
ggsave("04_figures/step08_TCGA/B1_QC/negative_control/XVII_permutation_results.pdf",
       p_perm, width = 8, height = 6)

cat("[XVII] Negative control permutation test complete\n\n")


# =========================================================================
# SECTION XVIII - Clinical Evidence Levels
# =========================================================================
cat("================================================================\n")
cat("SECTION XVIII - Clinical Evidence Levels\n")
cat("================================================================\n")

# Merge LUAD and LUSC recalculated results
luad_rec <- cox_recalc_results$LUAD
lusc_rec <- cox_recalc_results$LUSC
luad_rec$cohort <- "LUAD"
lusc_rec$cohort <- "LUSC"

# Also merge with original concordance data for rho >= 0.60 check
# Recalculate rho from ssGSEA/GSVA if needed
rho_luad <- data.frame()
rho_lusc <- data.frame()

for (pg in rownames(ssgsea_luad)) {
  if (pg %in% rownames(gsva_luad)) {
    x <- as.numeric(ssgsea_luad[pg, ])
    y <- as.numeric(gsva_luad[pg, ])
    complete <- !is.na(x) & !is.na(y)
    if (sum(complete) >= 3) {
      rho <- cor(x[complete], y[complete], method = "spearman")
    } else {
      rho <- NA_real_
    }
    rho_luad <- rbind(rho_luad, data.frame(program_id = pg, rho_ssgsea_gsva = rho,
                                            stringsAsFactors = FALSE))
  }
}
for (pg in rownames(ssgsea_lusc)) {
  if (pg %in% rownames(gsva_lusc)) {
    x <- as.numeric(ssgsea_lusc[pg, ])
    y <- as.numeric(gsva_lusc[pg, ])
    complete <- !is.na(x) & !is.na(y)
    if (sum(complete) >= 3) {
      rho <- cor(x[complete], y[complete], method = "spearman")
    } else {
      rho <- NA_real_
    }
    rho_lusc <- rbind(rho_lusc, data.frame(program_id = pg, rho_ssgsea_gsva = rho,
                                            stringsAsFactors = FALSE))
  }
}

# Build clinical evidence level table
all_programs <- unique(c(luad_rec$program_id, lusc_rec$program_id))

evidence_levels <- data.frame(
  program_id = all_programs,
  Level_A = FALSE,
  Level_B = FALSE,
  Level_C = FALSE,
  No_clinical_support = TRUE,
  evidence_reason = "",
  stringsAsFactors = FALSE
)

for (pg in all_programs) {
  lu_match <- luad_rec[luad_rec$program_id == pg, ]
  ls_match <- lusc_rec[lusc_rec$program_id == pg, ]
  
  lu_fdr <- ifelse(nrow(lu_match) > 0 && !is.na(lu_match$FDR), lu_match$FDR, NA)
  ls_fdr <- ifelse(nrow(ls_match) > 0 && !is.na(ls_match$FDR), ls_match$FDR, NA)
  
  lu_logHR <- ifelse(nrow(lu_match) > 0, lu_match$logHR, NA)
  ls_logHR <- ifelse(nrow(ls_match) > 0, ls_match$logHR, NA)
  lu_model <- ifelse(nrow(lu_match) > 0, lu_match$model_level, NA)
  ls_model <- ifelse(nrow(ls_match) > 0, ls_match$model_level, NA)
  
  # Concordance rho >= 0.60
  lu_rho <- rho_luad$rho_ssgsea_gsva[rho_luad$program_id == pg]
  ls_rho <- rho_lusc$rho_ssgsea_gsva[rho_lusc$program_id == pg]
  if (length(lu_rho) == 0) lu_rho <- NA
  if (length(ls_rho) == 0) ls_rho <- NA
  
  # Meta results
  meta_match <- meta_recalc[meta_recalc$program_id == pg, ]
  meta_fdr <- ifelse(nrow(meta_match) > 0, meta_match$meta_FDR_recalculated, NA)
  i2 <- ifelse(nrow(meta_match) > 0, meta_match$I2, NA)
  
  # PH results
  ph_lu <- ph_recalc_results$LUAD[ph_recalc_results$LUAD$program_id == pg, ]
  ph_ls <- ph_recalc_results$LUSC[ph_recalc_results$LUSC$program_id == pg, ]
  
  ph_lu_pass <- ifelse(nrow(ph_lu) > 0 && !is.na(ph_lu$GLOBAL_PH_p) && ph_lu$GLOBAL_PH_p > 0.05, TRUE, FALSE)
  ph_ls_pass <- ifelse(nrow(ph_ls) > 0 && !is.na(ph_ls$GLOBAL_PH_p) && ph_ls$GLOBAL_PH_p > 0.05, TRUE, FALSE)
  ph_pass <- ph_lu_pass || ph_ls_pass
  
  # Direction consistent
  direction_consistent <- (lu_logHR > 0 && ls_logHR > 0) || (lu_logHR < 0 && ls_logHR < 0)
  if (is.na(lu_logHR) || is.na(ls_logHR)) direction_consistent <- NA
  
  # Evidence level rules
  reason <- ""
  
  # Level_A: both cohorts rho>=0.60, direction consistent, meta FDR<0.05, I2<50%, PH pass
  lu_rho_pass <- !is.na(lu_rho) && lu_rho >= 0.60
  ls_rho_pass <- !is.na(ls_rho) && ls_rho >= 0.60
  
  if (lu_rho_pass && ls_rho_pass && !is.na(direction_consistent) && direction_consistent &&
      !is.na(meta_fdr) && meta_fdr < 0.05 && !is.na(i2) && i2 < 50 && ph_pass) {
    evidence_levels$Level_A[evidence_levels$program_id == pg] <- TRUE
    evidence_levels$No_clinical_support[evidence_levels$program_id == pg] <- FALSE
    reason <- "Level_A: both cohorts rho>=0.60, direction consistent, meta FDR<0.05, I2<50%, PH pass"
  }
  # Level_B: one cohort FDR<0.05 + other direction consistent, or meta FDR<0.10 + I2<50%
  else {
    one_fdr_low <- (!is.na(lu_fdr) && lu_fdr < 0.05 && !is.na(direction_consistent) && direction_consistent) ||
      (!is.na(ls_fdr) && ls_fdr < 0.05 && !is.na(direction_consistent) && direction_consistent)
    meta_moderate <- !is.na(meta_fdr) && meta_fdr < 0.10 && !is.na(i2) && i2 < 50
    
    if (one_fdr_low || meta_moderate) {
      evidence_levels$Level_B[evidence_levels$program_id == pg] <- TRUE
      evidence_levels$No_clinical_support[evidence_levels$program_id == pg] <- FALSE
      reason <- "Level_B: "
      if (one_fdr_low) reason <- paste0(reason, "one cohort FDR<0.05+direction consistent; ")
      if (meta_moderate) reason <- paste0(reason, "meta FDR<0.10+I2<50%; ")
    }
    # Level_C
    else if (!is.na(lu_fdr) && !is.na(ls_fdr)) {
      evidence_levels$Level_C[evidence_levels$program_id == pg] <- TRUE
      evidence_levels$No_clinical_support[evidence_levels$program_id == pg] <- FALSE
      reason <- "Level_C: stage-only or reduced model or PH warning or discordant"
    }
  }
  
  evidence_levels$evidence_reason[evidence_levels$program_id == pg] <- reason
}

cat("[XVIII] Evidence level counts:\n")
cat("    Level_A:", sum(evidence_levels$Level_A), "\n")
cat("    Level_B:", sum(evidence_levels$Level_B), "\n")
cat("    Level_C:", sum(evidence_levels$Level_C), "\n")
cat("    No_clinical_support:", sum(evidence_levels$No_clinical_support), "\n")

write.csv(evidence_levels,
          "03_results/step08_TCGA/B1_QC/clinical_evidence_levels.csv",
          row.names = FALSE)
cat("[XVIII] Saved clinical_evidence_levels.csv\n\n")


# =========================================================================
# SECTION XIX - Audit Status
# =========================================================================
cat("================================================================\n")
cat("SECTION XIX - Audit Status Determination\n")
cat("================================================================\n")

audit_status <- "PASS"
audit_flags <- character()

# Check Cox extraction
cox_pass_rate <- overall_pass_rate
if (cox_pass_rate < 0.95) {
  audit_status <- "FAIL_COX_EXTRACTION"
  audit_flags <- c(audit_flags, paste("Cox PASS rate:", round(cox_pass_rate * 100, 2), "%"))
}

# Check meta calculation
if (nrow(meta_recalc) > 0) {
  meta_max_diff <- max(abs(meta_recalc$meta_logHR_hand - meta_recalc$meta_logHR_rma), na.rm = TRUE)
  if (meta_max_diff > 1e-6) {
    audit_status <- "FAIL_META_CALCULATION"
    audit_flags <- c(audit_flags, paste("Meta hand vs rma diff:", meta_max_diff))
  }
}

# Check FDR calculation
if (exists("fdr_comparison_all") && nrow(fdr_comparison_all) > 0) {
  if ("FDR_orig" %in% names(fdr_comparison_all)) {
    fdr_n_diff <- sum(abs(fdr_comparison_all$FDR_recalculated - fdr_comparison_all$FDR_orig) > 1e-10,
                      na.rm = TRUE)
    if (fdr_n_diff > 0) {
      audit_status <- "FAIL_FDR_CALCULATION"
      audit_flags <- c(audit_flags, paste("FDR differences:", fdr_n_diff))
    }
  }
}

# Check program mapping
if (n_unique != 145 || n_validated != 145) {
  audit_status <- "FAIL_PROGRAM_MAPPING"
  audit_flags <- c(audit_flags, paste("Manifest:", n_unique, "Validated:", n_validated))
}

# Check redundancy
if (!is.null(redundancy_results$LUAD)) {
  n_high_rho_luad <- redundancy_results$LUAD$n_ge_080
  n_high_rho_lusc <- redundancy_results$LUSC$n_ge_080
  total_pairs_luad <- nrow(std_scores_luad) * (nrow(std_scores_luad) - 1) / 2
  pct_high <- (n_high_rho_luad + n_high_rho_lusc) / (2 * total_pairs_luad) * 100
  
  if (audit_status == "PASS" && pct_high > 20) {
    audit_status <- "WARN_REDUNDANT_PROGRAMS"
    audit_flags <- c(audit_flags, paste("High rho pairs:", round(pct_high, 2), "%"))
  }
}

cat("[XIX] Final Audit Status:", audit_status, "\n")
if (length(audit_flags) > 0) {
  cat("[XIX] Flags:\n")
  for (flag in audit_flags) cat("    -", flag, "\n")
}

# Save audit status
status_df <- data.frame(
  status = audit_status,
  flags = paste(audit_flags, collapse = "; "),
  timestamp = as.character(Sys.time()),
  cox_pass_rate = ifelse(exists("overall_pass_rate"), overall_pass_rate, NA),
  n_meta_programs = nrow(meta_recalc),
  n_level_a = sum(evidence_levels$Level_A),
  n_level_b = sum(evidence_levels$Level_B),
  stringsAsFactors = FALSE
)
write.csv(status_df, "03_results/step08_TCGA/B1_QC/audit_status.csv", row.names = FALSE)
cat("[XIX] Saved audit_status.csv\n\n")


# =========================================================================
# SECTION XX - B2 Candidate List
# =========================================================================
cat("================================================================\n")
cat("SECTION XX - B2 Candidate List\n")
cat("================================================================\n")

if (audit_status %in% c("PASS", "WARN_REDUNDANT_PROGRAMS")) {
  cat("[XX] Audit status allows B2 candidate generation:", audit_status, "\n")
  
  # Build candidate list
  candidates <- evidence_levels[evidence_levels$Level_A | evidence_levels$Level_B, ]
  
  # Add tier info
  if ("tier" %in% names(program_manifest)) {
    candidates <- merge(candidates, program_manifest[, c("program_id", "tier")],
                         by = "program_id", all.x = TRUE)
  } else {
    candidates$tier <- NA
  }
  
  # Add PH pass
  candidates$ph_pass <- FALSE
  for (i in seq_len(nrow(candidates))) {
    pg <- candidates$program_id[i]
    ph_lu <- ph_recalc_results$LUAD[ph_recalc_results$LUAD$program_id == pg, ]
    ph_ls <- ph_recalc_results$LUSC[ph_recalc_results$LUSC$program_id == pg, ]
    if (nrow(ph_lu) > 0 && !is.na(ph_lu$GLOBAL_PH_p) && ph_lu$GLOBAL_PH_p > 0.05) {
      candidates$ph_pass[i] <- TRUE
    }
    if (nrow(ph_ls) > 0 && !is.na(ph_ls$GLOBAL_PH_p) && ph_ls$GLOBAL_PH_p > 0.05) {
      candidates$ph_pass[i] <- TRUE
    }
  }
  
  # Direction concordant
  candidates$direction_concordant <- FALSE
  for (i in seq_len(nrow(candidates))) {
    pg <- candidates$program_id[i]
    lu_h <- cox_recalc_results$LUAD$logHR[cox_recalc_results$LUAD$program_id == pg]
    ls_h <- cox_recalc_results$LUSC$logHR[cox_recalc_results$LUSC$program_id == pg]
    if (length(lu_h) > 0 && length(ls_h) > 0 &&
        !is.na(lu_h) && !is.na(ls_h) &&
        ((lu_h > 0 && ls_h > 0) || (lu_h < 0 && ls_h < 0))) {
      candidates$direction_concordant[i] <- TRUE
    }
  }
  
  # Add cluster representative info
  candidates$cluster_rep <- FALSE
  for (cohort in c("LUAD", "LUSC")) {
    if (!is.null(redundancy_results[[cohort]]$representatives)) {
      reps <- redundancy_results[[cohort]]$representatives$representative
      candidates$cluster_rep[candidates$program_id %in% reps] <- TRUE
    }
  }
  
  # Priority sorting: Tier 1, Level_A, Level_B, PH pass, direction concordant, one per cluster
  candidates <- candidates[order(
    -candidates$Level_A,          # Level_A first
    candidates$tier,               # Lower tier number first
    -candidates$ph_pass,           # PH pass first
    -candidates$direction_concordant,  # Concordant first
    -candidates$cluster_rep        # Cluster representatives first
  ), ]
  
  # Target 20-50 programs
  if (nrow(candidates) > 50) {
    cat("[XX] More than", nrow(candidates), "candidates. Creating top 50 list.\n")
    candidates_top50 <- head(candidates, 50)
    write.csv(candidates_top50,
              "03_results/step08_TCGA/B1_QC/B2_candidate_list_top50.csv",
              row.names = FALSE)
    cat("[XX] Saved B2_candidate_list_top50.csv\n")
  }
  
  write.csv(candidates,
            "03_results/step08_TCGA/B1_QC/B2_candidate_list.csv",
            row.names = FALSE)
  cat("[XX] Saved B2_candidate_list.csv\n")
  cat("[XX] Total candidates:", nrow(candidates), "\n")
  cat("[XX]   Level_A:", sum(candidates$Level_A), "\n")
  cat("[XX]   Level_B:", sum(candidates$Level_B), "\n")
} else {
  cat("[XX] Audit status FAIL - B2 candidate list NOT generated\n")
  cat("[XX] Status:", audit_status, "\n")
}

cat("[XX] B2 candidate list complete\n\n")


# =========================================================================
# SECTION XXI - Completion Conditions
# =========================================================================
cat("================================================================\n")
cat("SECTION XXI - Completion Conditions Check\n")
cat("================================================================\n")

conditions <- list()
conditions$prerequisite_file <- file.exists("03_results/GSE243013_step08B1_COMPLETE.txt")
conditions$all_packages_loaded <- TRUE
conditions$program_id_integrity <- (n_unique == 145 && n_validated == 145 &&
                                     identical(sort(manifest_set), sort(validated_set)))
conditions$cox_refitting <- (nrow(cox_recalc_results$LUAD) > 0 && nrow(cox_recalc_results$LUSC) > 0)
conditions$cox_comparison <- (overall_pass_rate > 0.95)
conditions$fdr_recalculated <- TRUE
conditions$meta_recalculated <- (nrow(meta_recalc) > 0)
conditions$heterogeneity_computed <- all(!is.na(meta_recalc$I2))
conditions$meta_fdr_applied <- all(!is.na(meta_recalc$meta_FDR_recalculated))
conditions$ph_recalculated <- (nrow(ph_recalc_results$LUAD) > 0 && nrow(ph_recalc_results$LUSC) > 0)
conditions$redundancy_analyzed <- !is.null(redundancy_results$LUAD)
conditions$jaccard_computed <- (nrow(jaccard_pairs) > 0)
conditions$pca_completed <- !is.null(pca_results$LUAD)
conditions$negative_controls_run <- (nrow(perm_results) > 0)
conditions$evidence_levels_assigned <- (nrow(evidence_levels) > 0)
conditions$audit_status_determined <- TRUE
conditions$b2_candidate_list <- (audit_status %in% c("PASS", "WARN_REDUNDANT_PROGRAMS"))
conditions$no_fatal_failures <- !grepl("FAIL", audit_status)

n_met <- sum(unlist(conditions))
n_total <- length(conditions)

cat("[XXI] Completion conditions:", n_met, "/", n_total, " met\n\n")

for (nm in names(conditions)) {
  status <- ifelse(conditions[[nm]], "PASS", "FAIL")
  cat("    [", status, "] ", nm, "\n")
}

# Write completion file
completion_df <- data.frame(
  condition = names(conditions),
  met = unlist(conditions),
  status = ifelse(unlist(conditions), "PASS", "FAIL"),
  stringsAsFactors = FALSE
)
write.csv(completion_df,
          "03_results/step08_TCGA/B1_QC/completion_conditions.csv",
          row.names = FALSE)

# Write completion marker
completion_marker <- data.frame(
  step = "08B1_QC",
  status = ifelse(all(unlist(conditions)), "COMPLETE", "INCOMPLETE"),
  timestamp = as.character(Sys.time()),
  audit_status = audit_status,
  conditions_met = n_met,
  conditions_total = n_total,
  stringsAsFactors = FALSE
)
write.csv(completion_marker,
          "03_results/step08_TCGA/B1_QC/QC_COMPLETE.txt",
          row.names = FALSE)

cat("[XXI] Completion conditions check complete\n\n")


# =========================================================================
# SECTION XXII - Final Report
# =========================================================================
cat("================================================================\n")
cat("SECTION XXII - Final Report (23 Items)\n")
cat("================================================================\n\n")

cat("1.  Script: 08B1_QC_audit_TCGA_survival_meta_results.R\n")
cat("2.  Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("3.  Audit Status:", audit_status, "\n")
cat("4.  Completion: ", n_met, "/", n_total, " conditions met\n\n")

cat("5.  Program Manifest:", n_unique, "unique program_ids\n")
cat("6.  Validated Programs:", n_validated, "unique program_ids\n")
cat("7.  Manifest-Validated Match:", ifelse(identical(sort(manifest_set), sort(validated_set)), "YES", "NO"), "\n")
cat("8.  Identical Gene Sets (different IDs):", ifelse(nrow(gene_dup) > 0, "FOUND", "NONE"), "\n\n")

cat("9.  ssGSEA LUAD:", nrow(ssgsea_luad), "x", ncol(ssgsea_luad), "\n")
cat("10. ssGSEA LUSC:", nrow(ssgsea_lusc), "x", ncol(ssgsea_lusc), "\n")
cat("11. GSVA/ssGSEA Dimension Match:", ifelse(gsva_dims_match, "YES", "NO"), "\n")
cat("12. Theoretical Denominator:", theoretical_max, "\n")
cat("13. Programs with rho >= 0.60:", n_rho_ge_06, "/", nrow(recalc_concordance), "\n\n")

cat("14. Cox Refitting:\n")
cat("    LUAD FULL_MODEL:", sum(cox_recalc_results$LUAD$model_level == "FULL_MODEL", na.rm = TRUE), "/ 145\n")
cat("    LUSC FULL_MODEL:", sum(cox_recalc_results$LUSC$model_level == "FULL_MODEL", na.rm = TRUE), "/ 145\n")
cat("    Overall PASS rate:", round(overall_pass_rate * 100, 2), "%\n\n")

cat("15. Meta-Analysis:\n")
cat("    Programs in both cohorts:", length(both_full), "\n")
cat("    Hand vs rma.uni logHR max diff:", round(max(rma_diff, na.rm = TRUE), 12), "\n")
cat("    Cochran Q range:", round(range(meta_recalc$Cochran_Q, na.rm = TRUE), 3), "\n")
cat("    I2 range:", round(range(meta_recalc$I2, na.rm = TRUE), 2), "%\n\n")

cat("16. FDR Recalculation:\n")
cat("    Meta P < 0.05:", sum(meta_recalc$meta_P_hand < 0.05, na.rm = TRUE), "\n")
cat("    Meta FDR < 0.05:", sum(meta_recalc$meta_FDR_recalculated < 0.05, na.rm = TRUE), "\n")
cat("    Meta FDR < 0.10:", sum(meta_recalc$meta_FDR_recalculated < 0.10, na.rm = TRUE), "\n\n")

cat("17. Clinical Evidence Levels:\n")
cat("    Level_A:", sum(evidence_levels$Level_A), "\n")
cat("    Level_B:", sum(evidence_levels$Level_B), "\n")
cat("    Level_C:", sum(evidence_levels$Level_C), "\n")
cat("    No_clinical_support:", sum(evidence_levels$No_clinical_support), "\n\n")

cat("18. Redundancy Analysis:\n")
for (cohort in c("LUAD", "LUSC")) {
  cat("    ", cohort, "- Clusters (h=0.20):", redundancy_results[[cohort]]$n_clusters, "\n")
  cat("    ", cohort, "- |rho|>=0.80:", redundancy_results[[cohort]]$n_ge_080, "pairs\n")
}
cat("\n")

cat("19. Jaccard Analysis:\n")
cat("    Total pairs with Jaccard >= 0.25:", nrow(jaccard_pairs), "\n")
cat("    Jaccard >= 0.50:", count_at_threshold[">=0.5"], "\n")
cat("    Jaccard >= 0.75:", count_at_threshold[">=0.75"], "\n")
cat("    Jaccard = 1.00:", count_at_threshold[">=1"], "\n\n")

cat("20. PCA:\n")
for (cohort in c("LUAD", "LUSC")) {
  cat("    ", cohort, "- PC1 variance:", round(pca_results[[cohort]]$pc_var[1] * 100, 2), "%\n")
}
cat("\n")

cat("21. Negative Controls:\n")
cat("    Programs tested:", nrow(perm_results), "\n")
cat("    Mean empirical P:", round(mean(perm_results$empirical_P, na.rm = TRUE), 4), "\n")
cat("    Programs with empirical P < 0.05:", sum(perm_results$empirical_P < 0.05, na.rm = TRUE), "\n\n")

cat("22. B2 Candidates:\n")
if (exists("candidates")) {
  cat("    Total candidates:", nrow(candidates), "\n")
  if (nrow(candidates) > 50) cat("    Top 50 list created: YES\n")
} else {
  cat("    NOT GENERATED (audit status:", audit_status, ")\n")
}
cat("\n")

cat("23. Output Files:\n")
cat("    03_results/step08_TCGA/B1_QC/file_freeze_manifest.csv\n")
cat("    03_results/step08_TCGA/B1_QC/section_III_program_id_integrity.csv\n")
cat("    03_results/step08_TCGA/B1_QC/section_IV_concordance_denominator.csv\n")
cat("    03_results/step08_TCGA/B1_QC/cox/*.csv(.gz)\n")
cat("    03_results/step08_TCGA/B1_QC/meta/*.csv\n")
cat("    03_results/step08_TCGA/B1_QC/redundancy/*.csv\n")
cat("    03_results/step08_TCGA/B1_QC/negative_control/*.csv\n")
cat("    03_results/step08_TCGA/B1_QC/clinical_evidence_levels.csv\n")
cat("    03_results/step08_TCGA/B1_QC/B2_candidate_list.csv\n")
cat("    03_results/step08_TCGA/B1_QC/QC_COMPLETE.txt\n")
cat("    04_figures/step08_TCGA/B1_QC/*.pdf\n\n")

cat("================================================================\n")
cat("08B1 QC Audit COMPLETE\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")
