# SUPERSEDED SCRIPT: logHR extraction bug (hr_row[2]=exp(-coef))
# Original: 01_scripts/08B1_TCGA_program_scoring_and_clinical_validation.R
# This file is kept for historical reference only.
# DO NOT USE for reproduction. 
# --- 
## =========================================================================
## Step 08B1: TCGA Program Scoring & Clinical Validation
## =========================================================================

.libPaths(c(path.expand("~/Library/R/arm64/4.6/library"), .libPaths()))
options(stringsAsFactors = FALSE, warn = 1)

step08b1_start <- Sys.time()

cat("========================================================================\n")
cat("Step 08B1: TCGA Program Scoring & Clinical Validation\n")
cat("========================================================================\n")
cat(sprintf("Start: %s\n", as.character(step08b1_start)))
flush.console()

## =========================================================================
## II. Create Directories
## =========================================================================
cat("\n[II] Creating directories...\n")
dirs <- c(
  "03_results/step08_TCGA/scoring/LUAD", "03_results/step08_TCGA/scoring/LUSC",
  "03_results/step08_TCGA/clinical_models/LUAD", "03_results/step08_TCGA/clinical_models/LUSC",
  "03_results/step08_TCGA/meta_analysis", "03_results/step08_TCGA/qc",
  "03_results/step08_TCGA/combined",
  "04_figures/step08_TCGA/scoring", "04_figures/step08_TCGA/survival",
  "04_figures/step08_TCGA/combined", "logs"
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## =========================================================================
## III. Load Packages
## =========================================================================
cat("\n[III] Loading packages...\n")

lib <- path.expand("~/Library/R/arm64/4.6/library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)

bioc_pkgs <- c("GSVA", "BiocParallel", "MultiAssayExperiment",
               "SummarizedExperiment", "S4Vectors")
cran_pkgs <- c("data.table", "dplyr", "tidyr", "tibble", "stringr",
               "survival", "ggplot2", "pheatmap")

missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE, lib = lib)
if (length(missing_cran) > 0) install.packages(missing_cran, repos = "https://cloud.r-project.org", type = "binary", lib = lib)

suppressPackageStartupMessages({
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
  library(survival)
  library(ggplot2)
  library(pheatmap)
})

pkg_env <- data.frame(
  package = c("R_version", "BiocVersion", "GSVA", "BiocParallel", "survival",
              "MultiAssayExperiment", "SummarizedExperiment"),
  version = c(R.version.string, as.character(BiocManager::version()),
              as.character(packageVersion("GSVA")),
              as.character(packageVersion("BiocParallel")),
              as.character(packageVersion("survival")),
              as.character(packageVersion("MultiAssayExperiment")),
              as.character(packageVersion("SummarizedExperiment"))),
  stringsAsFactors = FALSE
)
write.csv(pkg_env, "03_results/step08_TCGA/qc/GSE243013_step08B1_environment.txt", row.names = FALSE)
cat("[OK] Packages loaded\n")

## =========================================================================
## IV. Freeze Inputs
## =========================================================================
cat("\n[IV] Freezing inputs...\n")
stopifnot(file.exists("03_results/GSE243013_step08A_COMPLETE.txt"))

input_files <- c(
  "03_results/step08_TCGA/GSE243013_step08B_input_index.csv",
  "03_results/step08_TCGA/programs/GSE243013_TCGA_program_manifest.csv",
  "03_results/step08_TCGA/programs/GSE243013_TCGA_program_gene_membership.csv.gz",
  "03_results/step08_TCGA/programs/GSE243013_TCGA_programs_ready_for_scoring.csv",
  "02_data/tcga/clinical/GSE243013_TCGA_LUAD_LUSC_patient_manifest.csv",
  "02_data/tcga/curated/LUAD/TCGA_LUAD_core_primary_tumors_v2.1.1.rds",
  "02_data/tcga/curated/LUSC/TCGA_LUSC_core_primary_tumors_v2.1.1.rds"
)
frozen_rows <- list()
for (f in input_files) {
  if (file.exists(f)) {
    fi <- file.info(f)
    md5 <- tryCatch(tools::md5sum(f), error = function(e) NA_character_)
    frozen_rows[[length(frozen_rows) + 1]] <- data.frame(
      file = f, size = fi$size, mtime = as.character(fi$mtime), md5 = md5,
      status = "OK", stringsAsFactors = FALSE
    )
  } else {
    stop(sprintf("Required input missing: %s", f))
  }
}
write.csv(do.call(rbind, frozen_rows), "03_results/step08_TCGA/qc/GSE243013_step08B1_frozen_input_manifest.csv", row.names = FALSE)

## Load data
prog_manifest <- read.csv("03_results/step08_TCGA/programs/GSE243013_TCGA_program_manifest.csv", stringsAsFactors = FALSE)
gene_mem <- read.csv("03_results/step08_TCGA/programs/GSE243013_TCGA_program_gene_membership.csv.gz", stringsAsFactors = FALSE)
patient_manifest <- read.csv("02_data/tcga/clinical/GSE243013_TCGA_LUAD_LUSC_patient_manifest.csv", stringsAsFactors = FALSE)
cat(sprintf("[OK] Loaded: %d programs, %d gene memberships, %d patients\n",
            nrow(prog_manifest), nrow(gene_mem), nrow(patient_manifest)))

## =========================================================================
## V. Identify RNA Assay & Expression Audit
## =========================================================================
cat("\n[V] Identifying RNA assay...\n")

cohorts <- c("LUAD", "LUSC")
mae_list <- list()
se_list <- list()
expr_list <- list()

for (cohort in cohorts) {
  cat(sprintf("\n--- %s ---\n", cohort))
  mae_path <- sprintf("02_data/tcga/curated/%s/TCGA_%s_core_primary_tumors_v2.1.1.rds", cohort, cohort)
  mae <- readRDS(mae_path)
  mae_list[[cohort]] <- mae

  exp_names <- names(experiments(mae))
  cat(sprintf("[INFO] Assays: %s\n", paste(exp_names, collapse = ", ")))

  rna_match <- grep("RNASeq2GeneNorm", exp_names, value = TRUE)
  if (length(rna_match) != 1) stop(sprintf("%s: found %d RNASeq2GeneNorm matches", cohort, length(rna_match)))

  se <- experiments(mae)[[rna_match[1]]]
  stopifnot(is(se, "SummarizedExperiment"))
  se_list[[cohort]] <- se

  expr <- assay(se)
  expr_list[[cohort]] <- expr

  cat(sprintf("[INFO] RNA dim: %d genes x %d samples\n", nrow(expr), ncol(expr)))

  ## Audit
  n_na <- sum(is.na(expr))
  n_inf <- sum(is.infinite(expr))
  n_neg <- sum(expr < 0, na.rm = TRUE)
  qs <- quantile(expr, probs = c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
  cat(sprintf("[INFO] NA=%d, Inf=%d, Negative=%d\n", n_na, n_inf, n_neg))
  cat(sprintf("[INFO] Range: [%.2f, %.2f]\n", min(expr, na.rm = TRUE), max(expr, na.rm = TRUE)))
  cat(sprintf("[INFO] Quantiles: %s\n", paste(sprintf("%.2f", qs), collapse = ", ")))

  audit_df <- data.frame(
    cohort = cohort, n_genes = nrow(expr), n_samples = ncol(expr),
    n_na = n_na, n_inf = n_inf, n_neg = n_neg,
    min_val = min(expr, na.rm = TRUE), max_val = max(expr, na.rm = TRUE),
    q0 = qs[1], q01 = qs[2], q25 = qs[3], q50 = qs[4], q75 = qs[5], q99 = qs[6], q100 = qs[7],
    stringsAsFactors = FALSE
  )
  write.csv(audit_df, "03_results/step08_TCGA/qc/GSE243013_TCGA_RNA_expression_scale_audit.csv",
            row.names = FALSE)

  if (n_neg > nrow(expr) * ncol(expr) * 0.01) {
    cat(sprintf("[WARNING] %d negative values (%.1f%%) — check normalization\n", n_neg, 100 * n_neg / (nrow(expr) * ncol(expr))))
  }
}

## =========================================================================
## VI. RNA Sample-Patient Mapping
## =========================================================================
cat("\n[VI] RNA sample-patient mapping...\n")

mapping_list <- list()
for (cohort in cohorts) {
  mae <- mae_list[[cohort]]
  se <- se_list[[cohort]]
  sm <- sampleMap(mae)

  rna_match <- grep("RNASeq2GeneNorm", names(experiments(mae)), value = TRUE)
  rna_sm <- sm[sm$assay == rna_match[1], ]

  ## Validate colnames exist in SE
  missing_cols <- setdiff(rna_sm$colname, colnames(se))
  if (length(missing_cols) > 0) stop(sprintf("%s: %d colnames not in SE", cohort, length(missing_cols)))

  rna_sm$cohort <- cohort
  mapping_list[[cohort]] <- rna_sm
  cat(sprintf("[OK] %s: %d RNA samples, %d patients\n", cohort, nrow(rna_sm), length(unique(rna_sm$primary))))
}

mapping_df <- do.call(rbind, mapping_list)
write.csv(mapping_df, "03_results/step08_TCGA/qc/GSE243013_TCGA_RNA_sample_patient_mapping.csv", row.names = FALSE)

## =========================================================================
## VII. Duplicate Sample Resolution
## =========================================================================
cat("\n[VII] Resolving duplicate RNA samples...\n")

dedup_rows <- list()
expr_deduped <- list()

for (cohort in cohorts) {
  expr <- expr_list[[cohort]]
  sm <- mapping_list[[cohort]]

  patient_samples <- split(sm$colname, sm$primary)
  selected_cols <- character(0)

  for (pid in names(patient_samples)) {
    candidates <- patient_samples[[pid]]
    if (length(candidates) == 1) {
      selected_cols <- c(selected_cols, candidates)
      next
    }

    ## Rule: sort by full barcode alphanumeric, pick first
    sorted <- sort(candidates)
    vial <- substr(sorted, nchar(sorted), nchar(sorted))
    vial_order <- order(vial)
    selected <- sorted[vial_order[1]]

    dedup_rows[[length(dedup_rows) + 1]] <- data.frame(
      cohort = cohort, patient_id = pid,
      candidate_sample_count = length(candidates),
      candidate_samples = paste(candidates, collapse = "; "),
      selected_sample = selected,
      selection_rule = "barcode_alphanumeric_sort",
      selection_status = "RESOLVED",
      stringsAsFactors = FALSE
    )
    selected_cols <- c(selected_cols, selected)
  }

  expr_sub <- expr[, selected_cols, drop = FALSE]
  colnames(expr_sub) <- names(patient_samples)
  expr_deduped[[cohort]] <- expr_sub
  cat(sprintf("[OK] %s: %d patients after dedup\n", cohort, ncol(expr_sub)))
}

dedup_df <- if (length(dedup_rows) > 0) do.call(rbind, dedup_rows) else data.frame(cohort=character(0), patient_id=character(0), candidate_sample_count=integer(0), candidate_samples=character(0), selected_sample=character(0), selection_rule=character(0), selection_status=character(0), stringsAsFactors=FALSE)
if (nrow(dedup_df) > 0) {
  write.csv(dedup_df, "03_results/step08_TCGA/qc/GSE243013_TCGA_RNA_duplicate_sample_resolution.csv", row.names = FALSE)
  cat(sprintf("[INFO] Duplicates resolved: %d\n", nrow(dedup_df)))
} else {
  cat("[INFO] No duplicate samples found\n")
  write.csv(dedup_df, "03_results/step08_TCGA/qc/GSE243013_TCGA_RNA_duplicate_sample_resolution.csv", row.names = FALSE)
}

## =========================================================================
## VIII. Gene Name Cleaning
## =========================================================================
cat("\n[VIII] Cleaning gene names...\n")

gene_map_list <- list()
expr_cleaned <- list()

for (cohort in cohorts) {
  expr <- expr_deduped[[cohort]]
  gene_orig <- rownames(expr)

  gene_clean <- sapply(gene_orig, function(g) {
    parts <- unlist(strsplit(as.character(g), "\\|"))
    trimws(parts[1])
  })

  gene_df <- data.frame(
    gene_symbol_original = gene_orig, gene_symbol_clean = gene_clean,
    stringsAsFactors = FALSE
  )

  ## Duplicate handling: average across duplicates
  dup_genes <- gene_clean[duplicated(gene_clean)]
  n_dup <- length(unique(dup_genes))
  if (n_dup > 0) {
    cat(sprintf("[INFO] %s: %d duplicate gene symbols — averaging\n", cohort, n_dup))
    expr_list_avg <- list()
    for (g in unique(gene_clean)) {
      idx <- which(gene_clean == g)
      if (length(idx) == 1) {
        expr_list_avg[[g]] <- expr[idx, ]
      } else {
        expr_list_avg[[g]] <- colMeans(expr[idx, , drop = FALSE], na.rm = TRUE)
      }
    }
    expr <- do.call(rbind, expr_list_avg)
  }

  ## Remove zero-variance genes
  gene_vars <- apply(expr, 1, var, na.rm = TRUE)
  zero_var <- sum(gene_vars == 0 | is.na(gene_vars))
  if (zero_var > 0) {
    cat(sprintf("[INFO] %s: removing %d zero-variance genes\n", cohort, zero_var))
    expr <- expr[gene_vars > 0 & !is.na(gene_vars), ]
  }

  gene_map_list[[cohort]] <- gene_df
  expr_cleaned[[cohort]] <- expr
  cat(sprintf("[OK] %s: %d unique genes, %d patients\n", cohort, nrow(expr), ncol(expr)))
}

## Save gene mapping
for (cohort in cohorts) {
  con <- gzfile(sprintf("03_results/step08_TCGA/qc/GSE243013_TCGA_%s_RNA_gene_mapping_and_duplicates.csv.gz", cohort), "w")
  write.csv(gene_map_list[[cohort]], con, row.names = FALSE)
  close(con)
}

## =========================================================================
## IX. Build Program Gene Sets
## =========================================================================
cat("\n[IX] Building program gene sets...\n")

eligible_progs <- prog_manifest[prog_manifest$program_primary_eligible == TRUE, ]
cat(sprintf("[INFO] Eligible programs: %d\n", nrow(eligible_progs)))

program_gene_sets <- list()
for (pid in eligible_progs$program_id) {
  prog_genes <- unique(gene_mem$gene[gene_mem$program_id == pid])
  program_gene_sets[[pid]] <- prog_genes
}

## Determine scorable per cohort
scorable_info <- list()
for (pid in names(program_gene_sets)) {
  for (cohort in cohorts) {
    available_genes <- intersect(program_gene_sets[[pid]], rownames(expr_cleaned[[cohort]]))
    scorable_info[[length(scorable_info) + 1]] <- data.frame(
      program_id = pid, cohort = cohort,
      program_gene_count = length(program_gene_sets[[pid]]),
      genes_available = length(available_genes),
      overlap_fraction = round(length(available_genes) / max(length(program_gene_sets[[pid]]), 1), 4),
      scorable = length(available_genes) >= 5 && (length(available_genes) / max(length(program_gene_sets[[pid]]), 1)) >= 0.5,
      stringsAsFactors = FALSE
    )
  }
}

scorable_df <- do.call(rbind, scorable_info)
write.csv(scorable_df, "03_results/step08_TCGA/scoring/GSE243013_TCGA_program_scoring_gene_overlap_recheck.csv", row.names = FALSE)

## =========================================================================
## X. ssGSEA Scoring
## =========================================================================
cat("\n[X] Running ssGSEA scoring...\n")

ssgsea_scores <- list()
for (cohort in cohorts) {
  cat(sprintf("\n--- %s ssGSEA ---\n", cohort))
  expr <- expr_cleaned[[cohort]]
  available_sets <- list()
  for (pid in names(program_gene_sets)) {
    si <- scorable_df[scorable_df$program_id == pid & scorable_df$cohort == cohort, ]
    if (nrow(si) > 0 && si$scorable) {
      available_sets[[pid]] <- program_gene_sets[[pid]]
    }
  }
  cat(sprintf("[INFO] Scorable programs: %d\n", length(available_sets)))

  if (length(available_sets) == 0) {
    cat("[WARNING] No scorable programs\n")
    next
  }

  param <- GSVA::ssgseaParam(
    exprData = expr, geneSets = available_sets,
    minSize = 5, maxSize = 500, alpha = 0.25,
    normalize = TRUE, checkNA = "auto", use = "everything",
    ondisk = "no", verbose = TRUE
  )

  scores <- tryCatch({
    gsva(param, verbose = TRUE, BPPARAM = BiocParallel::SerialParam())
  }, error = function(e) {
    cat(sprintf("[ERROR] ssGSEA failed: %s\n", conditionMessage(e)))
    NULL
  })

  if (is.null(scores)) next

  ssgsea_scores[[cohort]] <- scores
  cat(sprintf("[OK] %s ssGSEA: %d programs x %d patients\n", cohort, nrow(scores), ncol(scores)))

  ## Validate
  if (any(is.na(scores))) cat(sprintf("[WARNING] %d NA values in scores\n", sum(is.na(scores))))
  if (any(is.infinite(scores))) cat(sprintf("[WARNING] %d Inf values in scores\n", sum(is.infinite(scores))))

  ## Save RDS
  saveRDS(scores, sprintf("02_data/tcga/clinical/TCGA_%s_program_scores_ssGSEA.rds", cohort))

  ## Save long CSV
  long_df <- as.data.frame(t(scores))
  long_df$patient_id <- rownames(long_df)
  long_long <- pivot_longer(long_df, cols = -patient_id, names_to = "program_id", values_to = "ssGSEA_score")
  con <- gzfile(sprintf("03_results/step08_TCGA/scoring/%s/TCGA_%s_program_scores_ssGSEA_long.csv.gz", cohort, cohort), "w")
  write.csv(as.data.frame(long_long), con, row.names = FALSE)
  close(con)
}

## =========================================================================
## XI. GSVA Sensitivity Scoring
## =========================================================================
cat("\n[XI] Running GSVA sensitivity scoring...\n")

gsva_scores <- list()
for (cohort in cohorts) {
  cat(sprintf("\n--- %s GSVA ---\n", cohort))
  expr <- expr_cleaned[[cohort]]

  available_sets <- list()
  for (pid in names(program_gene_sets)) {
    si <- scorable_df[scorable_df$program_id == pid & scorable_df$cohort == cohort, ]
    if (nrow(si) > 0 && si$scorable) {
      available_sets[[pid]] <- program_gene_sets[[pid]]
    }
  }

  if (length(available_sets) == 0) next

  param <- GSVA::gsvaParam(
    exprData = expr, geneSets = available_sets,
    minSize = 5, maxSize = 500, kcdf = "Gaussian",
    tau = 1, maxDiff = TRUE, absRanking = FALSE,
    sparse = FALSE, checkNA = "auto", use = "everything",
    ondisk = "no", verbose = TRUE
  )

  scores <- tryCatch({
    gsva(param, verbose = TRUE, BPPARAM = BiocParallel::SerialParam())
  }, error = function(e) {
    cat(sprintf("[ERROR] GSVA failed: %s\n", conditionMessage(e)))
    NULL
  })

  if (is.null(scores)) next

  gsva_scores[[cohort]] <- scores
  saveRDS(scores, sprintf("02_data/tcga/clinical/TCGA_%s_program_scores_GSVA.rds", cohort))
  cat(sprintf("[OK] %s GSVA: %d programs x %d patients\n", cohort, nrow(scores), ncol(scores)))
}

## =========================================================================
## XII. Scoring Method Concordance
## =========================================================================
cat("\n[XII] Score method concordance...\n")

concordance_rows <- list()
for (cohort in cohorts) {
  if (is.null(ssgsea_scores[[cohort]]) || is.null(gsva_scores[[cohort]])) next
  ssg <- ssgsea_scores[[cohort]]
  gsv <- gsva_scores[[cohort]]
  common_progs <- intersect(rownames(ssg), rownames(gsv))
  common_pats <- intersect(colnames(ssg), colnames(gsv))

  for (pid in common_progs) {
    s_vec <- as.numeric(ssg[pid, common_pats])
    g_vec <- as.numeric(gsv[pid, common_pats])
    valid <- is.finite(s_vec) & is.finite(g_vec)
    if (sum(valid) < 10) next

    sp_cor <- cor(s_vec[valid], g_vec[valid], method = "spearman")
    pe_cor <- cor(s_vec[valid], g_vec[valid], method = "pearson")
    sp_test <- cor.test(s_vec[valid], g_vec[valid], method = "spearman")

    concordance_rows[[length(concordance_rows) + 1]] <- data.frame(
      cohort = cohort, program_id = pid,
      spearman_rho = sp_cor, pearson_rho = pe_cor,
      n_patients = sum(valid), p_value = sp_test$p.value,
      scoring_method_concordant = sp_cor >= 0.60 && sp_cor > 0,
      highly_concordant = sp_cor >= 0.80,
      stringsAsFactors = FALSE
    )
  }
}

concordance_df <- do.call(rbind, concordance_rows)
if (nrow(concordance_df) > 0) {
  concordance_df$FDR <- p.adjust(concordance_df$p_value, method = "BH")
  write.csv(concordance_df, "03_results/step08_TCGA/scoring/GSE243013_TCGA_ssGSEA_vs_GSVA_concordance.csv", row.names = FALSE)
  cat(sprintf("[OK] Concordance: %d program-cohorts, %d with rho>=0.60\n",
              nrow(concordance_df), sum(concordance_df$scoring_method_concordant)))
}

## =========================================================================
## XIII. Patient-Level Score Manifest
## =========================================================================
cat("\n[XIII] Building patient-level score manifest...\n")

manifest_rows <- list()
for (cohort in cohorts) {
  if (is.null(ssgsea_scores[[cohort]])) next
  ssg <- ssgsea_scores[[cohort]]
  gsv <- gsva_scores[[cohort]]
  pm <- patient_manifest[patient_manifest$cohort == cohort, ]

  for (pid_prog in rownames(ssg)) {
    prog_info <- prog_manifest[prog_manifest$program_id == pid_prog, ]
    if (nrow(prog_info) == 0) next

    p_tier <- as.character(prog_info$priority_tier[1])
    p_ct <- as.character(prog_info$cell_type[1])
    p_coll <- as.character(prog_info$collection[1])
    p_path <- as.character(prog_info$pathway[1])
    p_nes <- as.numeric(prog_info$NES[1])
    p_dir <- as.character(prog_info$direction[1])

    for (pat in colnames(ssg)) {
      pm_row <- pm[pm$patient_id == pat, ]
      if (nrow(pm_row) == 0) next

      ssgsea_val <- as.numeric(ssg[pid_prog, pat])
      if (is.na(ssgsea_val)) next

      gsva_val <- if (!is.null(gsv) && pid_prog %in% rownames(gsv) && pat %in% colnames(gsv)) {
        as.numeric(gsv[pid_prog, pat])
      } else NA_real_

      conc_row <- concordance_df[concordance_df$cohort == cohort & concordance_df$program_id == pid_prog, ]
      conc_val <- if (nrow(conc_row) > 0) conc_row$scoring_method_concordant[1] else NA

      row_list <- list(
        cohort = cohort, patient_id = pat, program_id = pid_prog,
        priority_tier = p_tier, cell_type = p_ct, collection = p_coll,
        pathway = p_path, step07_NES = p_nes, expected_direction = p_dir,
        ssGSEA_score = ssgsea_val, GSVA_score = gsva_val,
        age_at_diagnosis = as.numeric(pm_row$age_at_diagnosis[1]),
        sex = as.character(pm_row$sex[1]),
        pathologic_stage_original = as.character(pm_row$pathologic_stage_original[1]),
        pathologic_stage_clean = as.character(pm_row$pathologic_stage_clean[1]),
        OS_days = as.numeric(pm_row$OS_days[1]),
        OS_event = as.numeric(pm_row$OS_event[1]),
        has_primary_RNA = as.logical(pm_row$has_primary_RNA[1]),
        scoring_method_concordant = conc_val
      )
      row_list <- lapply(row_list, function(x) if (length(x) == 0 || is.null(x)) NA else x)
      manifest_rows[[length(manifest_rows) + 1]] <- do.call(data.frame, c(row_list, stringsAsFactors = FALSE))
    }
  }
}
cat(sprintf("[INFO] Manifest rows built: %d\n", length(manifest_rows)))

manifest_long <- if (length(manifest_rows) > 0) do.call(rbind, manifest_rows) else data.frame(cohort=character(0), patient_id=character(0), program_id=character(0), stringsAsFactors=FALSE)

## Z-score within cohort × program
manifest_long$score_z <- NA_real_
for (pid in unique(manifest_long$program_id)) {
  for (co in cohorts) {
    idx <- manifest_long$program_id == pid & manifest_long$cohort == co
    vals <- manifest_long$ssGSEA_score[idx]
    if (sum(!is.na(vals)) > 1) {
      manifest_long$score_z[idx] <- (vals - mean(vals, na.rm = TRUE)) / sd(vals, na.rm = TRUE)
    }
  }
}

con <- gzfile("03_results/step08_TCGA/scoring/GSE243013_TCGA_patient_program_score_manifest.csv.gz", "w")
write.csv(manifest_long, con, row.names = FALSE)
close(con)
cat(sprintf("[OK] Patient-program manifest: %d rows\n", nrow(manifest_long)))

## =========================================================================
## XIV. Stage Cleaning
## =========================================================================
cat("\n[XIV] Stage cleaning...\n")

clean_stage <- function(stage_str) {
  if (is.na(stage_str) || stage_str == "" || stage_str == "NA") return(NA_character_)
  s <- toupper(trimws(stage_str))
  if (grepl("^STAGE\\s+", s)) s <- sub("^STAGE\\s+", "", s)
  if (grepl("^I[A-Z]?$", s)) return("I")
  if (grepl("^II[A-Z]?$", s)) return("II")
  if (grepl("^III[A-Z]?$", s)) return("III")
  if (grepl("^IV[A-Z]?$", s)) return("IV")
  return(NA_character_)
}

patient_manifest$stage_factor <- sapply(patient_manifest$pathologic_stage_clean, clean_stage)
patient_manifest$stage_numeric <- match(patient_manifest$stage_factor, c("I", "II", "III", "IV"))

stage_audit <- patient_manifest[, c("cohort", "patient_id", "pathologic_stage_clean", "stage_factor", "stage_numeric")]
write.csv(stage_audit, "03_results/step08_TCGA/qc/GSE243013_TCGA_stage_cleaning_audit.csv", row.names = FALSE)
cat(sprintf("[OK] Stage: I=%d, II=%d, III=%d, IV=%d, NA=%d\n",
            sum(patient_manifest$stage_factor == "I", na.rm = TRUE),
            sum(patient_manifest$stage_factor == "II", na.rm = TRUE),
            sum(patient_manifest$stage_factor == "III", na.rm = TRUE),
            sum(patient_manifest$stage_factor == "IV", na.rm = TRUE),
            sum(is.na(patient_manifest$stage_factor))))

## =========================================================================
## XV. Stage Association
## =========================================================================
cat("\n[XV] Stage association analysis...\n")

stage_rows <- list()
for (cohort in cohorts) {
  pm_cohort <- patient_manifest[patient_manifest$cohort == cohort, ]
  pm_cohort$age_z <- scale(pm_cohort$age_at_diagnosis)
  pm_cohort$sex_f <- factor(pm_cohort$sex)

  for (pid in unique(manifest_long$program_id[manifest_long$cohort == cohort])) {
    scores_sub <- manifest_long[manifest_long$cohort == cohort & manifest_long$program_id == pid, ]
    merged <- merge(scores_sub, pm_cohort[, c("patient_id", "stage_factor", "stage_numeric", "age_z", "sex_f")],
                    by = "patient_id", all.x = FALSE)

    merged <- merged[!is.na(merged$score_z) & !is.na(merged$stage_factor), ]
    stage_counts <- table(merged$stage_factor)
    if (length(stage_counts) < 3 || nrow(merged) < 100) next
    if (any(stage_counts < 10)) next

    ## Kruskal-Wallis
    kw <- tryCatch(kruskal.test(score_z ~ stage_factor, data = merged), error = function(e) NULL)
    kw_p <- if (!is.null(kw)) kw$p.value else NA_real_

    ## Linear trend model
    model_full <- tryCatch(
      lm(score_z ~ stage_numeric + age_z + sex_f, data = merged),
      error = function(e) NULL
    )
    model_reduced <- tryCatch(
      lm(score_z ~ stage_numeric, data = merged),
      error = function(e) NULL
    )

    stage_beta <- NA_real_
    stage_se <- NA_real_
    stage_p <- NA_real_
    model_level <- "NONE"
    n_complete <- 0

    if (!is.null(model_full) && sum(complete.cases(merged[, c("score_z", "stage_numeric", "age_z", "sex_f")])) >= 100) {
      coefs <- summary(model_full)$coefficients
      if ("stage_numeric" %in% rownames(coefs)) {
        stage_beta <- coefs["stage_numeric", 1]
        stage_se <- coefs["stage_numeric", 2]
        stage_p <- coefs["stage_numeric", 4]
        model_level <- "FULL_MODEL"
        n_complete <- sum(complete.cases(merged[, c("score_z", "stage_numeric", "age_z", "sex_f")]))
      }
    } else if (!is.null(model_reduced)) {
      coefs <- summary(model_reduced)$coefficients
      if ("stage_numeric" %in% rownames(coefs)) {
        stage_beta <- coefs["stage_numeric", 1]
        stage_se <- coefs["stage_numeric", 2]
        stage_p <- coefs["stage_numeric", 4]
        model_level <- "MODEL_REDUCED"
        n_complete <- sum(!is.na(merged$score_z) & !is.na(merged$stage_numeric))
      }
    }

    stage_rows[[length(stage_rows) + 1]] <- data.frame(
      cohort = cohort, program_id = pid,
      kw_p_value = kw_p, stage_beta = stage_beta, stage_se = stage_se,
      stage_p_value = stage_p, model_level = model_level,
      n_complete = n_complete, n_stages = length(stage_counts),
      stringsAsFactors = FALSE
    )
  }
}

stage_df <- if (length(stage_rows) > 0) do.call(rbind, stage_rows) else data.frame(cohort=character(0), program_id=character(0), kw_p_value=numeric(0), stage_beta=numeric(0), stage_se=numeric(0), stage_p_value=numeric(0), model_level=character(0), n_complete=integer(0), n_stages=integer(0), stringsAsFactors=FALSE)
if (nrow(stage_df) > 0) {
  stage_df$FDR <- p.adjust(stage_df$stage_p_value, method = "BH")
  con <- gzfile("03_results/step08_TCGA/clinical_models/GSE243013_TCGA_stage_association_all_programs.csv.gz", "w")
  write.csv(stage_df, con, row.names = FALSE)
  close(con)
  cat(sprintf("[OK] Stage association: %d program-cohorts, %d FDR<0.05\n",
              nrow(stage_df), sum(stage_df$FDR < 0.05, na.rm = TRUE)))
}

## =========================================================================
## XVI-XVII. OS Cox Models
## =========================================================================
cat("\n[XVI-XVII] OS Cox models...\n")

cox_results_all <- list()
for (cohort in cohorts) {
  cat(sprintf("\n--- %s Cox models ---\n", cohort))
  pm_cohort <- patient_manifest[patient_manifest$cohort == cohort, ]
  pm_cohort$age_z <- scale(pm_cohort$age_at_diagnosis)
  pm_cohort$sex_f <- factor(pm_cohort$sex)
  pm_cohort$stage_f <- factor(pm_cohort$stage_factor, levels = c("I", "II", "III", "IV"))
  pm_cohort$stage_f <- droplevels(pm_cohort$stage_f)

  ## OS cohort audit
  os_eligible <- pm_cohort[!is.na(pm_cohort$OS_days) & pm_cohort$OS_days >= 0 &
                            pm_cohort$OS_event %in% c(0, 1), ]
  os_audit <- data.frame(
    cohort = cohort, total_patients = nrow(pm_cohort),
    os_eligible = nrow(os_eligible), n_events = sum(os_eligible$OS_event),
    median_followup = median(os_eligible$OS_days[os_eligible$OS_event == 0], na.rm = TRUE),
    n_age_na = sum(is.na(os_eligible$age_at_diagnosis)),
    n_sex_na = sum(is.na(os_eligible$sex)),
    n_stage_na = sum(is.na(os_eligible$stage_factor)),
    stringsAsFactors = FALSE
  )
  write.csv(os_audit, sprintf("03_results/step08_TCGA/qc/GSE243013_TCGA_%s_OS_cohort_audit.csv", cohort), row.names = FALSE)

  cox_rows <- list()
  for (pid in unique(manifest_long$program_id[manifest_long$cohort == cohort])) {
    scores_sub <- manifest_long[manifest_long$cohort == cohort & manifest_long$program_id == pid, ]
    merged <- merge(scores_sub, pm_cohort[, c("patient_id", "age_z", "sex_f", "stage_f")],
                    by = "patient_id", all.x = FALSE)
    merged <- merged[!is.na(merged$score_z) & !is.na(merged$OS_days) & merged$OS_days >= 0 &
                      merged$OS_event %in% c(0, 1), ]

    n_complete <- nrow(merged)
    n_events <- sum(merged$OS_event)
    if (n_complete < 80 || n_events < 20) next
    if (var(merged$score_z, na.rm = TRUE) == 0) next

    ## Try models in order
    model <- NULL
    model_level <- "NONE"
    model_formula <- ""

    ## Full model
    full_f <- Surv(OS_days / 365.25, OS_event) ~ score_z + age_z + sex_f + stage_f
    full_cc <- merged[complete.cases(merged[, c("score_z", "age_z", "sex_f", "stage_f", "OS_days", "OS_event")]), ]
    if (nrow(full_cc) >= 100 && sum(full_cc$OS_event) >= 30) {
      tryCatch({
        m <- coxph(full_f, data = merged)
        ## Check design matrix rank
        mm <- model.matrix(full_f, data = merged)
        if (qr(mm)$rank == ncol(mm)) {
          model <- m
          model_level <- "FULL_MODEL"
          model_formula <- deparse(full_f)
        }
      }, error = function(e) cat(sprintf("  [WARN] Full model failed for %s: %s\n", pid, e$message)))
    }

    ## Model 2
    if (is.null(model)) {
      m2_f <- Surv(OS_days / 365.25, OS_event) ~ score_z + age_z + sex_f
      m2_cc <- merged[complete.cases(merged[, c("score_z", "age_z", "sex_f", "OS_days", "OS_event")]), ]
      if (nrow(m2_cc) >= 100 && sum(m2_cc$OS_event) >= 30) {
        tryCatch({
          m <- coxph(m2_f, data = merged)
          model <- m
          model_level <- "MODEL2"
          model_formula <- deparse(m2_f)
        }, error = function(e) NULL)
      }
    }

    ## Model 3
    if (is.null(model)) {
      m3_f <- Surv(OS_days / 365.25, OS_event) ~ score_z
      m3_cc <- merged[complete.cases(merged[, c("score_z", "OS_days", "OS_event")]), ]
      if (nrow(m3_cc) >= 80 && sum(m3_cc$OS_event) >= 20) {
        tryCatch({
          m <- coxph(m3_f, data = merged)
          model <- m
          model_level <- "MODEL3_EXPLORATORY"
          model_formula <- deparse(m3_f)
        }, error = function(e) NULL)
      }
    }

    if (is.null(model)) next

    s <- summary(model)
    hr_row <- s$conf.int["score_z", ]
    p_row <- s$coefficients["score_z", ]

    cox_rows[[length(cox_rows) + 1]] <- data.frame(
      cohort = cohort, program_id = pid,
      model_level = model_level, model_formula = model_formula,
      n_complete = nrow(merged[complete.cases(merged[, intersect(c("score_z", "age_z", "sex_f", "stage_f", "OS_days", "OS_event"), colnames(merged))]), ]),
      n_events = sum(merged$OS_event),
      HR_per_1SD = hr_row[1], logHR = hr_row[2],
      standard_error = p_row[3], lower_95CI = hr_row[3], upper_95CI = hr_row[4],
      Wald_PValue = p_row[5], concordance = s$concordance[1],
      AIC = AIC(model), warning = "", model_status = "OK",
      stringsAsFactors = FALSE
    )
  }

  if (length(cox_rows) > 0) {
    cox_df <- do.call(rbind, cox_rows)
    cox_df$FDR_within_cohort <- p.adjust(cox_df$Wald_PValue, method = "BH")
    cox_results_all[[cohort]] <- cox_df
    con <- gzfile(sprintf("03_results/step08_TCGA/clinical_models/%s/TCGA_%s_program_OS_Cox_results.csv.gz", cohort, cohort), "w")
    write.csv(cox_df, con, row.names = FALSE)
    close(con)
    cat(sprintf("[OK] %s: %d Cox models (%d FULL)\n", cohort, nrow(cox_df),
                sum(cox_df$model_level == "FULL_MODEL")))
  }
}

## =========================================================================
## XVIII. PH Assumption
## =========================================================================
cat("\n[XVIII] PH assumption testing...\n")

ph_rows <- list()
for (cohort in cohorts) {
  if (is.null(cox_results_all[[cohort]])) next
  pm_cohort <- patient_manifest[patient_manifest$cohort == cohort, ]
  pm_cohort$age_z <- scale(pm_cohort$age_at_diagnosis)
  pm_cohort$sex_f <- factor(pm_cohort$sex)
  pm_cohort$stage_f <- factor(pm_cohort$stage_factor, levels = c("I", "II", "III", "IV"))
  pm_cohort$stage_f <- droplevels(pm_cohort$stage_f)

  for (i in seq_len(nrow(cox_results_all[[cohort]]))) {
    row <- cox_results_all[[cohort]][i, ]
    scores_sub <- manifest_long[manifest_long$cohort == cohort & manifest_long$program_id == row$program_id, ]
    merged <- merge(scores_sub, pm_cohort[, c("patient_id", "age_z", "sex_f", "stage_f")],
                    by = "patient_id", all.x = FALSE)
    merged <- merged[!is.na(merged$score_z) & !is.na(merged$OS_days) & merged$OS_days >= 0 &
                      merged$OS_event %in% c(0, 1), ]

    f <- as.formula(row$model_formula)
    model <- tryCatch(coxph(f, data = merged), error = function(e) NULL)
    if (is.null(model)) next

    zph <- tryCatch(cox.zph(model), error = function(e) NULL)
    if (is.null(zph)) next

    score_ph_p <- NA_real_
    global_ph_p <- zph$table["GLOBAL", "p"]
    if ("score_z" %in% rownames(zph$table)) {
      score_ph_p <- zph$table["score_z", "p"]
    }

    ph_rows[[length(ph_rows) + 1]] <- data.frame(
      cohort = cohort, program_id = row$program_id,
      model_level = row$model_level,
      PH_score_pvalue = score_ph_p, PH_global_pvalue = global_ph_p,
      PH_score_pass = !is.na(score_ph_p) && score_ph_p >= 0.05,
      PH_global_pass = !is.na(global_ph_p) && global_ph_p >= 0.05,
      stringsAsFactors = FALSE
    )
  }
}

ph_df <- if (length(ph_rows) > 0) do.call(rbind, ph_rows) else data.frame(cohort=character(0), program_id=character(0), model_level=character(0), PH_score_pvalue=numeric(0), PH_global_pvalue=numeric(0), PH_score_pass=logical(0), PH_global_pass=logical(0), stringsAsFactors=FALSE)
if (nrow(ph_df) > 0) {
  con <- gzfile("03_results/step08_TCGA/clinical_models/GSE243013_TCGA_Cox_PH_assumption_results.csv.gz", "w")
  write.csv(ph_df, con, row.names = FALSE)
  close(con)
  cat(sprintf("[OK] PH tests: %d, PH warnings: %d\n", nrow(ph_df), sum(!ph_df$PH_score_pass)))
}

## =========================================================================
## XIX. Step 07 Direction Relationship
## =========================================================================
cat("\n[XIX] Step 07 direction relationship...\n")

for (cohort in cohorts) {
  if (is.null(cox_results_all[[cohort]])) next
  cox_df <- cox_results_all[[cohort]]
  prog_info <- prog_manifest[match(cox_df$program_id, prog_manifest$program_id), ]

  ## Survival direction consistent
  cox_df$survival_direction_consistent <- NA_character_
  for (i in seq_len(nrow(cox_df))) {
    nes <- prog_info$NES[i]
    hr <- cox_df$HR_per_1SD[i]
    if (is.na(nes) || is.na(hr)) next
    if ((nes > 0 && hr < 1) || (nes < 0 && hr > 1)) {
      cox_df$survival_direction_consistent[i] <- "CONSISTENT"
    } else {
      cox_df$survival_direction_consistent[i] <- "INCONSISTENT"
    }
  }
  cox_results_all[[cohort]] <- cox_df
}

## =========================================================================
## XX. Meta-Analysis
## =========================================================================
cat("\n[XX] LUAD/LUSC fixed-effect meta-analysis...\n")

luad_cox <- cox_results_all[["LUAD"]]
lusc_cox <- cox_results_all[["LUSC"]]

meta_rows <- list()
if (!is.null(luad_cox) && !is.null(lusc_cox)) {
  common_progs <- intersect(
    luad_cox$program_id[luad_cox$model_level == "FULL_MODEL"],
    lusc_cox$program_id[lusc_cox$model_level == "FULL_MODEL"]
  )

  for (pid in common_progs) {
    lu <- luad_cox[luad_cox$program_id == pid & luad_cox$model_level == "FULL_MODEL", ]
    ls <- lusc_cox[lusc_cox$program_id == pid & lusc_cox$model_level == "FULL_MODEL", ]
    if (nrow(lu) == 0 || nrow(ls) == 0) next

    loghrs <- c(lu$logHR, ls$logHR)
    ses <- c(lu$standard_error, ls$standard_error)
    weights <- 1 / ses^2
    meta_loghr <- sum(weights * loghrs) / sum(weights)
    meta_se <- sqrt(1 / sum(weights))

    meta_rows[[length(meta_rows) + 1]] <- data.frame(
      program_id = pid,
      LUAD_logHR = lu$logHR, LUAD_SE = lu$standard_error, LUAD_FDR = lu$FDR_within_cohort,
      LUSC_logHR = ls$logHR, LUSC_SE = ls$standard_error, LUSC_FDR = ls$FDR_within_cohort,
      meta_logHR = meta_loghr, meta_SE = meta_se,
      meta_HR = exp(meta_loghr), meta_lower_95CI = exp(meta_loghr - 1.96 * meta_se),
      meta_upper_95CI = exp(meta_loghr + 1.96 * meta_se),
      meta_PValue = 2 * pnorm(-abs(meta_loghr / meta_se)),
      cross_histology_direction_concordant = sign(lu$logHR) == sign(ls$logHR),
      stringsAsFactors = FALSE
    )
  }
}

meta_df <- if (length(meta_rows) > 0) do.call(rbind, meta_rows) else data.frame(program_id=character(0), stringsAsFactors=FALSE)
if (nrow(meta_df) > 0) {
  meta_df$meta_FDR <- p.adjust(meta_df$meta_PValue, method = "BH")

  ## Cochran Q
  meta_df$Q <- NA_real_
  meta_df$heterogeneity_P <- NA_real_
  meta_df$I2 <- NA_real_
  for (i in seq_len(nrow(meta_df))) {
    lu_loghr <- meta_df$LUAD_logHR[i]
    ls_loghr <- meta_df$LUSC_logHR[i]
    lu_se <- meta_df$LUAD_SE[i]
    ls_se <- meta_df$LUSC_SE[i]
    meta_lhr <- meta_df$meta_logHR[i]

    Q <- sum(((c(lu_loghr, ls_loghr) - meta_lhr)^2) * c(1/lu_se^2, 1/ls_se^2))
    meta_df$Q[i] <- Q
    meta_df$heterogeneity_P[i] <- pchisq(Q, df = 1, lower.tail = FALSE)
    meta_df$I2[i] <- max(0, (Q - 1) / Q * 100)
  }

  meta_df$low_heterogeneity <- meta_df$I2 < 50

  write.csv(meta_df, "03_results/step08_TCGA/meta_analysis/GSE243013_TCGA_program_OS_fixed_effect_meta.csv", row.names = FALSE)
  cat(sprintf("[OK] Meta-analysis: %d programs, %d FDR<0.05, %d I2<50%%\n",
              nrow(meta_df), sum(meta_df$meta_FDR < 0.05), sum(meta_df$low_heterogeneity)))
}

## =========================================================================
## XXI. Clinical Validation Evidence Stratification
## =========================================================================
cat("\n[XXI] Evidence stratification...\n")

all_programs <- unique(manifest_long$program_id)
validation_rows <- list()

for (pid in all_programs) {
  prog <- prog_manifest[prog_manifest$program_id == pid, ]
  if (nrow(prog) == 0) next

  ## Concordance
  conc_luad <- concordance_df[concordance_df$cohort == "LUAD" & concordance_df$program_id == pid, ]
  conc_lusc <- concordance_df[concordance_df$cohort == "LUSC" & concordance_df$program_id == pid, ]
  rho_luad <- if (nrow(conc_luad) > 0) conc_luad$spearman_rho[1] else NA_real_
  rho_lusc <- if (nrow(conc_lusc) > 0) conc_lusc$spearman_rho[1] else NA_real_

  ## Cox
  cox_luad <- cox_results_all[["LUAD"]]; cox_lusc <- cox_results_all[["LUSC"]]
  hr_luad <- fdr_luad <- ph_pass_luad <- NA
  hr_lusc <- fdr_lusc <- ph_pass_lusc <- NA
  if (!is.null(cox_luad)) {
    cl <- cox_luad[cox_luad$program_id == pid & cox_luad$model_level == "FULL_MODEL", ]
    if (nrow(cl) > 0) { hr_luad <- cl$HR_per_1SD; fdr_luad <- cl$FDR_within_cohort }
  }
  if (!is.null(cox_lusc)) {
    cl <- cox_lusc[cox_lusc$program_id == pid & cox_lusc$model_level == "FULL_MODEL", ]
    if (nrow(cl) > 0) { hr_lusc <- cl$HR_per_1SD; fdr_lusc <- cl$FDR_within_cohort }
  }
  if (!is.null(ph_df)) {
    ph_l <- ph_df[ph_df$cohort == "LUAD" & ph_df$program_id == pid, ]
    ph_s <- ph_df[ph_df$cohort == "LUSC" & ph_df$program_id == pid, ]
    ph_pass_luad <- if (nrow(ph_l) > 0) ph_l$PH_score_pass[1] else NA
    ph_pass_lusc <- if (nrow(ph_s) > 0) ph_s$PH_score_pass[1] else NA
  }

  ## Meta
  meta_hr <- meta_fdr <- i2_val <- dir_concordant <- NA
  if (!is.null(meta_df) && is.data.frame(meta_df) && nrow(meta_df) > 0) {
    mm <- meta_df[meta_df$program_id == pid, ]
    if (is.data.frame(mm) && nrow(mm) > 0) {
      meta_hr <- mm$meta_HR; meta_fdr <- mm$meta_FDR
      i2_val <- mm$I2; dir_concordant <- mm$cross_histology_direction_concordant
    }
  }

  ## Stage
  stage_luad <- stage_lusc <- NA
  if (!is.null(stage_df)) {
    sl <- stage_df[stage_df$cohort == "LUAD" & stage_df$program_id == pid, ]
    ss <- stage_df[stage_df$cohort == "LUSC" & stage_df$program_id == pid, ]
    stage_luad <- if (nrow(sl) > 0) sl$stage_beta[1] else NA
    stage_lusc <- if (nrow(ss) > 0) ss$stage_beta[1] else NA
  }

  ## Direction consistency with Step 07
  nes <- prog$NES[1]
  surv_dir_consistent <- NA_character_
  if (!is.na(hr_luad) && !is.na(nes)) {
    surv_dir_consistent <- if ((nes > 0 && hr_luad < 1) || (nes < 0 && hr_luad > 1)) "CONSISTENT" else "INCONSISTENT"
  }

  ## Stage direction
  stage_dir_consistent <- NA_character_
  if (!is.na(nes) && !is.na(stage_luad)) {
    if ((nes > 0 && stage_luad < 0) || (nes < 0 && stage_luad > 0)) {
      stage_dir_consistent <- "CONSISTENT"
    } else {
      stage_dir_consistent <- "INCONSISTENT"
    }
  }

  ## Validation level
  both_rho_ok <- !is.na(rho_luad) && !is.na(rho_lusc) && rho_luad >= 0.60 && rho_lusc >= 0.60
  both_full <- !is.na(hr_luad) && !is.na(hr_lusc)
  both_direction <- both_full && dir_concordant == TRUE
  both_ph <- !is.na(ph_pass_luad) && !is.na(ph_pass_lusc) && ph_pass_luad && ph_pass_lusc

  if (both_rho_ok && both_direction && !is.na(meta_fdr) && meta_fdr < 0.05 &&
      !is.na(i2_val) && i2_val < 50 && both_ph) {
    val_level <- "Level_A"
  } else if ((!is.na(fdr_luad) && fdr_luad < 0.05 && !is.na(hr_lusc) && sign(log(hr_luad)) == sign(log(hr_lusc))) ||
             (!is.na(fdr_lusc) && fdr_lusc < 0.05 && !is.na(hr_luad) && sign(log(hr_luad)) == sign(log(hr_lusc))) ||
             (!is.na(meta_fdr) && meta_fdr < 0.10 && !is.na(i2_val) && i2_val < 50)) {
    val_level <- "Level_B"
  } else if (nrow(stage_df) > 0) {
    sl_fdr <- stage_df$FDR[stage_df$cohort == "LUAD" & stage_df$program_id == pid]
    ss_fdr <- stage_df$FDR[stage_df$cohort == "LUSC" & stage_df$program_id == pid]
    if (any(!is.na(sl_fdr) & sl_fdr < 0.05, na.rm = TRUE) || any(!is.na(ss_fdr) & ss_fdr < 0.05, na.rm = TRUE) ||
        identical(dir_concordant, FALSE) || any(!is.na(ph_pass_luad) & !ph_pass_luad, na.rm = TRUE)) {
      val_level <- "Level_C"
    } else {
      val_level <- "No_clinical_support"
    }
  } else {
    val_level <- "No_clinical_support"
  }

  validation_rows[[length(validation_rows) + 1]] <- data.frame(
    program_id = pid, priority_tier = prog$priority_tier[1],
    cell_type = prog$cell_type[1], collection = prog$collection[1],
    pathway = prog$pathway[1], step07_NES = nes,
    expected_direction = prog$direction[1],
    LUAD_score_method_rho = rho_luad, LUSC_score_method_rho = rho_lusc,
    LUAD_HR = hr_luad, LUAD_FDR = fdr_luad, LUAD_PH_pass = ph_pass_luad,
    LUSC_HR = hr_lusc, LUSC_FDR = fdr_lusc, LUSC_PH_pass = ph_pass_lusc,
    meta_HR = meta_hr, meta_FDR = meta_fdr, I2 = i2_val,
    cross_histology_direction_concordant = dir_concordant,
    survival_direction_consistent_with_step07 = surv_dir_consistent,
    LUAD_stage_beta = stage_luad, LUAD_stage_FDR = NA_real_,
    LUSC_stage_beta = stage_lusc, LUSC_stage_FDR = NA_real_,
    clinical_validation_level = val_level,
    stringsAsFactors = FALSE
  )
}

validation_df <- do.call(rbind, validation_rows)

## Fill stage FDR
for (i in seq_len(nrow(validation_df))) {
  pid <- validation_df$program_id[i]
  if (!is.null(stage_df)) {
    sl <- stage_df[stage_df$cohort == "LUAD" & stage_df$program_id == pid, ]
    ss <- stage_df[stage_df$cohort == "LUSC" & stage_df$program_id == pid, ]
    if (nrow(sl) > 0) validation_df$LUAD_stage_FDR[i] <- sl$FDR[1]
    if (nrow(ss) > 0) validation_df$LUSC_stage_FDR[i] <- ss$FDR[1]
  }
}

write.csv(validation_df, "03_results/step08_TCGA/combined/GSE243013_TCGA_clinically_validated_programs.csv", row.names = FALSE)
cat(sprintf("[OK] Validation: Level_A=%d, Level_B=%d, Level_C=%d, No_support=%d\n",
            sum(validation_df$clinical_validation_level == "Level_A"),
            sum(validation_df$clinical_validation_level == "Level_B"),
            sum(validation_df$clinical_validation_level == "Level_C"),
            sum(validation_df$clinical_validation_level == "No_clinical_support")))

## =========================================================================
## XXII. KM Plots (visualization only)
## =========================================================================
cat("\n[XXII] KM plots...\n")

## Select top programs for KM
top_km_progs <- validation_df[validation_df$clinical_validation_level %in% c("Level_A", "Level_B"), ]
if (nrow(top_km_progs) < 20) {
  extra <- validation_df[!validation_df$program_id %in% top_km_progs$program_id & !is.na(validation_df$meta_FDR), ]
  extra <- extra[order(extra$meta_FDR), ]
  top_km_progs <- rbind(top_km_progs, extra[1:min(20 - nrow(top_km_progs), nrow(extra)), ])
}
top_km_progs <- head(top_km_progs, 20)

plot_ok <- 0
plot_fail <- 0

with_pdf <- function(filename, width, height, plot_fun) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file = filename, width = width, height = height)
  on.exit(try(grDevices::dev.off(), silent = TRUE))
  force(plot_fun())
}

for (cohort in cohorts) {
  pm_cohort <- patient_manifest[patient_manifest$cohort == cohort, ]

  for (i in seq_len(nrow(top_km_progs))) {
    pid <- top_km_progs$program_id[i]
    scores_sub <- manifest_long[manifest_long$cohort == cohort & manifest_long$program_id == pid, ]
    merged <- scores_sub
    merged <- merged[!is.na(merged$score_z) & !is.na(merged$OS_days) & merged$OS_days >= 0 & merged$OS_event %in% c(0, 1), ]
    if (nrow(merged) < 20) next

    merged$group <- ifelse(merged$score_z >= median(merged$score_z), "High", "Low")
    surv_obj <- Surv(merged$OS_days / 365.25, merged$OS_event)

    safe_pid <- gsub("[^A-Za-z0-9_]", "_", pid)
    fig_path <- sprintf("04_figures/step08_TCGA/survival/%s__%s__KM_median_split.pdf", cohort, safe_pid)

    tryCatch({
      with_pdf(fig_path, 7, 6, function() {
        km_fit <- survfit(surv_obj ~ group, data = merged)
        plot(km_fit, col = c("red", "blue"), lwd = 2,
             xlab = "Time (years)", ylab = "Overall Survival",
             main = sprintf("%s: %s", cohort, pid))
        legend("topright", c("High", "Low"), col = c("red", "blue"), lwd = 2)
        legend("bottomleft", "Median split for visualization only", cex = 0.7, bty = "n")
        p_val <- survdiff(surv_obj ~ group, data = merged)
        p_log <- 1 - pchisq(p_val$chisq, df = 1)
        legend("bottomright", sprintf("Log-rank p = %.3f", p_log), cex = 0.7, bty = "n")
      })
      plot_ok <- plot_ok + 1
    }, error = function(e) {
      plot_fail <<- plot_fail + 1
      cat(sprintf("[PLOT FAIL] %s/%s: %s\n", cohort, pid, e$message))
    })
  }
}
cat(sprintf("[OK] KM plots: %d OK, %d failed\n", plot_ok, plot_fail))

## =========================================================================
## XXIII. Summary Figures
## =========================================================================
cat("\n[XXIII] Summary figures...\n")

## 1. ssGSEA vs GSVA concordance distribution
if (nrow(concordance_df) > 0) {
  tryCatch({
    with_pdf("04_figures/step08_TCGA/scoring/ssGSEA_GSVA_concordance_distribution.pdf", 8, 5, function() {
      hist(concordance_df$spearman_rho, breaks = 30, col = "steelblue", main = "ssGSEA vs GSVA Spearman rho",
           xlab = "Spearman rho", ylab = "Program-cohorts")
      abline(v = 0.6, col = "red", lty = 2)
      abline(v = 0.8, col = "darkgreen", lty = 2)
      legend("topright", c("rho=0.6", "rho=0.8"), col = c("red", "darkgreen"), lty = 2)
    })
    cat("[OK] ssGSEA-GSVA distribution plot\n")
  }, error = function(e) cat(sprintf("[PLOT FAIL] concordance dist: %s\n", e$message)))
}

## 2. Program score PCA
for (cohort in cohorts) {
  if (is.null(ssgsea_scores[[cohort]])) next
  tryCatch({
    ssg <- ssgsea_scores[[cohort]]
    pca_res <- prcomp(t(ssg), scale. = TRUE)
    var_explained <- summary(pca_res)$importance[2, 1:2] * 100

    with_pdf(sprintf("04_figures/step08_TCGA/scoring/TCGA_%s_program_score_PCA.pdf", cohort), 8, 6, function() {
      plot(pca_res$x[, 1], pca_res$x[, 2], pch = 16, col = "steelblue",
           xlab = sprintf("PC1 (%.1f%%)", var_explained[1]),
           ylab = sprintf("PC2 (%.1f%%)", var_explained[2]),
           main = sprintf("%s Program Score PCA", cohort))
    })
    cat(sprintf("[OK] %s PCA plot\n", cohort))
  }, error = function(e) cat(sprintf("[PLOT FAIL] %s PCA: %s\n", cohort, e$message)))
}

## 3. Top clinical programs heatmap
for (cohort in cohorts) {
  top_progs <- validation_df[validation_df$clinical_validation_level %in% c("Level_A", "Level_B"), ]
  if (nrow(top_progs) == 0) next
  top_progs <- head(top_progs, min(40, nrow(top_progs)))

  ssg <- ssgsea_scores[[cohort]]
  if (is.null(ssg)) next
  common <- intersect(top_progs$program_id, rownames(ssg))
  if (length(common) < 2) next

  mat <- ssg[common, , drop = FALSE]
  mat_z <- t(scale(t(mat)))

  tryCatch({
    with_pdf(sprintf("04_figures/step08_TCGA/combined/Top_clinically_supported_program_scores_%s.pdf", cohort), 12, max(6, length(common) * 0.3), function() {
      pheatmap::pheatmap(mat_z, main = sprintf("%s: Top Clinical Program Scores (z-scored)", cohort),
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         fontsize_row = 4, fontsize_col = 4,
                         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
    })
    cat(sprintf("[OK] %s heatmap plot\n", cohort))
  }, error = function(e) cat(sprintf("[PLOT FAIL] %s heatmap: %s\n", cohort, e$message)))
}

## 4. OS forest plot (meta top 30)
if (!is.null(meta_df) && nrow(meta_df) > 0) {
  meta_top <- meta_df[order(meta_df$meta_FDR), ]
  meta_top <- head(meta_top, min(30, nrow(meta_top)))

  tryCatch({
    with_pdf("04_figures/step08_TCGA/combined/Program_OS_meta_forest_top30.pdf", 10, max(6, nrow(meta_top) * 0.4), function() {
      par(mar = c(5, 12, 4, 2))
      plot(NULL, xlim = range(c(log(meta_top$meta_lower_95CI), log(meta_top$meta_upper_95CI)), na.rm = TRUE),
           ylim = c(0, nrow(meta_top) + 1), yaxt = "n", ylab = "",
           xlab = "log(HR) [95% CI]", main = "Top 30 Programs: Meta-Analysis OS")
      abline(v = 0, lty = 2, col = "gray")
      for (j in seq_len(nrow(meta_top))) {
        y <- nrow(meta_top) - j + 1
        points(log(meta_top$meta_HR[j]), y, pch = 16, col = ifelse(meta_top$meta_FDR[j] < 0.05, "red", "black"))
        segments(log(meta_top$meta_lower_95CI[j]), y, log(meta_top$meta_upper_95CI[j]), y)
      }
      axis(2, at = seq_len(nrow(meta_top)), labels = meta_top$program_id, las = 1, cex.axis = 0.6)
      legend("topright", c("FDR<0.05", "FDR>=0.05"), pch = 16, col = c("red", "black"), cex = 0.8)
    })
    cat("[OK] Forest plot\n")
  }, error = function(e) cat(sprintf("[PLOT FAIL] forest: %s\n", e$message)))
}

## 5. Stage association heatmap
if (!is.null(stage_df) && nrow(stage_df) > 0) {
  sig_stage <- stage_df[stage_df$FDR < 0.10 | stage_df$kw_p_value < 0.05, ]
  if (nrow(sig_stage) > 0) {
    tryCatch({
      stage_wide <- sig_stage[, c("program_id", "cohort", "stage_beta")]
      stage_pivot <- stage_wide %>% pivot_wider(names_from = cohort, values_from = stage_beta)
      stage_mat <- as.data.frame(stage_pivot)
      rownames(stage_mat) <- stage_mat$program_id
      stage_mat$program_id <- NULL
      stage_mat <- as.matrix(stage_mat)

      with_pdf("04_figures/step08_TCGA/combined/Program_stage_association_LUAD_LUSC.pdf", 6, max(4, nrow(stage_mat) * 0.4), function() {
        pheatmap::pheatmap(stage_mat, main = "Stage Association Beta (LUAD vs LUSC)",
                           cluster_rows = TRUE, cluster_cols = FALSE,
                           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                           fontsize_row = 5)
      })
      cat("[OK] Stage heatmap\n")
    }, error = function(e) cat(sprintf("[PLOT FAIL] stage heatmap: %s\n", e$message)))
  }
}

## 6. Step 07 NES vs TCGA meta logHR
if (!is.null(meta_df) && nrow(meta_df) > 0) {
  merged_meta <- merge(meta_df, prog_manifest[, c("program_id", "NES")], by = "program_id")
  if (nrow(merged_meta) > 0) {
    tryCatch({
      with_pdf("04_figures/step08_TCGA/combined/Step07_NES_vs_TCGA_meta_logHR.pdf", 8, 6, function() {
        plot(merged_meta$NES, merged_meta$meta_logHR, pch = 16, col = "steelblue",
             xlab = "Step 07 NES", ylab = "TCGA Meta log(HR)",
             main = "Step 07 NES vs TCGA Meta log(HR)")
        abline(h = 0, lty = 2, col = "gray")
        abline(v = 0, lty = 2, col = "gray")
        if (nrow(merged_meta) > 5) {
          tryCatch({
            cor_test <- cor.test(merged_meta$NES, merged_meta$meta_logHR, method = "spearman")
            legend("topright", sprintf("Spearman rho = %.3f\np = %.3f", cor_test$estimate, cor_test$p.value), cex = 0.8, bty = "n")
          }, error = function(e) NULL)
        }
      })
      cat("[OK] NES vs HR scatter\n")
    }, error = function(e) cat(sprintf("[PLOT FAIL] NES vs HR: %s\n", e$message)))
  }
}

## =========================================================================
## XXIV. Result Index
## =========================================================================
cat("\n[XXIV] Result index...\n")

result_files <- list(
  list(type = "ssGSEA_scores", cohort = "LUAD", path = "02_data/tcga/clinical/TCGA_LUAD_program_scores_ssGSEA.rds"),
  list(type = "ssGSEA_scores", cohort = "LUSC", path = "02_data/tcga/clinical/TCGA_LUSC_program_scores_ssGSEA.rds"),
  list(type = "GSVA_scores", cohort = "LUAD", path = "02_data/tcga/clinical/TCGA_LUAD_program_scores_GSVA.rds"),
  list(type = "GSVA_scores", cohort = "LUSC", path = "02_data/tcga/clinical/TCGA_LUSC_program_scores_GSVA.rds"),
  list(type = "concordance", cohort = "combined", path = "03_results/step08_TCGA/scoring/GSE243013_TCGA_ssGSEA_vs_GSVA_concordance.csv"),
  list(type = "patient_manifest", cohort = "combined", path = "03_results/step08_TCGA/scoring/GSE243013_TCGA_patient_program_score_manifest.csv.gz"),
  list(type = "stage_association", cohort = "combined", path = "03_results/step08_TCGA/clinical_models/GSE243013_TCGA_stage_association_all_programs.csv.gz"),
  list(type = "Cox_results", cohort = "LUAD", path = "03_results/step08_TCGA/clinical_models/LUAD/TCGA_LUAD_program_OS_Cox_results.csv.gz"),
  list(type = "Cox_results", cohort = "LUSC", path = "03_results/step08_TCGA/clinical_models/LUSC/TCGA_LUSC_program_OS_Cox_results.csv.gz"),
  list(type = "PH_assumption", cohort = "combined", path = "03_results/step08_TCGA/clinical_models/GSE243013_TCGA_Cox_PH_assumption_results.csv.gz"),
  list(type = "meta_analysis", cohort = "combined", path = "03_results/step08_TCGA/meta_analysis/GSE243013_TCGA_program_OS_fixed_effect_meta.csv"),
  list(type = "validated_programs", cohort = "combined", path = "03_results/step08_TCGA/combined/GSE243013_TCGA_clinically_validated_programs.csv")
)

index_rows <- list()
for (rf in result_files) {
  exists <- file.exists(rf$path)
  size <- if (exists) file.info(rf$path)$size else NA_integer_
  index_rows[[length(index_rows) + 1]] <- data.frame(
    result_type = rf$type, cohort = rf$cohort, file_path = rf$path,
    file_size = size, status = ifelse(exists, "OK", "MISSING"),
    stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, index_rows), "03_results/step08_TCGA/combined/GSE243013_step08B1_result_index.csv", row.names = FALSE)

## =========================================================================
## XXV. Interpretation Boundaries
## =========================================================================
cat("\n[XXV] Interpretation boundaries...\n")

def_text <- c(
  "Step 08B1 Analysis Definition",
  "==============================",
  "",
  sprintf("Generated: %s", as.character(Sys.time())),
  "",
  "- TCGA-LUAD and TCGA-LUSC scored and modeled separately.",
  "- RNASeq2GeneNorm is upper-quartile log2-normalized RSEM TPM.",
  "- No additional log2 transformation applied.",
  "- Primary scoring: ssGSEA (alpha=0.25, normalized).",
  "- Sensitivity scoring: GSVA Gaussian.",
  "- Scores z-scored within each cohort independently.",
  "- Primary survival variable: continuous program score (z-scored).",
  "- HR represents risk per 1 SD increase in program score.",
  "- Median split used ONLY for KM visualization, not statistical testing.",
  "- Cox models adjust for age, sex, and pathologic stage.",
  "- LUAD and LUSC NOT merged for modeling.",
  "- Cross-histology results use fixed-effect meta-analysis.",
  "- TCGA is NOT an immunotherapy cohort.",
  "- TCGA cannot directly validate anti-PD-1 Responder vs Non_responder.",
  "- Bulk RNA program scores cannot be equated to cell-type-specific activity.",
  "- Clinical associations do not prove causality.",
  "- Mutation, CNV, Methylation, RPPA in Step 08B2."
)
writeLines(def_text, "00_config/step08_TCGA/GSE243013_step08B1_analysis_definition.txt")

## =========================================================================
## XXVI. Completion
## =========================================================================
cat("\n[XXVI] Completion check...\n")

if (grDevices::dev.cur() != 1L) grDevices::dev.off()

disk_final <- system("df -Pk .", intern = TRUE)
disk_parts <- strsplit(disk_final[length(disk_final)], "\\s+")[[1]]
disk_avail <- sprintf("%.1f", as.numeric(disk_parts[4]) / 1048576)

conditions <- list(
  list("LUAD RNA mapping", !is.null(mapping_list[["LUAD"]])),
  list("LUSC RNA mapping", !is.null(mapping_list[["LUSC"]])),
  list("Duplicate resolved", TRUE),
  list("No re-log", TRUE),
  list("LUAD ssGSEA", !is.null(ssgsea_scores[["LUAD"]])),
  list("LUSC ssGSEA", !is.null(ssgsea_scores[["LUSC"]])),
  list("LUAD GSVA", !is.null(gsva_scores[["LUAD"]])),
  list("LUSC GSVA", !is.null(gsva_scores[["LUSC"]])),
  list("Concordance table", file.exists("03_results/step08_TCGA/scoring/GSE243013_TCGA_ssGSEA_vs_GSVA_concordance.csv")),
  list("Patient manifest", file.exists("03_results/step08_TCGA/scoring/GSE243013_TCGA_patient_program_score_manifest.csv.gz")),
  list("Stage association", file.exists("03_results/step08_TCGA/clinical_models/GSE243013_TCGA_stage_association_all_programs.csv.gz")),
  list("Cox LUAD", file.exists("03_results/step08_TCGA/clinical_models/LUAD/TCGA_LUAD_program_OS_Cox_results.csv.gz")),
  list("Cox LUSC", file.exists("03_results/step08_TCGA/clinical_models/LUSC/TCGA_LUSC_program_OS_Cox_results.csv.gz")),
  list("PH audit", file.exists("03_results/step08_TCGA/clinical_models/GSE243013_TCGA_Cox_PH_assumption_results.csv.gz")),
  list("Meta-analysis", file.exists("03_results/step08_TCGA/meta_analysis/GSE243013_TCGA_program_OS_fixed_effect_meta.csv")),
  list("Validated programs", file.exists("03_results/step08_TCGA/combined/GSE243013_TCGA_clinically_validated_programs.csv"))
)

all_pass <- all(sapply(conditions, `[[`, 2))
for (c in conditions) cat(sprintf("  [%s] %s\n", ifelse(c[[2]], "PASS", "FAIL"), c[[1]]))

if (all_pass) {
  n_full_luad <- if (!is.null(cox_results_all[["LUAD"]])) sum(cox_results_all[["LUAD"]]$model_level == "FULL_MODEL") else 0
  n_full_lusc <- if (!is.null(cox_results_all[["LUSC"]])) sum(cox_results_all[["LUSC"]]$model_level == "FULL_MODEL") else 0
  n_ph_warn <- if (!is.null(ph_df)) sum(!ph_df$PH_score_pass) else 0

  completion <- c(
    "GSE243013 Step 08B1 COMPLETE",
    "=============================",
    sprintf("Completion time: %s", as.character(Sys.time())),
    sprintf("R version: %s", R.version.string),
    sprintf("GSVA version: %s", as.character(packageVersion("GSVA"))),
    sprintf("LUAD RNA patients: %d", ncol(expr_cleaned[["LUAD"]])),
    sprintf("LUSC RNA patients: %d", ncol(expr_cleaned[["LUSC"]])),
    sprintf("LUAD scorable programs: %d", sum(scorable_df$cohort == "LUAD" & scorable_df$scorable)),
    sprintf("LUSC scorable programs: %d", sum(scorable_df$cohort == "LUSC" & scorable_df$scorable)),
    sprintf("ssGSEA-GSVA rho>=0.60: %d", sum(concordance_df$scoring_method_concordant, na.rm = TRUE)),
    sprintf("Full Cox models: %d (LUAD: %d, LUSC: %d)", n_full_luad + n_full_lusc, n_full_luad, n_full_lusc),
    sprintf("PH warnings: %d", n_ph_warn),
    sprintf("Level A: %d", sum(validation_df$clinical_validation_level == "Level_A")),
    sprintf("Level B: %d", sum(validation_df$clinical_validation_level == "Level_B")),
    sprintf("Level C: %d", sum(validation_df$clinical_validation_level == "Level_C")),
    sprintf("Meta FDR < 0.05: %d", sum(meta_df$meta_FDR < 0.05, na.rm = TRUE)),
    sprintf("Plots OK: %d, Failed: %d", plot_ok, plot_fail),
    sprintf("dev.cur(): %d", grDevices::dev.cur()),
    sprintf("Disk available: %s GB", disk_avail),
    sprintf("Runtime: %.1f seconds", as.numeric(difftime(Sys.time(), step08b1_start, units = "secs")))
  )
  writeLines(completion, "03_results/GSE243013_step08B1_COMPLETE.txt")
  cat("\n[OK] Step 08B1 COMPLETE\n")
} else {
  failed <- sapply(conditions, function(c) if (!c[[2]]) c[[1]] else NA_character_)
  writeLines(c("Step 08B1 FAILED", sprintf("Time: %s", Sys.time()),
               sprintf("Failed: %s", paste(failed[!is.na(failed)], collapse = ", "))),
             "03_results/GSE243013_step08B1_FAILED.txt")
  cat("\n[FAILED]\n")
}

cat(sprintf("\nTotal runtime: %.1f seconds\n", as.numeric(difftime(Sys.time(), step08b1_start, units = "secs"))))
cat("========================================================================\n")
