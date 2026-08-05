## =========================================================================
## Step 06: edgeR Patient-Level Differential Expression
## =========================================================================

.libPaths(c(path.expand("~/Library/R/arm64/4.6/library"), .libPaths()))
options(stringsAsFactors = FALSE)

cat("========================================================================\n")
cat("Step 06: edgeR Patient-Level Differential Expression\n")
cat("========================================================================\n\n")
flush.console()

step_start <- Sys.time()
failed_step <- NA_character_
error_message <- NA_character_

## =========================================================================
## III. Load Packages
## =========================================================================
cat("[III] Loading packages...\n")
flush.console()

required_pkgs <- c("edgeR", "limma", "data.table", "ggplot2", "Matrix", "statmod")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly=TRUE)]

if (length(missing_pkgs) > 0) {
  cat(sprintf("[INFO] Missing packages: %s\n", paste(missing_pkgs, collapse=", ")))
  lib <- path.expand("~/Library/R/arm64/4.6/library")
  bioc_miss <- intersect(missing_pkgs, c("edgeR","limma"))
  if (length(bioc_miss) > 0) {
    if (!requireNamespace("BiocManager", quietly=TRUE)) {
      install.packages("BiocManager", repos="https://cloud.r-project.org", type="binary", lib=lib)
    }
    BiocManager::install(bioc_miss, ask=FALSE, update=FALSE, lib=lib)
  }
  cran_miss <- setdiff(missing_pkgs, c("edgeR","limma"))
  if (length(cran_miss) > 0) {
    install.packages(cran_miss, repos="https://cloud.r-project.org", type="binary", lib=lib)
  }
}

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(statmod))

## Close any residual graphics devices from prior sessions
graphics.off()
cat("[INFO] graphics.off() called — residual devices cleared.\n")

## =========================================================================
## Safe PDF helper — guarantees device closure on error
## =========================================================================
with_pdf <- function(filename, width = 7, height = 6, plot_fun) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  grDevices::pdf(
    file = filename,
    width  = width,
    height = height,
    onefile = TRUE
  )

  opened_device <- grDevices::dev.cur()

  on.exit({
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && opened_device %in% open_devices) {
      try(grDevices::dev.off(which = opened_device), silent = TRUE)
    }
  }, add = TRUE)

  force(plot_fun())
  invisible(filename)
}

env_info <- c(
  sprintf("R version: %s", R.version.string),
  sprintf("R platform: %s", R.version$platform),
  sprintf("edgeR version: %s", packageVersion("edgeR")),
  sprintf("limma version: %s", packageVersion("limma")),
  sprintf("data.table version: %s", packageVersion("data.table")),
  sprintf("ggplot2 version: %s", packageVersion("ggplot2")),
  sprintf("Matrix version: %s", packageVersion("Matrix")),
  sprintf("statmod version: %s", packageVersion("statmod")),
  sprintf("libPaths: %s", paste(.libPaths(), collapse=", ")),
  sprintf("Date: %s", Sys.time())
)

dir.create("03_results/step06_edgeR", recursive=TRUE, showWarnings=FALSE)
dir.create("03_results/step06_edgeR/primary_anti_PD1", recursive=TRUE, showWarnings=FALSE)
dir.create("03_results/step06_edgeR/strict_chemoimmunotherapy", recursive=TRUE, showWarnings=FALSE)
dir.create("03_results/step06_edgeR/histology_stratified", recursive=TRUE, showWarnings=FALSE)
dir.create("03_results/step06_edgeR/combined", recursive=TRUE, showWarnings=FALSE)
dir.create("03_results/step06_edgeR/qc", recursive=TRUE, showWarnings=FALSE)
dir.create("04_figures/step06_edgeR", recursive=TRUE, showWarnings=FALSE)
dir.create("04_figures/step06_edgeR/primary", recursive=TRUE, showWarnings=FALSE)
dir.create("04_figures/step06_edgeR/sensitivity", recursive=TRUE, showWarnings=FALSE)

writeLines(env_info, "03_results/step06_edgeR/GSE243013_step06_environment.txt")
cat("[INFO] Saved: GSE243013_step06_environment.txt\n")
cat(sprintf("[INFO] edgeR v%s, limma v%s\n", packageVersion("edgeR"), packageVersion("limma")))
flush.console()

## =========================================================================
## IV. Confirm Step 05 Status
## =========================================================================
cat("\n[IV] Confirming Step 05 status...\n")
flush.console()

if (!file.exists("03_results/GSE243013_step05_COMPLETE.txt")) {
  cat("[FATAL] Step 05 not complete. Cannot proceed.\n")
  stop("Step 05 incomplete")
}
cat("[INFO] Step 05 COMPLETE confirmed.\n")

file_index <- read.csv("03_results/step05_pseudobulk/GSE243013_pseudobulk_file_index.csv",
                        stringsAsFactors=FALSE)
cat(sprintf("[INFO] Pseudobulk file index: %d rows\n", nrow(file_index)))

primary_eligible_ct <- tryCatch(
  read.csv("03_results/step05_pseudobulk/GSE243013_primary_DE_celltypes_min20.csv",
           stringsAsFactors=FALSE),
  error=function(e) NULL
)

exploratory_eligible_ct <- tryCatch(
  read.csv("03_results/step05_pseudobulk/GSE243013_exploratory_DE_celltypes_min20.csv",
           stringsAsFactors=FALSE),
  error=function(e) NULL
)

manifest <- read.csv("03_results/GSE243013_patient_manifest_revised.csv", stringsAsFactors=FALSE)
cat(sprintf("[INFO] Patient manifest: %d rows\n", nrow(manifest)))

cat(sprintf("[INFO] Primary eligible CT file: %s (%d rows)\n",
            ifelse(is.null(primary_eligible_ct), "NULL", "OK"),
            ifelse(is.null(primary_eligible_ct), 0, nrow(primary_eligible_ct))))
cat(sprintf("[INFO] Exploratory eligible CT file: %s (%d rows)\n",
            ifelse(is.null(exploratory_eligible_ct), "NULL", "OK"),
            ifelse(is.null(exploratory_eligible_ct), 0, nrow(exploratory_eligible_ct))))
flush.console()

## =========================================================================
## IVb. Input File Audit
## =========================================================================
cat("\n[IVb] Auditing input files...\n")

audit_rows <- list()
for (i in seq_len(nrow(file_index))) {
  fi <- file_index[i, ]
  rds_exists <- file.exists(fi$counts_rds_path)
  csv_exists <- file.exists(fi$sample_metadata_path)
  rds_valid <- FALSE
  if (rds_exists) {
    tryCatch({
      tmp <- readRDS(fi$counts_rds_path)
      rds_valid <- is.matrix(tmp) && nrow(tmp) > 0 && ncol(tmp) > 0
      rm(tmp)
    }, error=function(e) NULL)
  }
  audit_rows[[i]] <- data.frame(
    annotation_level = fi$annotation_level,
    cell_type = fi$cell_type_original,
    counts_rds_exists = rds_exists,
    sample_meta_csv_exists = csv_exists,
    rds_valid = rds_valid,
    validation_status = fi$validation_status,
    n_genes = fi$n_genes,
    n_pseudobulk_samples = fi$n_pseudobulk_samples,
    pass = rds_exists && csv_exists && rds_valid &&
           fi$validation_status == "complete" &&
           fi$n_genes > 0 && fi$n_pseudobulk_samples > 0,
    stringsAsFactors = FALSE
  )
}
input_audit <- do.call(rbind, audit_rows)
write.csv(input_audit, "03_results/step06_edgeR/qc/GSE243013_step06_input_file_audit.csv",
          row.names=FALSE)

n_pass <- sum(input_audit$pass)
cat(sprintf("[INFO] Input audit: %d/%d files pass\n", n_pass, nrow(input_audit)))
if (n_pass == 0) stop("No valid input files found")

## =========================================================================
## V. Define Analysis Objects
## =========================================================================
cat("\n[V] Defining analysis objects...\n")
flush.console()

col_names <- names(file_index)
ct_col <- if ("cell_type_original" %in% col_names) "cell_type_original" else "cell_type"
safe_col <- if ("cell_type_safe" %in% col_names) "cell_type_safe" else "cell_type"

make_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  make.unique(x, sep="_")
}

all_immune_rows <- file_index[file_index$annotation_level == "all_immune", ]
major_rows <- file_index[file_index$annotation_level == "major_cell_type", ]
sub_rows <- file_index[file_index$annotation_level == "sub_cell_type", ]

## The eligible CT files use "cell_type" column, file_index uses "cell_type_original"
ect_ct_col <- if ("cell_type_original" %in% names(primary_eligible_ct)) "cell_type_original" else "cell_type"

if (!is.null(primary_eligible_ct) && nrow(primary_eligible_ct) > 0) {
  eligible_sub <- intersect(primary_eligible_ct[[ect_ct_col]], sub_rows[[ct_col]])
  primary_sub_rows <- sub_rows[sub_rows[[ct_col]] %in% eligible_sub, ]
} else {
  primary_sub_rows <- sub_rows[0, ]
}

primary_analysis_objects <- list()

primary_analysis_objects[["all_immune"]] <- list(
  annotation_level = "all_immune",
  cell_type_original = all_immune_rows[[ct_col]],
  cell_type_safe = make_safe_name(all_immune_rows[[ct_col]]),
  counts_rds = all_immune_rows$counts_rds_path,
  sample_meta_csv = all_immune_rows$sample_metadata_path
)

for (i in seq_len(nrow(major_rows))) {
  mr <- major_rows[i, ]
  ct_name <- mr[[ct_col]]
  safe <- make_safe_name(ct_name)
  primary_analysis_objects[[safe]] <- list(
    annotation_level = "major_cell_type",
    cell_type_original = ct_name,
    cell_type_safe = safe,
    counts_rds = mr$counts_rds_path,
    sample_meta_csv = mr$sample_metadata_path
  )
}

for (i in seq_len(nrow(primary_sub_rows))) {
  sr <- primary_sub_rows[i, ]
  ct_name <- sr[[ct_col]]
  safe <- make_safe_name(ct_name)
  primary_analysis_objects[[safe]] <- list(
    annotation_level = "sub_cell_type",
    cell_type_original = ct_name,
    cell_type_safe = safe,
    counts_rds = sr$counts_rds_path,
    sample_meta_csv = sr$sample_metadata_path
  )
}

cat(sprintf("[INFO] Primary analysis objects: %d\n", length(primary_analysis_objects)))
cat(sprintf("  all_immune: 1\n"))
cat(sprintf("  major_cell_type: %d\n", nrow(major_rows)))
cat(sprintf("  sub_cell_type (primary eligible): %d\n", nrow(primary_sub_rows)))
flush.console()

ct_list <- data.frame(
  name = names(primary_analysis_objects),
  annotation_level = sapply(primary_analysis_objects, function(x) x$annotation_level),
  cell_type_original = sapply(primary_analysis_objects, function(x) x$cell_type_original),
  cell_type_safe = sapply(primary_analysis_objects, function(x) x$cell_type_safe),
  stringsAsFactors = FALSE
)
write.csv(ct_list, "03_results/step06_edgeR/GSE243013_primary_analysis_celltypes.csv",
          row.names=FALSE)

if (!is.null(exploratory_eligible_ct) && nrow(exploratory_eligible_ct) > 0) {
  only_exploratory <- setdiff(exploratory_eligible_ct[[ct_col]], primary_sub_rows[[ct_col]])
  if (length(only_exploratory) > 0) {
    expl_df <- data.frame(cell_type = only_exploratory, stringsAsFactors=FALSE)
    write.csv(expl_df, "03_results/step06_edgeR/GSE243013_exploratory_only_celltypes.csv",
              row.names=FALSE)
    cat(sprintf("[INFO] Exploratory-only cell types: %d\n", length(only_exploratory)))
  }
}

## =========================================================================
## VI. Predefined Analysis Cohorts
## =========================================================================
cat("\n[VI] Defining analysis cohorts...\n")

primary_sample_ids <- manifest$sampleID[manifest$primary_analysis_eligible == TRUE &
                                         !is.na(manifest$primary_analysis_eligible)]
strict_sample_ids <- manifest$sampleID[manifest$strict_sensitivity_analysis_eligible == TRUE &
                                        !is.na(manifest$strict_sensitivity_analysis_eligible)]

cat(sprintf("[INFO] Primary anti-PD1 eligible patients: %d\n", length(primary_sample_ids)))
cat(sprintf("[INFO] Strict chemoimmunotherapy eligible patients: %d\n", length(strict_sample_ids)))

## =========================================================================
## VII. Audit Patient Metadata for Each Pseudobulk
## =========================================================================
cat("\n[VII] Auditing patient metadata for pseudobulk matrices...\n")
flush.console()

audit_and_prep_pseudobulk <- function(counts_rds, sample_meta_csv, analysis_cohort_ids,
                                       model_formula_str, analysis_name, cell_type_safe) {
  counts <- readRDS(counts_rds)
  sample_meta <- read.csv(sample_meta_csv, stringsAsFactors=FALSE)

  stopifnot(is.matrix(counts))
  stopifnot(all(colnames(counts) == sample_meta$sampleID))
  stopifnot(!any(duplicated(sample_meta$sampleID)))
  stopifnot(all(counts >= 0))
  stopifnot(all(is.finite(counts)))

  if (max(abs(counts - round(counts))) > 1e-6) {
    cat(sprintf("[WARNING] %s/%s: counts not integer-like\n", analysis_name, cell_type_safe))
  }

  lib_sizes <- colSums(counts)
  stopifnot(all(lib_sizes > 0))

  if ("library_size" %in% names(sample_meta)) {
    lib_diff <- abs(sample_meta$library_size - lib_sizes)
    n_bad <- sum(lib_diff > 1)
    if (n_bad > 0) {
      cat(sprintf("[WARNING] %s/%s: %d samples with lib_size mismatch\n",
                  analysis_name, cell_type_safe, n_bad))
    }
  }

  valid_response <- sample_meta$response_binary %in% c("Responder", "Non_responder")
  valid_cancer <- sample_meta$cancer_type %in% c("LUAD", "LUSC") |
                  is.na(sample_meta$cancer_type)

  in_cohort <- sample_meta$sampleID %in% analysis_cohort_ids

  keep <- in_cohort & valid_response & (lib_sizes > 0)

  counts <- counts[, keep, drop=FALSE]
  sample_meta <- sample_meta[keep, ]

  sample_meta$response_binary <- factor(sample_meta$response_binary,
                                         levels=c("Non_responder", "Responder"))
  if ("cancer_type" %in% names(sample_meta)) {
    sample_meta$cancer_type <- factor(sample_meta$cancer_type, levels=c("LUSC", "LUAD"))
  }

  n_resp <- sum(sample_meta$response_binary == "Responder")
  n_nonresp <- sum(sample_meta$response_binary == "Non_responder")
  n_luad <- sum(sample_meta$cancer_type == "LUAD", na.rm=TRUE)
  n_lusc <- sum(sample_meta$cancer_type == "LUSC", na.rm=TRUE)

  sample_audit <- data.frame(
    sampleID = sample_meta$sampleID,
    response_binary = as.character(sample_meta$response_binary),
    cancer_type = as.character(sample_meta$cancer_type),
    library_size = as.integer(colSums(counts)),
    in_cohort = TRUE,
    pass = TRUE,
    stringsAsFactors = FALSE
  )

  audit_path <- file.path("03_results/step06_edgeR/qc",
                           sprintf("%s__%s__sample_audit.csv", analysis_name, cell_type_safe))
  write.csv(sample_audit, audit_path, row.names=FALSE)

  list(
    counts = counts,
    sample_meta = sample_meta,
    n_samples = ncol(counts),
    n_resp = n_resp,
    n_nonresp = n_nonresp,
    n_luad = n_luad,
    n_lusc = n_lusc,
    n_genes_before = nrow(counts)
  )
}

## =========================================================================
## VIII. Covariate Balance Audit
## =========================================================================
cat("\n[VIII] Covariate balance audit...\n")
flush.console()

all_balance_rows <- list()

## =========================================================================
## X. edgeR Analysis Function
## =========================================================================
cat("\n[X] Defining edgeR analysis function...\n")
flush.console()

run_one_edger_analysis <- function(counts, sample_metadata, analysis_name,
                                    cell_type_original, cell_type_safe,
                                    output_directory, figure_directory,
                                    model_type="primary", minimum_cells=20,
                                    treat_lfc=log2(1.2)) {

  cat(sprintf("\n--- edgeR: %s / %s (%s) ---\n", analysis_name, cell_type_original, model_type))
  flush.console()

  result <- list(
    status = "FAILED_MODEL_ERROR",
    n_samples = ncol(counts),
    n_resp = sum(sample_metadata$response_binary == "Responder"),
    n_nonresp = sum(sample_metadata$response_binary == "Non_responder"),
    n_luad = sum(sample_metadata$cancer_type == "LUAD", na.rm=TRUE),
    n_lusc = sum(sample_metadata$cancer_type == "LUSC", na.rm=TRUE),
    n_genes_before = nrow(counts),
    n_genes_after = 0,
    qlf_n_up_resp = 0,
    qlf_n_up_nonresp = 0,
    qlf_n_sig = 0,
    treat_n_up_resp = 0,
    treat_n_up_nonresp = 0,
    treat_n_sig = 0,
    min_fdr = NA_real_,
    max_abs_logfc = NA_real_,
    warnings = character(0),
    qlf_table = NULL,
    treat_table = NULL,
    fit_obj = NULL
  )

  tryCatch({
    ## Design matrix
    if (model_type == "histology") {
      design <- model.matrix(~ response_binary, data=sample_metadata)
    } else {
      design <- model.matrix(~ cancer_type + response_binary, data=sample_metadata)
    }

    target_coef <- "response_binaryResponder"
    if (!(target_coef %in% colnames(design))) {
      cat(sprintf("[FATAL] Target coefficient '%s' not found in design\n", target_coef))
      result$status <- "FAILED_MODEL_ERROR"
      return(result)
    }

    qr_check <- qr(design)
    if (qr_check$rank < ncol(design)) {
      cat("[FATAL] Design matrix not full rank\n")
      result$status <- "SKIPPED_DESIGN_NOT_FULL_RANK"

      design_df <- as.data.frame(design)
      write.csv(design_df, file.path(output_directory,
                sprintf("%s__%s__design_matrix.csv", analysis_name, cell_type_safe)),
                row.names=FALSE)
      writeLines(c("Design matrix not full rank",
                     sprintf("QR rank: %d, ncol: %d", qr_check$rank, ncol(design)),
                     sprintf("Columns: %s", paste(colnames(design), collapse=", "))),
                 file.path("03_results/step06_edgeR/qc",
                           sprintf("%s__%s__design_audit.txt", analysis_name, cell_type_safe)))
      return(result)
    }

    design_df <- as.data.frame(design)
    write.csv(design_df, file.path(output_directory,
              sprintf("%s__%s__design_matrix.csv", analysis_name, cell_type_safe)),
              row.names=FALSE)

    design_audit <- c(
      sprintf("nrow(design) = %d, ncol(counts) = %d", nrow(design), ncol(counts)),
      sprintf("QR rank = %d, ncol(design) = %d", qr_check$rank, ncol(design)),
      sprintf("Target coef: %s (present: %s)", target_coef, target_coef %in% colnames(design)),
      sprintf("Model formula: %s",
              if (model_type == "histology") "~ response_binary" else "~ cancer_type + response_binary"),
      sprintf("Responder: %d, Non_responder: %d",
              sum(sample_metadata$response_binary == "Responder"),
              sum(sample_metadata$response_binary == "Non_responder"))
    )
    writeLines(design_audit,
               file.path("03_results/step06_edgeR/qc",
                         sprintf("%s__%s__design_audit.txt", analysis_name, cell_type_safe)))

    ## DGEList
    y <- edgeR::DGEList(counts=counts, samples=sample_metadata,
                          group=sample_metadata$response_binary)

    ## Low expression filtering
    keep <- edgeR::filterByExpr(y, group=sample_metadata$response_binary)
    n_after <- sum(keep)
    cat(sprintf("[INFO] filterByExpr: %d -> %d genes\n", nrow(y), n_after))

    if (n_after < 500) {
      cat(sprintf("[WARNING] Too few expressed genes (%d < 500), skipping\n", n_after))
      result$status <- "SKIPPED_TOO_FEW_EXPRESSED_GENES"
      result$n_genes_after <- n_after
      return(result)
    }

    y <- y[keep, , keep.lib.sizes=FALSE]

    ## TMM normalization
    y <- edgeR::normLibSizes(y, method="TMM")
    cat(sprintf("[INFO] TMM normalization done. NF range: [%.4f, %.4f]\n",
                min(y$samples$norm.factors), max(y$samples$norm.factors)))

    extreme_nf <- any(y$samples$norm.factors < 0.5 | y$samples$norm.factors > 2)
    if (extreme_nf) cat("[WARNING] Extreme normalization factors detected\n")

    ## Dispersion estimation
    y <- edgeR::estimateDisp(y, design=design, robust=TRUE)
    cat(sprintf("[INFO] Common dispersion: %.6f\n", y$common.dispersion))

    ## QL fit
    fit <- edgeR::glmQLFit(y, design=design, robust=TRUE)
    cat(sprintf("[INFO] QL prior df: %.2f\n", fit$df.prior))

    ## QLF test
    qlf <- edgeR::glmQLFTest(fit, coef=target_coef)
    qlf_table <- edgeR::topTags(qlf, n=Inf, sort.by="PValue")$table
    qlf_table$gene <- rownames(qlf_table)

    qlf_table$FDR_within_celltype <- qlf_table$FDR
    qlf_table$analysis_name <- analysis_name
    qlf_table$cell_type <- cell_type_original
    qlf_table$annotation_level <- model_type
    qlf_table$n_samples <- ncol(counts)
    qlf_table$n_responder <- sum(sample_metadata$response_binary == "Responder")
    qlf_table$n_nonresponder <- sum(sample_metadata$response_binary == "Non_responder")
    qlf_table$n_LUAD <- sum(sample_metadata$cancer_type == "LUAD", na.rm=TRUE)
    qlf_table$n_LUSC <- sum(sample_metadata$cancer_type == "LUSC", na.rm=TRUE)
    qlf_table$n_genes_before_filter <- nrow(y) + sum(!keep)
    qlf_table$n_genes_after_filter <- nrow(y)
    qlf_table$minimum_cells_per_patient <- minimum_cells
    qlf_table$model_formula <- if (model_type == "histology") "~ response_binary" else "~ cancer_type + response_binary"
    qlf_table$positive_logFC_direction <- "Higher_in_Responder"
    qlf_table$edgeR_version <- as.character(packageVersion("edgeR"))

    qlf_n_sig <- sum(qlf_table$FDR_within_celltype < 0.05, na.rm=TRUE)
    qlf_n_up_resp <- sum(qlf_table$FDR_within_celltype < 0.05 & qlf_table$logFC > 0, na.rm=TRUE)
    qlf_n_up_nonresp <- sum(qlf_table$FDR_within_celltype < 0.05 & qlf_table$logFC < 0, na.rm=TRUE)

    cat(sprintf("[INFO] QLF: %d FDR<0.05 (%d up Responder, %d up Non_responder)\n",
                qlf_n_sig, qlf_n_up_resp, qlf_n_up_nonresp))

    ## TREAT test
    treat <- edgeR::glmTreat(fit, coef=target_coef, lfc=treat_lfc, null="interval")
    treat_table <- edgeR::topTags(treat, n=Inf, sort.by="PValue")$table
    treat_table$gene <- rownames(treat_table)

    treat_table$FDR_within_celltype <- treat_table$FDR
    treat_table$treat_log2FC_threshold <- treat_lfc
    treat_table$analysis_name <- analysis_name
    treat_table$cell_type <- cell_type_original
    treat_table$n_samples <- ncol(counts)
    treat_table$n_responder <- sum(sample_metadata$response_binary == "Responder")
    treat_table$n_nonresponder <- sum(sample_metadata$response_binary == "Non_responder")
    treat_table$positive_logFC_direction <- "Higher_in_Responder"

    treat_n_sig <- sum(treat_table$FDR_within_celltype < 0.05, na.rm=TRUE)
    treat_n_up_resp <- sum(treat_table$FDR_within_celltype < 0.05 & treat_table$logFC > 0, na.rm=TRUE)
    treat_n_up_nonresp <- sum(treat_table$FDR_within_celltype < 0.05 & treat_table$logFC < 0, na.rm=TRUE)

    cat(sprintf("[INFO] TREAT: %d FDR<0.05 (%d up Responder, %d up Non_responder)\n",
                treat_n_sig, treat_n_up_resp, treat_n_up_nonresp))

    ## Save results
    qlf_path <- file.path(output_directory, sprintf("%s__edgeR_QLF_all_genes.csv.gz", cell_type_safe))
    con <- gzfile(qlf_path, "w")
    write.csv(qlf_table, con, row.names=FALSE)
    close(con)

    treat_path <- file.path(output_directory, sprintf("%s__edgeR_TREAT_all_genes.csv.gz", cell_type_safe))
    con <- gzfile(treat_path, "w")
    write.csv(treat_table, con, row.names=FALSE)
    close(con)

    fit_path <- file.path(output_directory, sprintf("%s__edgeR_fit_objects.rds", cell_type_safe))
    saveRDS(list(
      design = design,
      target_coef = target_coef,
      n_genes_before = nrow(y) + sum(!keep),
      n_genes_after = nrow(y),
      common_dispersion = y$common.dispersion,
      qlf_prior_df = fit$df.prior,
      norm_factors_range = range(y$samples$norm.factors),
      model_type = model_type
    ), fit_path)

    ## Save full edgeR objects for plot recovery
    edger_obj_path <- file.path(output_directory, sprintf("%s__edgeR_objects.rds", cell_type_safe))
    saveRDS(list(
      y = y, fit = fit, qlf = qlf, treat = treat,
      qlf_table = qlf_table, treat_table = treat_table,
      sample_metadata = sample_metadata, counts = counts
    ), edger_obj_path)

    ## Track plot recovery status
    plot_recovery_rows <- list()
    log_plot_status <- function(plot_type, plot_path, status, error_msg = "") {
      dev_before <- grDevices::dev.cur()
      plot_recovery_rows[[length(plot_recovery_rows) + 1]] <<- data.frame(
        analysis_name = analysis_name,
        cell_type = cell_type_original,
        plot_type = plot_type,
        plot_path = plot_path,
        file_exists = file.exists(plot_path),
        file_size_bytes = if (file.exists(plot_path)) file.info(plot_path)$size else 0L,
        status = status,
        error_message = error_msg,
        device_before = dev_before,
        device_after = grDevices::dev.cur(),
        stringsAsFactors = FALSE
      )
      if (grDevices::dev.cur() != 1L) {
        cat(sprintf("[WARNING] Device leak detected after %s (cur=%d)\n",
                    plot_type, grDevices::dev.cur()))
        graphics.off()
      }
    }

    ## ---- Plot 1: Library size by response ----
    lib_path <- file.path(figure_directory, sprintf("%s__library_size_by_response.pdf", cell_type_safe))
    tryCatch({
      lib_df <- data.frame(
        sampleID = sample_metadata$sampleID,
        response = as.character(sample_metadata$response_binary),
        lib_size = colSums(counts)
      )
      with_pdf(lib_path, width = 8, height = 5, plot_fun = function() {
        print(ggplot(lib_df, aes(x = reorder(sampleID, lib_size), y = lib_size, fill = response)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          theme_minimal() +
          labs(title = sprintf("Library Size: %s", cell_type_original),
               x = "Patient", y = "Library Size (UMI count)") +
          scale_fill_manual(values = c("Non_responder" = "#E41A1C", "Responder" = "#377EB8")))
      })
      log_plot_status("library_size_by_response", lib_path, "OK")
    }, error = function(e) {
      log_plot_status("library_size_by_response", lib_path, "FAILED_PLOT", conditionMessage(e))
      cat(sprintf("[WARNING] Library size plot failed: %s\n", e$message))
    })

    ## ---- Plot 2: MDS by response ----
    mds_resp_path <- file.path(figure_directory, sprintf("%s__MDS_by_response.pdf", cell_type_safe))
    tryCatch({
      with_pdf(mds_resp_path, width = 8, height = 6, plot_fun = function() {
        plotMDS(y, col = ifelse(sample_metadata$response_binary == "Responder", "#377EB8", "#E41A1C"),
                main = sprintf("MDS: %s (by response)", cell_type_original))
        legend("topright", legend = c("Responder", "Non_responder"),
               col = c("#377EB8", "#E41A1C"), pch = 19)
      })
      log_plot_status("MDS_by_response", mds_resp_path, "OK")
    }, error = function(e) {
      log_plot_status("MDS_by_response", mds_resp_path, "FAILED_PLOT", conditionMessage(e))
      cat(sprintf("[WARNING] MDS response plot failed: %s\n", e$message))
    })

    ## ---- Plot 3: MDS by histology ----
    mds_hist_path <- file.path(figure_directory, sprintf("%s__MDS_by_histology.pdf", cell_type_safe))
    tryCatch({
      with_pdf(mds_hist_path, width = 8, height = 6, plot_fun = function() {
        hist_colors <- ifelse(sample_metadata$cancer_type == "LUAD", "#4DAF4A",
                       ifelse(sample_metadata$cancer_type == "LUSC", "#984EA3", "grey50"))
        plotMDS(y, col = hist_colors,
                main = sprintf("MDS: %s (by histology)", cell_type_original))
        legend("topright", legend = c("LUAD", "LUSC"),
               col = c("#4DAF4A", "#984EA3"), pch = 19)
      })
      log_plot_status("MDS_by_histology", mds_hist_path, "OK")
    }, error = function(e) {
      log_plot_status("MDS_by_histology", mds_hist_path, "FAILED_PLOT", conditionMessage(e))
      cat(sprintf("[WARNING] MDS histology plot failed: %s\n", e$message))
    })

    ## ---- Plot 4: BCV ----
    bcv_path <- file.path(figure_directory, sprintf("%s__BCV.pdf", cell_type_safe))
    tryCatch({
      with_pdf(bcv_path, width = 8, height = 6, plot_fun = function() {
        plotBCV(y, main = sprintf("BCV: %s", cell_type_original))
      })
      log_plot_status("BCV", bcv_path, "OK")
    }, error = function(e) {
      log_plot_status("BCV", bcv_path, "FAILED_PLOT", conditionMessage(e))
      cat(sprintf("[WARNING] BCV plot failed: %s\n", e$message))
    })

    ## ---- Plot 5: QL dispersion ----
    qldisp_path <- file.path(figure_directory, sprintf("%s__QL_dispersion.pdf", cell_type_safe))
    tryCatch({
      with_pdf(qldisp_path, width = 8, height = 6, plot_fun = function() {
        plotQLDisp(fit, main = sprintf("QL Dispersion: %s", cell_type_original))
      })
      log_plot_status("QL_dispersion", qldisp_path, "OK")
    }, error = function(e) {
      log_plot_status("QL_dispersion", qldisp_path, "FAILED_PLOT", conditionMessage(e))
      cat(sprintf("[WARNING] QL Dispersion plot failed: %s\n", e$message))
    })

    ## ---- Plot 6: MD plot (QLF) ----
    md_path <- file.path(figure_directory, sprintf("%s__MD_QLF.pdf", cell_type_safe))
    tryCatch({
      with_pdf(md_path, width = 8, height = 6, plot_fun = function() {
        plotMD(qlf, main = sprintf("MD Plot: %s (QLF)", cell_type_original))
        abline(h = c(-1, 1), lty = 2)
      })
      log_plot_status("MD_QLF", md_path, "OK")
    }, error = function(e) {
      log_plot_status("MD_QLF", md_path, "FAILED_PLOT", conditionMessage(e))
      cat(sprintf("[WARNING] MD plot failed: %s\n", e$message))
    })

    ## ---- Plot 7: Volcano ----
    volc_path <- file.path(figure_directory, sprintf("%s__volcano.pdf", cell_type_safe))
    tryCatch({
      if (!is.data.frame(qlf_table) || nrow(qlf_table) == 0) {
        cat("[INFO] qlf_table is empty or not a data.frame, skipping volcano plot\n")
        log_plot_status("volcano", volc_path, "SKIPPED_EMPTY_TABLE")
      } else {
        ## Coerce to plain data.frame (handles data.table / tibble)
        qlf_df <- as.data.frame(qlf_table[, c("gene", "logFC", "FDR_within_celltype")])
        treat_df <- as.data.frame(treat_table[, c("gene", "FDR_within_celltype")])
        n_qlf <- nrow(qlf_df)
        if (n_qlf == 0 || !all(sapply(qlf_df, length) == n_qlf)) {
          cat(sprintf("[INFO] qlf_table column length mismatch (n=%d), skipping volcano\n", n_qlf))
          log_plot_status("volcano", volc_path, "SKIPPED_COLUMN_MISMATCH")
        } else {
          treat_sig_genes <- treat_df$gene[treat_df$FDR_within_celltype < 0.05]
          qlf_only_genes  <- qlf_df$gene[qlf_df$FDR_within_celltype < 0.05 &
                                          !(qlf_df$gene %in% treat_sig_genes)]
          volc_df <- data.frame(
            gene          = qlf_df$gene,
            logFC         = qlf_df$logFC,
            neg_log10_fdr = -log10(pmax(qlf_df$FDR_within_celltype, 1e-300)),
            category      = "Not_significant",
            stringsAsFactors = FALSE
          )
          volc_df$category[volc_df$gene %in% treat_sig_genes] <- "TREAT_sig"
          volc_df$category[volc_df$gene %in% qlf_only_genes]  <- "QLF_only"
          volc_df$category <- factor(volc_df$category,
                                     levels = c("Not_significant", "QLF_only", "TREAT_sig"))
          volc_df$label <- ""
          top_genes <- volc_df$gene[volc_df$category == "TREAT_sig"][
            1:min(20, sum(volc_df$category == "TREAT_sig"))]
          volc_df$label[volc_df$gene %in% top_genes] <- volc_df$gene[volc_df$gene %in% top_genes]

          with_pdf(volc_path, width = 9, height = 7, plot_fun = function() {
            print(ggplot(volc_df, aes(x = logFC, y = neg_log10_fdr, color = category)) +
              geom_point(size = 0.8, alpha = 0.6) +
              geom_text(aes(label = label), size = 2.5, vjust = -0.5,
                        check_overlap = TRUE, show.legend = FALSE) +
              scale_color_manual(values = c("Not_significant" = "grey60",
                                            "QLF_only" = "#FF7F00",
                                            "TREAT_sig" = "#E41A1C")) +
              theme_minimal() +
              labs(title = sprintf("Volcano: %s", cell_type_original),
                   x = "log2(Fold Change)", y = "-log10(FDR within cell type)",
                   subtitle = sprintf("TREAT sig: %d, QLF-only: %d",
                                      sum(volc_df$category == "TREAT_sig"),
                                      sum(volc_df$category == "QLF_only"))))
          })
          log_plot_status("volcano", volc_path, "OK")
        }
      }
    }, error = function(e) {
      log_plot_status("volcano", volc_path, "FAILED_PLOT", conditionMessage(e))
      cat(sprintf("[WARNING] Volcano plot failed: %s\n", e$message))
    })

    ## ---- Plot 8: Heatmap ----
    heat_path <- file.path(figure_directory, sprintf("%s__top_gene_heatmap.pdf", cell_type_safe))
    tryCatch({
      qlf_df_ht <- as.data.frame(qlf_table[, c("gene", "FDR_within_celltype")])
      sig_genes <- qlf_df_ht$gene[qlf_df_ht$FDR_within_celltype < 0.05]
      if (length(sig_genes) < 2) {
        cat("[INFO] Fewer than 2 significant genes, skipping heatmap\n")
        log_plot_status("heatmap", heat_path, "SKIPPED_FEW_GENES")
      } else {
        top_n   <- min(30, length(sig_genes))
        top_genes <- sig_genes[1:top_n]
        cpm_mat <- edgeR::cpm(y, log = TRUE, prior.count = 2)
        matched <- top_genes[top_genes %in% rownames(cpm_mat)]
        if (length(matched) < 2) {
          cat("[INFO] Fewer than 2 significant genes in CPM matrix, skipping heatmap\n")
          log_plot_status("heatmap", heat_path, "SKIPPED_NO_MATCH")
        } else {
          plot_mat <- cpm_mat[matched, , drop = FALSE]
          plot_mat <- sweep(plot_mat, 1, rowMeans(plot_mat), "-")
          ha <- data.frame(
            Response  = as.character(sample_metadata$response_binary),
            Histology = as.character(sample_metadata$cancer_type),
            row.names = colnames(plot_mat)
          )
          with_pdf(heat_path,
                   width  = max(8, ncol(plot_mat) * 0.15),
                   height = max(6, length(matched) * 0.25),
                   plot_fun = function() {
                     tryCatch({
                       library(pheatmap)
                       pheatmap::pheatmap(plot_mat, annotation_col = ha, show_rownames = TRUE,
                                          show_colnames = FALSE,
                                          main = sprintf("Top DE Genes: %s", cell_type_original),
                                          cluster_cols = TRUE, cluster_rows = TRUE, fontsize_row = 6)
                     }, error = function(e2) {
                       image(t(plot_mat[nrow(plot_mat):1, ]), axes = FALSE,
                             main = sprintf("Top DE Genes: %s", cell_type_original))
                     })
                   })
          log_plot_status("heatmap", heat_path, "OK")
        }
      }
    }, error = function(e) {
      log_plot_status("heatmap", heat_path, "FAILED_PLOT", conditionMessage(e))
      cat(sprintf("[WARNING] Heatmap failed: %s\n", e$message))
    })

    ## Save plot recovery status for this cell type
    if (length(plot_recovery_rows) > 0) {
      pr_df <- do.call(rbind, plot_recovery_rows)
      pr_path <- file.path("03_results/step06_edgeR/qc",
                           sprintf("%s__%s__plot_recovery.csv", analysis_name, cell_type_safe))
      write.csv(pr_df, pr_path, row.names = FALSE)
    }

    result$status <- "COMPLETE"
    result$n_genes_after <- nrow(y)
    result$qlf_n_sig <- qlf_n_sig
    result$qlf_n_up_resp <- qlf_n_up_resp
    result$qlf_n_up_nonresp <- qlf_n_up_nonresp
    result$treat_n_sig <- treat_n_sig
    result$treat_n_up_resp <- treat_n_up_resp
    result$treat_n_up_nonresp <- treat_n_up_nonresp
    result$min_fdr <- min(qlf_table$FDR_within_celltype, na.rm=TRUE)
    result$max_abs_logfc <- max(abs(qlf_table$logFC), na.rm=TRUE)
    result$qlf_table <- qlf_table
    result$treat_table <- treat_table

    if (extreme_nf) result$warnings <- c(result$warnings, "extreme_normalization_factors")

    ## Model diagnostics
    diag_df <- data.frame(
      metric = c("common_dispersion", "qlf_prior_df", "min_norm_factor", "max_norm_factor",
                  "min_lib_size", "max_lib_size", "n_genes_before", "n_genes_after",
                  "extreme_nf", "residual_df_check"),
      value = c(y$common.dispersion, fit$df.prior,
                min(y$samples$norm.factors), max(y$samples$norm.factors),
                min(y$samples$lib.size), max(y$samples$lib.size),
                nrow(y) + sum(!keep), nrow(y),
                as.numeric(extreme_nf), nrow(design) - qr_check$rank),
      stringsAsFactors = FALSE
    )
    write.csv(diag_df, file.path("03_results/step06_edgeR/qc",
              sprintf("%s__%s__model_diagnostics.csv", analysis_name, cell_type_safe)),
              row.names=FALSE)

    rm(y, fit, qlf, treat)
    gc()

  }, error=function(e) {
    result$status <<- "FAILED_MODEL_ERROR"
    result$warnings <<- c(result$warnings, conditionMessage(e))
    cat(sprintf("[ERROR] %s\n", conditionMessage(e)))
  })

  return(result)
}

cat("[INFO] Function defined.\n")
flush.console()

## =========================================================================
## XI. Run Primary Analysis (all cell types)
## =========================================================================
cat("\n[XI] Running primary anti-PD1 analysis...\n")
flush.console()

primary_results <- list()
model_status_df_rows <- list()

for (obj_name in names(primary_analysis_objects)) {
  obj <- primary_analysis_objects[[obj_name]]
  cat(sprintf("\n[PRIMARY] %s (%s)\n", obj$cell_type_original, obj$annotation_level))
  flush.console()

  ap <- audit_and_prep_pseudobulk(
    counts_rds = obj$counts_rds,
    sample_meta_csv = obj$sample_meta_csv,
    analysis_cohort_ids = primary_sample_ids,
    model_formula_str = "~ cancer_type + response_binary",
    analysis_name = "primary_anti_PD1",
    cell_type_safe = obj$cell_type_safe
  )

  if (ap$n_resp < 10 || ap$n_nonresp < 10) {
    cat(sprintf("[SKIP] Insufficient replication: Resp=%d, NonResp=%d\n", ap$n_resp, ap$n_nonresp))
    primary_results[[obj_name]] <- list(status="SKIPPED_INSUFFICIENT_REPLICATION", obj=obj)
    model_status_df_rows[[length(model_status_df_rows)+1]] <- data.frame(
      analysis="primary_anti_PD1", cell_type=obj$cell_type_original,
      annotation_level=obj$annotation_level, status="SKIPPED_INSUFFICIENT_REPLICATION",
      n_samples=ap$n_samples, n_resp=ap$n_resp, n_nonresp=ap$n_nonresp,
      stringsAsFactors=FALSE)
    next
  }

  ## Covariate balance
  resp_ids <- ap$sample_meta$sampleID[ap$sample_meta$response_binary == "Responder"]
  nonresp_ids <- ap$sample_meta$sampleID[ap$sample_meta$response_binary == "Non_responder"]
  all_balance_rows[[length(all_balance_rows)+1]] <- data.frame(
    analysis="primary_anti_PD1", cell_type=obj$cell_type_original,
    group="Responder", n=length(resp_ids),
    n_LUAD=sum(ap$sample_meta$cancer_type[ap$sample_meta$response_binary == "Responder"] == "LUAD", na.rm=TRUE),
    n_LUSC=sum(ap$sample_meta$cancer_type[ap$sample_meta$response_binary == "Responder"] == "LUSC", na.rm=TRUE),
    median_cells=median(ap$sample_meta$n_cells[ap$sample_meta$response_binary == "Responder"]),
    median_lib=median(colSums(ap$counts)[ap$sample_meta$response_binary == "Responder"]),
    stringsAsFactors=FALSE)
  all_balance_rows[[length(all_balance_rows)+1]] <- data.frame(
    analysis="primary_anti_PD1", cell_type=obj$cell_type_original,
    group="Non_responder", n=length(nonresp_ids),
    n_LUAD=sum(ap$sample_meta$cancer_type[ap$sample_meta$response_binary == "Non_responder"] == "LUAD", na.rm=TRUE),
    n_LUSC=sum(ap$sample_meta$cancer_type[ap$sample_meta$response_binary == "Non_responder"] == "LUSC", na.rm=TRUE),
    median_cells=median(ap$sample_meta$n_cells[ap$sample_meta$response_binary == "Non_responder"]),
    median_lib=median(colSums(ap$counts)[ap$sample_meta$response_binary == "Non_responder"]),
    stringsAsFactors=FALSE)

  res <- run_one_edger_analysis(
    counts = ap$counts,
    sample_metadata = ap$sample_meta,
    analysis_name = "primary_anti_PD1",
    cell_type_original = obj$cell_type_original,
    cell_type_safe = obj$cell_type_safe,
    output_directory = "03_results/step06_edgeR/primary_anti_PD1",
    figure_directory = "04_figures/step06_edgeR/primary",
    model_type = "primary",
    minimum_cells = 20,
    treat_lfc = log2(1.2)
  )

  primary_results[[obj_name]] <- list(status=res$status, result=res, obj=obj)

  model_status_df_rows[[length(model_status_df_rows)+1]] <- data.frame(
    analysis="primary_anti_PD1", cell_type=obj$cell_type_original,
    annotation_level=obj$annotation_level, status=res$status,
    n_samples=res$n_samples, n_resp=res$n_resp, n_nonresp=res$n_nonresp,
    n_luad=res$n_luad, n_lusc=res$n_lusc,
    n_genes_before=res$n_genes_before, n_genes_after=res$n_genes_after,
    qlf_n_sig=res$qlf_n_sig, qlf_n_up_resp=res$qlf_n_up_resp,
    qlf_n_up_nonresp=res$qlf_n_up_nonresp,
    treat_n_sig=res$treat_n_sig, treat_n_up_resp=res$treat_n_up_resp,
    treat_n_up_nonresp=res$treat_n_up_nonresp,
    min_fdr=res$min_fdr, max_abs_logfc=res$max_abs_logfc,
    warnings=paste(res$warnings, collapse="; "),
    stringsAsFactors=FALSE)

  cat(sprintf("[DONE] %s: %s\n", obj$cell_type_original, res$status))
  flush.console()
}

## Save covariate balance
if (length(all_balance_rows) > 0) {
  balance_df <- do.call(rbind, all_balance_rows)
  write.csv(balance_df, "03_results/step06_edgeR/qc/GSE243013_covariate_balance_long.csv",
            row.names=FALSE)
  cat("[INFO] Saved: covariate balance\n")
}

## =========================================================================
## XII. Run Strict Sensitivity Analysis
## =========================================================================
cat("\n[XII] Running strict chemoimmunotherapy sensitivity analysis...\n")
flush.console()

strict_results <- list()

for (obj_name in names(primary_analysis_objects)) {
  obj <- primary_analysis_objects[[obj_name]]
  cat(sprintf("\n[STRICT] %s\n", obj$cell_type_original))
  flush.console()

  ap <- audit_and_prep_pseudobulk(
    counts_rds = obj$counts_rds,
    sample_meta_csv = obj$sample_meta_csv,
    analysis_cohort_ids = strict_sample_ids,
    model_formula_str = "~ cancer_type + response_binary",
    analysis_name = "strict_chemoimmunotherapy",
    cell_type_safe = obj$cell_type_safe
  )

  if (ap$n_resp < 10 || ap$n_nonresp < 10) {
    cat(sprintf("[SKIP] Insufficient: Resp=%d, NonResp=%d\n", ap$n_resp, ap$n_nonresp))
    strict_results[[obj_name]] <- list(status="SKIPPED_INSUFFICIENT_REPLICATION")
    model_status_df_rows[[length(model_status_df_rows)+1]] <- data.frame(
      analysis="strict_chemoimmunotherapy", cell_type=obj$cell_type_original,
      annotation_level=obj$annotation_level, status="SKIPPED_INSUFFICIENT_REPLICATION",
      n_samples=ap$n_samples, n_resp=ap$n_resp, n_nonresp=ap$n_nonresp,
      stringsAsFactors=FALSE)
    next
  }

  res <- run_one_edger_analysis(
    counts = ap$counts,
    sample_metadata = ap$sample_meta,
    analysis_name = "strict_chemoimmunotherapy",
    cell_type_original = obj$cell_type_original,
    cell_type_safe = obj$cell_type_safe,
    output_directory = "03_results/step06_edgeR/strict_chemoimmunotherapy",
    figure_directory = "04_figures/step06_edgeR/sensitivity",
    model_type = "primary",
    minimum_cells = 20,
    treat_lfc = log2(1.2)
  )

  strict_results[[obj_name]] <- list(status=res$status, result=res)
  model_status_df_rows[[length(model_status_df_rows)+1]] <- data.frame(
    analysis="strict_chemoimmunotherapy", cell_type=obj$cell_type_original,
    annotation_level=obj$annotation_level, status=res$status,
    n_samples=res$n_samples, n_resp=res$n_resp, n_nonresp=res$n_nonresp,
    n_genes_before=res$n_genes_before, n_genes_after=res$n_genes_after,
    qlf_n_sig=res$qlf_n_sig, treat_n_sig=res$treat_n_sig,
    stringsAsFactors=FALSE)
  cat(sprintf("[DONE] %s: %s\n", obj$cell_type_original, res$status))
  flush.console()
}

## =========================================================================
## XIII. Run Histology-Stratified Analysis (LUAD and LUSC)
## =========================================================================
cat("\n[XIII] Running histology-stratified analysis...\n")
flush.console()

hist_results <- list()

for (hist in c("LUAD", "LUSC")) {
  cat(sprintf("\n[HISTOLOGY] %s\n", hist))
  hist_ids <- manifest$sampleID[manifest$primary_analysis_eligible == TRUE &
                                 manifest$cancer_type == hist &
                                 !is.na(manifest$primary_analysis_eligible)]

  for (obj_name in names(primary_analysis_objects)) {
    obj <- primary_analysis_objects[[obj_name]]
    safe_name <- sprintf("%s_%s", hist, obj$cell_type_safe)

    ap <- tryCatch(
      audit_and_prep_pseudobulk(
        counts_rds = obj$counts_rds,
        sample_meta_csv = obj$sample_meta_csv,
        analysis_cohort_ids = hist_ids,
        model_formula_str = "~ response_binary",
        analysis_name = sprintf("primary_%s", hist),
        cell_type_safe = safe_name
      ),
      error=function(e) { cat(sprintf("[SKIP] %s/%s: %s\n", hist, obj$cell_type_original, e$message)); NULL }
    )

    if (is.null(ap) || ap$n_resp < 10 || ap$n_nonresp < 10) {
      cat(sprintf("[SKIP] %s/%s: Resp=%d, NonResp=%d (need >=10 each)\n",
                  hist, obj$cell_type_original,
                  ifelse(is.null(ap), 0, ap$n_resp),
                  ifelse(is.null(ap), 0, ap$n_nonresp)))
      hist_results[[safe_name]] <- list(status="SKIPPED_INSUFFICIENT_REPLICATION")
      model_status_df_rows[[length(model_status_df_rows)+1]] <- data.frame(
        analysis=sprintf("primary_%s", hist), cell_type=obj$cell_type_original,
        annotation_level=obj$annotation_level, status="SKIPPED_INSUFFICIENT_REPLICATION",
        n_samples=ifelse(is.null(ap), 0, ap$n_samples),
        n_resp=ifelse(is.null(ap), 0, ap$n_resp),
        n_nonresp=ifelse(is.null(ap), 0, ap$n_nonresp),
        stringsAsFactors=FALSE)
      next
    }

    res <- run_one_edger_analysis(
      counts = ap$counts,
      sample_metadata = ap$sample_meta,
      analysis_name = sprintf("primary_%s", hist),
      cell_type_original = obj$cell_type_original,
      cell_type_safe = safe_name,
      output_directory = "03_results/step06_edgeR/histology_stratified",
      figure_directory = "04_figures/step06_edgeR/sensitivity",
      model_type = "histology",
      minimum_cells = 20,
      treat_lfc = log2(1.2)
    )

    hist_results[[safe_name]] <- list(status=res$status, result=res)
    model_status_df_rows[[length(model_status_df_rows)+1]] <- data.frame(
      analysis=sprintf("primary_%s", hist), cell_type=obj$cell_type_original,
      annotation_level=obj$annotation_level, status=res$status,
      n_samples=res$n_samples, n_resp=res$n_resp, n_nonresp=res$n_nonresp,
      n_genes_before=res$n_genes_before, n_genes_after=res$n_genes_after,
      qlf_n_sig=res$qlf_n_sig, treat_n_sig=res$treat_n_sig,
      stringsAsFactors=FALSE)
    cat(sprintf("[DONE] %s/%s: %s\n", hist, obj$cell_type_original, res$status))
  }
}

## =========================================================================
## XIV. Save Model Status
## =========================================================================
cat("\n[XIV] Saving model status...\n")

## Ensure consistent columns
all_status_cols <- unique(unlist(lapply(model_status_df_rows, names)))
for (i in seq_along(model_status_df_rows)) {
  missing_cols <- setdiff(all_status_cols, names(model_status_df_rows[[i]]))
  for (mc in missing_cols) model_status_df_rows[[i]][[mc]] <- NA
}

model_status_df <- do.call(rbind, model_status_df_rows)
write.csv(model_status_df, "03_results/step06_edgeR/combined/GSE243013_edgeR_model_status.csv",
          row.names=FALSE)

n_complete_primary <- sum(model_status_df$analysis == "primary_anti_PD1" &
                           model_status_df$status == "COMPLETE")
n_complete_strict <- sum(model_status_df$analysis == "strict_chemoimmunotherapy" &
                          model_status_df$status == "COMPLETE")
cat(sprintf("[INFO] Primary COMPLETE: %d\n", n_complete_primary))
cat(sprintf("[INFO] Strict COMPLETE: %d\n", n_complete_strict))

## =========================================================================
## XV. Primary Analysis Summary
## =========================================================================
cat("\n[XV] Primary analysis summary...\n")
flush.console()

summary_rows <- list()
for (nm in names(primary_results)) {
  pr <- primary_results[[nm]]
  if (!is.null(pr$result) && !is.null(pr$result$qlf_table)) {
    qt <- pr$result$qlf_table
    tt <- pr$result$treat_table
    summary_rows[[length(summary_rows)+1]] <- data.frame(
      cell_type = pr$obj$cell_type_original,
      annotation_level = pr$obj$annotation_level,
      status = pr$result$status,
      n_samples = pr$result$n_samples,
      n_resp = pr$result$n_resp,
      n_nonresp = pr$result$n_nonresp,
      n_luad = pr$result$n_luad,
      n_lusc = pr$result$n_lusc,
      n_genes_before = pr$result$n_genes_before,
      n_genes_after = pr$result$n_genes_after,
      qlf_fdr05_n = pr$result$qlf_n_sig,
      qlf_up_responder = pr$result$qlf_n_up_resp,
      qlf_up_nonresponder = pr$result$qlf_n_up_nonresp,
      treat_fdr05_n = pr$result$treat_n_sig,
      treat_up_responder = pr$result$treat_n_up_resp,
      treat_up_nonresponder = pr$result$treat_n_up_nonresp,
      min_fdr = pr$result$min_fdr,
      max_abs_logfc = pr$result$max_abs_logfc,
      warnings = paste(pr$result$warnings, collapse="; "),
      stringsAsFactors = FALSE
    )
  } else {
    summary_rows[[length(summary_rows)+1]] <- data.frame(
      cell_type = pr$obj$cell_type_original,
      annotation_level = pr$obj$annotation_level,
      status = pr$status,
      n_samples = NA, n_resp = NA, n_nonresp = NA,
      n_luad = NA, n_lusc = NA,
      n_genes_before = NA, n_genes_after = NA,
      qlf_fdr05_n = NA, qlf_up_responder = NA, qlf_up_nonresponder = NA,
      treat_fdr05_n = NA, treat_up_responder = NA, treat_up_nonresponder = NA,
      min_fdr = NA, max_abs_logfc = NA, warnings = "",
      stringsAsFactors = FALSE
    )
  }
}
primary_summary <- do.call(rbind, summary_rows)
write.csv(primary_summary,
          "03_results/step06_edgeR/combined/GSE243013_primary_edgeR_summary.csv",
          row.names=FALSE)
cat("[INFO] Saved: primary summary\n")

## =========================================================================
## XVI. Cross-Cell-Type Global FDR
## =========================================================================
cat("\n[XVI] Cross-cell-type global FDR...\n")
flush.console()

## Collect all primary QLF results
all_qlf_list <- list()
for (nm in names(primary_results)) {
  pr <- primary_results[[nm]]
  if (!is.null(pr$result) && !is.null(pr$result$qlf_table) && pr$result$status == "COMPLETE") {
    tbl <- pr$result$qlf_table
    tbl$gene_celltype <- paste(tbl$gene, tbl$cell_type, sep="___")
    all_qlf_list[[nm]] <- tbl
  }
}

if (length(all_qlf_list) > 0) {
  all_qlf_combined <- do.call(rbind, all_qlf_list)
  rownames(all_qlf_combined) <- NULL

  all_qlf_combined$FDR_global_primary <- p.adjust(all_qlf_combined$PValue, method="BH")

  write.csv(all_qlf_combined,
            "03_results/step06_edgeR/combined/GSE243013_primary_edgeR_all_celltypes_all_genes.csv.gz",
            row.names=FALSE)
  cat(sprintf("[INFO] Saved combined QLF: %d gene x celltype rows\n", nrow(all_qlf_combined)))
} else {
  all_qlf_combined <- data.frame()
  cat("[WARNING] No primary QLF results to combine\n")
}

## Collect all primary TREAT results
all_treat_list <- list()
for (nm in names(primary_results)) {
  pr <- primary_results[[nm]]
  if (!is.null(pr$result) && !is.null(pr$result$treat_table) && pr$result$status == "COMPLETE") {
    tbl <- pr$result$treat_table
    tbl$gene_celltype <- paste(tbl$gene, tbl$cell_type, sep="___")
    all_treat_list[[nm]] <- tbl
  }
}

if (length(all_treat_list) > 0) {
  all_treat_combined <- do.call(rbind, all_treat_list)
  rownames(all_treat_combined) <- NULL

  all_treat_combined$FDR_global_TREAT <- p.adjust(all_treat_combined$PValue, method="BH")

  write.csv(all_treat_combined,
            "03_results/step06_edgeR/combined/GSE243013_primary_TREAT_all_celltypes_all_genes.csv.gz",
            row.names=FALSE)
  cat(sprintf("[INFO] Saved combined TREAT: %d gene x celltype rows\n", nrow(all_treat_combined)))
} else {
  all_treat_combined <- data.frame()
  cat("[WARNING] No primary TREAT results to combine\n")
}

## =========================================================================
## XVII. Priority Candidate Gene Lists
## =========================================================================
cat("\n[XVII] Priority candidate gene lists...\n")

## QLF significant within cell type
if (nrow(all_qlf_combined) > 0) {
  qlf_sig <- all_qlf_combined[all_qlf_combined$FDR_within_celltype < 0.05, ]
  write.csv(qlf_sig,
            "03_results/step06_edgeR/combined/GSE243013_QLF_significant_within_celltype.csv",
            row.names=FALSE)
  cat(sprintf("[INFO] QLF significant within celltype: %d rows\n", nrow(qlf_sig)))
} else {
  qlf_sig <- data.frame()
}

## TREAT significant within cell type
if (nrow(all_treat_combined) > 0) {
  treat_sig <- all_treat_combined[all_treat_combined$FDR_within_celltype < 0.05, ]
  write.csv(treat_sig,
            "03_results/step06_edgeR/combined/GSE243013_TREAT_significant_within_celltype.csv",
            row.names=FALSE)
  cat(sprintf("[INFO] TREAT significant within celltype: %d rows\n", nrow(treat_sig)))
} else {
  treat_sig <- data.frame()
}

## Global TREAT priority (both within and global FDR < 0.05)
if (nrow(all_treat_combined) > 0) {
  global_priority <- all_treat_combined[
    all_treat_combined$FDR_within_celltype < 0.05 &
    all_treat_combined$FDR_global_TREAT < 0.05, ]
  write.csv(global_priority,
            "03_results/step06_edgeR/combined/GSE243013_TREAT_global_priority_genes.csv",
            row.names=FALSE)
  cat(sprintf("[INFO] Global TREAT priority: %d rows\n", nrow(global_priority)))
} else {
  global_priority <- data.frame()
}

## =========================================================================
## XVIII. Gene Recurrence Across Cell Types
## =========================================================================
cat("\n[XVIII] Gene recurrence across cell types...\n")

if (nrow(all_qlf_combined) > 0) {
  all_genes <- unique(all_qlf_combined$gene)
  n_ct_tested <- length(unique(all_qlf_combined$cell_type))

  recurrence_rows <- list()
  for (g in all_genes) {
    g_qlf <- all_qlf_combined[all_qlf_combined$gene == g, ]
    g_treat <- if (nrow(all_treat_combined) > 0) {
      all_treat_combined[all_treat_combined$gene == g, ]
    } else data.frame()

    n_tested <- nrow(g_qlf)
    n_qlf_sig <- sum(g_qlf$FDR_within_celltype < 0.05, na.rm=TRUE)
    n_treat_sig <- if (nrow(g_treat) > 0) sum(g_treat$FDR_within_celltype < 0.05, na.rm=TRUE) else 0
    n_up_resp <- sum(g_qlf$FDR_within_celltype < 0.05 & g_qlf$logFC > 0, na.rm=TRUE)
    n_up_nonresp <- sum(g_qlf$FDR_within_celltype < 0.05 & g_qlf$logFC < 0, na.rm=TRUE)

    dirs <- g_qlf$logFC[g_qlf$FDR_within_celltype < 0.05]
    consistent_dir <- if (length(dirs) >= 2) {
      all(dirs > 0) || all(dirs < 0)
    } else if (length(dirs) == 1) TRUE else NA

    recurrence_rows[[length(recurrence_rows)+1]] <- data.frame(
      gene = g,
      n_celltypes_tested = n_tested,
      n_qlf_sig = n_qlf_sig,
      n_treat_sig = n_treat_sig,
      n_up_in_responder = n_up_resp,
      n_up_in_nonresponder = n_up_nonresp,
      direction_consistent = consistent_dir,
      min_fdr = min(g_qlf$FDR_within_celltype, na.rm=TRUE),
      median_logfc = median(g_qlf$logFC, na.rm=TRUE),
      max_abs_logfc = max(abs(g_qlf$logFC), na.rm=TRUE),
      celltypes_affected = paste(g_qlf$cell_type[g_qlf$FDR_within_celltype < 0.05], collapse="; "),
      stringsAsFactors = FALSE
    )
  }

  recurrence_df <- do.call(rbind, recurrence_rows)
  recurrence_df <- recurrence_df[order(-recurrence_df$n_treat_sig, -recurrence_df$n_qlf_sig,
                                        recurrence_df$min_fdr), ]

  write.csv(recurrence_df,
            "03_results/step06_edgeR/combined/GSE243013_gene_recurrence_across_celltypes.csv",
            row.names=FALSE)
  cat(sprintf("[INFO] Saved gene recurrence: %d genes\n", nrow(recurrence_df)))

  ## Cell-type specific TREAT candidates
  if (nrow(global_priority) > 0) {
    ct_specific <- global_priority[global_priority$FDR_within_celltype < 0.05, ]
    ct_specific_genes <- unique(ct_specific$gene)
    ct_specific_list <- list()
    for (g in ct_specific_genes) {
      g_all <- all_treat_combined[all_treat_combined$gene == g, ]
      if (nrow(g_all) == 1) {
        ct_specific_list[[length(ct_specific_list)+1]] <- g_all
      }
    }
    if (length(ct_specific_list) > 0) {
      ct_spec_df <- do.call(rbind, ct_specific_list)
      write.csv(ct_spec_df,
                "03_results/step06_edgeR/combined/GSE243013_celltype_specific_TREAT_candidates.csv",
                row.names=FALSE)
      cat(sprintf("[INFO] Cell-type specific TREAT candidates: %d rows\n", nrow(ct_spec_df)))
    }
  }
}

## =========================================================================
## XIX. Primary vs Strict Sensitivity Concordance
## =========================================================================
cat("\n[XIX] Primary vs strict sensitivity concordance...\n")
flush.console()

if (length(primary_results) > 0 && length(strict_results) > 0) {
  concordance_rows <- list()

  for (nm in names(primary_results)) {
    pr <- primary_results[[nm]]
    sr <- strict_results[[nm]]

    if (is.null(pr$result) || pr$result$status != "COMPLETE") next
    if (is.null(sr$result) || sr$result$status != "COMPLETE") next

    p_qlf <- pr$result$qlf_table
    s_qlf <- sr$result$qlf_table
    if (is.null(p_qlf) || nrow(p_qlf) == 0) next
    if (is.null(s_qlf) || nrow(s_qlf) == 0) next

    merged <- merge(
      data.frame(gene=p_qlf$gene, logFC_primary=p_qlf$logFC, FDR_primary=p_qlf$FDR_within_celltype,
                 stringsAsFactors=FALSE),
      data.frame(gene=s_qlf$gene, logFC_strict=s_qlf$logFC, FDR_strict=s_qlf$FDR_within_celltype,
                 stringsAsFactors=FALSE),
      by="gene", all=TRUE
    )
    merged$cell_type <- pr$obj$cell_type_original
    merged$direction_consistent <- sign(merged$logFC_primary) == sign(merged$logFC_strict)
    merged$both_sig <- merged$FDR_primary < 0.05 & merged$FDR_strict < 0.05

    if (sum(!is.na(merged$logFC_primary) & !is.na(merged$logFC_strict)) >= 3) {
      merged$correlation <- cor(merged$logFC_primary, merged$logFC_strict, use="complete.obs")
    } else {
      merged$correlation <- NA
    }

    concordance_rows[[length(concordance_rows)+1]] <- merged
  }

  if (length(concordance_rows) > 0) {
    concordance_df <- do.call(rbind, concordance_rows)
    write.csv(concordance_df,
              "03_results/step06_edgeR/combined/GSE243013_primary_vs_strict_concordance.csv",
              row.names=FALSE)
    cat(sprintf("[INFO] Saved primary vs strict concordance: %d rows\n", nrow(concordance_df)))
  }
}

## =========================================================================
## XX. Primary vs Histology-Stratified Concordance
## =========================================================================
cat("\n[XX] Primary vs histology-stratified concordance...\n")

for (hist in c("LUAD", "LUSC")) {
  hist_conc_rows <- list()
  for (nm in names(primary_results)) {
    pr <- primary_results[[nm]]
    hist_key <- sprintf("%s_%s", hist, pr$obj$cell_type_safe)

    if (is.null(pr$result) || pr$result$status != "COMPLETE") next
    if (is.null(hist_results[[hist_key]])) next
    hr <- hist_results[[hist_key]]
    if (is.null(hr$result) || hr$result$status != "COMPLETE") next

    p_qlf <- pr$result$qlf_table
    h_qlf <- hr$result$qlf_table
    if (is.null(p_qlf) || nrow(p_qlf) == 0) next
    if (is.null(h_qlf) || nrow(h_qlf) == 0) next

    merged <- merge(
      data.frame(gene=p_qlf$gene, logFC_primary=p_qlf$logFC, FDR_primary=p_qlf$FDR_within_celltype,
                 stringsAsFactors=FALSE),
      data.frame(gene=h_qlf$gene, logFC_hist=h_qlf$logFC, FDR_hist=h_qlf$FDR_within_celltype,
                 stringsAsFactors=FALSE),
      by="gene", all=TRUE
    )
    merged$cell_type <- pr$obj$cell_type_original
    merged$histology <- hist
    merged$direction_consistent <- sign(merged$logFC_primary) == sign(merged$logFC_hist)

    hist_conc_rows[[length(hist_conc_rows)+1]] <- merged
  }

  if (length(hist_conc_rows) > 0) {
    hist_conc_df <- do.call(rbind, hist_conc_rows)
    write.csv(hist_conc_df,
              sprintf("03_results/step06_edgeR/combined/GSE243013_primary_vs_%s_concordance.csv", hist),
              row.names=FALSE)
    cat(sprintf("[INFO] Saved primary vs %s concordance\n", hist))
  }
}

## =========================================================================
## XXI. Scientific Interpretation Boundary
## =========================================================================
cat("\n[XXI] Creating analysis definition...\n")

def_text <- c(
  "GSE243013 edgeR Analysis Definition",
  "====================================",
  "",
  "Statistical unit: patient sampleID (NOT individual cells)",
  "Input: patient-level raw UMI pseudobulk counts from Step 05",
  "Primary comparison: Responder vs Non_responder",
  "Responder definition: pCR + MPR",
  "Non_responder definition: non-MPR",
  "Primary cohort: anti-PD1 treated patients with known pathological response",
  "Minimum cells per patient per cell type: 20",
  "Primary model adjusts for LUAD vs LUSC histology",
  "Positive logFC: Responder expression higher",
  "Negative logFC: Non_responder expression higher",
  "QLF: general differential expression test",
  "glmTreat: formal test for >1.2x fold change threshold",
  "FDR_within_celltype: inference within single cell type",
  "FDR_global: cross-cell-type strict candidate screening",
  "Strict chemoimmunotherapy cohort: sensitivity analysis",
  "LUAD and LUSC stratified models: exploratory sensitivity analysis",
  "No per-cell statistical tests performed",
  "Results interpreted as post-treatment pathological response-associated immune expression programs",
  "Pathway enrichment and multi-omics integration not yet performed"
)
writeLines(def_text, "00_config/GSE243013_edgeR_analysis_definition.txt")
cat("[INFO] Saved: GSE243013_edgeR_analysis_definition.txt\n")

## =========================================================================
## XXI-b. Aggregate Plot Recovery Status
## =========================================================================
cat("\n[XXI-b] Aggregating plot recovery status...\n")

plot_recovery_files <- list.files("03_results/step06_edgeR/qc",
                                  pattern = "__plot_recovery\\.csv$",
                                  full.names = TRUE)
if (length(plot_recovery_files) > 0) {
  all_plot_recovery <- do.call(rbind, lapply(plot_recovery_files, read.csv, stringsAsFactors = FALSE))
  write.csv(all_plot_recovery,
            "03_results/step06_edgeR/qc/GSE243013_plot_recovery_status.csv",
            row.names = FALSE)

  n_total_plots    <- nrow(all_plot_recovery)
  n_ok             <- sum(all_plot_recovery$status == "OK", na.rm = TRUE)
  n_failed_plots   <- sum(all_plot_recovery$status == "FAILED_PLOT", na.rm = TRUE)
  n_skipped        <- sum(grepl("^SKIPPED", all_plot_recovery$status), na.rm = TRUE)
  n_device_leaks   <- sum(all_plot_recovery$device_after != 1L, na.rm = TRUE)

  cat(sprintf("[INFO] Total plots attempted: %d\n", n_total_plots))
  cat(sprintf("[INFO]   OK: %d\n", n_ok))
  cat(sprintf("[INFO]   FAILED_PLOT: %d\n", n_failed_plots))
  cat(sprintf("[INFO]   SKIPPED: %d\n", n_skipped))
  cat(sprintf("[INFO]   Device leaks remaining: %d\n", n_device_leaks))

  if (n_failed_plots > 0) {
    cat("\n--- Failed plots ---\n")
    failed <- all_plot_recovery[all_plot_recovery$status == "FAILED_PLOT", ]
    for (i in seq_len(nrow(failed))) {
      cat(sprintf("  %s / %s: %s\n", failed$cell_type[i], failed$plot_type[i],
                  failed$error_message[i]))
    }
  }
} else {
  cat("[WARNING] No per-cell-type plot recovery CSVs found\n")
}

## Final device check
remaining_devs <- grDevices::dev.list()
if (!is.null(remaining_devs) && length(remaining_devs) > 0) {
  cat(sprintf("[WARNING] %d graphics devices still open — calling graphics.off()\n",
              length(remaining_devs)))
  graphics.off()
} else {
  cat("[INFO] All graphics devices closed (dev.cur() == 1)\n")
}

## =========================================================================
## XXII. Completion Marker
## =========================================================================
cat("\n[XXII] Creating completion marker...\n")
flush.console()

total_runtime <- as.numeric(difftime(Sys.time(), step_start, units="secs"))

n_complete <- sum(model_status_df$status == "COMPLETE", na.rm=TRUE)
n_skipped <- sum(grepl("^SKIPPED", model_status_df$status), na.rm=TRUE)
n_failed <- sum(grepl("^FAILED", model_status_df$status), na.rm=TRUE)

qlf_sig_gene_ct <- if (nrow(all_qlf_combined) > 0) {
  sum(all_qlf_combined$FDR_within_celltype < 0.05, na.rm=TRUE)
} else 0

treat_sig_gene_ct <- if (nrow(all_treat_combined) > 0) {
  sum(all_treat_combined$FDR_within_celltype < 0.05, na.rm=TRUE)
} else 0

global_treat_gene_ct <- if (nrow(global_priority) > 0) {
  nrow(global_priority)
} else 0

all_checks_pass <- (
  n_complete_primary >= 1 &&
  file.exists("03_results/GSE243013_step05_COMPLETE.txt") &&
  file.exists("03_results/step06_edgeR/combined/GSE243013_edgeR_model_status.csv") &&
  file.exists("03_results/step06_edgeR/combined/GSE243013_primary_edgeR_summary.csv")
)

disk_after <- tryCatch({
  as.numeric(system("df -Pk . | tail -1 | awk '{print $4}'", intern=TRUE)) / (1024*1024)
}, error=function(e) NA_real_)

if (all_checks_pass) {
  complete_text <- c(
    "GSE243013 Step 06 COMPLETE",
    "===========================",
    "",
    sprintf("Completion time: %s", Sys.time()),
    sprintf("R version: %s", R.version.string),
    sprintf("edgeR version: %s", packageVersion("edgeR")),
    sprintf("limma version: %s", packageVersion("limma")),
    sprintf("Primary analysis cell types: %d", length(primary_analysis_objects)),
    sprintf("Primary COMPLETE models: %d", n_complete_primary),
    sprintf("Strict COMPLETE models: %d", n_complete_strict),
    sprintf("SKIPPED models: %d", n_skipped),
    sprintf("FAILED models: %d", n_failed),
    sprintf("QLF significant gene-celltype: %d", qlf_sig_gene_ct),
    sprintf("TREAT significant gene-celltype: %d", treat_sig_gene_ct),
    sprintf("Global TREAT priority gene-celltype: %d", global_treat_gene_ct),
    sprintf("Disk available: %.1f GB", disk_after),
    sprintf("Total runtime: %.1f seconds", total_runtime),
    "",
    "Statistical unit: patient sampleID",
    "No per-cell tests performed",
    "No Seurat objects created",
    "No pseudobulk RDS files modified"
  )
  writeLines(complete_text, "03_results/GSE243013_step06_COMPLETE.txt")
  cat("[INFO] Step 06 COMPLETE marker saved.\n")
} else {
  failed_text <- c(
    "GSE243013 Step 06 FAILED or INCOMPLETE",
    "=======================================",
    "",
    sprintf("Time: %s", Sys.time()),
    sprintf("Primary COMPLETE models: %d", n_complete_primary),
    sprintf("Total models attempted: %d", nrow(model_status_df)),
    sprintf("Check: all_checks_pass = %s", all_checks_pass),
    sprintf("R version: %s", R.version.string),
    sprintf("edgeR version: %s", packageVersion("edgeR"))
  )
  writeLines(failed_text, "03_results/GSE243013_step06_FAILED.txt")
  cat("[WARNING] Step 06 FAILED marker saved.\n")
}

## =========================================================================
## XXIII. Final Summary
## =========================================================================
cat("\n========================================================================\n")
cat("STEP 06 FINAL SUMMARY\n")
cat("========================================================================\n\n")

cat(sprintf("edgeR version: %s\n", packageVersion("edgeR")))
cat(sprintf("Primary analysis cell types: %d\n", length(primary_analysis_objects)))
cat(sprintf("Primary COMPLETE: %d\n", n_complete_primary))
cat(sprintf("Strict COMPLETE: %d\n", n_complete_strict))
cat(sprintf("SKIPPED: %d\n", n_skipped))
cat(sprintf("FAILED: %d\n", n_failed))
cat(sprintf("QLF FDR<0.05 gene-celltype: %d\n", qlf_sig_gene_ct))
cat(sprintf("TREAT FDR<0.05 gene-celltype: %d\n", treat_sig_gene_ct))
cat(sprintf("Global TREAT priority gene-celltype: %d\n", global_treat_gene_ct))

if (nrow(primary_summary) > 0) {
  cat("\n--- Top 10 cell types by QLF signal ---\n")
  qlf_sorted <- primary_summary[!is.na(primary_summary$qlf_fdr05_n), ]
  qlf_sorted <- qlf_sorted[order(-qlf_sorted$qlf_fdr05_n), ]
  for (i in seq_len(min(10, nrow(qlf_sorted)))) {
    cat(sprintf("  %s: %d QLF sig, %d TREAT sig\n",
                qlf_sorted$cell_type[i], qlf_sorted$qlf_fdr05_n[i],
                qlf_sorted$treat_fdr05_n[i]))
  }

  cat("\n--- Top 10 cell types by TREAT signal ---\n")
  treat_sorted <- primary_summary[!is.na(primary_summary$treat_fdr05_n), ]
  treat_sorted <- treat_sorted[order(-treat_sorted$treat_fdr05_n), ]
  for (i in seq_len(min(10, nrow(treat_sorted)))) {
    cat(sprintf("  %s: %d TREAT sig\n", treat_sorted$cell_type[i], treat_sorted$treat_fdr05_n[i]))
  }
}

if (exists("recurrence_df") && nrow(recurrence_df) > 0) {
  cat("\n--- Top 20 genes by cross-cell-type recurrence ---\n")
  for (i in seq_len(min(20, nrow(recurrence_df)))) {
    cat(sprintf("  %s: tested in %d CTs, %d QLF sig, %d TREAT sig, consistent=%s\n",
                recurrence_df$gene[i], recurrence_df$n_celltypes_tested[i],
                recurrence_df$n_qlf_sig[i], recurrence_df$n_treat_sig[i],
                recurrence_df$direction_consistent[i]))
  }
}

## Check design matrix issues
design_issues <- model_status_df$status[model_status_df$status == "SKIPPED_DESIGN_NOT_FULL_RANK"]
cat(sprintf("\nDesign matrix rank issues: %d\n", length(design_issues)))

## Check extreme NF
nf_warnings <- model_status_df$warnings[grepl("extreme", model_status_df$warnings)]
cat(sprintf("Extreme normalization factor warnings: %d\n", length(nf_warnings)))

cat(sprintf("\nStep 06 status: %s\n", ifelse(all_checks_pass, "COMPLETE", "INCOMPLETE")))
cat(sprintf("Ready for Step 07 pathway activity and TF analysis: %s\n",
            ifelse(all_checks_pass, "YES", "NO")))

cat("\n========================================================================\n")
cat("Step 06 completed.\n")
cat("========================================================================\n")
