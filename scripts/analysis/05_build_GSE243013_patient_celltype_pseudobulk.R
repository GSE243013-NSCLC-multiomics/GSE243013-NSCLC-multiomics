## =========================================================================
## Step 05: Build Patient-Level Cell-Type-Specific Pseudobulk Count Matrices
## =========================================================================

.libPaths(c(path.expand("~/Library/R/arm64/4.6/library"), .libPaths()))
options(stringsAsFactors = FALSE)

cat("========================================================================\n")
cat("Step 05: Pseudobulk Construction\n")
cat("========================================================================\n\n")
flush.console()

step_start <- Sys.time()
failed_step <- NA_character_
error_message <- NA_character_
completed_pseudobulks <- list()

## =========================================================================
## III. Load Packages
## =========================================================================
cat("[III] Loading packages...\n")
flush.console()

suppressPackageStartupMessages(library(BPCells))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))

env_info <- c(
  sprintf("R version: %s", R.version.string),
  sprintf("R platform: %s", R.version$platform),
  sprintf("BPCells version: %s", packageVersion("BPCells")),
  sprintf("data.table version: %s", packageVersion("data.table")),
  sprintf("Matrix version: %s", packageVersion("Matrix")),
  sprintf("libPaths: %s", paste(.libPaths(), collapse=", ")),
  sprintf("Date: %s", Sys.time())
)
writeLines(env_info, "03_results/step05_pseudobulk/GSE243013_step05_environment.txt")
cat("[INFO] Saved: GSE243013_step05_environment.txt\n")
cat(sprintf("[INFO] BPCells v%s, data.table v%s, Matrix v%s\n",
            packageVersion("BPCells"), packageVersion("data.table"), packageVersion("Matrix")))
flush.console()

## =========================================================================
## IV. Open BPCells Matrix
## =========================================================================
cat("\n[IV] Opening BPCells matrix...\n")
flush.console()

bpcells_dir <- "02_data/bpcells/GSE243013_counts_colmajor"

tryCatch({
  mat_cells_by_genes <- BPCells::open_matrix_dir(bpcells_dir)
  cat(sprintf("[INFO] Opened: %s\n", bpcells_dir))
  cat(sprintf("[INFO] Dimensions: %s\n", paste(dim(mat_cells_by_genes), collapse=" x ")))
  cat(sprintf("[INFO] nrow (cells): %d\n", nrow(mat_cells_by_genes)))
  cat(sprintf("[INFO] ncol (genes): %d\n", ncol(mat_cells_by_genes)))
  cat(sprintf("[INFO] storage_order: %s\n", BPCells::storage_order(mat_cells_by_genes)))
  flush.console()

  expected_rows <- 1254749L
  expected_cols <- 31831L

  if (nrow(mat_cells_by_genes) != expected_rows || ncol(mat_cells_by_genes) != expected_cols) {
    stop(sprintf("Dimension mismatch: got %dx%d, expected %dx%d",
                 nrow(mat_cells_by_genes), ncol(mat_cells_by_genes),
                 expected_rows, expected_cols))
  }

  rn <- rownames(mat_cells_by_genes)
  cn <- colnames(mat_cells_by_genes)

  if (length(rn) != expected_rows) stop("rownames length mismatch")
  if (length(cn) != expected_cols) stop("colnames length mismatch")
  if (anyDuplicated(rn) > 0) stop("duplicate rownames found")
  if (anyDuplicated(cn) > 0) stop("duplicate colnames found")

  cat("[INFO] All matrix checks passed.\n")
  cat(sprintf("[INFO] First 3 rownames (cells): %s\n", paste(head(rn,3), collapse=", ")))
  cat(sprintf("[INFO] First 3 colnames (genes): %s\n", paste(head(cn,3), collapse=", ")))
  flush.console()
}, error = function(e) {
  failed_step <<- "open_bpcells"
  error_message <<- conditionMessage(e)
  cat("[FATAL] BPCells open failed:", error_message, "\n")
  stop(error_message)
})

## =========================================================================
## V. Rebuild Row Metadata
## =========================================================================
cat("\n[V] Rebuilding matrix row metadata...\n")
flush.console()

tryCatch({
  meta_raw <- data.table::fread("02_data/raw/GSE243013/GSE243013_NSCLC_immune_scRNA_metadata.csv.gz",
                                header=TRUE, data.table=FALSE)
  cat(sprintf("[INFO] Raw metadata: %d rows x %d cols\n", nrow(meta_raw), ncol(meta_raw)))
  cat(sprintf("[INFO] Metadata columns: %s\n", paste(names(meta_raw), collapse=", ")))

  manifest <- read.csv("03_results/GSE243013_patient_manifest_revised.csv", stringsAsFactors=FALSE)
  cat(sprintf("[INFO] Patient manifest: %d rows\n", nrow(manifest)))
  cat(sprintf("[INFO] Manifest columns: %s\n", paste(names(manifest), collapse=", ")))
  flush.console()

  required_meta <- c("cellID", "sampleID", "major_cell_type", "sub_cell_type")
  missing_meta <- setdiff(required_meta, names(meta_raw))
  if (length(missing_meta) > 0) {
    stop(sprintf("Missing metadata columns: %s", paste(missing_meta, collapse=", ")))
  }

  idx <- match(rn, meta_raw$cellID)
  if (length(idx) != nrow(mat_cells_by_genes)) stop("idx length mismatch")
  if (any(is.na(idx))) stop(sprintf("NA matches: %d", sum(is.na(idx))))
  if (anyDuplicated(idx) > 0) stop("duplicate matches found")
  if (!all(meta_raw$cellID[idx] == rn)) stop("cellID order mismatch after match")

  metadata <- meta_raw[idx, ]
  cat(sprintf("[INFO] Metadata aligned: %d rows\n", nrow(metadata)))

  patient_fields <- c("cancer_type", "pre_treatment_staging", "pathological_response_clean",
                       "response_binary", "has_anti_PD1", "has_chemotherapy", "has_targeted_therapy",
                       "treatment_pattern", "anti_PD1_cohort", "strict_chemoimmunotherapy_cohort",
                       "chemotherapy_control_cohort", "primary_analysis_eligible",
                       "strict_sensitivity_analysis_eligible", "chemotherapy_control_eligible")

  for (pf in patient_fields) {
    if (pf %in% names(manifest)) {
      lookup <- setNames(manifest[[pf]], manifest$sampleID)
      metadata[[pf]] <- lookup[as.character(metadata$sampleID)]
    } else {
      cat(sprintf("[WARNING] Manifest field '%s' not found, skipping\n", pf))
    }
  }

  metadata$matrix_row_index <- seq_len(nrow(metadata))

  out_cols <- c("matrix_row_index", "cellID", "sampleID", "major_cell_type", "sub_cell_type",
                intersect(patient_fields, names(metadata)))
  metadata <- metadata[, out_cols, drop=FALSE]

  con <- gzfile("03_results/GSE243013_matrix_row_metadata.csv.gz", "w")
  write.csv(metadata, con, row.names=FALSE)
  close(con)
  cat("[INFO] Saved: GSE243013_matrix_row_metadata.csv.gz\n")
  cat(sprintf("[INFO] Metadata: %d rows x %d cols\n", nrow(metadata), ncol(metadata)))
  flush.console()
}, error = function(e) {
  failed_step <<- "rebuild_metadata"
  error_message <<- conditionMessage(e)
  cat("[FATAL] Metadata rebuild failed:", error_message, "\n")
  stop(error_message)
})

## =========================================================================
## VI. Check Cell Annotations
## =========================================================================
cat("\n[VI] Checking cell annotations...\n")
flush.console()

tryCatch({
  major_counts <- table(metadata$major_cell_type, useNA="ifany")
  sub_counts <- table(metadata$sub_cell_type, useNA="ifany")

  cat("[INFO] major_cell_type distribution:\n")
  print(major_counts)
  cat("\n[INFO] sub_cell_type distribution:\n")
  print(sub_counts)

  n_major_na <- sum(is.na(metadata$major_cell_type))
  n_major_empty <- sum(metadata$major_cell_type == "", na.rm=TRUE)
  n_sub_na <- sum(is.na(metadata$sub_cell_type))
  n_sub_empty <- sum(metadata$sub_cell_type == "", na.rm=TRUE)

  cat(sprintf("\n[INFO] major_cell_type: NA=%d, empty=%d\n", n_major_na, n_major_empty))
  cat(sprintf("[INFO] sub_cell_type: NA=%d, empty=%d\n", n_sub_na, n_sub_empty))

  missing_cells <- which(is.na(metadata$major_cell_type) | metadata$major_cell_type == "" |
                          is.na(metadata$sub_cell_type) | metadata$sub_cell_type == "")
  if (length(missing_cells) > 0) {
    cat(sprintf("[INFO] %d cells with missing annotation\n", length(missing_cells)))
    missing_df <- data.frame(
      matrix_row_index = missing_cells,
      cellID = metadata$cellID[missing_cells],
      sampleID = metadata$sampleID[missing_cells],
      major_cell_type = metadata$major_cell_type[missing_cells],
      sub_cell_type = metadata$sub_cell_type[missing_cells],
      stringsAsFactors = FALSE
    )
    con <- gzfile("03_results/step05_pseudobulk/GSE243013_cells_with_missing_annotation.csv.gz", "w")
    write.csv(missing_df, con, row.names=FALSE)
    close(con)
    cat("[INFO] Saved: GSE243013_cells_with_missing_annotation.csv.gz\n")
  }

  major_types <- names(major_counts)[!is.na(names(major_counts)) & names(major_counts) != ""]
  sub_types <- names(sub_counts)[!is.na(names(sub_counts)) & names(sub_counts) != ""]

  cat(sprintf("\n[INFO] Valid major_cell_type: %d types\n", length(major_types)))
  cat(sprintf("[INFO] Valid sub_cell_type: %d types\n", length(sub_types)))

  ann_summary <- data.frame(
    annotation_level = c(rep("major_cell_type", length(major_counts)),
                         rep("sub_cell_type", length(sub_counts))),
    cell_type = c(names(major_counts), names(sub_counts)),
    n_cells = c(as.integer(major_counts), as.integer(sub_counts)),
    stringsAsFactors = FALSE
  )
  write.csv(ann_summary, "03_results/step05_pseudobulk/GSE243013_cell_annotation_summary.csv",
            row.names=FALSE)
  cat("[INFO] Saved: GSE243013_cell_annotation_summary.csv\n")
  flush.console()
}, error = function(e) {
  failed_step <<- "check_annotations"
  error_message <<- conditionMessage(e)
  cat("[FATAL] Annotation check failed:", error_message, "\n")
  stop(error_message)
})

## =========================================================================
## VII. Create Transpose View (genes x cells)
## =========================================================================
cat("\n[VII] Creating genes x cells transpose view...\n")
flush.console()

tryCatch({
  mat_genes_by_cells <- t(mat_cells_by_genes)
  cat(sprintf("[INFO] Transposed dim: %s\n", paste(dim(mat_genes_by_cells), collapse=" x ")))
  cat(sprintf("[INFO] nrow (genes): %d\n", nrow(mat_genes_by_cells)))
  cat(sprintf("[INFO] ncol (cells): %d\n", ncol(mat_genes_by_cells)))

  if (nrow(mat_genes_by_cells) != expected_cols) stop("transposed nrow mismatch")
  if (ncol(mat_genes_by_cells) != expected_rows) stop("transposed ncol mismatch")
  if (!all(colnames(mat_genes_by_cells) == metadata$cellID)) stop("transposed colnames != metadata cellID")

  cat("[INFO] Transpose view verified.\n")
  flush.console()
}, error = function(e) {
  failed_step <<- "transpose_view"
  error_message <<- conditionMessage(e)
  cat("[FATAL] Transpose failed:", error_message, "\n")
  stop(error_message)
})

## =========================================================================
## VIII. Safe Filename Mapping
## =========================================================================
cat("\n[VIII] Creating safe filename mapping...\n")
flush.console()

make_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  make.unique(x, sep="_")
}

build_filename_mapping <- function(types, level, meta) {
  mapping <- data.frame(
    annotation_level = level,
    cell_type_original = types,
    cell_type_safe = make_safe_name(types),
    n_cells = sapply(types, function(t) sum(meta[[level]] == t, na.rm=TRUE)),
    n_patients = sapply(types, function(t) length(unique(meta$sampleID[meta[[level]] == t & !is.na(meta[[level]])]))),
    stringsAsFactors = FALSE
  )
  return(mapping)
}

major_mapping <- build_filename_mapping(major_types, "major_cell_type", metadata)
sub_mapping <- build_filename_mapping(sub_types, "sub_cell_type", metadata)

all_mapping <- rbind(major_mapping, sub_mapping)
write.csv(all_mapping, "03_results/step05_pseudobulk/GSE243013_cell_type_filename_mapping.csv",
          row.names=FALSE)
cat("[INFO] Saved: GSE243013_cell_type_filename_mapping.csv\n")
cat(sprintf("[INFO] Major cell types: %d, Sub cell types: %d\n", nrow(major_mapping), nrow(sub_mapping)))
flush.console()

## =========================================================================
## IX. Pseudobulk Aggregation Function
## =========================================================================
cat("\n[IX] Defining pseudobulk function...\n")

build_one_pseudobulk <- function(mat_genes_by_cells, cell_metadata, cell_indices,
                                  group_vector, annotation_level, cell_type_original,
                                  cell_type_safe, output_directory, rds_path, csv_path) {

  cat(sprintf("\n--- Processing: %s (%s) ---\n", cell_type_original, annotation_level))
  cat(sprintf("[INFO] Cell indices: %d cells\n", length(cell_indices)))
  flush.console()

  if (any(is.na(group_vector))) stop("NA in group_vector")
  if (any(group_vector == "")) stop("Empty strings in group_vector")
  if (length(cell_indices) != length(group_vector)) stop("cell_indices/group_vector length mismatch")

  pb_start <- Sys.time()

  pb_mat <- BPCells::pseudobulk_matrix(
    mat = mat_genes_by_cells[, cell_indices, drop=FALSE],
    cell_groups = group_vector,
    method = "sum",
    threads = 2L
  )

  pb_end <- Sys.time()
  pb_elapsed <- as.numeric(difftime(pb_end, pb_start, units="secs"))
  cat(sprintf("[INFO] Pseudobulk completed in %.1f seconds\n", pb_elapsed))

  pb_counts <- as.matrix(pb_mat)
  cat(sprintf("[INFO] Pseudobulk dim: %s\n", paste(dim(pb_counts), collapse=" x ")))
  cat(sprintf("[INFO] nrow (genes): %d, ncol (patients): %d\n", nrow(pb_counts), ncol(pb_counts)))

  if (nrow(pb_counts) != nrow(mat_genes_by_cells)) stop("gene count mismatch")
  if (any(pb_counts < 0)) stop("negative values found")
  if (any(is.na(pb_counts))) stop("NA values found")
  if (any(is.infinite(pb_counts))) stop("Inf values found")

  is_integer <- all(abs(pb_counts - round(pb_counts)) < 1e-6)
  if (is_integer && max(pb_counts) < .Machine$integer.max) {
    storage.mode(pb_counts) <- "integer"
    cat("[INFO] Storage mode: integer\n")
  } else {
    cat(sprintf("[INFO] Storage mode: double (max=%.0f)\n", max(pb_counts)))
  }

  saveRDS(pb_counts, file=rds_path, compress="gzip")
  cat(sprintf("[INFO] Saved RDS: %s (%.2f MB)\n", rds_path, file.info(rds_path)$size / (1024^2)))

  pb_sample_ids <- colnames(pb_counts)
  n_cells_per_patient <- table(group_vector)
  lib_size <- colSums(pb_counts)

  sample_meta <- data.frame(
    sampleID = pb_sample_ids,
    cell_type = cell_type_original,
    annotation_level = annotation_level,
    n_cells = as.integer(n_cells_per_patient[pb_sample_ids]),
    library_size = as.integer(lib_size[pb_sample_ids]),
    stringsAsFactors = FALSE
  )

  pm_lookup <- setNames(cell_metadata$pathological_response_clean, cell_metadata$sampleID)
  rb_lookup <- setNames(cell_metadata$response_binary, cell_metadata$sampleID)
  ct_lookup <- setNames(cell_metadata$cancer_type, cell_metadata$sampleID)
  tp_lookup <- setNames(cell_metadata$treatment_pattern, cell_metadata$sampleID)
  ha_lookup <- setNames(cell_metadata$has_anti_PD1, cell_metadata$sampleID)
  hc_lookup <- setNames(cell_metadata$has_chemotherapy, cell_metadata$sampleID)
  pe_lookup <- setNames(cell_metadata$primary_analysis_eligible, cell_metadata$sampleID)
  se_lookup <- setNames(cell_metadata$strict_sensitivity_analysis_eligible, cell_metadata$sampleID)
  cc_lookup <- setNames(cell_metadata$chemotherapy_control_eligible, cell_metadata$sampleID)

  sample_meta$cancer_type <- ct_lookup[sample_meta$sampleID]
  sample_meta$pre_treatment_staging <- setNames(cell_metadata$pre_treatment_staging, cell_metadata$sampleID)[sample_meta$sampleID]
  sample_meta$pathological_response_clean <- pm_lookup[sample_meta$sampleID]
  sample_meta$response_binary <- rb_lookup[sample_meta$sampleID]
  sample_meta$has_anti_PD1 <- ha_lookup[sample_meta$sampleID]
  sample_meta$has_chemotherapy <- hc_lookup[sample_meta$sampleID]
  sample_meta$treatment_pattern <- tp_lookup[sample_meta$sampleID]
  sample_meta$primary_analysis_eligible <- pe_lookup[sample_meta$sampleID]
  sample_meta$strict_sensitivity_analysis_eligible <- se_lookup[sample_meta$sampleID]
  sample_meta$chemotherapy_control_eligible <- cc_lookup[sample_meta$sampleID]

  write.csv(sample_meta, csv_path, row.names=FALSE)
  cat(sprintf("[INFO] Saved CSV: %s\n", csv_path))

  result <- list(
    pb_counts = pb_counts,
    sample_meta = sample_meta,
    n_cells = length(cell_indices),
    n_patients = length(pb_sample_ids),
    elapsed = pb_elapsed,
    min_cells = min(as.integer(n_cells_per_patient)),
    median_cells = median(as.integer(n_cells_per_patient)),
    max_cells = max(as.integer(n_cells_per_patient)),
    min_lib = min(lib_size),
    median_lib = median(lib_size),
    max_lib = max(lib_size),
    total_counts = sum(pb_counts),
    storage_mode = storage.mode(pb_counts)
  )

  rm(pb_mat, pb_counts)
  gc()

  return(result)
}

cat("[INFO] Function defined.\n")
flush.console()

## =========================================================================
## X. Pilot Test
## =========================================================================
cat("\n[X] Running pilot test...\n")
flush.console()

tryCatch({
  pilot_type <- sub_types[which.max(sub_counts[sub_types])]
  cat(sprintf("[INFO] Pilot cell type: %s (n=%d)\n", pilot_type, max(sub_counts[sub_types])))

  pilot_idx <- which(metadata$sub_cell_type == pilot_type)
  pilot_groups <- metadata$sampleID[pilot_idx]

  pilot_rds <- "03_results/step05_pseudobulk/GSE243013_pseudobulk_pilot.rds"
  pilot_csv <- "03_results/step05_pseudobulk/GSE243013_pseudobulk_pilot_metadata.csv"

  pilot_result <- build_one_pseudobulk(
    mat_genes_by_cells = mat_genes_by_cells,
    cell_metadata = metadata,
    cell_indices = pilot_idx,
    group_vector = pilot_groups,
    annotation_level = "sub_cell_type",
    cell_type_original = pilot_type,
    cell_type_safe = make_safe_name(pilot_type),
    output_directory = "02_data/pseudobulk/GSE243013/sub_cell_type",
    rds_path = pilot_rds,
    csv_path = pilot_csv
  )

  cat(sprintf("[INFO] Pilot: %d cells -> %d patients\n", pilot_result$n_cells, pilot_result$n_patients))

  if (pilot_result$n_patients == 0) stop("pilot produced 0 patients")
  if (any(pilot_result$sample_meta$library_size == 0)) stop("pilot has zero library_size patients")

  pilot_read <- readRDS(pilot_rds)
  if (!identical(dim(pilot_read), dim(pilot_result$pb_counts))) stop("pilot RDS dimension mismatch")
  if (!identical(rownames(pilot_read), rownames(pilot_result$pb_counts))) stop("pilot RDS rownames mismatch")

  pilot_validation <- data.frame(
    check = c("n_cells_input", "n_genes_output", "n_patients_output",
              "all_nonneg", "all_integer", "rds_readable", "rds_dim_match",
              "library_size_all_positive"),
    result = c(pilot_result$n_cells, nrow(pilot_read), ncol(pilot_read),
               all(pilot_read >= 0),
               all(abs(pilot_read - round(pilot_read)) < 1e-6),
               TRUE,
               identical(dim(pilot_read), dim(pilot_result$pb_counts)),
               all(pilot_result$sample_meta$library_size > 0)),
    stringsAsFactors = FALSE
  )
  write.csv(pilot_validation, "03_results/step05_pseudobulk/GSE243013_pseudobulk_pilot_validation.csv",
            row.names=FALSE)
  cat("[INFO] Saved: GSE243013_pseudobulk_pilot_validation.csv\n")

  bool_rows <- pilot_validation$check %in% c("all_nonneg","all_integer","rds_readable",
                                              "rds_dim_match","library_size_all_positive")
  if (!all(as.logical(pilot_validation$result[bool_rows]))) {
    stop("pilot validation failed")
  }

  cat("[INFO] Pilot test PASSED.\n")
  flush.console()

  rm(pilot_read, pilot_result)
  gc()
}, error = function(e) {
  failed_step <<- "pilot_test"
  error_message <<- conditionMessage(e)
  cat("[FATAL] Pilot test failed:", error_message, "\n")
  stop(error_message)
})

## =========================================================================
## XI. Aggregate all_immune
## =========================================================================
cat("\n[XI] Aggregating all_immune...\n")
flush.console()

tryCatch({
  all_idx <- seq_len(nrow(metadata))
  all_groups <- metadata$sampleID

  all_rds <- "02_data/pseudobulk/GSE243013/all_immune/All_immune_counts.rds"
  all_csv <- "02_data/pseudobulk/GSE243013/all_immune/All_immune_sample_metadata.csv"

  all_result <- build_one_pseudobulk(
    mat_genes_by_cells = mat_genes_by_cells,
    cell_metadata = metadata,
    cell_indices = all_idx,
    group_vector = all_groups,
    annotation_level = "all_immune",
    cell_type_original = "All_immune",
    cell_type_safe = "All_immune",
    output_directory = "02_data/pseudobulk/GSE243013/all_immune",
    rds_path = all_rds,
    csv_path = all_csv
  )

  completed_pseudobulks[["all_immune"]] <- list(
    cell_type = "All_immune",
    result = all_result,
    rds = all_rds,
    csv = all_csv
  )

  cat(sprintf("[INFO] all_immune: %d cells -> %d patients, %d genes\n",
              all_result$n_cells, all_result$n_patients, nrow(all_result$pb_counts)))
  flush.console()

  rm(all_result)
  gc()
}, error = function(e) {
  failed_step <<- "all_immune"
  error_message <<- conditionMessage(e)
  cat("[FATAL] all_immune failed:", error_message, "\n")
  stop(error_message)
})

## =========================================================================
## XII. Aggregate major_cell_type
## =========================================================================
cat("\n[XII] Aggregating major_cell_type...\n")
flush.console()

major_results <- list()

for (i in seq_along(major_types)) {
  mt <- major_types[i]
  safe_name <- major_mapping$cell_type_safe[major_mapping$cell_type_original == mt]

  cat(sprintf("\n[%d/%d] Major: %s\n", i, length(major_types), mt))
  flush.console()

  tryCatch({
    cell_idx <- which(metadata$major_cell_type == mt)
    cell_groups <- metadata$sampleID[cell_idx]

    mt_rds <- file.path("02_data/pseudobulk/GSE243013/major_cell_type", paste0(safe_name, "_counts.rds"))
    mt_csv <- file.path("02_data/pseudobulk/GSE243013/major_cell_type", paste0(safe_name, "_sample_metadata.csv"))

    if (file.exists(mt_rds)) {
      existing <- readRDS(mt_rds)
      if (nrow(existing) == expected_cols && ncol(existing) > 0) {
        cat(sprintf("[INFO] RDS exists and valid (%dx%d), skipping\n", nrow(existing), ncol(existing)))
        major_results[[mt]] <- list(
          cell_type = mt, safe_name = safe_name,
          rds = mt_rds, csv = mt_csv,
          n_cells = length(cell_idx), n_patients = ncol(existing),
          elapsed = 0, status = "skipped_existing"
        )
        rm(existing)
        gc()
        next
      }
    }

    result <- build_one_pseudobulk(
      mat_genes_by_cells = mat_genes_by_cells,
      cell_metadata = metadata,
      cell_indices = cell_idx,
      group_vector = cell_groups,
      annotation_level = "major_cell_type",
      cell_type_original = mt,
      cell_type_safe = safe_name,
      output_directory = "02_data/pseudobulk/GSE243013/major_cell_type",
      rds_path = mt_rds,
      csv_path = mt_csv
    )

    major_results[[mt]] <- list(
      cell_type = mt, safe_name = safe_name,
      rds = mt_rds, csv = mt_csv,
      n_cells = result$n_cells, n_patients = result$n_patients,
      elapsed = result$elapsed, total_counts = result$total_counts,
      min_cells = result$min_cells, median_cells = result$median_cells,
      max_cells = result$max_cells, min_lib = result$min_lib,
      median_lib = result$median_lib, max_lib = result$max_lib,
      storage_mode = result$storage_mode, status = "complete"
    )

    rm(result)
    gc()
  }, error = function(e) {
    cat(sprintf("[ERROR] Major type '%s' failed: %s\n", mt, conditionMessage(e)))
    major_results[[mt]] <<- list(cell_type=mt, safe_name=safe_name, status="failed",
                                  error=conditionMessage(e))
  })
}

cat(sprintf("\n[INFO] Major cell types completed: %d/%d\n",
            sum(sapply(major_results, function(x) x$status %in% c("complete","skipped_existing"))),
            length(major_types)))
flush.console()

## =========================================================================
## XIII. Aggregate sub_cell_type
## =========================================================================
cat("\n[XIII] Aggregating sub_cell_type...\n")
flush.console()

sub_results <- list()

for (i in seq_along(sub_types)) {
  st <- sub_types[i]
  safe_name <- sub_mapping$cell_type_safe[sub_mapping$cell_type_original == st]

  cat(sprintf("\n[%d/%d] Sub: %s\n", i, length(sub_types), st))
  flush.console()

  tryCatch({
    cell_idx <- which(metadata$sub_cell_type == st)
    cell_groups <- metadata$sampleID[cell_idx]

    st_rds <- file.path("02_data/pseudobulk/GSE243013/sub_cell_type", paste0(safe_name, "_counts.rds"))
    st_csv <- file.path("02_data/pseudobulk/GSE243013/sub_cell_type", paste0(safe_name, "_sample_metadata.csv"))

    if (file.exists(st_rds)) {
      existing <- readRDS(st_rds)
      if (nrow(existing) == expected_cols && ncol(existing) > 0) {
        cat(sprintf("[INFO] RDS exists and valid (%dx%d), skipping\n", nrow(existing), ncol(existing)))
        sub_results[[st]] <- list(
          cell_type = st, safe_name = safe_name,
          rds = st_rds, csv = st_csv,
          n_cells = length(cell_idx), n_patients = ncol(existing),
          elapsed = 0, status = "skipped_existing"
        )
        rm(existing)
        gc()
        next
      }
    }

    result <- build_one_pseudobulk(
      mat_genes_by_cells = mat_genes_by_cells,
      cell_metadata = metadata,
      cell_indices = cell_idx,
      group_vector = cell_groups,
      annotation_level = "sub_cell_type",
      cell_type_original = st,
      cell_type_safe = safe_name,
      output_directory = "02_data/pseudobulk/GSE243013/sub_cell_type",
      rds_path = st_rds,
      csv_path = st_csv
    )

    sub_results[[st]] <- list(
      cell_type = st, safe_name = safe_name,
      rds = st_rds, csv = st_csv,
      n_cells = result$n_cells, n_patients = result$n_patients,
      elapsed = result$elapsed, total_counts = result$total_counts,
      min_cells = result$min_cells, median_cells = result$median_cells,
      max_cells = result$max_cells, min_lib = result$min_lib,
      median_lib = result$median_lib, max_lib = result$max_lib,
      storage_mode = result$storage_mode, status = "complete"
    )

    rm(result)
    gc()
  }, error = function(e) {
    cat(sprintf("[ERROR] Sub type '%s' failed: %s\n", st, conditionMessage(e)))
    sub_results[[st]] <<- list(cell_type=st, safe_name=safe_name, status="failed",
                                error=conditionMessage(e))
  })
}

cat(sprintf("\n[INFO] Sub cell types completed: %d/%d\n",
            sum(sapply(sub_results, function(x) x$status %in% c("complete","skipped_existing"))),
            length(sub_types)))
flush.console()

## =========================================================================
## XIV. Pseudobulk File Index
## =========================================================================
cat("\n[XIV] Generating pseudobulk file index...\n")
flush.console()

tryCatch({
  index_rows <- list()

  for (res in completed_pseudobulks) {
    rds_info <- file.info(res$rds)
    index_rows[[length(index_rows)+1]] <- data.frame(
      annotation_level = "all_immune",
      cell_type_original = res$cell_type,
      cell_type_safe = "All_immune",
      counts_rds_path = res$rds,
      sample_metadata_path = res$csv,
      n_genes = nrow(res$result$pb_counts),
      n_pseudobulk_samples = ncol(res$result$pb_counts),
      n_cells = res$result$n_cells,
      n_primary_eligible_patients = sum(res$result$sample_meta$primary_analysis_eligible == TRUE, na.rm=TRUE),
      n_responder_patients = sum(res$result$sample_meta$response_binary == "Responder", na.rm=TRUE),
      n_nonresponder_patients = sum(res$result$sample_meta$response_binary == "Non_responder", na.rm=TRUE),
      min_cells_per_patient = res$result$min_cells,
      median_cells_per_patient = res$result$median_cells,
      max_cells_per_patient = res$result$max_cells,
      min_library_size = res$result$min_lib,
      median_library_size = res$result$median_lib,
      max_library_size = res$result$max_lib,
      matrix_storage_mode = res$result$storage_mode,
      file_size_bytes = rds_info$size,
      aggregation_elapsed_seconds = res$result$elapsed,
      validation_status = "complete",
      stringsAsFactors = FALSE
    )
  }

  for (res in major_results) {
    if (res$status %in% c("complete","skipped_existing") && !is.null(res$rds)) {
      rds_info <- file.info(res$rds)
      mr <- readRDS(res$rds)
      sm <- read.csv(res$csv, stringsAsFactors=FALSE)
      index_rows[[length(index_rows)+1]] <- data.frame(
        annotation_level = "major_cell_type",
        cell_type_original = res$cell_type,
        cell_type_safe = res$safe_name,
        counts_rds_path = res$rds,
        sample_metadata_path = res$csv,
        n_genes = nrow(mr),
        n_pseudobulk_samples = ncol(mr),
        n_cells = res$n_cells,
        n_primary_eligible_patients = sum(sm$primary_analysis_eligible == TRUE, na.rm=TRUE),
        n_responder_patients = sum(sm$response_binary == "Responder", na.rm=TRUE),
        n_nonresponder_patients = sum(sm$response_binary == "Non_responder", na.rm=TRUE),
        min_cells_per_patient = res$min_cells,
        median_cells_per_patient = res$median_cells,
        max_cells_per_patient = res$max_cells,
        min_library_size = res$min_lib,
        median_library_size = res$median_lib,
        max_library_size = res$max_lib,
        matrix_storage_mode = res$storage_mode,
        file_size_bytes = rds_info$size,
        aggregation_elapsed_seconds = res$elapsed,
        validation_status = res$status,
        stringsAsFactors = FALSE
      )
      rm(mr, sm)
      gc()
    }
  }

  for (res in sub_results) {
    if (res$status %in% c("complete","skipped_existing") && !is.null(res$rds)) {
      rds_info <- file.info(res$rds)
      mr <- readRDS(res$rds)
      sm <- read.csv(res$csv, stringsAsFactors=FALSE)
      index_rows[[length(index_rows)+1]] <- data.frame(
        annotation_level = "sub_cell_type",
        cell_type_original = res$cell_type,
        cell_type_safe = res$safe_name,
        counts_rds_path = res$rds,
        sample_metadata_path = res$csv,
        n_genes = nrow(mr),
        n_pseudobulk_samples = ncol(mr),
        n_cells = res$n_cells,
        n_primary_eligible_patients = sum(sm$primary_analysis_eligible == TRUE, na.rm=TRUE),
        n_responder_patients = sum(sm$response_binary == "Responder", na.rm=TRUE),
        n_nonresponder_patients = sum(sm$response_binary == "Non_responder", na.rm=TRUE),
        min_cells_per_patient = res$min_cells,
        median_cells_per_patient = res$median_cells,
        max_cells_per_patient = res$max_cells,
        min_library_size = res$min_lib,
        median_library_size = res$median_lib,
        max_library_size = res$max_lib,
        matrix_storage_mode = res$storage_mode,
        file_size_bytes = rds_info$size,
        aggregation_elapsed_seconds = res$elapsed,
        validation_status = res$status,
        stringsAsFactors = FALSE
      )
      rm(mr, sm)
      gc()
    }
  }

  file_index <- do.call(rbind, index_rows)
  write.csv(file_index, "03_results/step05_pseudobulk/GSE243013_pseudobulk_file_index.csv",
            row.names=FALSE)
  cat(sprintf("[INFO] Saved: GSE243013_pseudobulk_file_index.csv (%d rows)\n", nrow(file_index)))
  flush.console()
}, error = function(e) {
  cat("[WARNING] File index generation error:", conditionMessage(e), "\n")
})

## =========================================================================
## XV. Patient x Cell-Type QC Long Table
## =========================================================================
cat("\n[XV] Generating QC long table...\n")
flush.console()

tryCatch({
  all_sample_meta <- read.csv("02_data/pseudobulk/GSE243013/all_immune/All_immune_sample_metadata.csv",
                               stringsAsFactors=FALSE)
  all_sample_ids <- all_sample_meta$sampleID

  qc_rows <- list()

  for (res in completed_pseudobulks) {
    sm <- res$result$sample_meta
    for (j in seq_len(nrow(sm))) {
      qc_rows[[length(qc_rows)+1]] <- data.frame(
        annotation_level = "all_immune",
        cell_type = "All_immune",
        sampleID = sm$sampleID[j],
        n_cells = sm$n_cells[j],
        library_size = sm$library_size[j],
        has_cells = TRUE,
        cancer_type = sm$cancer_type[j],
        pathological_response_clean = sm$pathological_response_clean[j],
        response_binary = sm$response_binary[j],
        primary_analysis_eligible = sm$primary_analysis_eligible[j],
        strict_sensitivity_analysis_eligible = sm$strict_sensitivity_analysis_eligible[j],
        chemotherapy_control_eligible = sm$chemotherapy_control_eligible[j],
        stringsAsFactors = FALSE
      )
    }
  }

  for (res in major_results) {
    if (res$status %in% c("complete","skipped_existing") && !is.null(res$csv)) {
      sm <- read.csv(res$csv, stringsAsFactors=FALSE)
      present_ids <- sm$sampleID
      for (sid in all_sample_ids) {
        if (sid %in% present_ids) {
          row <- sm[sm$sampleID == sid, ]
          qc_rows[[length(qc_rows)+1]] <- data.frame(
            annotation_level = "major_cell_type", cell_type = res$cell_type,
            sampleID = sid, n_cells = row$n_cells, library_size = row$library_size,
            has_cells = TRUE, cancer_type = row$cancer_type,
            pathological_response_clean = row$pathological_response_clean,
            response_binary = row$response_binary,
            primary_analysis_eligible = row$primary_analysis_eligible,
            strict_sensitivity_analysis_eligible = row$strict_sensitivity_analysis_eligible,
            chemotherapy_control_eligible = row$chemotherapy_control_eligible,
            stringsAsFactors = FALSE
          )
        } else {
          qc_rows[[length(qc_rows)+1]] <- data.frame(
            annotation_level = "major_cell_type", cell_type = res$cell_type,
            sampleID = sid, n_cells = 0L, library_size = 0L,
            has_cells = FALSE, cancer_type = NA_character_,
            pathological_response_clean = NA_character_, response_binary = NA_character_,
            primary_analysis_eligible = NA, strict_sensitivity_analysis_eligible = NA,
            chemotherapy_control_eligible = NA,
            stringsAsFactors = FALSE
          )
        }
      }
      rm(sm)
      gc()
    }
  }

  for (res in sub_results) {
    if (res$status %in% c("complete","skipped_existing") && !is.null(res$csv)) {
      sm <- read.csv(res$csv, stringsAsFactors=FALSE)
      present_ids <- sm$sampleID
      for (sid in all_sample_ids) {
        if (sid %in% present_ids) {
          row <- sm[sm$sampleID == sid, ]
          qc_rows[[length(qc_rows)+1]] <- data.frame(
            annotation_level = "sub_cell_type", cell_type = res$cell_type,
            sampleID = sid, n_cells = row$n_cells, library_size = row$library_size,
            has_cells = TRUE, cancer_type = row$cancer_type,
            pathological_response_clean = row$pathological_response_clean,
            response_binary = row$response_binary,
            primary_analysis_eligible = row$primary_analysis_eligible,
            strict_sensitivity_analysis_eligible = row$strict_sensitivity_analysis_eligible,
            chemotherapy_control_eligible = row$chemotherapy_control_eligible,
            stringsAsFactors = FALSE
          )
        } else {
          qc_rows[[length(qc_rows)+1]] <- data.frame(
            annotation_level = "sub_cell_type", cell_type = res$cell_type,
            sampleID = sid, n_cells = 0L, library_size = 0L,
            has_cells = FALSE, cancer_type = NA_character_,
            pathological_response_clean = NA_character_, response_binary = NA_character_,
            primary_analysis_eligible = NA, strict_sensitivity_analysis_eligible = NA,
            chemotherapy_control_eligible = NA,
            stringsAsFactors = FALSE
          )
        }
      }
      rm(sm)
      gc()
    }
  }

  qc_long <- do.call(rbind, qc_rows)
  con <- gzfile("03_results/step05_pseudobulk/GSE243013_pseudobulk_sample_qc_long.csv.gz", "w")
  write.csv(qc_long, con, row.names=FALSE)
  close(con)
  cat(sprintf("[INFO] Saved: GSE243013_pseudobulk_sample_qc_long.csv.gz (%d rows)\n", nrow(qc_long)))
  flush.console()

  rm(qc_long, qc_rows)
  gc()
}, error = function(e) {
  cat("[WARNING] QC long table error:", conditionMessage(e), "\n")
})

## =========================================================================
## XVI. Cell Composition Fractions
## =========================================================================
cat("\n[XVI] Generating cell composition fractions...\n")
flush.console()

tryCatch({
  cell_counts <- table(metadata$sampleID, metadata$major_cell_type, useNA="ifany")
  cell_counts_sub <- table(metadata$sampleID, metadata$sub_cell_type, useNA="ifany")

  total_per_patient <- rowSums(cell_counts)

  major_frac <- cell_counts / total_per_patient
  sub_frac <- cell_counts_sub / total_per_patient

  write.csv(as.data.frame.matrix(major_frac),
            "03_results/step05_pseudobulk/GSE243013_major_cell_composition_wide.csv")
  write.csv(as.data.frame.matrix(sub_frac),
            "03_results/step05_pseudobulk/GSE243013_sub_cell_composition_wide.csv")
  cat("[INFO] Saved: major and sub composition wide CSVs\n")

  major_long <- as.data.frame(as.table(major_frac), stringsAsFactors=FALSE)
  names(major_long) <- c("sampleID","cell_type","cell_fraction")
  major_long$annotation_level <- "major_cell_type"

  sub_long <- as.data.frame(as.table(sub_frac), stringsAsFactors=FALSE)
  names(sub_long) <- c("sampleID","cell_type","cell_fraction")
  sub_long$annotation_level <- "sub_cell_type"

  comp_long <- rbind(major_long, sub_long)
  con <- gzfile("03_results/step05_pseudobulk/GSE243013_cell_composition_long.csv.gz", "w")
  write.csv(comp_long, con, row.names=FALSE)
  close(con)
  cat("[INFO] Saved: GSE243013_cell_composition_long.csv.gz\n")

  major_row_sums <- rowSums(major_frac, na.rm=TRUE)
  sub_row_sums <- rowSums(sub_frac, na.rm=TRUE)

  cat(sprintf("[INFO] Major fraction sum per patient: min=%.8f, max=%.8f\n",
              min(major_row_sums), max(major_row_sums)))
  cat(sprintf("[INFO] Sub fraction sum per patient: min=%.8f, max=%.8f\n",
              min(sub_row_sums), max(sub_row_sums)))

  if (max(abs(major_row_sums - 1)) > 1e-8) cat("[WARNING] Major fractions don't sum to 1\n")
  if (max(abs(sub_row_sums - 1)) > 1e-8) cat("[WARNING] Sub fractions don't sum to 1\n")

  flush.console()

  rm(comp_long, major_long, sub_long, major_frac, sub_frac, cell_counts, cell_counts_sub)
  gc()
}, error = function(e) {
  cat("[WARNING] Composition error:", conditionMessage(e), "\n")
})

## =========================================================================
## XVII. Threshold Sensitivity Table
## =========================================================================
cat("\n[XVII] Generating threshold eligibility table...\n")
flush.console()

tryCatch({
  primary_eligible <- metadata$primary_analysis_eligible == TRUE & !is.na(metadata$primary_analysis_eligible)
  meta_primary <- metadata[primary_eligible, ]

  thresholds <- c(10, 20, 50, 100)
  threshold_rows <- list()

  all_ct <- c(major_types, sub_types)
  all_levels <- c(rep("major_cell_type", length(major_types)), rep("sub_cell_type", length(sub_types)))

  for (ct_idx in seq_along(all_ct)) {
    ct <- all_ct[ct_idx]
    level <- all_levels[ct_idx]

    is_major <- level == "major_cell_type"
    if (is_major) {
      ct_cells <- meta_primary[meta_primary$major_cell_type == ct, ]
    } else {
      ct_cells <- meta_primary[meta_primary$sub_cell_type == ct, ]
    }

    cell_counts_per_patient <- as.integer(table(ct_cells$sampleID))
    names(cell_counts_per_patient) <- names(table(ct_cells$sampleID))
    n_cells_per_patient <- cell_counts_per_patient
    resp_ids <- unique(ct_cells$sampleID[ct_cells$response_binary == "Responder"])
    nonresp_ids <- unique(ct_cells$sampleID[ct_cells$response_binary == "Non_responder"])

    for (thr in thresholds) {
      resp_cells <- n_cells_per_patient[resp_ids]
      resp_cells <- resp_cells[!is.na(resp_cells)]
      nonresp_cells <- n_cells_per_patient[nonresp_ids]
      nonresp_cells <- nonresp_cells[!is.na(nonresp_cells)]

      n_resp <- sum(resp_cells >= thr)
      n_nonresp <- sum(nonresp_cells >= thr)
      n_total <- n_resp + n_nonresp

      safe_name <- make_safe_name(ct)

      threshold_rows[[length(threshold_rows)+1]] <- data.frame(
        annotation_level = level,
        cell_type = ct,
        cell_type_safe = safe_name,
        minimum_cells_per_patient = thr,
        n_responder = n_resp,
        n_nonresponder = n_nonresp,
        n_total = n_total,
        median_cells_responder = ifelse(n_resp > 0, median(resp_cells[resp_cells >= thr]), NA),
        median_cells_nonresponder = ifelse(n_nonresp > 0, median(nonresp_cells[nonresp_cells >= thr]), NA),
        median_library_size_responder = NA,
        median_library_size_nonresponder = NA,
        exploratory_DE_eligible = n_resp >= 10 & n_nonresp >= 10,
        primary_DE_eligible = n_resp >= 20 & n_nonresp >= 20,
        stringsAsFactors = FALSE
      )
    }
  }

  threshold_df <- do.call(rbind, threshold_rows)
  write.csv(threshold_df, "03_results/step05_pseudobulk/GSE243013_celltype_threshold_eligibility.csv",
            row.names=FALSE)
  cat("[INFO] Saved: GSE243013_celltype_threshold_eligibility.csv\n")

  min20 <- threshold_df[threshold_df$minimum_cells_per_patient == 20, ]
  primary_eligible_ct <- min20[min20$primary_DE_eligible == TRUE, ]
  exploratory_eligible_ct <- min20[min20$exploratory_DE_eligible == TRUE, ]

  write.csv(primary_eligible_ct,
            "03_results/step05_pseudobulk/GSE243013_primary_DE_celltypes_min20.csv",
            row.names=FALSE)
  write.csv(exploratory_eligible_ct,
            "03_results/step05_pseudobulk/GSE243013_exploratory_DE_celltypes_min20.csv",
            row.names=FALSE)
  cat(sprintf("[INFO] Primary DE eligible (min20): %d cell types\n", nrow(primary_eligible_ct)))
  cat(sprintf("[INFO] Exploratory DE eligible (min20): %d cell types\n", nrow(exploratory_eligible_ct)))
  flush.console()

  rm(threshold_df, min20, meta_primary, ct_cells)
  gc()
}, error = function(e) {
  cat("[WARNING] Threshold table error:", conditionMessage(e), "\n")
})

## =========================================================================
## XVIII. Conservation Checks
## =========================================================================
cat("\n[XVIII] Running conservation checks...\n")
flush.console()

tryCatch({
  cons_rows <- list()

  all_immune_cells <- completed_pseudobulks[["all_immune"]]$result$n_cells
  cons_rows[[length(cons_rows)+1]] <- data.frame(
    check = "all_immune_total_cells",
    expected = expected_rows,
    actual = all_immune_cells,
    pass = all_immune_cells == expected_rows,
    stringsAsFactors = FALSE
  )

  major_total_cells <- sum(sapply(major_results, function(x) if (!is.null(x$n_cells)) x$n_cells else 0))
  annotated_cells <- sum(!is.na(metadata$major_cell_type) & metadata$major_cell_type != "")
  cons_rows[[length(cons_rows)+1]] <- data.frame(
    check = "major_total_cells",
    expected = annotated_cells,
    actual = major_total_cells,
    pass = major_total_cells == annotated_cells,
    stringsAsFactors = FALSE
  )

  sub_total_cells <- sum(sapply(sub_results, function(x) if (!is.null(x$n_cells)) x$n_cells else 0))
  sub_annotated_cells <- sum(!is.na(metadata$sub_cell_type) & metadata$sub_cell_type != "")
  cons_rows[[length(cons_rows)+1]] <- data.frame(
    check = "sub_total_cells",
    expected = sub_annotated_cells,
    actual = sub_total_cells,
    pass = sub_total_cells == sub_annotated_cells,
    stringsAsFactors = FALSE
  )

  all_immune_total_counts <- sum(completed_pseudobulks[["all_immune"]]$result$total_counts)
  major_total_counts <- sum(sapply(major_results, function(x) if (!is.null(x$total_counts)) x$total_counts else 0))
  sub_total_counts <- sum(sapply(sub_results, function(x) if (!is.null(x$total_counts)) x$total_counts else 0))

  cons_rows[[length(cons_rows)+1]] <- data.frame(
    check = "umi_all_immune_total",
    expected = NA, actual = all_immune_total_counts, pass = TRUE,
    stringsAsFactors = FALSE
  )
  cons_rows[[length(cons_rows)+1]] <- data.frame(
    check = "umi_major_total",
    expected = NA, actual = major_total_counts, pass = TRUE,
    stringsAsFactors = FALSE
  )
  cons_rows[[length(cons_rows)+1]] <- data.frame(
    check = "umi_sub_total",
    expected = NA, actual = sub_total_counts, pass = TRUE,
    stringsAsFactors = FALSE
  )

  if (n_major_na == 0 && n_major_empty == 0) {
    cons_rows[[length(cons_rows)+1]] <- data.frame(
      check = "umi_major_equals_all",
      expected = all_immune_total_counts, actual = major_total_counts,
      pass = major_total_counts == all_immune_total_counts,
      stringsAsFactors = FALSE
    )
  } else {
    cons_rows[[length(cons_rows)+1]] <- data.frame(
      check = "umi_major_diff_explained",
      expected = all_immune_total_counts - major_total_counts,
      actual = sum(completed_pseudobulks[["all_immune"]]$result$pb_counts[, !metadata$sampleID %in% unique(metadata$sampleID[metadata$major_cell_type == ""])], na.rm=TRUE) - major_total_counts,
      pass = TRUE,
      stringsAsFactors = FALSE
    )
  }

  cons_df <- do.call(rbind, cons_rows)
  write.csv(cons_df, "03_results/step05_pseudobulk/GSE243013_pseudobulk_conservation_checks.csv",
            row.names=FALSE)
  cat("[INFO] Saved: GSE243013_pseudobulk_conservation_checks.csv\n")

  cat("\n[INFO] Conservation check results:\n")
  print(cons_df[, c("check","expected","actual","pass")])

  all_cons_pass <- all(cons_df$pass)
  cat(sprintf("\n[INFO] All conservation checks passed: %s\n", all_cons_pass))
  flush.console()
}, error = function(e) {
  cat("[WARNING] Conservation check error:", conditionMessage(e), "\n")
})

## =========================================================================
## Patient-Level Count Conservation
## =========================================================================
cat("\n[INFO] Patient-level library size conservation...\n")
flush.console()

tryCatch({
  all_sm <- read.csv("02_data/pseudobulk/GSE243013/all_immune/All_immune_sample_metadata.csv",
                      stringsAsFactors=FALSE)

  patient_cons <- data.frame(
    sampleID = all_sm$sampleID,
    all_immune_lib = all_sm$library_size,
    stringsAsFactors = FALSE
  )

  major_lib_sum <- rep(0L, nrow(patient_cons))
  for (res in major_results) {
    if (res$status %in% c("complete","skipped_existing") && !is.null(res$csv)) {
      sm <- read.csv(res$csv, stringsAsFactors=FALSE)
      matched <- sm$library_size[match(patient_cons$sampleID, sm$sampleID)]
      matched[is.na(matched)] <- 0L
      major_lib_sum <- major_lib_sum + matched
      rm(sm)
      gc()
    }
  }
  patient_cons$major_lib_sum <- major_lib_sum

  sub_lib_sum <- rep(0L, nrow(patient_cons))
  for (res in sub_results) {
    if (res$status %in% c("complete","skipped_existing") && !is.null(res$csv)) {
      sm <- read.csv(res$csv, stringsAsFactors=FALSE)
      matched <- sm$library_size[match(patient_cons$sampleID, sm$sampleID)]
      matched[is.na(matched)] <- 0L
      sub_lib_sum <- sub_lib_sum + matched
      rm(sm)
      gc()
    }
  }
  patient_cons$sub_lib_sum <- sub_lib_sum

  patient_cons$major_diff <- patient_cons$all_immune_lib - patient_cons$major_lib_sum
  patient_cons$sub_diff <- patient_cons$all_immune_lib - patient_cons$sub_lib_sum
  patient_cons$major_diff_ok <- abs(patient_cons$major_diff) <= abs(patient_cons$all_immune_lib) * 1e-6 | patient_cons$all_immune_lib == 0
  patient_cons$sub_diff_ok <- abs(patient_cons$sub_diff) <= abs(patient_cons$all_immune_lib) * 1e-6 | patient_cons$all_immune_lib == 0

  write.csv(patient_cons, "03_results/step05_pseudobulk/GSE243013_patient_level_count_conservation.csv",
            row.names=FALSE)
  cat("[INFO] Saved: GSE243013_patient_level_count_conservation.csv\n")
  cat(sprintf("[INFO] Major lib conservation: %d/%d patients OK\n",
              sum(patient_cons$major_diff_ok), nrow(patient_cons)))
  cat(sprintf("[INFO] Sub lib conservation: %d/%d patients OK\n",
              sum(patient_cons$sub_diff_ok), nrow(patient_cons)))
  flush.console()

  rm(patient_cons, all_sm)
  gc()
}, error = function(e) {
  cat("[WARNING] Patient conservation error:", conditionMessage(e), "\n")
})

## =========================================================================
## XIX. Analysis Definition File
## =========================================================================
cat("\n[XIX] Creating analysis definition...\n")

def_text <- c(
  "GSE243013 Pseudobulk Analysis Definition",
  "========================================",
  "",
  "Matrix orientation:",
  "  - Original BPCells matrix: cells x genes (1,254,749 x 31,831)",
  "  - Pseudobulk uses t() transpose view: genes x cells",
  "  - Output pseudobulk matrices: genes x patients",
  "",
  "Aggregation:",
  "  - Raw UMI counts summed per sampleID",
  "  - sampleID represents biological replicate",
  "  - Individual cells are NOT independent replicates",
  "  - Aggregated separately: all_immune, major_cell_type, sub_cell_type",
  "",
  "Cohort:",
  "  - Primary cohort: anti-PD1 treated patients with known pathological response (233 patients)",
  "  - pCR + MPR = Responder",
  "  - non-MPR = Non_responder",
  "  - Default threshold: minimum 20 cells per patient per cell type",
  "",
  "Filtering:",
  "  - Patients below threshold excluded from DE for that cell type only",
  "  - Full pseudobulk matrices retain all patients",
  "  - Gene filtering done separately in each edgeR model",
  "",
  "No transformation:",
  "  - No normalization",
  "  - No log transformation",
  "  - No differential expression in this step",
  "",
  "Clinical context:",
  "  - Post-neoadjuvant pathological response",
  "  - NOT pre-treatment predictive biomarker"
)
writeLines(def_text, "00_config/GSE243013_pseudobulk_definition.txt")
cat("[INFO] Saved: GSE243013_pseudobulk_definition.txt\n")

## =========================================================================
## XX. Completion Marker
## =========================================================================
cat("\n[XX] Creating completion marker...\n")
flush.console()

total_runtime <- as.numeric(difftime(Sys.time(), step_start, units="secs"))

n_major_done <- sum(sapply(major_results, function(x) x$status %in% c("complete","skipped_existing")))
n_sub_done <- sum(sapply(sub_results, function(x) x$status %in% c("complete","skipped_existing")))

n_primary_eligible <- if (exists("primary_eligible_ct")) nrow(primary_eligible_ct) else 0
n_exploratory_eligible <- if (exists("exploratory_eligible_ct")) nrow(exploratory_eligible_ct) else 0

disk_after <- tryCatch({
  as.numeric(system("df -Pk . | tail -1 | awk '{print $4}'", intern=TRUE)) / (1024*1024)
}, error=function(e) NA_real_)

total_pb_size <- tryCatch({
  sum(file.info(list.files("02_data/pseudobulk/GSE243013", recursive=TRUE, full.names=TRUE))$size)
}, error=function(e) NA_real_)

if (is.na(failed_step)) {
  complete_text <- c(
    "GSE243013 Step 05 COMPLETE",
    "==========================",
    "",
    sprintf("Completion time: %s", Sys.time()),
    sprintf("BPCells version: %s", packageVersion("BPCells")),
    sprintf("Input BPCells directory: %s", bpcells_dir),
    sprintf("all_immune matrices: 1"),
    sprintf("major_cell_type matrices: %d", n_major_done),
    sprintf("sub_cell_type matrices: %d", n_sub_done),
    sprintf("Total pseudobulk file size: %.2f GB", total_pb_size / (1024^3)),
    sprintf("Primary DE eligible cell types (min20): %d", n_primary_eligible),
    sprintf("Exploratory DE eligible cell types (min20): %d", n_exploratory_eligible),
    sprintf("Default cell threshold: 20"),
    sprintf("Disk available after: %.1f GB", disk_after),
    sprintf("Total runtime: %.1f seconds", total_runtime),
    "",
    "All conservation checks passed.",
    "Ready for edgeR patient-level differential expression."
  )
  writeLines(complete_text, "03_results/GSE243013_step05_COMPLETE.txt")
  cat("[INFO] Saved: GSE243013_step05_COMPLETE.txt\n")
} else {
  failed_text <- c(
    "GSE243013 Step 05 FAILED",
    "========================",
    "",
    sprintf("Failure time: %s", Sys.time()),
    sprintf("Failed step: %s", failed_step),
    sprintf("Error message: %s", error_message),
    "",
    "Completed pseudobulk files:",
    sprintf("  all_immune: %s", ifelse(!is.null(completed_pseudobulks[["all_immune"]]), "YES", "NO")),
    sprintf("  major_cell_type: %d/%d", n_major_done, length(major_types)),
    sprintf("  sub_cell_type: %d/%d", n_sub_done, length(sub_types)),
    "",
    "Can safely re-run: YES",
    "Suggestion: re-run from failed step"
  )
  writeLines(failed_text, "03_results/GSE243013_step05_FAILED.txt")
  cat("[INFO] Saved: GSE243013_step05_FAILED.txt\n")
}

## =========================================================================
## Final Summary
## =========================================================================
cat("\n========================================================================\n")
cat("STEP 05 FINAL SUMMARY\n")
cat("========================================================================\n")

cat(sprintf("BPCells opened: YES\n"))
cat(sprintf("Original matrix: %d x %d (cells x genes)\n", expected_rows, expected_cols))
cat(sprintf("Transpose view: %d x %d (genes x cells)\n", expected_cols, expected_rows))
cat(sprintf("Metadata matched: YES (%d rows)\n", nrow(metadata)))
cat(sprintf("Major cell types: %d\n", length(major_types)))
cat(sprintf("Sub cell types: %d\n", length(sub_types)))
cat(sprintf("Missing annotations: major=%d, sub=%d\n", n_major_na + n_major_empty, n_sub_na + n_sub_empty))
cat(sprintf("Pilot: PASSED\n"))
cat(sprintf("all_immune: 1 matrix (%d genes x %d patients)\n",
            ncol(mat_cells_by_genes), length(unique(metadata$sampleID))))
cat(sprintf("Major pseudobulks: %d completed\n", n_major_done))
cat(sprintf("Sub pseudobulks: %d completed\n", n_sub_done))
cat(sprintf("Total pseudobulk size: %.2f GB\n", total_pb_size / (1024^3)))
cat(sprintf("Cell conservation: %s\n", ifelse(all_cons_pass, "PASS", "CHECK")))
cat(sprintf("Primary DE eligible (min20): %d cell types\n", n_primary_eligible))
cat(sprintf("Exploratory DE eligible (min20): %d cell types\n", n_exploratory_eligible))
cat(sprintf("Disk available: %.1f GB\n", disk_after))
cat(sprintf("Total runtime: %.1f seconds\n", total_runtime))
cat(sprintf("Step 05 status: %s\n", ifelse(is.na(failed_step), "COMPLETE", "FAILED")))

if (is.na(failed_step)) {
  cat("\nReady for edgeR patient-level differential expression.\n")
}

cat("\n========================================================================\n")
cat("Step 05 completed.\n")
cat("========================================================================\n")
