## =========================================================================
## Step 01: Download & Inspect GSE243013 Metadata
## =========================================================================

options(timeout = 3600)

cat("========================================================================\n")
cat("Step 01: Download & Inspect GSE243013 Metadata\n")
cat("========================================================================\n\n")

RAW_DIR   <- "02_data/raw/GSE243013"
MANIFEST  <- "02_data/manifest/GSE243013_supplementary_files.csv"
RESULTS   <- "03_results"

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)

BASE_URL <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE243013&format=file&file="

ALLOWED_FILES <- c(
  "GSE243013_NMF_all_group_5.csv.gz",
  "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz",
  "GSE243013_T_with_TCR_annotation.csv.gz",
  "GSE243013_UMAP_info.tar.gz",
  "GSE243013_barcodes.csv.gz",
  "GSE243013_genes.csv.gz"
)

FORBIDDEN_FILES <- c(
  "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz",
  "GSE243013_RAW.tar"
)

## =========================================================================
## I. Download Files
## =========================================================================
cat("[I] Downloading allowed supplementary files...\n\n")

download_status <- data.frame(
  file_name            = character(),
  url                  = character(),
  local_path           = character(),
  expected_or_remote_size = character(),
  downloaded_size      = character(),
  status               = character(),
  error_message        = character(),
  stringsAsFactors     = FALSE
)

read_manifest_url <- function(fname) {
  tryCatch({
    if (file.exists(MANIFEST)) {
      m <- read.csv(MANIFEST, stringsAsFactors = FALSE)
      row <- m[m$file_name == fname, ]
      if (nrow(row) > 0 && !is.na(row$url) && nchar(row$url) > 0) {
        return(list(url = row$url, size = row$size_bytes))
      }
    }
    return(NULL)
  }, error = function(e) return(NULL))
}

for (fname in ALLOWED_FILES) {
  cat(sprintf("  Processing: %s\n", fname))
  local_path <- file.path(RAW_DIR, fname)

  ## Check if already downloaded
  if (file.exists(local_path) && file.info(local_path)$size > 0) {
    existing_size <- file.info(local_path)$size
    cat(sprintf("    Already exists (%d bytes). Skipping download.\n", existing_size))
    download_status <- rbind(download_status, data.frame(
      file_name            = fname,
      url                  = "(already present)",
      local_path           = local_path,
      expected_or_remote_size = as.character(existing_size),
      downloaded_size      = as.character(existing_size),
      status               = "skipped_exists",
      error_message        = "",
      stringsAsFactors     = FALSE
    ))
    next
  }

  ## Determine URL — always use HTTPS download API (FTP URLs return 403)
  manifest_info <- read_manifest_url(fname)
  dl_url <- paste0(BASE_URL, utils::URLencode(fname, reserved = TRUE))
  expected_size <- ifelse(!is.null(manifest_info), manifest_info$size, NA_character_)

  ## Download with tryCatch
  dl_ok <- FALSE
  dl_err <- ""
  tryCatch({
    download.file(
      url      = dl_url,
      destfile = local_path,
      method   = "libcurl",
      mode     = "wb",
      quiet    = FALSE
    )
    dl_ok <- TRUE
  }, error = function(e) {
    dl_err <<- conditionMessage(e)
    cat(sprintf("    [ERROR] Download failed: %s\n", dl_err))
  })

  ## Verify download
  dl_size <- NA_character_
  status  <- "failed"
  if (dl_ok && file.exists(local_path)) {
    dl_size <- as.character(file.info(local_path)$size)
    if (as.numeric(dl_size) > 0) {
      status <- "success"
    } else {
      status <- "empty_file"
      dl_err <- "Downloaded file is 0 bytes"
    }
  }

  download_status <- rbind(download_status, data.frame(
    file_name            = fname,
    url                  = dl_url,
    local_path           = local_path,
    expected_or_remote_size = ifelse(is.na(expected_size), "", expected_size),
    downloaded_size      = ifelse(is.na(dl_size), "", dl_size),
    status               = status,
    error_message        = dl_err,
    stringsAsFactors     = FALSE
  ))
  cat(sprintf("    Status: %s\n\n", status))
}

## Check for forbidden files (should NOT exist)
cat("[I] Verifying forbidden files are NOT present...\n")
for (fname in FORBIDDEN_FILES) {
  fp <- file.path(RAW_DIR, fname)
  if (file.exists(fp)) {
    cat(sprintf("  [WARNING] Forbidden file exists: %s (%d bytes). NOT deleting.\n",
                fp, file.info(fp)$size))
  } else {
    cat(sprintf("  [OK] %s not present (correct).\n", fname))
  }
}

## Save download status
cat("\n[I] Saving download status...\n")
tryCatch({
  write.csv(download_status, file.path(RESULTS, "GSE243013_download_status.csv"),
            row.names = FALSE)
  cat("[I] Saved: 03_results/GSE243013_download_status.csv\n")
}, error = function(e) {
  cat("[ERROR] Failed to save download_status:", conditionMessage(e), "\n")
})

## =========================================================================
## Verify downloaded files
## =========================================================================
cat("\n[II] Verifying downloaded files...\n")
for (fname in ALLOWED_FILES) {
  fp <- file.path(RAW_DIR, fname)
  cat(sprintf("  %s: ", fname))
  if (!file.exists(fp)) {
    cat("NOT FOUND\n")
    next
  }
  fsize <- file.info(fp)$size
  if (fsize == 0) { cat("EMPTY\n"); next }

  if (grepl("\\.csv\\.gz$", fname)) {
    tryCatch({
      con <- gzfile(fp, "r")
      header <- readLines(con, n = 1)
      close(con)
      cat(sprintf("OK (%d bytes, CSV.GZ readable, header: %s)\n", fsize,
                  substr(header, 1, 80)))
    }, error = function(e) {
      cat(sprintf("CSV.GZ READ ERROR: %s\n", conditionMessage(e)))
    })
  } else if (grepl("UMAP_info\\.tar\\.gz$", fname)) {
    tryCatch({
      contents <- untar(fp, list = TRUE)
      cat(sprintf("OK (%d bytes, %d files inside)\n", fsize, length(contents)))
      cat(sprintf("    Contents: %s\n", paste(head(contents, 10), collapse = ", ")))
    }, error = function(e) {
      cat(sprintf("TAR READ ERROR: %s\n", conditionMessage(e)))
    })
  } else {
    cat(sprintf("OK (%d bytes)\n", fsize))
  }
}

## =========================================================================
## II. Read Main Metadata
## =========================================================================
cat("\n[III] Reading main metadata...\n")

META_PATH <- file.path(RAW_DIR, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")
meta <- NULL

## Try data.table::fread first
tryCatch({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cat("[INFO] data.table not installed. Installing from CRAN...\n")
    install.packages("data.table", repos = "https://cloud.r-project.org")
  }
  meta <- data.table::fread(META_PATH, header = TRUE, data.table = FALSE)
  cat(sprintf("[INFO] Read with data.table::fread: %d rows x %d columns\n",
              nrow(meta), ncol(meta)))
}, error = function(e) {
  cat("[WARNING] data.table::fread failed:", conditionMessage(e), "\n")
  cat("[INFO] Falling back to read.csv...\n")
  tryCatch({
    meta <- read.csv(gzfile(META_PATH), check.names = FALSE)
    cat(sprintf("[INFO] Read with read.csv: %d rows x %d columns\n",
                nrow(meta), ncol(meta)))
  }, error = function(e2) {
    cat("[ERROR] read.csv also failed:", conditionMessage(e2), "\n")
  })
})

if (is.null(meta)) {
  cat("[FATAL] Cannot read metadata. Stopping.\n")
  quit(save = "no", status = 1)
}

## =========================================================================
## III. Column Profile
## =========================================================================
cat("\n[IV] Generating column profile...\n")

col_profile <- data.frame(
  column_name          = character(),
  column_class         = character(),
  number_of_rows       = integer(),
  number_of_unique_values = integer(),
  number_of_missing_values = integer(),
  missing_percentage   = numeric(),
  top10_examples       = character(),
  stringsAsFactors     = FALSE
)

for (cn in names(meta)) {
  col <- meta[[cn]]
  n_miss <- sum(is.na(col))
  n_unique <- length(unique(col[!is.na(col)]))
  top10 <- head(unique(col[!is.na(col)]), 10)
  top10_str <- paste(top10, collapse = "; ")
  col_profile <- rbind(col_profile, data.frame(
    column_name          = cn,
    column_class         = class(col)[1],
    number_of_rows       = nrow(meta),
    number_of_unique_values = n_unique,
    number_of_missing_values = n_miss,
    missing_percentage   = round(100 * n_miss / nrow(meta), 2),
    top10_examples       = top10_str,
    stringsAsFactors     = FALSE
  ))
}

tryCatch({
  write.csv(col_profile, file.path(RESULTS, "GSE243013_metadata_column_profile.csv"),
            row.names = FALSE)
  cat("[INFO] Saved: GSE243013_metadata_column_profile.csv\n")
}, error = function(e) {
  cat("[ERROR] Failed to save column profile:", conditionMessage(e), "\n")
})

## Dimensions summary
tryCatch({
  col_names <- paste(names(meta), collapse = ", ")
  obj_size <- format(object.size(meta), units = "MB")
  dup_rows <- sum(duplicated(meta))
  dims_text <- c(
    sprintf("Total rows: %d", nrow(meta)),
    sprintf("Total columns: %d", ncol(meta)),
    sprintf("Column names: %s", col_names),
    sprintf("Memory usage: %s", obj_size),
    sprintf("Duplicate rows: %d", dup_rows)
  )
  writeLines(dims_text, file.path(RESULTS, "GSE243013_metadata_dimensions.txt"))
  cat("[INFO] Saved: GSE243013_metadata_dimensions.txt\n")
  cat(paste(dims_text, collapse = "\n"), "\n\n")
}, error = function(e) {
  cat("[ERROR] Failed to save dimensions:", conditionMessage(e), "\n")
})

## =========================================================================
## IV. Candidate Column Selection
## =========================================================================
cat("\n[V] Selecting candidate columns...\n")

keyword_groups <- list(
  patient        = c("patient", "patient_id", "subject", "subject_id", "donor", "case", "individual"),
  sample         = c("sample", "sample_id", "specimen", "tissue"),
  response       = c("response", "pathological", "pathologic", "MPR", "pCR", "NMPR",
                     "non-MPR", "TRG", "residual", "regression"),
  treatment      = c("treatment", "therapy", "regimen", "chemo", "immunotherapy",
                     "neoadjuvant", "ICI", "PD.1", "PD.L1", "PD-1", "PD-L1"),
  histology      = c("histology", "subtype", "LUAD", "LUSC", "adenocarcinoma", "squamous"),
  cell_type      = c("cell_type", "celltype", "cell.type", "annotation", "cluster",
                     "major", "minor", "lineage"),
  barcode        = c("barcode", "cell", "cell_id")
)

candidate_df <- data.frame(
  category            = character(),
  column_name         = character(),
  column_class        = character(),
  unique_count        = integer(),
  missing_count       = integer(),
  example_values      = character(),
  interpretation_status = character(),
  stringsAsFactors    = FALSE
)

for (cat_name in names(keyword_groups)) {
  keywords <- keyword_groups[[cat_name]]
  for (cn in names(meta)) {
    cn_lower <- tolower(cn)
    matched <- any(sapply(keywords, function(kw) grepl(tolower(kw), cn_lower)))
    if (matched) {
      col <- meta[[cn]]
      n_miss <- sum(is.na(col))
      n_unique <- length(unique(col[!is.na(col)]))
      ex <- paste(head(unique(col[!is.na(col)]), 5), collapse = "; ")
      candidate_df <- rbind(candidate_df, data.frame(
        category            = cat_name,
        column_name         = cn,
        column_class        = class(col)[1],
        unique_count        = n_unique,
        missing_count       = n_miss,
        example_values      = ex,
        interpretation_status = "需要人工确认",
        stringsAsFactors    = FALSE
      ))
    }
  }
}

cat(sprintf("[INFO] Found %d candidate columns across %d categories\n",
            nrow(candidate_df), length(unique(candidate_df$category))))

tryCatch({
  write.csv(candidate_df, file.path(RESULTS, "GSE243013_candidate_metadata_columns.csv"),
            row.names = FALSE)
  cat("[INFO] Saved: GSE243013_candidate_metadata_columns.csv\n")
}, error = function(e) {
  cat("[ERROR] Failed to save candidates:", conditionMessage(e), "\n")
})

## =========================================================================
## V. Value Counts for Candidates
## =========================================================================
cat("\n[VI] Computing value counts for candidate columns...\n")

value_counts_df <- data.frame(
  column_name    = character(),
  value          = character(),
  count          = integer(),
  stringsAsFactors = FALSE
)

for (cn in unique(candidate_df$column_name)) {
  col <- meta[[cn]]
  n_unique <- length(unique(col[!is.na(col)]))
  tbl <- sort(table(col), decreasing = TRUE)
  if (n_unique <= 100) {
    vals <- names(tbl)
    cnts <- as.integer(tbl)
  } else {
    vals <- head(names(tbl), 20)
    cnts <- head(as.integer(tbl), 20)
  }
  for (i in seq_along(vals)) {
    value_counts_df <- rbind(value_counts_df, data.frame(
      column_name = cn,
      value       = vals[i],
      count       = cnts[i],
      stringsAsFactors = FALSE
    ))
  }
  n_miss <- sum(is.na(col))
  if (n_miss > 0) {
    value_counts_df <- rbind(value_counts_df, data.frame(
      column_name = cn,
      value       = "<NA>",
      count       = n_miss,
      stringsAsFactors = FALSE
    ))
  }
}

tryCatch({
  write.csv(value_counts_df, file.path(RESULTS, "GSE243013_candidate_value_counts.csv"),
            row.names = FALSE)
  cat("[INFO] Saved: GSE243013_candidate_value_counts.csv\n")
}, error = function(e) {
  cat("[ERROR] Failed to save value counts:", conditionMessage(e), "\n")
})

## =========================================================================
## VI. Patient Candidate Summary
## =========================================================================
cat("\n[VII] Analyzing patient candidates...\n")

patient_candidates <- candidate_df[candidate_df$category == "patient", ]

if (nrow(patient_candidates) == 0) {
  cat("[INFO] No patient candidate columns found by keyword.\n")
  ## Try broader search
  for (cn in names(meta)) {
    col <- meta[[cn]]
    n_miss <- sum(is.na(col))
    n_unique <- length(unique(col[!is.na(col)]))
    if (n_miss < nrow(meta) * 0.5 && n_unique >= 5 && n_unique <= 500) {
      cat(sprintf("  Potential patient column (by heuristic): %s (%d unique)\n", cn, n_unique))
      patient_candidates <- rbind(patient_candidates, data.frame(
        category = "patient_heuristic", column_name = cn,
        column_class = class(col)[1], unique_count = n_unique,
        missing_count = n_miss, example_values = paste(head(unique(col[!is.na(col)]), 5), collapse = "; "),
        interpretation_status = "需要人工确认", stringsAsFactors = FALSE
      ))
    }
  }
}

patient_summary <- data.frame(
  patient_column = character(),
  n_patients     = integer(),
  min_cells      = integer(),
  median_cells   = integer(),
  max_cells      = integer(),
  multi_sample   = logical(),
  has_missing    = logical(),
  stringsAsFactors = FALSE
)

if (nrow(patient_candidates) > 0) {
  for (pc in patient_candidates$column_name) {
    col <- meta[[pc]]
    n_miss <- sum(is.na(col))
    non_miss <- col[!is.na(col)]
    n_pat <- length(unique(non_miss))
    cell_counts <- table(non_miss)
    patient_summary <- rbind(patient_summary, data.frame(
      patient_column = pc,
      n_patients     = n_pat,
      min_cells      = as.integer(min(cell_counts)),
      median_cells   = as.integer(median(cell_counts)),
      max_cells      = as.integer(max(cell_counts)),
      multi_sample   = FALSE,
      has_missing    = n_miss > 0,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  Column '%s': %d patients, cells/patient: min=%s median=%s max=%s, missing=%d\n",
                pc, n_pat, min(cell_counts), median(cell_counts), max(cell_counts), n_miss))
  }
} else {
  cat("[INFO] No patient candidate columns found.\n")
  cat("[INFO] 未能自动确定患者编号字段，需要结合论文或元数据说明人工确认。\n")
}

tryCatch({
  write.csv(patient_summary, file.path(RESULTS, "GSE243013_patient_candidate_summary.csv"),
            row.names = FALSE)
  cat("[INFO] Saved: GSE243013_patient_candidate_summary.csv\n")
}, error = function(e) {
  cat("[ERROR] Failed to save patient summary:", conditionMessage(e), "\n")
})

## =========================================================================
## VII. Response Candidate Summary
## =========================================================================
cat("\n[VIII] Analyzing response candidates...\n")

response_keywords <- c("response", "pathological", "pathologic", "MPR", "pCR",
                        "NMPR", "non-MPR", "TRG", "residual", "regression",
                        "treatment.response", "trg")

response_candidates <- data.frame(
  column_name = character(), stringsAsFactors = FALSE
)
for (cn in names(meta)) {
  cn_lower <- tolower(cn)
  if (any(sapply(response_keywords, function(kw) grepl(tolower(kw), cn_lower)))) {
    response_candidates <- rbind(response_candidates,
                                 data.frame(column_name = cn, stringsAsFactors = FALSE))
  }
}

response_summary <- data.frame(
  column_name    = character(),
  pattern        = character(),
  n_cells        = integer(),
  n_samples      = integer(),
  n_patients     = character(),
  stringsAsFactors = FALSE
)

target_patterns <- c("pCR", "MPR", "NMPR", "non-MPR", "major pathological response",
                     "non-major pathological response", "residual viable tumor", "RVT")

if (nrow(response_candidates) > 0) {
  ## Find patient column for sample/patient counting
  pat_col <- NULL
  if (nrow(patient_candidates) > 0) {
    pat_col <- patient_candidates$column_name[1]
  }
  ## Find sample column
  samp_candidates <- candidate_df[candidate_df$category == "sample", ]
  samp_col <- if (nrow(samp_candidates) > 0) samp_candidates$column_name[1] else NULL

  for (rc in response_candidates$column_name) {
    col <- meta[[rc]]
    col_str <- as.character(col)
    for (pat in target_patterns) {
      hits <- grepl(pat, col_str, ignore.case = TRUE)
      n_cells <- sum(hits, na.rm = TRUE)
      if (n_cells > 0) {
        n_samp <- if (!is.null(samp_col)) length(unique(meta[[samp_col]][hits])) else NA_integer_
        n_pat  <- if (!is.null(pat_col)) length(unique(meta[[pat_col]][hits])) else "patient column unknown"
        response_summary <- rbind(response_summary, data.frame(
          column_name = rc, pattern = pat, n_cells = n_cells,
          n_samples = n_samp, n_patients = as.character(n_pat),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
} else {
  cat("[INFO] No response candidate columns found.\n")
}

if (nrow(response_summary) > 0) {
  cat("\n  Response pattern matches:\n")
  for (i in seq_len(nrow(response_summary))) {
    cat(sprintf("    [%s] pattern='%s': %d cells, %s samples, %s patients\n",
                response_summary$column_name[i], response_summary$pattern[i],
                response_summary$n_cells[i],
                response_summary$n_samples[i],
                response_summary$n_patients[i]))
  }
} else {
  cat("[INFO] No response pattern matches found.\n")
}

tryCatch({
  write.csv(response_summary, file.path(RESULTS, "GSE243013_response_candidate_summary.csv"),
            row.names = FALSE)
  cat("[INFO] Saved: GSE243013_response_candidate_summary.csv\n")
}, error = function(e) {
  cat("[ERROR] Failed to save response summary:", conditionMessage(e), "\n")
})

## =========================================================================
## VIII. Barcode Check
## =========================================================================
cat("\n[IX] Checking cell barcodes...\n")

barcode_candidates <- candidate_df[candidate_df$category == "barcode", ]

barcode_check <- data.frame(
  column_name    = character(),
  n_missing      = integer(),
  n_distinct     = integer(),
  n_total        = integer(),
  globally_unique = logical(),
  sample_barcode_unique = logical(),
  stringsAsFactors = FALSE
)

if (nrow(barcode_candidates) > 0) {
  samp_col <- if (nrow(samp_candidates) > 0) samp_candidates$column_name[1] else NULL

  for (bc in barcode_candidates$column_name) {
    col <- meta[[bc]]
    n_miss <- sum(is.na(col))
    n_distinct <- length(unique(col[!is.na(col)]))
    n_total <- nrow(meta)
    globally_unique <- n_distinct == (n_total - n_miss)

    sb_unique <- NA
    if (!is.null(samp_col) && n_miss == 0) {
      sb <- paste(meta[[samp_col]], meta[[bc]], sep = "_")
      sb_unique <- length(unique(sb)) == n_total
    }

    barcode_check <- rbind(barcode_check, data.frame(
      column_name = bc,
      n_missing = n_miss,
      n_distinct = n_distinct,
      n_total = n_total,
      globally_unique = globally_unique,
      sample_barcode_unique = sb_unique,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  Column '%s': missing=%d, distinct=%d/%d, globally_unique=%s",
                bc, n_miss, n_distinct, n_total, globally_unique))
    if (!is.na(sb_unique)) {
      cat(sprintf(", sample+barcode_unique=%s", sb_unique))
    }
    cat("\n")
  }
} else {
  cat("[INFO] No barcode candidate columns found.\n")
}

tryCatch({
  write.csv(barcode_check, file.path(RESULTS, "GSE243013_barcode_check.csv"),
            row.names = FALSE)
  cat("[INFO] Saved: GSE243013_barcode_check.csv\n")
}, error = function(e) {
  cat("[ERROR] Failed to save barcode check:", conditionMessage(e), "\n")
})

## =========================================================================
## IX. Inspect Other Small Files
## =========================================================================
cat("\n[X] Inspecting other small files...\n")

small_files <- c(
  "GSE243013_barcodes.csv.gz",
  "GSE243013_genes.csv.gz",
  "GSE243013_NMF_all_group_5.csv.gz",
  "GSE243013_T_with_TCR_annotation.csv.gz"
)

small_summary <- data.frame(
  file_name       = character(),
  n_rows          = integer(),
  n_cols          = integer(),
  column_names    = character(),
  duplicate_rows  = integer(),
  total_na_cells  = integer(),
  stringsAsFactors = FALSE
)

for (sf in small_files) {
  fp <- file.path(RAW_DIR, sf)
  cat(sprintf("\n  --- %s ---\n", sf))
  if (!file.exists(fp)) {
    cat("  NOT FOUND\n")
    small_summary <- rbind(small_summary, data.frame(
      file_name = sf, n_rows = NA_integer_, n_cols = NA_integer_,
      column_names = "", duplicate_rows = NA_integer_, total_na_cells = NA_integer_,
      stringsAsFactors = FALSE
    ))
    next
  }
  tryCatch({
    if (requireNamespace("data.table", quietly = TRUE)) {
      df <- data.table::fread(fp, header = TRUE, data.table = FALSE)
    } else {
      df <- read.csv(gzfile(fp), check.names = FALSE)
    }
    cat(sprintf("  Rows: %d, Columns: %d\n", nrow(df), ncol(df)))
    cat(sprintf("  Column names: %s\n", paste(names(df), collapse = ", ")))
    cat("  Structure:\n")
    print(head(df, 3))
    n_dup <- sum(duplicated(df))
    n_na  <- sum(is.na(df))
    cat(sprintf("  Duplicate rows: %d, Total NA cells: %d\n", n_dup, n_na))
    small_summary <- rbind(small_summary, data.frame(
      file_name = sf, n_rows = nrow(df), n_cols = ncol(df),
      column_names = paste(names(df), collapse = "; "),
      duplicate_rows = n_dup, total_na_cells = n_na,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    cat(sprintf("  ERROR reading %s: %s\n", sf, conditionMessage(e)))
    small_summary <- rbind(small_summary, data.frame(
      file_name = sf, n_rows = NA_integer_, n_cols = NA_integer_,
      column_names = paste("ERROR:", conditionMessage(e)),
      duplicate_rows = NA_integer_, total_na_cells = NA_integer_,
      stringsAsFactors = FALSE
    ))
  })
}

tryCatch({
  write.csv(small_summary, file.path(RESULTS, "GSE243013_small_file_summary.csv"),
            row.names = FALSE)
  cat("\n[INFO] Saved: GSE243013_small_file_summary.csv\n")
}, error = function(e) {
  cat("[ERROR] Failed to save small file summary:", conditionMessage(e), "\n")
})

## UMAP archive contents
cat("\n[XI] Listing UMAP_info.tar.gz contents...\n")
umap_file <- file.path(RAW_DIR, "GSE243013_UMAP_info.tar.gz")
if (file.exists(umap_file)) {
  tryCatch({
    umap_contents <- untar(umap_file, list = TRUE)
    umap_df <- data.frame(
      archive_name = "GSE243013_UMAP_info.tar.gz",
      internal_file = umap_contents,
      stringsAsFactors = FALSE
    )
    write.csv(umap_df, file.path(RESULTS, "GSE243013_UMAP_archive_contents.csv"),
              row.names = FALSE)
    cat(sprintf("[INFO] UMAP archive contains %d files:\n", length(umap_contents)))
    cat(paste("  ", head(umap_contents, 20), collapse = "\n"), "\n")
    cat("[INFO] Saved: GSE243013_UMAP_archive_contents.csv\n")
  }, error = function(e) {
    cat(sprintf("[ERROR] Cannot list UMAP archive: %s\n", conditionMessage(e)))
  })
} else {
  cat("[INFO] UMAP_info.tar.gz not found.\n")
}

## =========================================================================
## X. Final Summary
## =========================================================================
cat("\n========================================================================\n")
cat("SUMMARY\n")
cat("========================================================================\n")

cat("\nDownload results:\n")
for (i in seq_len(nrow(download_status))) {
  cat(sprintf("  %-55s  %s\n", download_status$file_name[i], download_status$status[i]))
}

cat(sprintf("\nMetadata dimensions: %d rows x %d columns\n", nrow(meta), ncol(meta)))
cat(sprintf("Column names: %s\n", paste(names(meta), collapse = ", ")))

cat("\nCandidate columns by category:\n")
for (cat_name in unique(candidate_df$category)) {
  cols <- candidate_df$column_name[candidate_df$category == cat_name]
  cat(sprintf("  %s: %s\n", cat_name, paste(cols, collapse = ", ")))
}

if (nrow(patient_candidates) > 0) {
  cat(sprintf("\nMost likely patient column: %s\n", patient_candidates$column_name[1]))
  cat(sprintf("Patient count: %d\n", patient_summary$n_patients[1]))
}

if (nrow(response_candidates) > 0) {
  cat(sprintf("Most likely response column: %s\n", response_candidates$column_name[1]))
}

if (nrow(barcode_candidates) > 0) {
  cat(sprintf("Barcode column: %s\n", barcode_candidates$column_name[1]))
}

cat("\n========================================================================\n")
cat("Step 01 completed.\n")
cat("========================================================================\n")
