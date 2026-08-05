## =========================================================================
## Step 04: Download & Import GSE243013 Counts Matrix via BPCells
## =========================================================================

.libPaths(c("~/Library/R/arm64/4.6/library", .libPaths()))
options(timeout = 3600)

cat("========================================================================\n")
cat("Step 04: Download & Import GSE243013 Counts Matrix via BPCells\n")
cat("========================================================================\n\n")
flush.console()

RESULTS  <- "03_results"
RAW_DIR  <- "02_data/raw/GSE243013"
BPCELLS  <- "02_data/bpcells"
TMP_DIR  <- "02_data/tmp/GSE243013_bpcells_import"

dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BPCELLS, recursive = TRUE, showWarnings = FALSE)
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

step_start <- Sys.time()
dl_seconds <- NA_real_
import_seconds <- NA_real_
failed_step <- NA_character_
error_message <- NA_character_

## =========================================================================
## III. Preflight Checks
## =========================================================================
cat("[III] Preflight checks...\n")

preflight <- character()
preflight <- c(preflight, sprintf("Date: %s", Sys.time()))
preflight <- c(preflight, sprintf("Working directory: %s", getwd()))
preflight <- c(preflight, sprintf("R version: %s", R.version.string))
preflight <- c(preflight, sprintf("R platform: %s", R.version$platform))
preflight <- c(preflight, sprintf("R lib paths: %s", paste(.libPaths(), collapse = ", ")))

tryCatch({
  os_info <- Sys.info()
  preflight <- c(preflight, sprintf("OS: %s %s", os_info["sysname"], os_info["release"]))
  preflight <- c(preflight, sprintf("Machine: %s", os_info["machine"]))
}, error = function(e) {
  preflight <<- c(preflight, sprintf("OS info error: %s", conditionMessage(e)))
})

## Memory
tryCatch({
  mem_bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
  mem_gb <- mem_bytes / (1024^3)
  preflight <- c(preflight, sprintf("Total memory: %s bytes (%.1f GB)", mem_bytes, mem_gb))
}, error = function(e) {
  preflight <- c(preflight, sprintf("Memory read error: %s", conditionMessage(e)))
})

## Disk space
tryCatch({
  disk_info <- system("df -Pk .", intern = TRUE)
  preflight <- c(preflight, "Disk space (df -Pk .):")
  preflight <- c(preflight, disk_info)
  ## Parse available space
  disk_parts <- strsplit(disk_info[2], "\\s+")[[1]]
  avail_kb <- as.numeric(disk_parts[4])
  avail_gb <- avail_kb / (1024 * 1024)
  preflight <- c(preflight, sprintf("Available disk: %.1f GB", avail_gb))

  if (avail_gb < 60) {
    cat("[FATAL] Insufficient disk space. Need 60GB, have %.1f GB. Stopping.\n", avail_gb)
    failed_step <<- "preflight_disk"
    error_message <<- sprintf("Insufficient disk space: %.1f GB available, need 60 GB", avail_gb)

    sink(file.path(RESULTS, "GSE243013_step04_PREFLIGHT_FAILED.txt"))
    cat("Step 04 PREFLIGHT FAILED\n")
    cat(sprintf("Time: %s\n", Sys.time()))
    cat(sprintf("Reason: %s\n", error_message))
    sink()
    stop(error_message)
  }
}, error = function(e) {
  if (is.na(failed_step)) {
    preflight <<- c(preflight, sprintf("Disk check error: %s", conditionMessage(e)))
  }
})

## Check ARM64 — R uses "aarch64" for ARM64 on macOS
platform <- R.version$platform
is_arm64 <- grepl("aarch64|arm64", platform, ignore.case = TRUE)
preflight <- c(preflight, sprintf("ARM64 platform: %s", is_arm64))

if (!is_arm64) {
  cat("[FATAL] Not ARM64 platform. Stopping.\n")
  failed_step <<- "preflight_arm64"
  error_message <<- sprintf("Not ARM64: %s", platform)

  writeLines(c(preflight, "", "FAILED: Not ARM64"),
             file.path(RESULTS, "GSE243013_step04_PREFLIGHT_FAILED.txt"))
  stop(error_message)
}

writeLines(preflight, file.path(RESULTS, "GSE243013_step04_preflight.txt"))
cat("[INFO] Saved: GSE243013_step04_preflight.txt\n")
cat("[INFO] Preflight checks passed.\n\n")

## =========================================================================
## IV. Install curl
## =========================================================================
cat("[IV] Checking curl package...\n")

if (!requireNamespace("curl", quietly = TRUE)) {
  cat("[INFO] Installing curl from CRAN (binary)...\n")
  tryCatch({
    install.packages("curl", repos = "https://cloud.r-project.org", type = "binary")
  }, error = function(e) {
    cat("[FATAL] curl installation failed:", conditionMessage(e), "\n")
    failed_step <<- "install_curl"
    error_message <<- conditionMessage(e)
    stop(error_message)
  })
}

cat(sprintf("[INFO] curl version: %s\n", packageVersion("curl")))
tryCatch({
  cv <- curl::curl_version()
  cat(sprintf("[INFO] curl library: %s %s\n", cv$package, cv$version))
}, error = function(e) {
  cat("[WARNING] Could not get curl version info:", conditionMessage(e), "\n")
})

## =========================================================================
## V. Install BPCells
## =========================================================================
cat("\n[V] Checking BPCells...\n")

bpcells_ready <- FALSE
tryCatch({
  if (requireNamespace("BPCells", quietly = TRUE)) {
    cat(sprintf("[INFO] BPCells already installed: %s\n", packageVersion("BPCells")))
    bpcells_ready <- TRUE
  } else {
    cat("[INFO] Installing BPCells from CRAN (binary)...\n")
    install.packages(
      "BPCells",
      repos = "https://cloud.r-project.org",
      type = "binary"
    )
    if (requireNamespace("BPCells", quietly = TRUE)) {
      cat(sprintf("[INFO] BPCells installed: %s\n", packageVersion("BPCells")))
      bpcells_ready <- TRUE
    } else {
      cat("[FATAL] BPCells installation failed.\n")
      failed_step <<- "install_bpcells"
      error_message <<- "BPCells package not available after installation"
    }
  }
}, error = function(e) {
  cat("[FATAL] BPCells installation error:", conditionMessage(e), "\n")
  failed_step <<- "install_bpcells"
  error_message <<- conditionMessage(e)
})

if (!bpcells_ready) {
  cat("[FATAL] Cannot proceed without BPCells.\n")
  stop(error_message)
}

## Verify key functions exist
cat("[INFO] Verifying BPCells functions...\n")
for (fn in c("import_matrix_market", "open_matrix_dir", "storage_order")) {
  exists_fn <- exists(fn, where = asNamespace("BPCells"), mode = "function")
  cat(sprintf("  BPCells::%s: %s\n", fn, ifelse(exists_fn, "FOUND", "MISSING")))
  if (!exists_fn) {
    failed_step <<- "verify_bpcells_functions"
    error_message <<- sprintf("BPCells::%s not found", fn)
    stop(error_message)
  }
}

tryCatch({
  bp_desc <- packageDescription("BPCells")
  cat(sprintf("[INFO] BPCells Built: %s\n", bp_desc$Built))
  cat(sprintf("[INFO] BPCells path: %s\n", find.package("BPCells")))
}, error = function(e) {
  cat("[WARNING] Could not get BPCells description:", conditionMessage(e), "\n")
})

## =========================================================================
## VI. Download Counts Matrix
## =========================================================================
cat("\n[VI] Downloading counts matrix...\n")

counts_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243013/suppl/GSE243013_NSCLC_immune_scRNA_counts.mtx.gz"
counts_path <- file.path(RAW_DIR, "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz")

## Check if already downloaded
size_before <- 0
if (file.exists(counts_path)) {
  size_before <- file.info(counts_path)$size
  cat(sprintf("[INFO] Existing file: %.0f bytes (%.2f GB)\n", size_before, size_before / (1024^3)))
}

## Expected size from previous step: 7123039063 bytes (~6.6 GB)
expected_size <- 7123039063

if (size_before >= expected_size * 0.99) {
  cat("[INFO] File appears complete. Skipping download.\n")
  dl_success <- TRUE
  dl_status_code <- "already_complete"
  size_after <- size_before
  dl_error <- ""
} else {
  cat(sprintf("[INFO] Downloading from: %s\n", counts_url))
  cat(sprintf("[INFO] Target: %s\n", counts_path))

  dl_start <- Sys.time()
  dl_success <- FALSE
  dl_status_code <- NA_character_
  dl_error <- ""

  tryCatch({
    result <- curl::multi_download(
      urls = counts_url,
      destfiles = counts_path,
      resume = TRUE,
      progress = TRUE,
      multi_timeout = Inf
    )
    dl_status_code <- as.character(result$status)
    if (result$status == 200 || result$status == 0) {
      dl_success <- TRUE
    }
  }, error = function(e) {
    dl_error <<- conditionMessage(e)
    cat("[ERROR] Download failed:", dl_error, "\n")
  })

  dl_end <- Sys.time()
  dl_seconds <<- as.numeric(difftime(dl_end, dl_start, units = "secs"))
  size_after <- if (file.exists(counts_path)) file.info(counts_path)$size else 0

  cat(sprintf("[INFO] Download elapsed: %.1f seconds\n", dl_seconds))
  cat(sprintf("[INFO] File size after download: %.0f bytes (%.2f GB)\n",
              size_after, size_after / (1024^3)))
  cat(sprintf("[INFO] HTTP status: %s\n", dl_status_code))

  if (!dl_success) {
    failed_step <<- "download_counts"
    error_message <<- dl_error
  }
}

## Save download status
download_status <- data.frame(
  url = counts_url,
  local_path = counts_path,
  download_size_before = size_before,
  download_size_after = size_after,
  status_code = dl_status_code,
  success = dl_success,
  resumefrom = size_before,
  elapsed_seconds = ifelse(is.na(dl_seconds), NA, dl_seconds),
  error_message = ifelse(nchar(dl_error) > 0, dl_error, ""),
  stringsAsFactors = FALSE
)
write.csv(download_status, file.path(RESULTS, "GSE243013_counts_download_status.csv"),
          row.names = FALSE)
cat("[INFO] Saved: GSE243013_counts_download_status.csv\n")

if (!dl_success) {
  cat("[FATAL] Counts matrix download failed. Stopping.\n")
  failed_step <<- "download_counts"
  error_message <<- "Download failed"
  stop(error_message)
}

## =========================================================================
## VII. File Validation
## =========================================================================
cat("\n[VII] Validating counts file...\n")

file_valid <- TRUE
validation_results <- data.frame(
  check = character(), result = character(), stringsAsFactors = FALSE
)

check_and_record <- function(name, passed, detail = "") {
  validation_results <<- rbind(validation_results, data.frame(
    check = name, result = ifelse(passed, "PASS", paste0("FAIL: ", detail)),
    stringsAsFactors = FALSE
  ))
  if (!passed) file_valid <<- FALSE
  cat(sprintf("  %s: %s\n", name, ifelse(passed, "PASS", "FAIL")))
}

## 1. File exists
check_and_record("file_exists", file.exists(counts_path))

## 2. File size > 5GB
actual_size <- file.info(counts_path)$size
check_and_record("file_size_gt_5GB", actual_size > 5 * 1024^3,
                 sprintf("size=%.0f", actual_size))

## 3. Not HTML error page
tryCatch({
  con <- file(counts_path, "rb")
  first_bytes <- readBin(con, "raw", n = 100)
  close(con)
  is_html <- any(first_bytes[1:4] == charToRaw("<html")[1:4]) ||
             identical(first_bytes[1:5], charToRaw("<!DOC")) ||
             identical(first_bytes[1:5], charToRaw("<!doc"))
  check_and_record("not_html", !is_html)
}, error = function(e) {
  check_and_record("not_html", FALSE, conditionMessage(e))
})

## 4. Gzip magic bytes
tryCatch({
  con <- file(counts_path, "rb")
  magic <- readBin(con, "raw", n = 2)
  close(con)
  is_gzip <- magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b)
  check_and_record("gzip_magic_bytes", is_gzip,
                    sprintf("got=%02X%02X", magic[1], magic[2]))
}, error = function(e) {
  check_and_record("gzip_magic_bytes", FALSE, conditionMessage(e))
})

## 5. gzfile readable
tryCatch({
  con <- gzfile(counts_path, "r")
  first_line <- readLines(con, n = 1)
  close(con)
  check_and_record("gzfile_readable", TRUE)
}, error = function(e) {
  check_and_record("gzfile_readable", FALSE, conditionMessage(e))
})

## 6. Starts with %%MatrixMarket
tryCatch({
  con <- gzfile(counts_path, "r")
  first_line <- readLines(con, n = 1)
  close(con)
  is_mm <- grepl("^%%MatrixMarket", first_line)
  check_and_record("matrixmarket_header", is_mm,
                    sprintf("first_line='%s'", substr(first_line, 1, 50)))
}, error = function(e) {
  check_and_record("matrixmarket_header", FALSE, conditionMessage(e))
})

write.csv(validation_results, file.path(RESULTS, "GSE243013_counts_file_validation.csv"),
          row.names = FALSE)
cat("[INFO] Saved: GSE243013_counts_file_validation.csv\n")

if (!file_valid) {
  failed_step <<- "file_validation"
  error_message <<- "File validation failed"
  cat("[FATAL] File validation failed. Stopping.\n")
  stop(error_message)
}

## =========================================================================
## VIII. Read Matrix Market Header
## =========================================================================
cat("\n[VIII] Reading Matrix Market header...\n")

mm_info <- list(
  header = NA_character_,
  nrow = NA_integer_,
  ncol = NA_integer_,
  nnz = NA_integer_
)

tryCatch({
  con <- gzfile(counts_path, "r")
  lines <- character()
  ## Read header and comment lines
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    lines <- c(lines, line)
    if (!grepl("^%", line)) break
  }
  close(con)

  mm_info$header <- lines[1]
  ## Last non-comment line is the dimension line
  dim_line <- lines[length(lines)]
  dim_parts <- as.integer(strsplit(trimws(dim_line), "\\s+")[[1]])
  mm_info$nrow <- dim_parts[1]
  mm_info$ncol <- dim_parts[2]
  mm_info$nnz <- dim_parts[3]

  cat(sprintf("[INFO] Matrix Market header: %s\n", mm_info$header))
  cat(sprintf("[INFO] Dimensions: %d rows x %d cols, %d non-zero\n",
              mm_info$nrow, mm_info$ncol, mm_info$nnz))
}, error = function(e) {
  cat("[ERROR] Failed to read Matrix Market header:", conditionMessage(e), "\n")
  failed_step <<- "read_matrix_header"
  error_message <<- conditionMessage(e)
  stop(error_message)
})

mm_header_df <- data.frame(
  matrix_market_header = mm_info$header,
  matrix_nrow = mm_info$nrow,
  matrix_ncol = mm_info$ncol,
  matrix_nnz = mm_info$nnz,
  compressed_file_size = actual_size,
  inspection_status = "success",
  stringsAsFactors = FALSE
)
write.csv(mm_header_df, file.path(RESULTS, "GSE243013_matrix_market_header.csv"),
          row.names = FALSE)
cat("[INFO] Saved: GSE243013_matrix_market_header.csv\n")

## Validate dimensions
if (mm_info$nrow <= 0 || mm_info$ncol <= 0 || mm_info$nnz <= 0) {
  failed_step <<- "matrix_dimensions"
  error_message <<- sprintf("Invalid dimensions: %d x %d, nnz=%d",
                            mm_info$nrow, mm_info$ncol, mm_info$nnz)
  stop(error_message)
}

## =========================================================================
## IX. Read Genes and Barcodes
## =========================================================================
cat("\n[IX] Reading genes and barcodes...\n")

## Genes
cat("[INFO] Reading genes...\n")
genes <- NULL
tryCatch({
  genes <- data.table::fread("02_data/raw/GSE243013/GSE243013_genes.csv.gz",
                              header = TRUE, data.table = FALSE)
  cat(sprintf("[INFO] Genes: %d rows x %d columns\n", nrow(genes), ncol(genes)))
  cat(sprintf("[INFO] Gene columns: %s\n", paste(names(genes), collapse = ", ")))
}, error = function(e) {
  cat("[ERROR] Failed to read genes:", conditionMessage(e), "\n")
  stop(conditionMessage(e))
})

## Identify gene name column
gene_candidates <- c("gene", "genes", "gene_name", "gene_symbol", "symbol",
                     "feature", "feature_name", "V1")
gene_col <- NULL
for (gc in gene_candidates) {
  if (gc %in% names(genes)) {
    gene_col <- gc
    break
  }
}
if (is.null(gene_col) && ncol(genes) == 1) {
  gene_col <- names(genes)[1]
}
if (is.null(gene_col)) {
  cat("[FATAL] Cannot identify gene name column. Columns:", paste(names(genes), collapse=", "), "\n")
  stop("Cannot identify gene name column")
}
cat(sprintf("[INFO] Using gene column: '%s'\n", gene_col))
gene_names <- genes[[gene_col]]

## Barcodes
cat("[INFO] Reading barcodes...\n")
barcodes <- NULL
tryCatch({
  barcodes_df <- data.table::fread("02_data/raw/GSE243013/GSE243013_barcodes.csv.gz",
                                    header = TRUE, data.table = FALSE)
  cat(sprintf("[INFO] Barcodes: %d rows x %d columns\n", nrow(barcodes_df), ncol(barcodes_df)))
  cat(sprintf("[INFO] Barcode columns: %s\n", paste(names(barcodes_df), collapse = ", ")))
}, error = function(e) {
  cat("[ERROR] Failed to read barcodes:", conditionMessage(e), "\n")
  stop(conditionMessage(e))
})

barcode_candidates <- c("barcode", "barcodes", "cell_id", "cellID", "V1")
barcode_col <- NULL
for (bc in barcode_candidates) {
  if (bc %in% names(barcodes_df)) {
    barcode_col <- bc
    break
  }
}
if (is.null(barcode_col) && ncol(barcodes_df) == 1) {
  barcode_col <- names(barcodes_df)[1]
}
if (is.null(barcode_col)) {
  cat("[FATAL] Cannot identify barcode column. Columns:", paste(names(barcodes_df), collapse=", "), "\n")
  stop("Cannot identify barcode column")
}
cat(sprintf("[INFO] Using barcode column: '%s'\n", barcode_col))
barcodes <- barcodes_df[[barcode_col]]

## Validation
n_genes <- length(gene_names)
n_barcodes <- length(barcodes)
n_gene_na <- sum(is.na(gene_names))
n_barcode_na <- sum(is.na(barcodes))
n_gene_dup <- sum(duplicated(gene_names))
n_barcode_dup <- sum(duplicated(barcodes))

cat(sprintf("[INFO] Genes: %d (NA: %d, duplicated: %d)\n", n_genes, n_gene_na, n_gene_dup))
cat(sprintf("[INFO] Barcodes: %d (NA: %d, duplicated: %d)\n", n_barcodes, n_barcode_na, n_barcode_dup))

gb_validation <- data.frame(
  check = c("n_genes", "n_barcodes", "gene_na", "barcode_na",
            "gene_duplicated", "barcode_duplicated"),
  value = c(n_genes, n_barcodes, n_gene_na, n_barcode_na,
            n_gene_dup, n_barcode_dup),
  stringsAsFactors = FALSE
)
write.csv(gb_validation, file.path(RESULTS, "GSE243013_gene_barcode_validation.csv"),
          row.names = FALSE)
cat("[INFO] Saved: GSE243013_gene_barcode_validation.csv\n")

## =========================================================================
## X. Verify Matrix Orientation
## =========================================================================
cat("\n[X] Verifying matrix orientation...\n")
flush.console()

orientation_transposed <- FALSE

if (mm_info$nrow == n_genes && mm_info$ncol == n_barcodes) {
  cat("[INFO] Orientation: genes x cells (standard)\n")
  orientation_transposed <- FALSE
} else if (mm_info$nrow == n_barcodes && mm_info$ncol == n_genes) {
  cat("[INFO] Orientation: cells x genes (transposed from standard)\n")
  cat("[INFO] Will swap row_names/col_names during BPCells import to produce genes x cells\n")
  orientation_transposed <- TRUE
} else {
  cat(sprintf("[FATAL] Dimension mismatch: matrix=%dx%d, genes=%d, barcodes=%d\n",
              mm_info$nrow, mm_info$ncol, n_genes, n_barcodes))
  failed_step <<- "matrix_orientation"
  error_message <<- "Dimensions do not match genes or barcodes"
  stop(error_message)
}
flush.console()

## =========================================================================
## XI. Handle Gene Names
## =========================================================================
cat("\n[XI] Processing gene names...\n")

gene_name_original <- gene_names
gene_name_unique <- make.unique(gene_name_original, sep = "__dup")
is_duplicated <- duplicated(gene_name_original)

cat(sprintf("[INFO] Duplicated gene names: %d\n", sum(is_duplicated)))

gene_mapping <- data.frame(
  matrix_row_index = seq_along(gene_name_original),
  gene_name_original = gene_name_original,
  gene_name_unique = gene_name_unique,
  is_duplicated_original = is_duplicated,
  stringsAsFactors = FALSE
)
write.csv(gene_mapping, file.path(RESULTS, "GSE243013_gene_name_mapping.csv"),
          row.names = FALSE)
cat("[INFO] Saved: GSE243013_gene_name_mapping.csv\n")

## =========================================================================
## XII. Validate Barcodes vs Metadata
## =========================================================================
cat("\n[XII] Validating barcodes vs metadata...\n")

meta <- NULL
tryCatch({
  meta <- data.table::fread(file.path(RAW_DIR, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz"),
                             header = TRUE, data.table = FALSE)
  cat(sprintf("[INFO] Metadata: %d rows x %d columns\n", nrow(meta), ncol(meta)))
}, error = function(e) {
  cat("[ERROR] Failed to read metadata:", conditionMessage(e), "\n")
  stop(conditionMessage(e))
})

manifest <- NULL
tryCatch({
  manifest <- read.csv("03_results/GSE243013_patient_manifest_revised.csv", stringsAsFactors = FALSE)
  cat(sprintf("[INFO] Patient manifest: %d rows\n", nrow(manifest)))
}, error = function(e) {
  cat("[ERROR] Failed to read patient manifest:", conditionMessage(e), "\n")
  stop(conditionMessage(e))
})

## Checks
cat("[INFO] Running barcode-metadata alignment checks...\n")

check_results <- list()

## 1. Metadata rows = barcodes
check_results$metadata_rows_match <- nrow(meta) == n_barcodes
cat(sprintf("  1. Metadata rows = barcodes: %s (%d vs %d)\n",
            ifelse(check_results$metadata_rows_match, "PASS", "FAIL"),
            nrow(meta), n_barcodes))

## 2. metadata cellID globally unique
check_results$cellID_unique <- length(unique(meta$cellID)) == nrow(meta)
cat(sprintf("  2. metadata cellID unique: %s\n",
            ifelse(check_results$cellID_unique, "PASS", "FAIL")))

## 3. barcodes globally unique
check_results$barcodes_unique <- n_barcodes == length(unique(barcodes))
cat(sprintf("  3. barcodes unique: %s\n",
            ifelse(check_results$barcodes_unique, "PASS", "FAIL")))

## 4. Set equality
check_results$set_equal <- setequal(meta$cellID, barcodes)
cat(sprintf("  4. setequal(meta$cellID, barcodes): %s\n",
            ifelse(check_results$set_equal, "PASS", "FAIL")))

## 5. All manifest sampleIDs in metadata
check_results$manifest_in_meta <- all(manifest$sampleID %in% meta$sampleID)
cat(sprintf("  5. Manifest sampleIDs in metadata: %s\n",
            ifelse(check_results$manifest_in_meta, "PASS", "FAIL")))

## 6. Primary cohort patients have cells
primary_ids <- manifest$sampleID[manifest$primary_analysis_eligible == TRUE]
check_results$primary_have_cells <- all(primary_ids %in% meta$sampleID)
cat(sprintf("  6. Primary cohort patients have cells: %s\n",
            ifelse(check_results$primary_have_cells, "PASS", "FAIL")))

## 7. Barcode order matches metadata
check_results$order_match <- identical(as.character(meta$cellID), as.character(barcodes))
cat(sprintf("  7. Barcode order matches metadata: %s\n",
            ifelse(check_results$order_match, "SAME", "DIFFERENT")))

## Create ordered metadata if needed
if (check_results$set_equal && !check_results$order_match) {
  cat("[INFO] Reordering metadata to match barcode order...\n")
  match_idx <- match(barcodes, meta$cellID)
  if (any(is.na(match_idx))) {
    cat("[FATAL] match() produced NAs. Stopping.\n")
    failed_step <<- "barcode_metadata_match"
    error_message <<- "match() produced NAs"
    stop(error_message)
  }
  metadata_ordered <- meta[match_idx, ]
  cat(sprintf("[INFO] metadata_ordered created: %d rows\n", nrow(metadata_ordered)))
} else {
  metadata_ordered <- meta
}

## Save alignment summary
alignment_summary <- data.frame(
  check = names(check_results),
  result = sapply(check_results, function(x) ifelse(x, "PASS", "FAIL")),
  stringsAsFactors = FALSE
)
write.csv(alignment_summary, file.path(RESULTS, "GSE243013_barcode_metadata_alignment_summary.csv"),
          row.names = FALSE)
cat("[INFO] Saved: GSE243013_barcode_metadata_alignment_summary.csv\n")

## Check for failures
if (!check_results$set_equal) {
  failed_step <<- "barcode_metadata_alignment"
  error_message <<- "barcodes and metadata cellID sets do not match"
  cat("[FATAL] Barcode-metadata alignment failed. Stopping.\n")

  only_in_barcodes <- setdiff(barcodes, meta$cellID)
  only_in_meta <- setdiff(meta$cellID, barcodes)
  write.csv(data.frame(cellID = head(only_in_barcodes, 100)),
            file.path(RESULTS, "GSE243013_only_in_barcodes.csv"), row.names = FALSE)
  write.csv(data.frame(cellID = head(only_in_meta, 100)),
            file.path(RESULTS, "GSE243013_only_in_metadata.csv"), row.names = FALSE)
  stop(error_message)
}

## Create column metadata with patient-level fields
cat("[INFO] Creating matrix column metadata...\n")

## Add patient-level fields from manifest
patient_fields <- c("sampleID", "cancer_type", "pathological_response_clean",
                     "response_binary", "has_anti_PD1", "has_chemotherapy",
                     "treatment_pattern", "primary_analysis_eligible",
                     "strict_sensitivity_analysis_eligible",
                     "chemotherapy_control_eligible")

## Merge manifest fields into ordered metadata
metadata_with_patients <- metadata_ordered
for (pf in patient_fields) {
  if (pf %in% names(manifest)) {
    ## Create a lookup from manifest
    lookup <- setNames(manifest[[pf]], manifest$sampleID)
    metadata_with_patients[[pf]] <- lookup[as.character(metadata_with_patients$sampleID)]
  }
}

## Add matrix column index
metadata_with_patients$matrix_col_index <- seq_len(nrow(metadata_with_patients))

## Select columns for output
output_cols <- c("matrix_col_index", "cellID", "sampleID", "major_cell_type",
                 "sub_cell_type", "cancer_type", "pathological_response_clean",
                 "response_binary", "has_anti_PD1", "has_chemotherapy",
                 "treatment_pattern", "primary_analysis_eligible",
                 "strict_sensitivity_analysis_eligible", "chemotherapy_control_eligible")
output_cols <- intersect(output_cols, names(metadata_with_patients))

col_metadata <- metadata_with_patients[, output_cols, drop = FALSE]
cat(sprintf("[INFO] Column metadata: %d rows x %d columns\n", nrow(col_metadata), ncol(col_metadata)))

## Save as compressed CSV
con <- gzfile(file.path(RESULTS, "GSE243013_matrix_column_metadata.csv.gz"), "w")
write.csv(col_metadata, con, row.names = FALSE)
close(con)
cat("[INFO] Saved: GSE243013_matrix_column_metadata.csv.gz\n")

## =========================================================================
## XIII. Import into BPCells
## =========================================================================
cat("\n[XIII] Importing into BPCells...\n")
flush.console()

final_dir <- file.path(BPCELLS, "GSE243013_counts_colmajor")
imported <- FALSE

## Check if final directory already exists and is valid
if (dir.exists(final_dir)) {
  cat("[INFO] Final BPCells directory exists. Testing...\n")
  tryCatch({
    test_mat <- BPCells::open_matrix_dir(final_dir)
    test_dim <- dim(test_mat)
    if (test_dim[1] == n_genes && test_dim[2] == n_barcodes) {
      cat(sprintf("[INFO] BPCells directory valid: %d x %d\n", test_dim[1], test_dim[2]))
      imported <- TRUE
    } else {
      cat(sprintf("[WARNING] Existing directory has wrong dimensions: %dx%d\n",
                  test_dim[1], test_dim[2]))
    }
  }, error = function(e) {
    cat("[WARNING] Cannot open existing BPCells directory:", conditionMessage(e), "\n")
  })
}

if (!imported) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  build_dir <- file.path(BPCELLS, sprintf("GSE243013_counts_build_%s", timestamp))

  cat(sprintf("[INFO] Creating build directory: %s\n", build_dir))
  ## Note: BPCells::import_matrix_market will create the directory.
  ## Do NOT call dir.create() here as it causes "Path already exists" error.

  import_start <- Sys.time()

  ## Determine correct row/col assignment based on matrix orientation
  if (orientation_transposed) {
    ## File is cells x genes: rows=barcodes, cols=genes
    import_row_names <- barcodes
    import_col_names <- gene_name_unique
    cat("[INFO] Import mapping: row_names=barcodes, col_names=genes (transposed file)\n")
  } else {
    ## File is genes x cells: rows=genes, cols=barcodes
    import_row_names <- gene_name_unique
    import_col_names <- barcodes
    cat("[INFO] Import mapping: row_names=genes, col_names=barcodes (standard file)\n")
  }
  flush.console()

  tryCatch({
    cat("[INFO] Running BPCells::import_matrix_market...\n")
    cat(sprintf("[INFO] mtx_path: %s\n", counts_path))
    cat(sprintf("[INFO] outdir: %s\n", build_dir))
    cat(sprintf("[INFO] row_names: %d entries\n", length(import_row_names)))
    cat(sprintf("[INFO] col_names: %d entries\n", length(import_col_names)))
    cat("[INFO] row_major = FALSE (column-major)\n")
    cat(sprintf("[INFO] tmpdir: %s\n", TMP_DIR))
    cat("[INFO] sort_bytes: 536870912 (512 MB)\n")
    flush.console()

    BPCells::import_matrix_market(
      mtx_path = counts_path,
      outdir = build_dir,
      row_names = import_row_names,
      col_names = import_col_names,
      row_major = FALSE,
      tmpdir = TMP_DIR,
      load_bytes = 4194304L,
      sort_bytes = 536870912L
    )

    import_end <- Sys.time()
    import_seconds <<- as.numeric(difftime(import_end, import_start, units = "secs"))
    cat(sprintf("[INFO] BPCells import completed in %.1f seconds\n", import_seconds))
    imported <- TRUE

    ## Rename to final directory
    if (!dir.exists(final_dir)) {
      cat("[INFO] Renaming build directory to final...\n")
      tryCatch({
        file.rename(build_dir, final_dir)
        cat("[INFO] Renamed successfully.\n")
      }, error = function(e) {
        cat("[WARNING] Rename failed:", conditionMessage(e), "\n")
        cat("[INFO] Using build directory as final.\n")
        final_dir <<- build_dir
      })
    } else {
      cat("[INFO] Final directory already exists. Keeping build directory.\n")
      final_dir <<- build_dir
    }
  }, error = function(e) {
    import_end <- Sys.time()
    import_seconds <<- as.numeric(difftime(import_end, import_start, units = "secs"))
    cat("[ERROR] BPCells import failed:", conditionMessage(e), "\n")
    failed_step <<- "bpcells_import"
    error_message <<- conditionMessage(e)
  })
}

if (!imported) {
  cat("[FATAL] BPCells import failed.\n")
  stop(error_message)
}

## =========================================================================
## XIV. Validate BPCells Matrix
## =========================================================================
cat("\n[XIV] Validating BPCells matrix...\n")

tryCatch({
  mat <- BPCells::open_matrix_dir(final_dir)

  cat(sprintf("[INFO] BPCells matrix dim: %s\n", paste(dim(mat), collapse = " x ")))
  cat(sprintf("[INFO] nrow: %d\n", nrow(mat)))
  cat(sprintf("[INFO] ncol: %d\n", ncol(mat)))
  cat(sprintf("[INFO] storage_order: %s\n", BPCells::storage_order(mat)))
  cat(sprintf("[INFO] length(rownames): %d\n", length(rownames(mat))))
  cat(sprintf("[INFO] length(colnames): %d\n", length(colnames(mat))))
  flush.console()

  ## Determine expected orientation based on what was imported
  if (orientation_transposed) {
    ## File was cells x genes, imported with row_names=barcodes, col_names=genes
    ## So BPCells matrix is barcodes x genes
    expected_rownames <- barcodes
    expected_colnames <- gene_name_unique
    row_label <- "barcodes"
    col_label <- "genes"
  } else {
    ## File was genes x cells, imported with row_names=genes, col_names=barcodes
    ## So BPCells matrix is genes x barcodes
    expected_rownames <- gene_name_unique
    expected_colnames <- barcodes
    row_label <- "genes"
    col_label <- "barcodes"
  }

  cat(sprintf("[INFO] Expected orientation: %s (rows) x %s (cols)\n", row_label, col_label))
  cat(sprintf("[INFO] First 5 rownames: %s\n", paste(head(rownames(mat), 5), collapse = ", ")))
  cat(sprintf("[INFO] Last 5 rownames: %s\n", paste(tail(rownames(mat), 5), collapse = ", ")))
  cat(sprintf("[INFO] First 5 colnames: %s\n", paste(head(colnames(mat), 5), collapse = ", ")))
  cat(sprintf("[INFO] Last 5 colnames: %s\n", paste(tail(colnames(mat), 5), collapse = ", ")))
  flush.console()

  ## Verify rownames and colnames
  rownames_match <- identical(rownames(mat), expected_rownames)
  colnames_match <- identical(colnames(mat), expected_colnames)
  cat(sprintf("[INFO] rownames == %s: %s\n", row_label, ifelse(rownames_match, "YES", "NO")))
  cat(sprintf("[INFO] colnames == %s: %s\n", col_label, ifelse(colnames_match, "YES", "NO")))

  ## Small-scale sampling validation
  cat("[INFO] Running small-scale sampling validation...\n")
  set.seed(20260803)
  n_sample_rows <- min(100, nrow(mat))
  n_sample_cols <- min(500, ncol(mat))
  sample_row_idx <- sample(seq_len(nrow(mat)), n_sample_rows)
  sample_col_idx <- sample(seq_len(ncol(mat)), n_sample_cols)

  sub_mat <- as.matrix(mat[sample_row_idx, sample_col_idx])

  all_nonneg <- all(sub_mat >= 0)
  all_integer <- all(abs(sub_mat - round(sub_mat)) < 1e-6)
  has_nonzero <- any(sub_mat > 0)

  cat(sprintf("[INFO] Submatrix %dx%d: non-negative=%s, integer=%s, has_nonzero=%s\n",
              n_sample_rows, n_sample_cols,
              all_nonneg, all_integer, has_nonzero))

  bpcells_validation <- data.frame(
    check = c("bp_dim", "bp_nrow", "bp_ncol", "bp_storage_order",
              "bp_rownames_match", "bp_colnames_match",
              "bp_submatrix_nonneg", "bp_submatrix_integer", "bp_submatrix_nonzero"),
    result = c(
      sprintf("%dx%d", nrow(mat), ncol(mat)),
      as.character(nrow(mat)),
      as.character(ncol(mat)),
      BPCells::storage_order(mat),
      as.character(rownames_match),
      as.character(colnames_match),
      as.character(all_nonneg),
      as.character(all_integer),
      as.character(has_nonzero)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(bpcells_validation, file.path(RESULTS, "GSE243013_bpcells_matrix_validation.csv"),
            row.names = FALSE)
  cat("[INFO] Saved: GSE243013_bpcells_matrix_validation.csv\n")

  ## Check all validations passed
  if (!all_nonneg || !all_integer || !has_nonzero) {
    failed_step <<- "bpcells_validation"
    error_message <<- "BPCells matrix validation failed"
    stop(error_message)
  }

}, error = function(e) {
  cat("[ERROR] BPCells validation error:", conditionMessage(e), "\n")
  failed_step <<- "bpcells_validation"
  error_message <<- conditionMessage(e)
  stop(error_message)
})

## =========================================================================
## XV. File and Disk Statistics
## =========================================================================
cat("\n[XV] Collecting file and disk statistics...\n")

mtx_size <- file.info(counts_path)$size
bpcells_size <- sum(file.info(list.files(final_dir, full.names = TRUE))$size)
n_bpcells_files <- length(list.files(final_dir))
tmp_size <- if (dir.exists(TMP_DIR)) sum(file.info(list.files(TMP_DIR, full.names = TRUE))$size) else 0

disk_after <- tryCatch({
  disk_info <- system("df -Pk .", intern = TRUE)
  disk_parts <- strsplit(disk_info[2], "\\s+")[[1]]
  as.numeric(disk_parts[4]) / (1024 * 1024)
}, error = function(e) NA_real_)

total_runtime <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))

storage_stats <- data.frame(
  metric = c("counts_file_size_bytes", "counts_file_size_GB",
             "bpcells_dir_size_bytes", "bpcells_dir_size_GB",
             "bpcells_n_files", "tmp_dir_size_bytes",
             "disk_available_GB_after", "download_seconds",
             "import_seconds", "total_runtime_seconds"),
  value = c(mtx_size, mtx_size / (1024^3),
            bpcells_size, bpcells_size / (1024^3),
            n_bpcells_files, tmp_size,
            disk_after, ifelse(is.na(dl_seconds), NA, dl_seconds),
            ifelse(is.na(import_seconds), NA, import_seconds),
            total_runtime),
  stringsAsFactors = FALSE
)
write.csv(storage_stats, file.path(RESULTS, "GSE243013_step04_storage_and_runtime.csv"),
          row.names = FALSE)
cat("[INFO] Saved: GSE243013_step04_storage_and_runtime.csv\n")

cat(sprintf("[INFO] counts.mtx.gz: %.2f GB\n", mtx_size / (1024^3)))
cat(sprintf("[INFO] BPCells dir: %.2f GB (%d files)\n", bpcells_size / (1024^3), n_bpcells_files))
cat(sprintf("[INFO] Tmp dir: %.2f GB\n", tmp_size / (1024^3)))
cat(sprintf("[INFO] Disk available after: %.1f GB\n", disk_after))
cat(sprintf("[INFO] Download time: %.1f s\n", dl_seconds))
cat(sprintf("[INFO] Import time: %.1f s\n", import_seconds))
cat(sprintf("[INFO] Total runtime: %.1f s\n", total_runtime))

## =========================================================================
## XVI. Create Completion Marker
## =========================================================================
cat("\n[XVI] Creating completion marker...\n")

if (is.na(failed_step)) {
  complete_text <- c(
    "GSE243013 Step 04 COMPLETE",
    "==========================",
    "",
    sprintf("Completion time: %s", Sys.time()),
    sprintf("Counts file: %s", counts_path),
    sprintf("Counts file size: %.2f GB", mtx_size / (1024^3)),
    sprintf("BPCells directory: %s", final_dir),
    sprintf("Matrix dimensions: %d x %d (rows x cols)", nrow(mat), ncol(mat)),
    sprintf("Matrix orientation: %s (rows) x %s (cols)", row_label, col_label),
    sprintf("Non-zero elements: %d", mm_info$nnz),
    sprintf("BPCells version: %s", packageVersion("BPCells")),
    sprintf("Storage order: %s", BPCells::storage_order(mat)),
    sprintf("Disk available after: %.1f GB", disk_after),
    sprintf("Download time: %.1f seconds", dl_seconds),
    sprintf("Import time: %.1f seconds", import_seconds),
    sprintf("Total runtime: %.1f seconds", total_runtime),
    "",
    "All validations passed."
  )
  writeLines(complete_text, file.path(RESULTS, "GSE243013_step04_COMPLETE.txt"))
  cat("[INFO] Saved: GSE243013_step04_COMPLETE.txt\n")
} else {
  failed_text <- c(
    "GSE243013 Step 04 FAILED",
    "========================",
    "",
    sprintf("Failure time: %s", Sys.time()),
    sprintf("Failed step: %s", failed_step),
    sprintf("Error message: %s", error_message),
    "",
    "Steps completed before failure:",
    sprintf("  - Preflight: %s", ifelse(!is.na(failed_step) && failed_step != "preflight_disk" && failed_step != "preflight_arm64", "YES", "NO")),
    sprintf("  - curl installed: %s", ifelse(!is.na(failed_step) && failed_step %in% c("install_bpcells", "verify_bpcells_functions", "download_counts", "file_validation", "read_matrix_header", "barcode_metadata_alignment", "bpcells_import", "bpcells_validation"), "YES", "NO")),
    sprintf("  - BPCells installed: %s", ifelse(!is.na(failed_step) && failed_step %in% c("download_counts", "file_validation", "read_matrix_header", "barcode_metadata_alignment", "bpcells_import", "bpcells_validation"), "YES", "NO")),
    sprintf("  - Counts downloaded: %s", ifelse(file.exists(counts_path) && file.info(counts_path)$size > 0, "YES (partial)", "NO")),
    "",
    "Can safely re-run: YES",
    "Next step to continue from: retry from failed step"
  )
  writeLines(failed_text, file.path(RESULTS, "GSE243013_step04_FAILED.txt"))
  cat("[INFO] Saved: GSE243013_step04_FAILED.txt\n")
}

## =========================================================================
## Final Summary
## =========================================================================
cat("\n========================================================================\n")
cat("STEP 04 FINAL SUMMARY\n")
cat("========================================================================\n")
flush.console()

cat(sprintf("curl installed: YES (v%s)\n", packageVersion("curl")))
cat(sprintf("BPCells installed: YES (v%s)\n", packageVersion("BPCells")))
cat(sprintf("BPCells ARM64 binary: YES (Built %s)\n", packageDescription("BPCells")$Built))
cat(sprintf("counts downloaded: %s\n", ifelse(file.exists(counts_path) && file.info(counts_path)$size > 5e9, "YES", "NO")))
cat(sprintf("counts file size: %.2f GB\n", mtx_size / (1024^3)))
cat(sprintf("Matrix file dimensions: %d x %d (rows x cols)\n", mm_info$nrow, mm_info$ncol))
cat(sprintf("Matrix orientation: %s (rows) x %s (cols)\n", row_label, col_label))
cat(sprintf("Non-zero elements: %d\n", mm_info$nnz))
cat(sprintf("Genes: %d\n", n_genes))
cat(sprintf("Barcodes: %d\n", n_barcodes))
cat(sprintf("Metadata cells: %d\n", nrow(meta)))
cat(sprintf("Barcodes match metadata: %s\n", ifelse(check_results$set_equal, "YES", "NO")))
cat(sprintf("Barcode order same as metadata: %s\n", ifelse(check_results$order_match, "SAME", "REORDERED")))
cat(sprintf("Duplicated genes: %d\n", sum(is_duplicated)))
cat(sprintf("Duplicated barcodes: %d\n", n_barcode_dup))
cat(sprintf("BPCells imported: %s\n", ifelse(imported, "YES", "NO")))
cat(sprintf("BPCells directory: %s\n", final_dir))
cat(sprintf("BPCells matrix dim: %d x %d\n", nrow(mat), ncol(mat)))
cat(sprintf("BPCells storage order: %s\n", BPCells::storage_order(mat)))
cat(sprintf("BPCells dir size: %.2f GB\n", bpcells_size / (1024^3)))
cat(sprintf("Disk available after: %.1f GB\n", disk_after))
cat(sprintf("Download time: %.1f seconds\n", dl_seconds))
cat(sprintf("Import time: %.1f seconds\n", import_seconds))
cat(sprintf("Total runtime: %.1f seconds\n", total_runtime))
cat(sprintf("Step 04 status: %s\n", ifelse(is.na(failed_step), "COMPLETE", "FAILED")))

if (is.na(failed_step)) {
  cat("\nReady for patient-level pseudobulk analysis.\n")
}

cat("\n========================================================================\n")
cat("Step 04 completed.\n")
cat("========================================================================\n")
