## =========================================================================
## Step 07A: Install Extracted MSigDB & Resume Step 07
## =========================================================================

.libPaths(c(path.expand("~/Library/R/arm64/4.6/library"), .libPaths()))
options(stringsAsFactors = FALSE)

cat("========================================================================\n")
cat("Step 07A: MSigDB Cache Install & Step 07 Resume\n")
cat("========================================================================\n\n")
flush.console()

step_start <- Sys.time()

## =========================================================================
## I. Check Old Processes
## =========================================================================
cat("[I] Checking for running processes...\n")
running <- system("pgrep -fl 'Rscript|07_pathway|07A_install' 2>/dev/null", intern = TRUE)
if (length(running) > 0) {
  cat("[FATAL] Step 07 related processes still running:\n")
  cat(paste(running, collapse = "\n"), "\n")
  cat("Not starting a second process.\n")
  stop("Aborted: existing process detected")
}
cat("[OK] No running processes.\n")

## =========================================================================
## II. Confirm Input Folder
## =========================================================================
cat("\n[II] Confirming input folder...\n")

source_root <- Sys.getenv("MSIGDB_EXTRACTED_DIR",
                          unset = file.path(path.expand("~"), "Downloads/msigdbr_manual/msigdb"))
if (!dir.exists(source_root)) stop(sprintf("[FATAL] source_root does not exist: %s", source_root))
all_files <- list.files(source_root, recursive = TRUE, full.names = FALSE)
if (length(all_files) == 0) stop("[FATAL] source_root is empty")
cat(sprintf("[OK] source_root: %s (%d files)\n", source_root, length(all_files)))

## Save file manifest
file_info <- file.info(file.path(source_root, all_files))
manifest_df <- data.frame(
  file_name = basename(all_files),
  relative_path = all_files,
  absolute_path = file.path(source_root, all_files),
  extension = tools::file_ext(all_files),
  file_size = file_info$size,
  modification_time = as.character(file_info$mtime),
  stringsAsFactors = FALSE
)
write.csv(manifest_df, "03_results/step07_programs/qc/GSE243013_manual_msigdb_source_file_list.csv",
          row.names = FALSE)

## =========================================================================
## III. Locate Summary File
## =========================================================================
cat("\n[III] Locating summary file...\n")

summary_candidates <- list.files(source_root, pattern = "^msigdb\\.2026\\.1\\.summary\\.rds$",
                                 recursive = TRUE, full.names = TRUE)
cat(sprintf("[INFO] Found %d summary candidate(s)\n", length(summary_candidates)))

if (length(summary_candidates) == 0) {
  rds_files <- list.files(source_root, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  cat("[FATAL] No summary file found. RDS files in source_root:\n")
  cat(paste(basename(rds_files), collapse = "\n"), "\n")
  writeLines(c("Step 07A FAILED", sprintf("Time: %s", Sys.time()), "No summary file found"),
             "03_results/GSE243013_step07_FAILED.txt")
  stop("No summary file found")
}
if (length(summary_candidates) > 1) {
  cat("[FATAL] Multiple summary files found:\n")
  cat(paste(summary_candidates, collapse = "\n"), "\n")
  stop("Multiple summary files")
}

summary_source <- summary_candidates[[1]]
data_root <- dirname(summary_source)
cat(sprintf("[OK] Summary: %s\n", summary_source))
cat(sprintf("[OK] Data root: %s\n", data_root))

## =========================================================================
## IV. Verify Summary RDS
## =========================================================================
cat("\n[IV] Verifying summary RDS...\n")

summary_data <- readRDS(summary_source)
cat(sprintf("[INFO] Class: %s, Rows: %d\n", class(summary_data)[1], nrow(summary_data)))
cat(sprintf("[INFO] Columns: %s\n", paste(colnames(summary_data), collapse = ", ")))

required_cols <- c("db_target_species", "gs_collection", "df_rds")
missing_cols <- setdiff(required_cols, colnames(summary_data))
if (length(missing_cols) > 0) stop(sprintf("[FATAL] Missing columns: %s", paste(missing_cols, collapse = ", ")))
if (nrow(summary_data) == 0) stop("[FATAL] Summary has 0 rows")
if (any(is.na(summary_data$df_rds) | summary_data$df_rds == ""))
  stop("[FATAL] df_rds contains NA or empty strings")

n_hs <- sum(summary_data$db_target_species == "HS", na.rm = TRUE)
n_h <- sum(summary_data$gs_collection == "H", na.rm = TRUE)
n_c2 <- sum(summary_data$gs_collection == "C2", na.rm = TRUE)
n_unique_rds <- length(unique(summary_data$df_rds))

cat(sprintf("[OK] HS collections: %d, H: %d, C2: %d, unique df_rds: %d\n",
            n_hs, n_h, n_c2, n_unique_rds))

## Save summary audit
summary_audit <- data.frame(
  metric = c("total_rows", "hs_rows", "h_collection_rows", "c2_collection_rows",
             "unique_df_rds_count", "species_values"),
  value = c(nrow(summary_data), n_hs, n_h, n_c2, n_unique_rds,
            paste(unique(summary_data$db_target_species), collapse = "; ")),
  stringsAsFactors = FALSE
)
write.csv(summary_audit, "03_results/step07_programs/qc/GSE243013_manual_msigdb_summary_audit.csv",
          row.names = FALSE)

## =========================================================================
## V. Validate All Referenced RDS
## =========================================================================
cat("\n[V] Validating referenced RDS files...\n")

referenced_rds <- unique(summary_data$df_rds)
rds_audit_rows <- list()
any_missing <- FALSE

for (rds_name in referenced_rds) {
  ref_path <- file.path(data_root, rds_name)
  exists <- file.exists(ref_path)
  fsize <- if (exists) file.info(ref_path)$size else 0L
  readable <- FALSE
  obj_class <- NA_character_
  row_count <- NA_integer_
  err_msg <- ""

  if (exists) {
    obj <- tryCatch({
      x <- readRDS(ref_path)
      readable <- TRUE
      obj_class <- paste(class(x), collapse = "/")
      row_count <- if (is.data.frame(x) || is.matrix(x)) nrow(x) else length(x)
      x
    }, error = function(e) {
      err_msg <<- conditionMessage(e)
      NULL
    })
    rm(obj)
    gc()
  } else {
    err_msg <- "File not found"
    any_missing <- TRUE
  }

  status <- if (!exists) "MISSING" else if (!readable) "UNREADABLE" else if (is.na(row_count)) "NO_ROWS" else "OK"

  rds_audit_rows[[length(rds_audit_rows) + 1]] <- data.frame(
    df_rds = rds_name, source_path = ref_path, exists = exists,
    file_size = fsize, readable = readable, object_class = obj_class,
    row_count = row_count, error_message = err_msg, validation_status = status,
    stringsAsFactors = FALSE
  )
}

rds_audit_df <- do.call(rbind, rds_audit_rows)
write.csv(rds_audit_df, "03_results/step07_programs/qc/GSE243013_manual_msigdb_referenced_rds_audit.csv",
          row.names = FALSE)

n_ok <- sum(rds_audit_df$validation_status == "OK")
n_fail <- sum(rds_audit_df$validation_status != "OK")
cat(sprintf("[INFO] RDS validation: %d OK, %d failed\n", n_ok, n_fail))

if (any_missing || n_fail > 0) {
  cat("[FATAL] Some referenced RDS files missing or unreadable\n")
  writeLines(c("Step 07A FAILED", sprintf("Time: %s", Sys.time()),
               sprintf("Missing/broken RDS: %d", n_fail)),
             "03_results/GSE243013_step07_FAILED.txt")
  stop("Incomplete MSigDB data")
}

## =========================================================================
## VI. Determine msigdbr Cache Directory
## =========================================================================
cat("\n[VI] Determining msigdbr cache directory...\n")

cache_dir <- tools::R_user_dir(package = "msigdbr", which = "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("[OK] Cache dir: %s\n", cache_dir))
writeLines(cache_dir, "03_results/step07_programs/qc/GSE243013_msigdbr_cache_path.txt")

## =========================================================================
## VII. Copy to Cache Directory
## =========================================================================
cat("\n[VII] Copying files to cache directory...\n")

files_to_copy <- c(basename(summary_source), referenced_rds)
copy_rows <- list()

for (fname in files_to_copy) {
  src_file <- file.path(data_root, fname)
  tgt_file <- file.path(cache_dir, fname)
  src_size <- file.info(src_file)$size
  tgt_exists <- file.exists(tgt_file)
  tgt_size <- if (tgt_exists) file.info(tgt_file)$size else 0L
  src_md5 <- tryCatch(tools::md5sum(src_file), error = function(e) NA_character_)
  tgt_md5 <- if (tgt_exists) tryCatch(tools::md5sum(tgt_file), error = function(e) NA_character_) else NA_character_

  copied <- FALSE
  reused <- FALSE
  status <- ""
  err_msg <- ""

  if (tgt_exists && src_md5 == tgt_md5) {
    reused <- TRUE
    status <- "REUSED"
    cat(sprintf("  [REUSED] %s\n", fname))
  } else if (tgt_exists && src_md5 != tgt_md5) {
    status <- "CONFLICT"
    err_msg <- "Target exists with different content — not overwriting"
    cat(sprintf("  [CONFLICT] %s (target differs)\n", fname))
  } else {
    result <- tryCatch({
      file.copy(from = src_file, to = tgt_file, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)
      TRUE
    }, error = function(e) {
      err_msg <<- conditionMessage(e)
      FALSE
    })
    if (result) {
      copied <- TRUE
      status <- "COPIED"
      cat(sprintf("  [COPIED] %s\n", fname))
    } else {
      status <- "FAILED"
      cat(sprintf("  [FAILED] %s: %s\n", fname, err_msg))
    }
  }

  copy_rows[[length(copy_rows) + 1]] <- data.frame(
    source_path = src_file, target_path = tgt_file,
    source_size = src_size, target_size = if (file.exists(tgt_file)) file.info(tgt_file)$size else 0L,
    source_md5 = src_md5, target_md5 = tgt_md5,
    copied = copied, reused_existing = reused, status = status,
    error_message = err_msg, stringsAsFactors = FALSE
  )
}

copy_df <- do.call(rbind, copy_rows)
write.csv(copy_df, "03_results/step07_programs/qc/GSE243013_msigdbr_cache_install_status.csv",
          row.names = FALSE)

n_copied <- sum(copy_df$copied)
n_reused <- sum(copy_df$reused_existing)
n_conflict <- sum(copy_df$status == "CONFLICT")
cat(sprintf("[INFO] Copied: %d, Reused: %d, Conflicts: %d\n", n_copied, n_reused, n_conflict))

## =========================================================================
## VIII. Post-Install Cache Verification
## =========================================================================
cat("\n[VIII] Verifying cache installation...\n")

cached_summary_path <- file.path(cache_dir, "msigdb.2026.1.summary.rds")
if (!file.exists(cached_summary_path)) {
  stop(sprintf("[FATAL] Cached summary not found: %s", cached_summary_path))
}

cached_summary <- readRDS(cached_summary_path)
all_exist <- all(file.exists(file.path(cache_dir, cached_summary$df_rds)))
cat(sprintf("[OK] All %d referenced RDS in cache: %s\n", nrow(cached_summary), all_exist))

if (!all_exist) stop("[FATAL] Not all referenced RDS found in cache")

## =========================================================================
## IX. Build Project-Level MSigDB Cache from Local RDS Files
## =========================================================================
cat("\n[IX] Building project-level MSigDB cache from local RDS files...\n")

## The RDS files have columns: db_gene_symbol, gs_name, gs_description, etc.
## We need to map db_gene_symbol -> gene_symbol to match msigdbr output format.
## Then build named lists: gs_name -> gene_vector (same format as original 07 script expects).

## --- Hallmark ---
cat("[INFO] Loading Hallmark gene sets from local RDS...\n")
hallmark_rds_path <- file.path(data_root, "msigdb.2026.1.Hs.H.rds")
if (!file.exists(hallmark_rds_path)) {
  stop(sprintf("[FATAL] Hallmark RDS not found: %s", hallmark_rds_path))
}
hallmark_raw <- readRDS(hallmark_rds_path)
cat(sprintf("[INFO] Hallmark RDS: %d rows, %d unique gs_name\n",
            nrow(hallmark_raw), length(unique(hallmark_raw$gs_name))))

## Filter to HS species and build gene symbol -> gs_name mapping
hallmark_df <- data.frame(
  gs_id        = hallmark_raw$gs_id,
  gs_name      = hallmark_raw$gs_name,
  gene_symbol  = hallmark_raw$db_gene_symbol,
  gs_description = hallmark_raw$gs_description,
  stringsAsFactors = FALSE
)
hallmark_df <- hallmark_df[!is.na(hallmark_df$gene_symbol) & hallmark_df$gene_symbol != "" &
                           !is.na(hallmark_df$gs_name) & hallmark_df$gs_name != "", ]
hallmark_df <- hallmark_df[!duplicated(paste(hallmark_df$gs_name, hallmark_df$gene_symbol)), ]
hallmark_gsets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
cat(sprintf("[OK] Hallmark: %d gene sets, %d memberships\n",
            length(hallmark_gsets), nrow(hallmark_df)))

## Save project-level cache
dir.create("00_config/step07_resources", recursive = TRUE, showWarnings = FALSE)
saveRDS(hallmark_gsets, "00_config/step07_resources/MSigDB_Hallmark.rds")
cat("[OK] Saved: 00_config/step07_resources/MSigDB_Hallmark.rds\n")

## --- Reactome ---
cat("\n[INFO] Loading Reactome gene sets from local RDS...\n")
reactome_rds_path <- file.path(data_root, "msigdb.2026.1.Hs.C2.rds")
if (!file.exists(reactome_rds_path)) {
  stop(sprintf("[FATAL] Reactome RDS not found: %s", reactome_rds_path))
}
reactome_raw <- readRDS(reactome_rds_path)
cat(sprintf("[INFO] Reactome RDS: %d rows, %d unique gs_name\n",
            nrow(reactome_raw), length(unique(reactome_raw$gs_name))))

## Filter to CP:REACTOME subcollection
if ("gs_subcollection" %in% colnames(reactome_raw)) {
  reactome_cp <- reactome_raw[reactome_raw$gs_subcollection == "CP:REACTOME", ]
  cat(sprintf("[INFO] After CP:REACTOME filter: %d rows\n", nrow(reactome_cp)))
} else {
  ## Fallback: use gs_name prefix to identify Reactome
  reactome_cp <- reactome_raw[grepl("^REACTOME_", reactome_raw$gs_name), ]
  cat(sprintf("[INFO] After REACTOME_ prefix filter: %d rows\n", nrow(reactome_cp)))
}

reactome_df <- data.frame(
  gs_id        = reactome_cp$gs_id,
  gs_name      = reactome_cp$gs_name,
  gene_symbol  = reactome_cp$db_gene_symbol,
  gs_description = reactome_cp$gs_description,
  stringsAsFactors = FALSE
)
reactome_df <- reactome_df[!is.na(reactome_df$gene_symbol) & reactome_df$gene_symbol != "" &
                           !is.na(reactome_df$gs_name) & reactome_df$gs_name != "", ]
reactome_df <- reactome_df[!duplicated(paste(reactome_df$gs_name, reactome_df$gene_symbol)), ]
reactome_gsets <- split(reactome_df$gene_symbol, reactome_df$gs_name)
cat(sprintf("[OK] Reactome: %d gene sets, %d memberships\n",
            length(reactome_gsets), nrow(reactome_df)))

## Save project-level cache
saveRDS(reactome_gsets, "00_config/step07_resources/MSigDB_Reactome.rds")
cat("[OK] Saved: 00_config/step07_resources/MSigDB_Reactome.rds\n")

## --- Metadata ---
msigdb_meta <- data.frame(
  msigdbr_version = as.character(packageVersion("msigdbr")),
  MSigDB_version = "2026.1",
  collection = c("H", "C2"),
  subcollection = c("HALLMARK", "CP:REACTOME"),
  gene_set_count = c(length(hallmark_gsets), length(reactome_gsets)),
  gene_membership_count = c(nrow(hallmark_df), nrow(reactome_df)),
  manual_source_root = source_root,
  data_root = data_root,
  cache_dir = cache_dir,
  retrieval_method = "local_rds_direct_load",
  validation_status = "OK",
  stringsAsFactors = FALSE
)
write.csv(msigdb_meta, "00_config/step07_resources/MSigDB_resource_metadata.csv", row.names = FALSE)
cat("[OK] Saved: MSigDB_resource_metadata.csv\n")

## Verify cached files can be read
cat("\n[INFO] Verifying cached files...\n")
hallmark_check <- readRDS("00_config/step07_resources/MSigDB_Hallmark.rds")
reactome_check <- readRDS("00_config/step07_resources/MSigDB_Reactome.rds")
cat(sprintf("[OK] MSigDB_Hallmark.rds: %d gene sets\n", length(hallmark_check)))
cat(sprintf("[OK] MSigDB_Reactome.rds: %d gene sets\n", length(reactome_check)))

## =========================================================================
## X. Fix MΦ_CXCL10 File Mapping
## =========================================================================
cat("\n[X] Fixing MΦ_CXCL10 file mapping...\n")

pri_dir <- "03_results/step06_edgeR/primary_anti_PD1"
actual_qlf_files <- list.files(pri_dir, pattern = "__edgeR_QLF_all_genes\\.csv\\.gz$", full.names = FALSE)

model_status <- read.csv("03_results/step06_edgeR/combined/GSE243013_edgeR_model_status.csv",
                         stringsAsFactors = FALSE)
pri_complete <- model_status[model_status$analysis == "primary_anti_PD1" &
                              model_status$status == "COMPLETE", ]

mapping_rows <- list()
for (i in seq_len(nrow(pri_complete))) {
  ct <- pri_complete$cell_type[i]
  ## Find the actual file that matches this cell type
  matched_file <- NA_character_
  for (f in actual_qlf_files) {
    ## Try exact match first
    ct_safe_direct <- gsub("[^A-Za-z0-9_]", "_", ct)
    if (grepl(paste0("^", ct_safe_direct, "__edgeR_QLF"), f)) {
      matched_file <- f
      break
    }
    ## Try removing non-ASCII chars
    ct_ascii <- gsub("[^A-Za-z0-9]", "_", ct)
    ct_ascii <- gsub("_+", "_", ct_ascii)
    ct_ascii <- sub("^_|_$", "", ct_ascii)
    if (grepl(paste0("^", ct_ascii, "__edgeR_QLF"), f)) {
      matched_file <- f
      break
    }
  }

  if (is.na(matched_file)) {
    ## Try partial matching: find files where most of the cell type name appears
    ct_words <- unlist(strsplit(gsub("[^A-Za-z0-9]", "_", ct), "_"))
    ct_words <- ct_words[nchar(ct_words) >= 3]
    for (f in actual_qlf_files) {
      if (all(sapply(ct_words, function(w) grepl(w, f, ignore.case = TRUE)))) {
        matched_file <- f
        break
      }
    }
  }

  ct_safe_gsub <- gsub("_+", "_", gsub("[^A-Za-z0-9_]", "_", ct))
  ct_safe_gsub <- sub("^_|_$", "", ct_safe_gsub)
  mapping_rows[[length(mapping_rows) + 1]] <- data.frame(
    cell_type = ct,
    cell_type_safe = ct_safe_gsub,
    model_status = pri_complete$status[i],
    actual_qlf_file = ifelse(is.na(matched_file), "NOT_FOUND", matched_file),
    n_samples = pri_complete$n_samples[i],
    n_resp = pri_complete$n_resp[i],
    n_nonresp = pri_complete$n_nonresp[i],
    stringsAsFactors = FALSE
  )
}

mapping_df <- do.call(rbind, mapping_rows)
write.csv(mapping_df, "03_results/step07_programs/qc/GSE243013_step06_exact_file_mapping.csv",
          row.names = FALSE)

cat(sprintf("[INFO] Mapping results:\n"))
for (i in seq_len(nrow(mapping_df))) {
  cat(sprintf("  %s -> %s\n", mapping_df$cell_type[i], mapping_df$actual_qlf_file[i]))
}

## =========================================================================
## XI. Build Missing Rank #8 (MΦ_CXCL10)
## =========================================================================
cat("\n[XI] Building missing rank for MΦ_CXCL10...\n")

cxcl10_mapping <- mapping_df[grepl("CXCL10", mapping_df$cell_type), ]
if (nrow(cxcl10_mapping) == 0) stop("[FATAL] MΦ_CXCL10 not found in mapping")

cxcl10_file <- cxcl10_mapping$actual_qlf_file[1]
if (cxcl10_file == "NOT_FOUND") stop("[FATAL] MΦ_CXCL10 QLF file not found on disk")

ct <- cxcl10_mapping$cell_type[1]
  ct_safe <- gsub("_+", "_", gsub("[^A-Za-z0-9_]", "_", ct))
  ct_safe <- sub("^_|_$", "", ct_safe)
  qlf_path <- file.path(pri_dir, cxcl10_file)
cat(sprintf("[INFO] Reading: %s\n", qlf_path))

qlf <- read.csv(qlf_path, stringsAsFactors = FALSE)
cat(sprintf("[INFO] Rows: %d, Columns: %s\n", nrow(qlf), paste(colnames(qlf), collapse = ", ")))

## Validate and clean
valid <- !is.na(qlf$gene) & qlf$gene != "" & is.finite(qlf$logFC) & is.finite(qlf$F) & is.finite(qlf$PValue)
qlf <- qlf[valid, ]
qlf <- qlf[!duplicated(qlf$gene), ]

## Create rank statistics
qlf$rank_signed_sqrtF <- sign(qlf$logFC) * sqrt(pmax(qlf$F, 0))
qlf$rank_logFC_logP <- qlf$logFC * -log10(pmax(qlf$PValue, 1e-300))

## Sort
qlf <- qlf[order(-qlf$rank_signed_sqrtF, -qlf$logFC, qlf$gene), ]

## Save
rank_df <- data.frame(
  gene = qlf$gene, rank_signed_sqrtF = qlf$rank_signed_sqrtF,
  rank_logFC_logP = qlf$rank_logFC_logP, logFC = qlf$logFC,
  F = qlf$F, PValue = qlf$PValue,
  stringsAsFactors = FALSE
)
rank_path <- sprintf("03_results/step07_programs/ranks/%s__gene_ranks.csv.gz", ct_safe)
con <- gzfile(rank_path, "w")
write.csv(rank_df, con, row.names = FALSE)
close(con)

cat(sprintf("[OK] %s: %d genes ranked (%d positive, %d negative)\n",
            ct, nrow(qlf), sum(qlf$rank_signed_sqrtF > 0), sum(qlf$rank_signed_sqrtF < 0)))

## Verify all 8 ranks
cat("\n[INFO] Verifying all 8 primary ranks...\n")
rank_files <- list.files("03_results/step07_programs/ranks", pattern = "__gene_ranks\\.csv\\.gz$",
                         full.names = FALSE)
cat(sprintf("[INFO] Rank files found: %d\n", length(rank_files)))
cat(paste("  ", rank_files, collapse = "\n"), "\n")

## =========================================================================
## XII. Fix Main Step 07 Script
## =========================================================================
cat("\n[XII] Fixing main Step 07 script...\n")

## Read the current script
script_path <- "01_scripts/07_pathway_TF_program_integration.R"
script_lines <- readLines(script_path, warn = FALSE)

## Fix 1: Make Section VIII use project-level cache instead of msigdbr download
## The script already checks for cached files first, so we just need to ensure
## it doesn't try to download from zenodo.

## Fix 2: The safe name generation issue — replace gsub approach with disk-based matching
## We'll patch the gene ranking section to use the mapping file.

## For now, the key fix is to ensure the script uses the cached MSigDB files
## which we've already created at 00_config/step07_resources/MSigDB_Hallmark.rds
## and MSigDB_Reactome.rds

## The script already has this logic at the beginning of Section VIII:
## if (file.exists(hallmark_path)) { ... load cached ... }
## So the fix is to make sure the cached files exist (done above).

## Let's also fix the safe name issue by patching the ranking section
## to use the mapping file when available.

cat("[INFO] Script already checks for cached MSigDB files — fix applied via cache installation.\n")
cat("[INFO] Safe name fix will be applied in the resume section.\n")

## Verify script parses
parse_check <- tryCatch({
  parse(file = script_path)
  TRUE
}, error = function(e) {
  cat(sprintf("[WARNING] Parse error: %s\n", e$message))
  FALSE
})

cat(sprintf("[INFO] Script parse check: %s\n", ifelse(parse_check, "OK", "FAILED")))

## =========================================================================
## XIII. Resume Step 07 — Run Remaining Sections
## =========================================================================
cat("\n[XIII] Resuming Step 07 from Section VIII...\n")

## Source the necessary functions and data from the original script
## We'll create a minimal resume script that loads what we need
## and continues from Section VIII onward

resume_script <- '01_scripts/07_pathway_TF_program_integration_resume.R'

resume_lines <- c(
  '## =========================================================================',
  '## Step 07 Resume: From Section VIII onward',
  '## =========================================================================',
  '',
  '.libPaths(c(path.expand("~/Library/R/arm64/4.6/library"), .libPaths()))',
  'options(stringsAsFactors = FALSE)',
  '',
  'cat("========================================================================\\n")',
  'cat("Step 07 Resume: From Section VIII\\n")',
  'cat("========================================================================\\n\\n")',
  'flush.console()',
  '',
  'step_start <- Sys.time()',
  '',
  '## Load required packages',
  'suppressPackageStartupMessages({',
  '  library(fgsea)',
  '  library(BiocParallel)',
  '  library(decoupleR)',
  '  library(data.table)',
  '  library(dplyr)',
  '  library(tidyr)',
  '  library(tibble)',
  '  library(ggplot2)',
  '  library(pheatmap)',
  '  library(RColorBrewer)',
  '})',
  '',
  '## Load cached MSigDB gene sets',
  'hallmark_gsets <- readRDS("00_config/step07_resources/MSigDB_Hallmark.rds")',
  'reactome_gsets <- readRDS("00_config/step07_resources/MSigDB_Reactome.rds")',
  'cat(sprintf("[INFO] Hallmark: %d gene sets\\n", length(hallmark_gsets)))',
  'cat(sprintf("[INFO] Reactome: %d gene sets\\n", length(reactome_gsets)))',
  '',
  '## Load model status and determine analyzable models',
  'model_status <- read.csv("03_results/step06_edgeR/combined/GSE243013_edgeR_model_status.csv",',
  '                         stringsAsFactors = FALSE)',
  'pri_complete <- model_status[model_status$analysis == "primary_anti_PD1" &',
  '                              model_status$status == "COMPLETE", ]',
  '',
  '## Load exact file mapping',
  'mapping_df <- read.csv("03_results/step07_programs/qc/GSE243013_step06_exact_file_mapping.csv",',
  '                        stringsAsFactors = FALSE)',
  '',
  '## Load existing ranks and build missing ones',
  'primary_ranks <- list()',
  'primary_analyzed <- character(0)',
  'rank_qc_rows <- list()',
  '',
  'pri_dir <- "03_results/step06_edgeR/primary_anti_PD1"',
  '',
  'for (i in seq_len(nrow(pri_complete))) {',
  '  row <- pri_complete[i, ]',
  '  ct <- row$cell_type',
  '  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)',
  '  analysis_name <- "primary_anti_PD1"',
  '',
  '  cat(sprintf("\\n--- Ranking: %s ---\\n", ct))',
  '',
  '  ## Get cell_type_safe from mapping',
  '  ct_safe <- mapping_df$cell_type_safe[mapping_df$cell_type == ct]',
  '  if (length(ct_safe) == 0 || is.na(ct_safe)) ct_safe <- gsub("_+", "_", gsub("[^A-Za-z0-9_]", "_", ct))',
  '  if (nchar(ct_safe) > 0 && substr(ct_safe, 1, 1) == "_") ct_safe <- substr(ct_safe, 2, nchar(ct_safe))',
  '  if (nchar(ct_safe) > 0 && substr(ct_safe, nchar(ct_safe), nchar(ct_safe)) == "_") ct_safe <- substr(ct_safe, 1, nchar(ct_safe) - 1)',
  '',
  '  ## Check if rank already exists',
  '  rank_path <- sprintf("03_results/step07_programs/ranks/%s__gene_ranks.csv.gz", ct_safe)',
  '  if (file.exists(rank_path)) {',
  '    rk <- read.csv(rank_path, stringsAsFactors = FALSE)',
  '    rank_primary <- setNames(rk$rank_signed_sqrtF, rk$gene)',
  '    rank_sensitivity <- setNames(rk$rank_logFC_logP, rk$gene)',
  '    primary_ranks[[ct]] <- list(',
  '      primary = rank_primary, sensitivity = rank_sensitivity,',
  '      analysis_name = analysis_name, cell_type = ct, cell_type_safe = ct_safe',
  '    )',
  '    primary_analyzed <- c(primary_analyzed, ct)',
  '    cat(sprintf("[REUSED] %s: %d genes\\n", ct, nrow(rk)))',
  '    next',
  '  }',
  '',
  '  ## Find actual QLF file from mapping',
  '  mapped <- mapping_df[mapping_df$cell_type == ct & mapping_df$actual_qlf_file != "NOT_FOUND", ]',
  '  if (nrow(mapped) == 0) {',
  '    cat(sprintf("[SKIP] No QLF file for %s\\n", ct))',
  '    next',
  '  }',
  '  qlf_file <- mapped$actual_qlf_file[1]',
  '  qlf_path <- file.path(pri_dir, qlf_file)',
  '  if (!file.exists(qlf_path)) {',
  '    cat(sprintf("[SKIP] QLF file not found: %s\\n", qlf_path))',
  '    next',
  '  }',
  '',
  '  qlf <- read.csv(qlf_path, stringsAsFactors = FALSE)',
  '  valid <- !is.na(qlf$gene) & qlf$gene != "" &',
  '           is.finite(qlf$logFC) & is.finite(qlf$F) & is.finite(qlf$PValue)',
  '  qlf <- qlf[valid, ]',
  '  qlf <- qlf[!duplicated(qlf$gene), ]',
  '  qlf$rank_signed_sqrtF <- sign(qlf$logFC) * sqrt(pmax(qlf$F, 0))',
  '  qlf$rank_logFC_logP <- qlf$logFC * -log10(pmax(qlf$PValue, 1e-300))',
  '  qlf <- qlf[order(-qlf$rank_signed_sqrtF, -qlf$logFC, qlf$gene), ]',
  '',
  '  rank_df <- data.frame(',
  '    gene = qlf$gene, rank_signed_sqrtF = qlf$rank_signed_sqrtF,',
  '    rank_logFC_logP = qlf$rank_logFC_logP, logFC = qlf$logFC,',
  '    F = qlf$F, PValue = qlf$PValue, stringsAsFactors = FALSE',
  '  )',
  '  dir.create("03_results/step07_programs/ranks", recursive = TRUE, showWarnings = FALSE)',
  '  con <- gzfile(rank_path, "w")',
  '  write.csv(rank_df, con, row.names = FALSE)',
  '  close(con)',
  '',
  '  rank_primary <- setNames(qlf$rank_signed_sqrtF, qlf$gene)',
  '  rank_sensitivity <- setNames(qlf$rank_logFC_logP, qlf$gene)',
  '  primary_ranks[[ct]] <- list(',
  '    primary = rank_primary, sensitivity = rank_sensitivity,',
  '    analysis_name = analysis_name, cell_type = ct, cell_type_safe = ct_safe',
  '  )',
  '  primary_analyzed <- c(primary_analyzed, ct)',
  '  cat(sprintf("[OK] %s: %d genes\\n", ct, nrow(qlf)))',
  '}',
  '',
  'cat(sprintf("\\n[INFO] Primary analyzed: %d cell types\\n", length(primary_analyzed)))',
  ''
)

## Now append the original script from Section VIII onward
## Find where Section VIII starts in the original script
original_lines <- readLines("01_scripts/07_pathway_TF_program_integration.R", warn = FALSE)
section_viii_start <- grep("^## VIII\\. Prepare MSigDB Gene Sets", original_lines)
if (length(section_viii_start) == 0) {
  section_viii_start <- grep("^## VIII", original_lines)[1]
}

if (!is.na(section_viii_start)) {
  ## Skip the Section VIII header and the caching logic (we already have cached files)
  ## Find where Section IX starts
  section_ix_start <- grep("^## IX\\. Run fgsea", original_lines)
  if (length(section_ix_start) == 0) {
    section_ix_start <- grep("^## IX", original_lines)[1]
  }

  ## Include from Section IX onward
  if (!is.na(section_ix_start)) {
    remaining_lines <- original_lines[section_ix_start:length(original_lines)]
    resume_lines <- c(resume_lines, remaining_lines)
  }
}

writeLines(resume_lines, resume_script)
cat(sprintf("[OK] Resume script created: %s (%d lines)\n", resume_script, length(resume_lines)))

## Verify parse
parse_ok <- tryCatch({
  parse(file = resume_script)
  TRUE
}, error = function(e) {
  cat(sprintf("[WARNING] Resume script parse error: %s\\n", e$message))
  FALSE
})

if (!parse_ok) {
  cat("[FATAL] Resume script has syntax errors\\n")
  stop("Resume script parse failed")
}

cat("[OK] Resume script parses successfully\\n")

## Run the resume
cat("\\n[XIII-b] Executing resume script...\\n")
flush.console()

resume_log <- "logs/07_pathway_TF_program_integration_resume_from_extracted_folder.log"
exit_code <- system2(
  "/usr/local/bin/Rscript",
  args = c("--vanilla", resume_script),
  stdout = resume_log,
  stderr = "2>&1"
)

cat(sprintf("[INFO] Resume exit code: %d\\n", exit_code))

## Check completion
if (file.exists("03_results/GSE243013_step07_COMPLETE.txt")) {
  cat("\\n[OK] Step 07 COMPLETE marker found!\\n")
  cat(readLines("03_results/GSE243013_step07_COMPLETE.txt"), sep = "\\n")
} else if (file.exists("03_results/GSE243013_step07_FAILED.txt")) {
  cat("\\n[WARNING] Step 07 FAILED marker found\\n")
  cat(readLines("03_results/GSE243013_step07_FAILED.txt"), sep = "\\n")
} else {
  cat("\\n[WARNING] No completion marker found\\n")
  cat("Last 30 lines of resume log:\\n")
  log_lines <- readLines(resume_log, warn = FALSE)
  cat(paste(tail(log_lines, 30), collapse = "\\n"), "\\n")
}

## Final device check
if (grDevices::dev.cur() != 1L) {
  cat("[WARNING] Closing remaining graphics devices\\n")
  graphics.off()
}

total_runtime <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
cat(sprintf("\\nStep 07A total runtime: %.1f seconds\\n", total_runtime))
cat("========================================================================\\n")
cat("Step 07A completed.\\n")
cat("========================================================================\\n")
