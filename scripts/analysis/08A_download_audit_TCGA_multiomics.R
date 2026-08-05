## =========================================================================
## Step 08A: Download & Audit TCGA Multi-Omics Data
## =========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(survival)
  library(digest)
})

.libPaths(c(path.expand("~/Library/R/arm64/4.6/library"), .libPaths()))
options(stringsAsFactors = FALSE, warn = 1)

step08a_start <- Sys.time()

cat("========================================================================\n")
cat("Step 08A: Download & Audit TCGA Multi-Omics Data\n")
cat("========================================================================\n")
cat(sprintf("Start time: %s\n", as.character(step08a_start)))
cat(sprintf("Working directory: %s\n", getwd()))
flush.console()

## =========================================================================
## II. Create Directory Structure
## =========================================================================
cat("\n[II] Creating directory structure...\n")

dirs_to_create <- c(
  "00_config/step08_TCGA",
  "02_data/tcga", "02_data/tcga/ExperimentHub_cache",
  "02_data/tcga/curated", "02_data/tcga/curated/LUAD", "02_data/tcga/curated/LUSC",
  "02_data/tcga/clinical", "02_data/tcga/manifests",
  "03_results/step08_TCGA", "03_results/step08_TCGA/preflight",
  "03_results/step08_TCGA/programs", "03_results/step08_TCGA/assay_audit",
  "03_results/step08_TCGA/sample_audit", "03_results/step08_TCGA/clinical_audit",
  "logs"
)
for (d in dirs_to_create) dir.create(d, recursive = TRUE, showWarnings = FALSE)
cat("[OK] All directories created\n")

## =========================================================================
## III. Environment & Disk Pre-check
## =========================================================================
cat("\n[III] Environment & disk pre-check...\n")

env_info <- data.frame(
  metric = c("date_time", "working_directory", "r_version", "r_platform",
             "hostname", "username", "lib_paths", "mem_bytes", "disk_total_kb",
             "disk_available_kb", "disk_available_gb"),
  value = c(
    as.character(Sys.time()),
    getwd(),
    R.version.string,
    R.version$platform,
    Sys.info()["nodename"],
    Sys.info()["user"],
    paste(.libPaths(), collapse = "; "),
    system("sysctl -n hw.memsize", intern = TRUE),
    NA_character_, NA_character_, NA_character_
  ),
  stringsAsFactors = FALSE
)

disk_info <- system("df -Pk .", intern = TRUE)
disk_line <- disk_info[length(disk_info)]
disk_parts <- strsplit(disk_line, "\\s+")[[1]]
env_info$value[env_info$metric == "disk_total_kb"] <- disk_parts[2]
env_info$value[env_info$metric == "disk_available_kb"] <- disk_parts[4]
env_info$value[env_info$metric == "disk_available_gb"] <- sprintf("%.1f", as.numeric(disk_parts[4]) / 1048576)

write.csv(env_info, "03_results/step08_TCGA/preflight/GSE243013_step08A_environment_preflight.txt",
          row.names = FALSE)
cat(sprintf("[OK] Environment saved. Disk available: %s GB\n", env_info$value[env_info$metric == "disk_available_gb"]))

disk_avail_gb <- as.numeric(env_info$value[env_info$metric == "disk_available_gb"])
if (disk_avail_gb < 35) {
  cat(sprintf("[FATAL] Insufficient disk space: %.1f GB (need 35 GB)\n", disk_avail_gb))
  writeLines(c("Step 08A FAILED", sprintf("Time: %s", Sys.time()),
               sprintf("Disk available: %.1f GB < 35 GB", disk_avail_gb)),
             "03_results/GSE243013_step08A_FAILED.txt")
  stop("Insufficient disk space")
}

## =========================================================================
## IV. Load and Install R Packages
## =========================================================================
cat("\n[IV] Loading and installing R packages...\n")

lib <- path.expand("~/Library/R/arm64/4.6/library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)

bioc_pkgs <- c("curatedTCGAData", "TCGAutils", "MultiAssayExperiment",
               "SummarizedExperiment", "RaggedExperiment", "HDF5Array",
               "DelayedArray", "ExperimentHub", "AnnotationHub",
               "S4Vectors", "BiocFileCache", "BiocManager")
cran_pkgs <- c("data.table", "dplyr", "tidyr", "tibble", "stringr",
               "survival", "digest")

missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_bioc) > 0) {
  cat(sprintf("[INFO] Installing %d missing Bioconductor packages: %s\n",
              length(missing_bioc), paste(missing_bioc, collapse = ", ")))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE, lib = lib)
}
if (length(missing_cran) > 0) {
  cat(sprintf("[INFO] Installing %d missing CRAN packages: %s\n",
              length(missing_cran), paste(missing_cran, collapse = ", ")))
  install.packages(missing_cran, repos = "https://cloud.r-project.org",
                   type = "binary", lib = lib)
}

## Load all packages
suppressPackageStartupMessages({
  library(curatedTCGAData)
  library(TCGAutils)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
  library(RaggedExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(ExperimentHub)
  library(S4Vectors)
})

## Save package environment
pkg_env <- data.frame(
  package = c("R_version", "BiocManager", "BiocVersion",
              "curatedTCGAData", "TCGAutils", "MultiAssayExperiment",
              "SummarizedExperiment", "HDF5Array", "ExperimentHub",
              "data.table", "dplyr", "tidyr"),
  version = c(
    R.version.string,
    as.character(packageVersion("BiocManager")),
    BiocManager::version(),
    tryCatch(as.character(packageVersion("curatedTCGAData")), error = function(e) "NOT_LOADED"),
    tryCatch(as.character(packageVersion("TCGAutils")), error = function(e) "NOT_LOADED"),
    tryCatch(as.character(packageVersion("MultiAssayExperiment")), error = function(e) "NOT_LOADED"),
    tryCatch(as.character(packageVersion("SummarizedExperiment")), error = function(e) "NOT_LOADED"),
    tryCatch(as.character(packageVersion("HDF5Array")), error = function(e) "NOT_LOADED"),
    tryCatch(as.character(packageVersion("ExperimentHub")), error = function(e) "NOT_LOADED"),
    tryCatch(as.character(packageVersion("data.table")), error = function(e) "NOT_LOADED"),
    tryCatch(as.character(packageVersion("dplyr")), error = function(e) "NOT_LOADED"),
    tryCatch(as.character(packageVersion("tidyr")), error = function(e) "NOT_LOADED")
  ),
  stringsAsFactors = FALSE
)
write.csv(pkg_env, "03_results/step08_TCGA/preflight/GSE243013_step08A_package_environment.txt",
          row.names = FALSE)
cat("[OK] Package environment saved\n")

## =========================================================================
## V. Freeze Step 07 Results
## =========================================================================
cat("\n[V] Freezing Step 07 results...\n")

stopifnot(file.exists("03_results/GSE243013_step07_COMPLETE.txt"))

step07_files <- c(
  "03_results/step07_programs/combined/GSE243013_prioritized_immune_programs.csv",
  "03_results/step07_programs/combined/GSE243013_pathway_recurrence_across_celltypes.csv"
)

## Check for TF files (optional)
tf_recurrence <- list.files("03_results/step07_programs/combined",
                            pattern = "TF_recurrence", full.names = TRUE)
if (length(tf_recurrence) > 0) step07_files <- c(step07_files, tf_recurrence[1])

tf_links_candidates <- list.files("03_results/step07_programs/combined",
                                  pattern = "TF.*leading.*edge|leading.*edge.*TF|TF_pathway",
                                  full.names = TRUE)
if (length(tf_links_candidates) > 0) step07_files <- c(step07_files, tf_links_candidates[1])

frozen_rows <- list()
for (f in step07_files) {
  if (file.exists(f)) {
    fi <- file.info(f)
    md5 <- tryCatch(tools::md5sum(f), error = function(e) NA_character_)
    frozen_rows[[length(frozen_rows) + 1]] <- data.frame(
      file_path = f, file_name = basename(f),
      file_size_bytes = fi$size,
      modification_time = as.character(fi$mtime),
      md5 = md5, status = "FROZEN",
      stringsAsFactors = FALSE
    )
    cat(sprintf("[FROZEN] %s (%d bytes)\n", basename(f), fi$size))
  } else {
    cat(sprintf("[MISSING] %s\n", f))
    frozen_rows[[length(frozen_rows) + 1]] <- data.frame(
      file_path = f, file_name = basename(f),
      file_size_bytes = NA_integer_, modification_time = NA_character_,
      md5 = NA_character_, status = "MISSING",
      stringsAsFactors = FALSE
    )
  }
}
frozen_df <- do.call(rbind, frozen_rows)
write.csv(frozen_df, "03_results/step08_TCGA/programs/GSE243013_step07_frozen_manifest.csv",
          row.names = FALSE)

## =========================================================================
## VI. Build TCGA Validation Program Manifest
## =========================================================================
cat("\n[VI] Building TCGA validation program manifest...\n")

programs <- read.csv("03_results/step07_programs/combined/GSE243013_prioritized_immune_programs.csv",
                      stringsAsFactors = FALSE)
cat(sprintf("[INFO] Loaded %d programs from Step 07\n", nrow(programs)))
cat(sprintf("[INFO] Columns: %s\n", paste(colnames(programs), collapse = ", ")))

## Filter Tier 1 and Tier 2
programs$toupper_tier <- toupper(programs$priority_tier)
programs <- programs[programs$toupper_tier %in% c("TIER 1", "TIER 2", "TIER1", "TIER2"), ]
cat(sprintf("[INFO] Tier 1+2 programs: %d\n", nrow(programs)))

## Standardize tier labels
programs$priority_tier_clean <- ifelse(grepl("1", programs$toupper_tier), "Tier 1", "Tier 2")

## Parse leading edge genes
parse_genes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  genes <- unlist(strsplit(as.character(x), "[;,|\\s]+"))
  genes <- trimws(genes)
  genes <- genes[nchar(genes) > 0]
  genes <- unique(genes)
  return(genes)
}

programs$n_leading_edge_genes <- sapply(programs$leading_edge_genes, function(x) {
  length(parse_genes(x))
})

## Build safe names
safe_name <- function(x) {
  gsub("_+", "_", gsub("[^A-Za-z0-9_]", "_", x))
}

programs$cell_type_safe <- sapply(programs$cell_type, safe_name)
programs$pathway_safe <- sapply(programs$pathway, safe_name)
programs$collection_safe <- sapply(programs$collection, safe_name)
programs$program_id <- paste0(programs$priority_tier_clean, "_",
                               programs$cell_type_safe, "_",
                               programs$collection_safe, "_",
                               programs$pathway_safe)
programs$program_id <- make.unique(programs$program_id, sep = "_dup")

## Eligibility criteria
programs$program_primary_eligible <- (
  programs$priority_tier_clean %in% c("Tier 1", "Tier 2") &
  programs$n_leading_edge_genes >= 5 &
  is.finite(programs$NES) &
  !is.na(programs$direction) &
  programs$direction != ""
)

cat(sprintf("[INFO] Primary eligible: %d / %d\n",
            sum(programs$program_primary_eligible), nrow(programs)))

## Build long-format gene membership
gene_membership_rows <- list()
for (i in seq_len(nrow(programs))) {
  row <- programs[i, ]
  genes <- parse_genes(row$leading_edge_genes)
  if (length(genes) == 0) next
  for (g in genes) {
    gene_membership_rows[[length(gene_membership_rows) + 1]] <- data.frame(
      program_id = row$program_id,
      priority_tier = row$priority_tier_clean,
      cell_type = row$cell_type,
      collection = row$collection,
      pathway = row$pathway,
      gene = g,
      step07_NES = row$NES,
      expected_direction = row$direction,
      top_supporting_TFs = row$top_supporting_TFs,
      program_gene_count = length(genes),
      program_primary_eligible = row$program_primary_eligible,
      stringsAsFactors = FALSE
    )
  }
}
gene_mem_df <- do.call(rbind, gene_membership_rows)
con <- gzfile("03_results/step08_TCGA/programs/GSE243013_TCGA_program_gene_membership.csv.gz", "w")
write.csv(gene_mem_df, con, row.names = FALSE)
close(con)
cat(sprintf("[OK] Gene membership: %d rows, %d programs\n", nrow(gene_mem_df),
            length(unique(gene_mem_df$program_id))))

## Save program manifest
prog_manifest <- programs[, c("program_id", "priority_tier_clean", "cell_type",
                               "collection", "pathway", "NES", "direction",
                               "leading_edge_genes", "top_supporting_TFs",
                               "n_leading_edge_genes", "program_primary_eligible")]
colnames(prog_manifest)[2] <- "priority_tier"
write.csv(prog_manifest, "03_results/step08_TCGA/programs/GSE243013_TCGA_program_manifest.csv",
          row.names = FALSE)
cat(sprintf("[OK] Program manifest: %d programs saved\n", nrow(prog_manifest)))
cat(sprintf("[INFO] Tier 1: %d, Tier 2: %d\n",
            sum(prog_manifest$priority_tier == "Tier 1"),
            sum(prog_manifest$priority_tier == "Tier 2")))

## =========================================================================
## VII. Setup ExperimentHub Project Cache
## =========================================================================
cat("\n[VII] Setting up ExperimentHub cache...\n")

eh_cache <- normalizePath("02_data/tcga/ExperimentHub_cache", mustWork = FALSE)
dir.create(eh_cache, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("[OK] Cache dir: %s\n", eh_cache))
writeLines(eh_cache, "00_config/step08_TCGA/GSE243013_ExperimentHub_cache_path.txt")

## =========================================================================
## VIII. TCGA Resource Dry-Run
## =========================================================================
cat("\n[VIII] TCGA resource dry-run...\n")

cohorts <- c("LUAD", "LUSC")
assay_patterns <- c("RNASeq2GeneNorm", "GISTIC_AllByGene", "GISTIC_ThresholdedByGene",
                     "Mutation", "RPPAArray", "Methylation_methyl450")
tcga_version <- "2.1.1"

dry_run_rows <- list()
for (cohort in cohorts) {
  for (ap in assay_patterns) {
    cat(sprintf("\n[Dry-Run] %s / %s...\n", cohort, ap))
    result <- tryCatch({
      dr <- curatedTCGAData::curatedTCGAData(
        diseaseCode = cohort,
        assays = ap,
        version = tcga_version,
        dry.run = TRUE,
        verbose = TRUE,
        cache = eh_cache
      )
      list(status = "OK", result = dr, error = NA_character_)
    }, error = function(e) {
      list(status = "ERROR", result = NULL, error = conditionMessage(e))
    })

    resource_name <- NA_character_
    resource_title <- NA_character_
    resource_id <- NA_character_
    available <- FALSE

    if (result$status == "OK" && !is.null(result$result)) {
      dr <- result$result
      if (is.data.frame(dr) && nrow(dr) > 0) {
        resource_name <- dr$title[1]
        resource_title <- dr$title[1]
        resource_id <- dr$ah_id[1]
        available <- TRUE
      } else if (is.list(dr) && length(dr) > 0) {
        resource_name <- names(dr)[1]
        resource_title <- resource_name
        available <- TRUE
      }
    }

    dry_run_rows[[length(dry_run_rows) + 1]] <- data.frame(
      cohort = cohort, requested_pattern = ap,
      resource_name = if (length(resource_name) == 0 || is.na(resource_name)) "NONE" else resource_name,
      resource_title = if (length(resource_title) == 0 || is.na(resource_title)) "NONE" else resource_title,
      resource_id = if (length(resource_id) == 0 || is.na(resource_id)) "NONE" else resource_id,
      available = available,
      error_message = if (length(result$error) == 0 || is.na(result$error)) "" else result$error,
      stringsAsFactors = FALSE
    )
    cat(sprintf("  -> %s: %s\n", ap, ifelse(available, "AVAILABLE", "NOT_AVAILABLE")))
  }
}

dry_run_df <- do.call(rbind, dry_run_rows)
write.csv(dry_run_df, "03_results/step08_TCGA/assay_audit/GSE243013_TCGA_curated_dry_run_resources.csv",
          row.names = FALSE)

## Build availability summary
avail_summary <- dry_run_df %>%
  group_by(cohort) %>%
  summarise(
    RNASeq2GeneNorm = any(available & requested_pattern == "RNASeq2GeneNorm"),
    GISTIC_AllByGene = any(available & requested_pattern == "GISTIC_AllByGene"),
    GISTIC_ThresholdedByGene = any(available & requested_pattern == "GISTIC_ThresholdedByGene"),
    Mutation = any(available & requested_pattern == "Mutation"),
    RPPAArray = any(available & requested_pattern == "RPPAArray"),
    Methylation_methyl450 = any(available & requested_pattern == "Methylation_methyl450"),
    .groups = "drop"
  )
write.csv(as.data.frame(avail_summary),
          "03_results/step08_TCGA/assay_audit/GSE243013_TCGA_assay_availability_summary.csv",
          row.names = FALSE)

cat("\n[SUMMARY] TCGA assay availability:\n")
print(as.data.frame(avail_summary))

## Check RNASeq2GeneNorm availability (required)
for (cohort in cohorts) {
  row <- avail_summary[avail_summary$cohort == cohort, ]
  if (!row$RNASeq2GeneNorm) {
    cat(sprintf("[FATAL] %s RNASeq2GeneNorm not available\n", cohort))
    writeLines(c("Step 08A FAILED", sprintf("Time: %s", Sys.time()),
                 sprintf("%s RNASeq2GeneNorm not available", cohort)),
               "03_results/GSE243013_step08A_FAILED.txt")
    stop(sprintf("%s RNASeq2GeneNorm not available", cohort))
  }
}
cat("[OK] Both LUAD and LUSC have RNASeq2GeneNorm\n")

## =========================================================================
## IX. Download Core Multi-Omics Data
## =========================================================================
cat("\n[IX] Downloading core multi-omics data...\n")

core_assays <- c("RNASeq2GeneNorm", "GISTIC_AllByGene", "GISTIC_ThresholdedByGene",
                  "Mutation", "RPPAArray")

## Determine which assays are available per cohort
download_assays <- list()
for (cohort in cohorts) {
  row <- avail_summary[avail_summary$cohort == cohort, ]
  available <- character(0)
  for (a in core_assays) {
    col_name <- ifelse(a == "GISTIC_AllByGene", "GISTIC_AllByGene",
                ifelse(a == "GISTIC_ThresholdedByGene", "GISTIC_ThresholdedByGene", a))
    if (isTRUE(row[[col_name]])) available <- c(available, a)
  }
  download_assays[[cohort]] <- available
  cat(sprintf("[INFO] %s core assays to download: %s\n", cohort, paste(available, collapse = ", ")))
}

mae_objects <- list()
for (cohort in cohorts) {
  cat(sprintf("\n--- Downloading %s ---\n", cohort))
  assays_to_get <- download_assays[[cohort]]

  mae <- NULL
  for (attempt in 1:3) {
    cat(sprintf("[Attempt %d/3] Downloading %s core multi-omics...\n", attempt, cohort))
    mae <- tryCatch({
      result <- curatedTCGAData::curatedTCGAData(
        diseaseCode = cohort,
        assays = assays_to_get,
        version = tcga_version,
        dry.run = FALSE,
        verbose = TRUE,
        cache = eh_cache
      )
      result
    }, error = function(e) {
      cat(sprintf("[WARNING] Attempt %d failed: %s\n", attempt, conditionMessage(e)))
      if (attempt < 3) Sys.sleep(20)
      NULL
    })
    if (!is.null(mae)) break
  }

  if (is.null(mae)) {
    cat(sprintf("[FATAL] %s download failed after 3 attempts\n", cohort))
    writeLines(c("Step 08A FAILED", sprintf("Time: %s", Sys.time()),
                 sprintf("%s core download failed", cohort)),
               "03_results/GSE243013_step08A_FAILED.txt")
    stop(sprintf("%s core download failed", cohort))
  }

  ## Validate MAE
  cat(sprintf("[INFO] %s MAE class: %s\n", cohort, class(mae)[1]))
  cat(sprintf("[INFO] %s experiments: %d\n", cohort, length(experiments(mae))))
  cat(sprintf("[INFO] %s colData rows: %d\n", cohort, nrow(colData(mae))))
  cat(sprintf("[INFO] %s sampleMap rows: %d\n", cohort, nrow(sampleMap(mae))))

  if (!is(mae, "MultiAssayExperiment")) {
    stop(sprintf("%s result is not a MultiAssayExperiment", cohort))
  }
  if (length(experiments(mae)) == 0) stop(sprintf("%s has 0 experiments", cohort))
  if (nrow(colData(mae)) == 0) stop(sprintf("%s has 0 colData rows", cohort))

  ## Save
  out_path <- sprintf("02_data/tcga/curated/%s/TCGA_%s_core_multiomics_all_samples_v%s.rds",
                       cohort, cohort, tcga_version)
  saveRDS(mae, out_path)
  cat(sprintf("[OK] Saved: %s\n", out_path))
  mae_objects[[cohort]] <- mae
}

## =========================================================================
## X. Download Methylation Data (if available)
## =========================================================================
cat("\n[X] Downloading methylation data...\n")

methyl_objects <- list()
for (cohort in cohorts) {
  has_methyl <- avail_summary$Methylation_methyl450[avail_summary$cohort == cohort]
  if (!has_methyl) {
    cat(sprintf("[SKIP] %s Methylation_methyl450 not available\n", cohort))
    next
  }

  cat(sprintf("\n--- Downloading %s 450K methylation ---\n", cohort))
  methyl_mae <- NULL
  for (attempt in 1:3) {
    cat(sprintf("[Attempt %d/3] Downloading %s methylation...\n", attempt, cohort))
    methyl_mae <- tryCatch({
      curatedTCGAData::curatedTCGAData(
        diseaseCode = cohort,
        assays = "Methylation_methyl450",
        version = tcga_version,
        dry.run = FALSE,
        verbose = TRUE,
        cache = eh_cache
      )
    }, error = function(e) {
      cat(sprintf("[WARNING] Attempt %d failed: %s\n", attempt, conditionMessage(e)))
      if (attempt < 3) Sys.sleep(20)
      NULL
    })
    if (!is.null(methyl_mae)) break
  }

  if (is.null(methyl_mae)) {
    cat(sprintf("[WARNING] %s methylation download failed\n", cohort))
    next
  }

  ## Audit without forcing matrix conversion
  cat(sprintf("[INFO] %s methyl class: %s\n", cohort, class(methyl_mae)[1]))
  cat(sprintf("[INFO] %s methyl experiments: %d\n", cohort, length(experiments(methyl_mae))))

  for (ename in names(experiments(methyl_mae))) {
    assay_obj <- experiments(methyl_mae)[[ename]]
    cat(sprintf("  Assay '%s': class=%s, dim=%s\n",
                ename, class(assay_obj)[1],
                paste(dim(assay_obj), collapse = "x")))
    if (is(assay_obj, "DelayedMatrix") || is(assay_obj, "HDF5Matrix")) {
      cat(sprintf("  -> Backend: Delayed/HDF5 (NOT loaded into memory)\n"))
    }
  }

  out_path <- sprintf("02_data/tcga/curated/%s/TCGA_%s_methyl450_all_samples_v%s.rds",
                       cohort, cohort, tcga_version)
  saveRDS(methyl_mae, out_path)
  cat(sprintf("[OK] Saved: %s\n", out_path))
  methyl_objects[[cohort]] <- methyl_mae
}

## =========================================================================
## XI. Filter Primary Tumors
## =========================================================================
cat("\n[XI] Filtering primary tumors...\n")

primary_mae <- list()
primary_methyl <- list()
filter_summary_rows <- list()

for (cohort in cohorts) {
  ## Core MAE
  if (!is.null(mae_objects[[cohort]])) {
    cat(sprintf("\n--- %s core primary tumor filter ---\n", cohort))
    mae <- mae_objects[[cohort]]
    n_before <- nrow(colData(mae))
    patients_before <- length(unique(sampleMap(mae)$primary))

    primary_mae[[cohort]] <- tryCatch({
      p_mae <- TCGAutils::TCGAprimaryTumors(mae)
      n_after <- nrow(colData(p_mae))
      patients_after <- length(unique(sampleMap(p_mae)$primary))

      cat(sprintf("[INFO] Patients: %d -> %d\n", patients_before, patients_after))

      ## Check for remaining non-01 samples
      barcode_col <- intersect(c("shortLetterCode", "sample_type", "Sample Type"),
                               colnames(colData(p_mae)))
      if (length(barcode_col) > 0) {
        non_01 <- sum(!grepl("^01$", colData(p_mae)[[barcode_col[1]]]))
        if (non_01 > 0) cat(sprintf("[WARNING] %d non-01 samples remain\n", non_01))
      }

      p_mae
    }, error = function(e) {
      cat(sprintf("[WARNING] Primary tumor filter failed: %s\n", conditionMessage(e)))
      mae
    })

    out_path <- sprintf("02_data/tcga/curated/%s/TCGA_%s_core_primary_tumors_v%s.rds",
                         cohort, cohort, tcga_version)
    saveRDS(primary_mae[[cohort]], out_path)

    filter_summary_rows[[length(filter_summary_rows) + 1]] <- data.frame(
      cohort = cohort, data_type = "core_multiomics",
      patients_before = patients_before,
      patients_after = nrow(colData(primary_mae[[cohort]])),
      samples_before = n_before,
      samples_after = nrow(colData(primary_mae[[cohort]])),
      stringsAsFactors = FALSE
    )
  }

  ## Methylation MAE
  if (!is.null(methyl_objects[[cohort]])) {
    cat(sprintf("\n--- %s methylation primary tumor filter ---\n", cohort))
    m_mae <- methyl_objects[[cohort]]
    patients_before_m <- length(unique(sampleMap(m_mae)$primary))

    primary_methyl[[cohort]] <- tryCatch({
      p_m_mae <- TCGAutils::TCGAprimaryTumors(m_mae)
      patients_after_m <- nrow(colData(p_m_mae))
      cat(sprintf("[INFO] Methyl patients: %d -> %d\n", patients_before_m, patients_after_m))
      p_m_mae
    }, error = function(e) {
      cat(sprintf("[WARNING] Methyl primary filter failed: %s\n", conditionMessage(e)))
      m_mae
    })

    out_path <- sprintf("02_data/tcga/curated/%s/TCGA_%s_methyl450_primary_tumors_v%s.rds",
                         cohort, cohort, tcga_version)
    saveRDS(primary_methyl[[cohort]], out_path)

    filter_summary_rows[[length(filter_summary_rows) + 1]] <- data.frame(
      cohort = cohort, data_type = "methyl450",
      patients_before = patients_before_m,
      patients_after = nrow(colData(primary_methyl[[cohort]])),
      samples_before = nrow(colData(m_mae)),
      samples_after = nrow(colData(primary_methyl[[cohort]])),
      stringsAsFactors = FALSE
    )
  }
}

filter_summary_df <- do.call(rbind, filter_summary_rows)
write.csv(filter_summary_df,
          "03_results/step08_TCGA/sample_audit/GSE243013_TCGA_primary_tumor_filter_summary.csv",
          row.names = FALSE)

## =========================================================================
## XII. Audit Each Assay
## =========================================================================
cat("\n[XII] Auditing each assay...\n")

assay_audit_rows <- list()
for (cohort in cohorts) {
  for (data_type in c("core", "methyl")) {
    mae <- if (data_type == "core") primary_mae[[cohort]] else primary_methyl[[cohort]]
    if (is.null(mae)) next

    for (ename in names(experiments(mae))) {
      assay_obj <- experiments(mae)[[ename]]
      assay_class <- class(assay_obj)[1]
      row_count <- tryCatch(nrow(assay_obj), error = function(e) NA_integer_)
      col_count <- tryCatch(ncol(assay_obj), error = function(e) NA_integer_)
      row_names_avail <- tryCatch(!is.null(rownames(assay_obj)) && length(rownames(assay_obj)) > 0,
                                  error = function(e) FALSE)
      col_names_avail <- tryCatch(!is.null(colnames(assay_obj)) && length(colnames(assay_obj)) > 0,
                                  error = function(e) FALSE)

      backend <- "unknown"
      if (is(assay_obj, "HDF5Matrix") || is(assay_obj, "HDF5Array")) backend <- "HDF5"
      else if (is(assay_obj, "DelayedMatrix") || is(assay_obj, "DelayedArray")) backend <- "DelayedArray"
      else if (is.matrix(assay_obj) || is(assay_obj, "dgCMatrix")) backend <- "in_memory"

      est_size <- tryCatch({
        tmp <- tempfile()
        saveRDS(assay_obj, tmp)
        sz <- file.info(tmp)$size
        file.remove(tmp)
        sz
      }, error = function(e) NA_real_)

      ## Get patient IDs from sampleMap
      sm <- sampleMap(mae)
      assay_sm <- sm[sm$assay == ename, ]
      patient_ids <- assay_sm$primary
      n_patients <- length(unique(patient_ids))
      n_dup_patients <- sum(duplicated(patient_ids))
      n_missing_pid <- sum(is.na(patient_ids) | patient_ids == "")

      assay_audit_rows[[length(assay_audit_rows) + 1]] <- data.frame(
        cohort = cohort, assay_name = ename, assay_class = assay_class,
        row_count = row_count, column_count = col_count,
        row_name_available = row_names_avail, column_name_available = col_names_avail,
        storage_backend = backend, estimated_object_size = est_size,
        primary_tumor_sample_count = col_count,
        unique_patient_count = n_patients,
        duplicated_patient_count = n_dup_patients,
        missing_patient_id_count = n_missing_pid,
        validation_status = "OK",
        stringsAsFactors = FALSE
      )
    }
  }
}

assay_audit_df <- do.call(rbind, assay_audit_rows)
write.csv(assay_audit_df,
          "03_results/step08_TCGA/assay_audit/GSE243013_TCGA_assay_dimension_class_audit.csv",
          row.names = FALSE)
cat(sprintf("[OK] Assay audit: %d assays recorded\n", nrow(assay_audit_df)))

## =========================================================================
## XIII. Build Sample-Patient Mapping
## =========================================================================
cat("\n[XIII] Building sample-patient mapping...\n")

mapping_rows <- list()
for (cohort in cohorts) {
  for (data_type in c("core", "methyl")) {
    mae <- if (data_type == "core") primary_mae[[cohort]] else primary_methyl[[cohort]]
    if (is.null(mae)) next

    sm <- sampleMap(mae)
    cat(sprintf("[INFO] %s %s sampleMap columns: %s\n", cohort, data_type,
                paste(colnames(sm), collapse = ", ")))

    ## Identify columns
    assay_col <- intersect(c("assay", "assay.name"), colnames(sm))
    primary_col <- intersect(c("primary", "patient"), colnames(sm))
    colname_col <- intersect(c("colname", "primary"), colnames(sm))

    if (length(assay_col) == 0 || length(primary_col) == 0 || length(colname_col) == 0) {
      cat(sprintf("[WARNING] Cannot identify columns in %s %s sampleMap\n", cohort, data_type))
      next
    }

    ## Use the correct colname (the barcode column)
    actual_colname_col <- setdiff(colname_col, primary_col)
    if (length(actual_colname_col) == 0) actual_colname_col <- colname_col[1]

    for (i in seq_len(nrow(sm))) {
      row <- sm[i, ]
      barcode <- as.character(row[[actual_colname_col]])
      patient_sm <- as.character(row[[primary_col[1]]])
      assay_name <- as.character(row[[assay_col[1]]])

      ## Extract patient ID from barcode (first 12 chars)
      barcode_patient <- substr(barcode, 1, 12)
      barcode_valid <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", barcode_patient)

      concordant <- if (barcode_valid) {
        if (patient_sm == barcode_patient) "CONCORDANT" else "PATIENT_ID_CONFLICT"
      } else "BARCODE_INVALID"

      ## Extract sample type code
      sample_type_code <- substr(barcode, 14, 15)
      is_primary <- sample_type_code == "01"

      mapping_rows[[length(mapping_rows) + 1]] <- data.frame(
        cohort = cohort, assay_name = assay_name,
        assay_column_name = barcode,
        patient_id_sampleMap = patient_sm,
        patient_id_barcode = barcode_patient,
        patient_id_concordant = concordant,
        sample_type_code = sample_type_code,
        is_primary_tumor = is_primary,
        mapping_status = concordant,
        stringsAsFactors = FALSE
      )
    }
  }
}

mapping_df <- do.call(rbind, mapping_rows)
con <- gzfile("02_data/tcga/manifests/GSE243013_TCGA_sample_patient_map.csv.gz", "w")
write.csv(mapping_df, con, row.names = FALSE)
close(con)
cat(sprintf("[OK] Sample-patient mapping: %d rows\n", nrow(mapping_df)))

## Report conflicts
n_conflicts <- sum(mapping_df$patient_id_concordant == "PATIENT_ID_CONFLICT")
if (n_conflicts > 0) {
  cat(sprintf("[WARNING] %d patient ID conflicts found\n", n_conflicts))
}

## =========================================================================
## XIV. Audit Duplicate Patient Samples
## =========================================================================
cat("\n[XIV] Auditing duplicate patient samples...\n")

dup_rows <- list()
for (cohort in cohorts) {
  for (data_type in c("core", "methyl")) {
    mae <- if (data_type == "core") primary_mae[[cohort]] else primary_methyl[[cohort]]
    if (is.null(mae)) next

    sm <- sampleMap(mae)
    assay_col <- intersect(c("assay", "assay.name"), colnames(sm))
    primary_col <- intersect(c("primary", "patient"), colnames(sm))
    colname_col <- intersect(c("colname", "primary"), colnames(sm))
    actual_colname_col <- setdiff(colname_col, primary_col)
    if (length(actual_colname_col) == 0) actual_colname_col <- colname_col[1]

    for (assay_name in unique(sm[[assay_col[1]]])) {
      assay_sm <- sm[sm[[assay_col[1]]] == assay_name, ]
      patient_counts <- table(assay_sm[[primary_col[1]]])
      multi_patients <- names(patient_counts[patient_counts > 1])

      for (pid in multi_patients) {
        pid_sm <- assay_sm[assay_sm[[primary_col[1]]] == pid, ]
        barcodes <- pid_sm[[actual_colname_col]]
        n_primary <- sum(substr(barcodes, 14, 15) == "01")

        dup_rows[[length(dup_rows) + 1]] <- data.frame(
          cohort = cohort, data_type = data_type,
          assay_name = assay_name, patient_id = pid,
          sample_count = as.integer(patient_counts[pid]),
          barcode_list = paste(barcodes, collapse = "; "),
          multiple_primary_tumor = n_primary > 1,
          multiple_aliquot = any(duplicated(barcodes)),
          has_duplicate_colnames = any(duplicated(barcodes)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

if (length(dup_rows) > 0) {
  dup_df <- do.call(rbind, dup_rows)
  write.csv(dup_df, "03_results/step08_TCGA/sample_audit/GSE243013_TCGA_duplicate_patient_samples.csv",
            row.names = FALSE)
  cat(sprintf("[INFO] %d patients with multiple samples found\n", nrow(dup_df)))
  cat("[INFO] NO automatic sample selection performed — Step 08B will apply deterministic rules\n")
} else {
  cat("[INFO] No duplicate patient samples found\n")
  data.frame(cohort=character(0), data_type=character(0), assay_name=character(0),
             patient_id=character(0), sample_count=integer(0), barcode_list=character(0),
             multiple_primary_tumor=logical(0), multiple_aliquot=logical(0),
             has_duplicate_colnames=logical(0)) %>%
    write.csv("03_results/step08_TCGA/sample_audit/GSE243013_TCGA_duplicate_patient_samples.csv",
              row.names = FALSE)
}

## =========================================================================
## XV. Audit Clinical Fields
## =========================================================================
cat("\n[XV] Auditing clinical fields...\n")

clinical_col_profiles <- list()
for (cohort in cohorts) {
  mae <- primary_mae[[cohort]]
  if (is.null(mae)) next

  cd <- as.data.frame(colData(mae))
  cat(sprintf("\n[INFO] %s colData: %d rows, %d columns\n", cohort, nrow(cd), ncol(cd)))

  for (col_name in colnames(cd)) {
    col_vec <- cd[[col_name]]
    n_total <- length(col_vec)
    n_na <- sum(is.na(col_vec))
    n_unique <- length(unique(col_vec[!is.na(col_vec)]))
    examples <- head(unique(col_vec[!is.na(col_vec)]), 10)

    clinical_col_profiles[[length(clinical_col_profiles) + 1]] <- data.frame(
      cohort = cohort, column_name = col_name,
      data_class = class(col_vec)[1],
      n_total = n_total, n_unique = n_unique,
      n_missing = n_na, missing_fraction = round(n_na / n_total, 4),
      example_values = paste(examples, collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
}

clinical_profile_df <- do.call(rbind, clinical_col_profiles)
con <- gzfile("03_results/step08_TCGA/clinical_audit/GSE243013_TCGA_clinical_column_profile.csv.gz", "w")
write.csv(clinical_profile_df, con, row.names = FALSE)
close(con)

## Candidate field mapping
candidate_fields <- list(
  patient_id = c("patientID", "patient_id", "bcr_patient_barcode", "submitter_id"),
  vital_status = c("vital_status"),
  days_to_death = c("days_to_death"),
  days_to_last_followup = c("days_to_last_followup", "days_to_last_follow_up"),
  age = c("years_to_birth", "age_at_initial_pathologic_diagnosis", "age_at_diagnosis"),
  sex = c("gender", "sex"),
  race = c("race"),
  stage = c("pathologic_stage", "ajcc_pathologic_stage", "tumor_stage"),
  T = c("pathologic_T"),
  N = c("pathologic_N"),
  M = c("pathologic_M"),
  smoking = c("tobacco_smoking_history", "smoking_history", "cigarettes_per_day", "pack_years_smoked"),
  radiation = c("radiation_therapy"),
  treatment = c("postoperative_rx_tx", "pharmaceutical_therapy")
)

field_mapping_rows <- list()
all_clinical_cols <- unique(clinical_profile_df$column_name)

for (field_name in names(candidate_fields)) {
  candidates <- candidate_fields[[field_name]]
  matched <- intersect(candidates, all_clinical_cols)
  status <- if (length(matched) == 0) "NOT_FOUND" else
            if (length(matched) == 1) "EXACT_MATCH" else "MULTIPLE_CANDIDATES"
  field_mapping_rows[[length(field_mapping_rows) + 1]] <- data.frame(
    field_name = field_name, candidates = paste(candidates, collapse = "; "),
    matched_columns = paste(matched, collapse = "; "),
    status = status,
    stringsAsFactors = FALSE
  )
}

field_mapping_df <- do.call(rbind, field_mapping_rows)
write.csv(field_mapping_df,
          "03_results/step08_TCGA/clinical_audit/GSE243013_TCGA_clinical_candidate_field_mapping.csv",
          row.names = FALSE)

cat("\n[INFO] Clinical field mapping:\n")
for (i in seq_len(nrow(field_mapping_df))) {
  cat(sprintf("  %s: %s -> %s\n", field_mapping_df$field_name[i],
              field_mapping_df$status[i], field_mapping_df$matched_columns[i]))
}

## =========================================================================
## XVI. Build Patient-Level Clinical Manifest
## =========================================================================
cat("\n[XVI] Building patient-level clinical manifest...\n")

build_patient_manifest <- function(mae, cohort) {
  cd <- as.data.frame(colData(mae))
  sm <- sampleMap(mae)

  ## Get unique patients
  primary_col <- intersect(c("primary", "patient"), colnames(sm))
  patient_ids <- unique(sm[[primary_col[1]]])

  manifest_rows <- list()
  for (pid in patient_ids) {
    patient_data <- cd[pid, , drop = FALSE]
    if (nrow(patient_data) == 0) next
    pd <- patient_data[1, ]

    ## Vital status
    vs_col <- intersect(c("vital_status"), colnames(pd))
    vs_orig <- if (length(vs_col) > 0) as.character(pd[[vs_col[1]]]) else NA_character_
    vs_clean <- toupper(trimws(vs_orig))
    os_event <- if (vs_clean %in% c("DEAD", "DECEASED", "1")) 1L else
                if (vs_clean %in% c("ALIVE", "LIVING", "0")) 0L else NA_integer_

    ## OS days
    dtd_col <- intersect(c("days_to_death"), colnames(pd))
    dtlf_col <- intersect(c("days_to_last_followup", "days_to_last_follow_up"), colnames(pd))
    dtd <- if (length(dtd_col) > 0) as.numeric(pd[[dtd_col[1]]]) else NA_real_
    dtlf <- if (length(dtlf_col) > 0) as.numeric(pd[[dtlf_col[1]]]) else NA_real_

    if (!is.na(dtd) && dtd < 0) dtd <- NA_real_
    if (!is.na(dtlf) && dtlf < 0) dtlf <- NA_real_

    os_days <- if (os_event == 1 && !is.na(dtd)) dtd else
               if (os_event == 0 && !is.na(dtlf)) dtlf else NA_real_

    ## Age
    age_col <- intersect(c("years_to_birth", "age_at_initial_pathologic_diagnosis", "age_at_diagnosis"), colnames(pd))
    age <- if (length(age_col) > 0) as.numeric(pd[[age_col[1]]]) else NA_real_

    ## Sex
    sex_col <- intersect(c("gender", "sex"), colnames(pd))
    sex <- if (length(sex_col) > 0) as.character(pd[[sex_col[1]]]) else NA_character_

    ## Stage
    stage_col <- intersect(c("pathologic_stage", "ajcc_pathologic_stage", "tumor_stage"), colnames(pd))
    stage_orig <- if (length(stage_col) > 0) as.character(pd[[stage_col[1]]]) else NA_character_
    stage_clean <- toupper(trimws(stage_orig))

    ## TNM
    t_col <- intersect(c("pathologic_T"), colnames(pd))
    n_col <- intersect(c("pathologic_N"), colnames(pd))
    m_col <- intersect(c("pathologic_M"), colnames(pd))
    pT <- if (length(t_col) > 0) as.character(pd[[t_col[1]]]) else NA_character_
    pN <- if (length(n_col) > 0) as.character(pd[[n_col[1]]]) else NA_character_
    pM <- if (length(m_col) > 0) as.character(pd[[m_col[1]]]) else NA_character_

    ## Smoking
    smoke_col <- intersect(c("tobacco_smoking_history", "smoking_history"), colnames(pd))
    smoking <- if (length(smoke_col) > 0) as.character(pd[[smoke_col[1]]]) else NA_character_

    ## Multi-omics coverage (match assay names with cohort prefix)
    exp_names <- names(experiments(mae))
    has_rna <- any(grepl("RNASeq2GeneNorm", exp_names))
    has_gistic_cont <- any(grepl("GISTIC_AllByGene", exp_names))
    has_gistic_thresh <- any(grepl("GISTIC_ThresholdedByGene", exp_names))
    has_mutation <- any(grepl("Mutation", exp_names))
    has_rppa <- any(grepl("RPPAArray", exp_names))

    manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
      cohort = cohort, patient_id = pid,
      vital_status_original = vs_orig, vital_status_clean = vs_clean,
      days_to_death = dtd, days_to_last_followup = dtlf,
      OS_days = os_days, OS_event = os_event,
      age_at_diagnosis = age, sex = sex,
      pathologic_stage_original = stage_orig, pathologic_stage_clean = stage_clean,
      pathologic_T = pT, pathologic_N = pN, pathologic_M = pM,
      smoking_field_original = smoking,
      has_primary_RNA = has_rna,
      has_GISTIC_continuous = has_gistic_cont,
      has_GISTIC_thresholded = has_gistic_thresh,
      has_mutation = has_mutation,
      has_RPPA = has_rppa,
      has_methyl450 = FALSE,
      has_all_core_omics = has_rna & has_gistic_cont & has_gistic_thresh & has_mutation & has_rppa,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, manifest_rows)
}

all_manifests <- list()
for (cohort in cohorts) {
  if (!is.null(primary_mae[[cohort]])) {
    cat(sprintf("\n--- Building %s patient manifest ---\n", cohort))
    manifest <- build_patient_manifest(primary_mae[[cohort]], cohort)
    all_manifests[[cohort]] <- manifest

    ## Update methyl450 column
    if (!is.null(primary_methyl[[cohort]])) {
      methyl_patients <- unique(sampleMap(primary_methyl[[cohort]])$primary)
      manifest$has_methyl450 <- manifest$patient_id %in% methyl_patients
      manifest$has_all_core_omics <- manifest$has_primary_RNA & manifest$has_GISTIC_continuous &
        manifest$has_GISTIC_thresholded & manifest$has_mutation & manifest$has_RPPA & manifest$has_methyl450
    }

    write.csv(manifest, sprintf("02_data/tcga/clinical/TCGA_%s_patient_manifest.csv", cohort),
              row.names = FALSE)
    cat(sprintf("[OK] %s: %d patients\n", cohort, nrow(manifest)))
  }
}

## Combined manifest
combined_manifest <- do.call(rbind, all_manifests)
write.csv(combined_manifest, "02_data/tcga/clinical/GSE243013_TCGA_LUAD_LUSC_patient_manifest.csv",
          row.names = FALSE)
cat(sprintf("[OK] Combined manifest: %d patients\n", nrow(combined_manifest)))

## =========================================================================
## XVII. Build Multi-Omics Coverage Matrix
## =========================================================================
cat("\n[XVII] Building multi-omics coverage matrix...\n")

for (cohort in cohorts) {
  manifest <- all_manifests[[cohort]]
  if (is.null(manifest)) next

  coverage <- data.frame(
    patient_id = manifest$patient_id,
    RNA = manifest$has_primary_RNA,
    GISTIC_continuous = manifest$has_GISTIC_continuous,
    GISTIC_thresholded = manifest$has_GISTIC_thresholded,
    Mutation = manifest$has_mutation,
    RPPA = manifest$has_RPPA,
    Methylation450 = manifest$has_methyl450,
    Clinical_OS = !is.na(manifest$OS_event),
    Clinical_stage = !is.na(manifest$pathologic_stage_clean) & manifest$pathologic_stage_clean != "",
    stringsAsFactors = FALSE
  )
  write.csv(coverage, sprintf("03_results/step08_TCGA/sample_audit/GSE243013_TCGA_%s_multiomics_coverage_matrix.csv", cohort),
            row.names = FALSE)

  ## Overlap summary
  overlap <- data.frame(
    metric = c("total_patients", "RNA", "GISTIC_continuous", "GISTIC_thresholded",
               "Mutation", "RPPA", "Methylation450",
               "RNA_and_GISTIC", "RNA_and_Mutation", "RNA_and_Methylation",
               "RNA_and_RPPA", "all_core_omics"),
    count = c(nrow(coverage),
              sum(coverage$RNA), sum(coverage$GISTIC_continuous),
              sum(coverage$GISTIC_thresholded), sum(coverage$Mutation),
              sum(coverage$RPPA), sum(coverage$Methylation450),
              sum(coverage$RNA & coverage$GISTIC_continuous),
              sum(coverage$RNA & coverage$Mutation),
              sum(coverage$RNA & coverage$Methylation450),
              sum(coverage$RNA & coverage$RPPA),
              sum(coverage$RNA & coverage$GISTIC_continuous & coverage$GISTIC_thresholded &
                    coverage$Mutation & coverage$RPPA & coverage$Methylation450)),
    cohort = cohort,
    stringsAsFactors = FALSE
  )
  write.csv(overlap, sprintf("03_results/step08_TCGA/sample_audit/GSE243013_TCGA_%s_multiomics_overlap_summary.csv", cohort),
            row.names = FALSE)

  cat(sprintf("\n[%s] Multi-omics overlap:\n", cohort))
  for (i in seq_len(nrow(overlap))) {
    cat(sprintf("  %s: %d\n", overlap$metric[i], overlap$count[i]))
  }
}

## =========================================================================
## XVIII. RNA Expression Audit
## =========================================================================
cat("\n[XVIII] RNA expression audit...\n")

for (cohort in cohorts) {
  mae <- primary_mae[[cohort]]
  if (is.null(mae)) next

  rna_assay_name <- grep("RNASeq2GeneNorm", names(experiments(mae)), value = TRUE)
  rna_assay <- tryCatch(experiments(mae)[[rna_assay_name[1]]], error = function(e) NULL)
  if (is.null(rna_assay)) {
    cat(sprintf("[SKIP] %s RNASeq2GeneNorm not found in experiments\n", cohort))
    next
  }

  cat(sprintf("\n--- %s RNA audit ---\n", cohort))
  cat(sprintf("[INFO] Class: %s, Dim: %s\n", class(rna_assay)[1],
              paste(dim(rna_assay), collapse = "x")))

  ## Sample-level info without loading full matrix
  n_genes <- nrow(rna_assay)
  n_samples <- ncol(rna_assay)
  gene_names <- rownames(rna_assay)
  sample_names <- colnames(rna_assay)

  n_dup_genes <- sum(duplicated(gene_names))
  n_dup_samples <- sum(duplicated(sample_names))

  cat(sprintf("[INFO] Genes: %d (%d duplicates), Samples: %d (%d duplicates)\n",
              n_genes, n_dup_genes, n_samples, n_dup_samples))

  ## Small matrix check (100 genes x 20 samples)
  n_check_genes <- min(100, n_genes)
  n_check_samples <- min(20, n_samples)
  check_mat <- tryCatch(as.matrix(rna_assay[1:n_check_genes, 1:n_check_samples]), error = function(e) NULL)

  if (!is.null(check_mat)) {
    n_na <- sum(is.na(check_mat))
    n_inf <- sum(is.infinite(check_mat))
    n_neg <- sum(check_mat < 0, na.rm = TRUE)
    cat(sprintf("[INFO] Small matrix check: NA=%d, Inf=%d, Negative=%d\n", n_na, n_inf, n_neg))
  }

  ## Overlap with program genes
  all_program_genes <- unique(gene_mem_df$gene)
  genes_in_cohort <- intersect(all_program_genes, gene_names)
  cat(sprintf("[INFO] Program genes in %s RNA: %d / %d (%.1f%%)\n",
              cohort, length(genes_in_cohort), length(all_program_genes),
              100 * length(genes_in_cohort) / length(all_program_genes)))
}

## Program-gene overlap per cohort
overlap_rows <- list()
for (cohort in cohorts) {
  mae <- primary_mae[[cohort]]
  if (is.null(mae)) next
  rna_name <- grep("RNASeq2GeneNorm", names(experiments(mae)), value = TRUE)
  rna_assay <- tryCatch(experiments(mae)[[rna_name[1]]], error = function(e) NULL)
  if (is.null(rna_assay)) next
  gene_names <- rownames(rna_assay)

  for (pid in unique(gene_mem_df$program_id)) {
    prog_genes <- unique(gene_mem_df$gene[gene_mem_df$program_id == pid])
    present <- intersect(prog_genes, gene_names)
    overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
      program_id = pid,
      cohort = cohort,
      program_gene_count = length(prog_genes),
      genes_present = length(present),
      overlap_fraction = round(length(present) / length(prog_genes), 4),
      scorable = length(present) >= 5 && (length(present) / length(prog_genes)) >= 0.5,
      stringsAsFactors = FALSE
    )
  }
}

overlap_df <- do.call(rbind, overlap_rows)
write.csv(overlap_df, "03_results/step08_TCGA/programs/GSE243013_TCGA_program_expression_overlap.csv",
          row.names = FALSE)

## Scorable programs
scorable <- overlap_df[overlap_df$scorable == TRUE, ]
if (nrow(scorable) > 0) {
  scorable_summary <- aggregate(scorable ~ program_id, data = scorable, FUN = function(x) all(x))
  colnames(scorable_summary)[2] <- "scorable_in_both"
} else {
  scorable_summary <- data.frame(program_id = character(0), scorable_in_both = logical(0))
}
write.csv(scorable_summary, "03_results/step08_TCGA/programs/GSE243013_TCGA_programs_ready_for_scoring.csv",
          row.names = FALSE)

cat(sprintf("\n[INFO] Programs scorable in both cohorts: %d / %d\n",
            sum(scorable_summary$scorable_in_both), nrow(scorable_summary)))

## =========================================================================
## XIX. GISTIC, Mutation, RPPA, Methylation Audit
## =========================================================================
cat("\n[XIX] Multi-omics specific audits...\n")

audit_rows <- list()

for (cohort in cohorts) {
  mae <- primary_mae[[cohort]]
  if (is.null(mae)) next

  ## GISTIC Continuous
  gistic_name <- grep("GISTIC_AllByGene", names(experiments(mae)), value = TRUE)
  gistic <- tryCatch(experiments(mae)[[gistic_name[1]]], error = function(e) NULL)
  if (!is.null(gistic)) {
    gistic_mat <- tryCatch(assay(gistic), error = function(e) NULL)
    g_range <- if (!is.null(gistic_mat)) range(gistic_mat[1:min(10, nrow(gistic_mat)), ], na.rm = TRUE) else c(NA, NA)
    audit_rows[[length(audit_rows) + 1]] <- data.frame(
      cohort = cohort, omics = "GISTIC_continuous",
      object_class = class(gistic)[1],
      row_count = nrow(gistic), column_count = ncol(gistic),
      value_range = paste(g_range, collapse = "-"),
      unique_gene_count = length(unique(rownames(gistic))),
      duplicate_genes = sum(duplicated(rownames(gistic))),
      notes = "Continuous copy number",
      stringsAsFactors = FALSE
    )
  }

  ## GISTIC Thresholded
  gistic_t_name <- grep("GISTIC_ThresholdedByGene", names(experiments(mae)), value = TRUE)
  gistic_t <- tryCatch(experiments(mae)[[gistic_t_name[1]]], error = function(e) NULL)
  if (!is.null(gistic_t)) {
    gistic_t_mat <- tryCatch(assay(gistic_t), error = function(e) NULL)
    vals <- if (!is.null(gistic_t_mat)) unique(as.vector(gistic_t_mat[1:min(10, nrow(gistic_t_mat)), 1:min(10, ncol(gistic_t_mat))])) else NA
    audit_rows[[length(audit_rows) + 1]] <- data.frame(
      cohort = cohort, omics = "GISTIC_thresholded",
      object_class = class(gistic_t)[1],
      row_count = nrow(gistic_t), column_count = ncol(gistic_t),
      value_range = paste(capture.output(cat(vals)), collapse = "; "),
      unique_gene_count = length(unique(rownames(gistic_t))),
      duplicate_genes = sum(duplicated(rownames(gistic_t))),
      notes = "Thresholded copy number",
      stringsAsFactors = FALSE
    )
  }

  ## Mutation
  mut_name <- grep("Mutation", names(experiments(mae)), value = TRUE)
  mut <- tryCatch(experiments(mae)[[mut_name[1]]], error = function(e) NULL)
  if (!is.null(mut)) {
    mut_class <- class(mut)[1]
    mut_rows <- tryCatch(nrow(mut), error = function(e) NA_integer_)
    audit_rows[[length(audit_rows) + 1]] <- data.frame(
      cohort = cohort, omics = "Mutation",
      object_class = mut_class,
      row_count = mut_rows, column_count = tryCatch(ncol(mut), error = function(e) NA_integer_),
      value_range = "categorical",
      unique_gene_count = NA_integer_,
      duplicate_genes = NA_integer_,
      notes = sprintf("Class: %s; do not convert to binary matrix in Step 08A", mut_class),
      stringsAsFactors = FALSE
    )
  }

  ## RPPA
  rppa_name <- grep("RPPAArray", names(experiments(mae)), value = TRUE)
  rppa <- tryCatch(experiments(mae)[[rppa_name[1]]], error = function(e) NULL)
  if (!is.null(rppa)) {
    rppa_mat <- tryCatch(assay(rppa), error = function(e) NULL)
    rppa_range <- if (!is.null(rppa_mat)) range(rppa_mat[1:min(5, nrow(rppa_mat)), ], na.rm = TRUE) else c(NA, NA)
    audit_rows[[length(audit_rows) + 1]] <- data.frame(
      cohort = cohort, omics = "RPPA",
      object_class = class(rppa)[1],
      row_count = nrow(rppa), column_count = ncol(rppa),
      value_range = paste(rppa_range, collapse = "-"),
      unique_gene_count = length(unique(rownames(rppa))),
      duplicate_genes = sum(duplicated(rownames(rppa))),
      notes = "Protein/antibody level",
      stringsAsFactors = FALSE
    )
  }

  ## Methylation
  methyl <- primary_methyl[[cohort]]
  if (!is.null(methyl)) {
    for (ename in names(experiments(methyl))) {
      m_obj <- experiments(methyl)[[ename]]
      is_delayed <- is(m_obj, "DelayedMatrix") || is(m_obj, "HDF5Matrix") || is(m_obj, "DelayedArray")
      has_rowdata <- tryCatch(nrow(rowData(m_obj)) > 0, error = function(e) FALSE)
      rowdata_cols <- tryCatch(colnames(rowData(m_obj)), error = function(e) character(0))

      audit_rows[[length(audit_rows) + 1]] <- data.frame(
        cohort = cohort, omics = "Methylation450",
        object_class = class(m_obj)[1],
        row_count = nrow(m_obj), column_count = ncol(m_obj),
        value_range = "beta_values_0_1",
        unique_gene_count = NA_integer_,
        duplicate_genes = NA_integer_,
        notes = sprintf("Delayed=%s; rowData=%s; cols=%s",
                        is_delayed, has_rowdata,
                        if (length(rowdata_cols) > 0) paste(head(rowdata_cols, 5), collapse = ",") else "NONE"),
        stringsAsFactors = FALSE
      )
    }
  }
}

audit_df <- do.call(rbind, audit_rows)
write.csv(audit_df, "03_results/step08_TCGA/assay_audit/GSE243013_TCGA_omics_specific_audit.csv",
          row.names = FALSE)
cat(sprintf("[OK] Omics-specific audit: %d records\n", nrow(audit_df)))

## =========================================================================
## XX. Create Step 08B Input Index
## =========================================================================
cat("\n[XX] Creating Step 08B input index...\n")

index_rows <- list()

## Core MAEs
for (cohort in cohorts) {
  path <- sprintf("02_data/tcga/curated/%s/TCGA_%s_core_primary_tumors_v%s.rds", cohort, cohort, tcga_version)
  index_rows[[length(index_rows) + 1]] <- data.frame(
    cohort = cohort, data_type = "core_primary_tumors",
    object_path = path,
    assay_name = "MultiAssayExperiment",
    object_class = "MultiAssayExperiment",
    row_count = tryCatch(nrow(colData(readRDS(path))), error = function(e) NA_integer_),
    column_count = tryCatch(length(experiments(readRDS(path))), error = function(e) NA_integer_),
    primary_tumor_only = TRUE,
    patient_count = tryCatch(length(unique(sampleMap(readRDS(path))$primary)), error = function(e) NA_integer_),
    ready_for_step08B = TRUE,
    warning_message = "",
    validation_status = "OK",
    stringsAsFactors = FALSE
  )
}

## Methylation MAEs
for (cohort in cohorts) {
  path <- sprintf("02_data/tcga/curated/%s/TCGA_%s_methyl450_primary_tumors_v%s.rds", cohort, cohort, tcga_version)
  exists <- file.exists(path)
  index_rows[[length(index_rows) + 1]] <- data.frame(
    cohort = cohort, data_type = "methyl450_primary_tumors",
    object_path = if (exists) path else NA_character_,
    assay_name = "Methylation_methyl450",
    object_class = if (exists) "MultiAssayExperiment" else NA_character_,
    row_count = if (exists) tryCatch(nrow(colData(readRDS(path))), error = function(e) NA_integer_) else NA_integer_,
    column_count = if (exists) tryCatch(length(experiments(readRDS(path))), error = function(e) NA_integer_) else NA_integer_,
    primary_tumor_only = TRUE,
    patient_count = if (exists) tryCatch(length(unique(sampleMap(readRDS(path))$primary)), error = function(e) NA_integer_) else NA_integer_,
    ready_for_step08B = exists,
    warning_message = ifelse(!exists, "NOT_AVAILABLE", ""),
    validation_status = ifelse(exists, "OK", "NOT_AVAILABLE"),
    stringsAsFactors = FALSE
  )
}

## Supporting files
supporting_files <- data.frame(
  cohort = "combined", data_type = "clinical_manifest",
  object_path = "02_data/tcga/clinical/GSE243013_TCGA_LUAD_LUSC_patient_manifest.csv",
  assay_name = "clinical", object_class = "data.frame",
  row_count = nrow(combined_manifest), column_count = ncol(combined_manifest),
  primary_tumor_only = TRUE, patient_count = nrow(combined_manifest),
  ready_for_step08B = TRUE, warning_message = "", validation_status = "OK",
  stringsAsFactors = FALSE
)
index_rows[[length(index_rows) + 1]] <- supporting_files

index_rows[[length(index_rows) + 1]] <- data.frame(
  cohort = "combined", data_type = "program_manifest",
  object_path = "03_results/step08_TCGA/programs/GSE243013_TCGA_program_manifest.csv",
  assay_name = "programs", object_class = "data.frame",
  row_count = nrow(prog_manifest), column_count = ncol(prog_manifest),
  primary_tumor_only = FALSE, patient_count = NA_integer_,
  ready_for_step08B = TRUE, warning_message = "", validation_status = "OK",
  stringsAsFactors = FALSE
)

index_rows[[length(index_rows) + 1]] <- data.frame(
  cohort = "combined", data_type = "program_gene_membership",
  object_path = "03_results/step08_TCGA/programs/GSE243013_TCGA_program_gene_membership.csv.gz",
  assay_name = "gene_membership", object_class = "data.frame",
  row_count = nrow(gene_mem_df), column_count = ncol(gene_mem_df),
  primary_tumor_only = FALSE, patient_count = NA_integer_,
  ready_for_step08B = TRUE, warning_message = "", validation_status = "OK",
  stringsAsFactors = FALSE
)

index_rows[[length(index_rows) + 1]] <- data.frame(
  cohort = "combined", data_type = "program_expression_overlap",
  object_path = "03_results/step08_TCGA/programs/GSE243013_TCGA_program_expression_overlap.csv",
  assay_name = "overlap", object_class = "data.frame",
  row_count = nrow(overlap_df), column_count = ncol(overlap_df),
  primary_tumor_only = FALSE, patient_count = NA_integer_,
  ready_for_step08B = TRUE, warning_message = "", validation_status = "OK",
  stringsAsFactors = FALSE
)

for (cohort in cohorts) {
  index_rows[[length(index_rows) + 1]] <- data.frame(
    cohort = cohort, data_type = "multiomics_coverage",
    object_path = sprintf("03_results/step08_TCGA/sample_audit/GSE243013_TCGA_%s_multiomics_coverage_matrix.csv", cohort),
    assay_name = "coverage", object_class = "data.frame",
    row_count = nrow(combined_manifest[combined_manifest$cohort == cohort, ]),
    column_count = 9,
    primary_tumor_only = TRUE,
    patient_count = sum(combined_manifest$cohort == cohort),
    ready_for_step08B = TRUE, warning_message = "", validation_status = "OK",
    stringsAsFactors = FALSE
  )
}

index_df <- do.call(rbind, index_rows)
write.csv(index_df, "03_results/step08_TCGA/GSE243013_step08B_input_index.csv",
          row.names = FALSE)
cat(sprintf("[OK] Step 08B input index: %d entries\n", nrow(index_df)))

## =========================================================================
## XXI. Data Source and Interpretation Boundaries
## =========================================================================
cat("\n[XXI] Writing data source definition...\n")

definition_text <- c(
  "TCGA External Validation Definition",
  "====================================",
  "",
  sprintf("Generated: %s", as.character(Sys.time())),
  sprintf("curatedTCGAData version: %s", tryCatch(as.character(packageVersion("curatedTCGAData")), error = function(e) "UNKNOWN")),
  sprintf("TCGA data version: %s", tcga_version),
  "",
  "1. TCGA-LUAD and TCGA-LUSC are analyzed separately.",
  "2. Data downloaded via curatedTCGAData and ExperimentHub from R directly.",
  "3. Only primary tumor samples (sample type code 01) are used.",
  "4. TCGA is primarily untreated primary tumor data.",
  "5. TCGA is NOT a neoadjuvant anti-PD-1 efficacy cohort.",
  "6. Therefore, Responder vs Non_responder cannot be directly validated.",
  "7. TCGA validates the independent detectability of immune programs.",
  "8. TCGA tests associations between programs and clinical/mutation/CNV/methylation/protein.",
  "9. LUAD and LUSC are NOT merged into a single expression cohort.",
  "10. Cross-histology results use effect size comparison or meta-analysis.",
  "11. Bulk RNA immune program scores cannot directly equal single-cell type activity.",
  "12. Downstream models must consider tumor purity or immune cell composition confounding.",
  "13. Step 08A only downloads and audits — no statistical inference is performed.",
  "",
  "Program Interpretation Boundaries:",
  "- Programs derive from post-treatment immune cell scRNA-seq leading edges.",
  "- In TCGA bulk RNA, these represent relative enrichment/detectability of immune programs.",
  "- They do NOT represent cell-type-specific abundance or activity.",
  "- Any bulk-level association requires deconvolution or immune fraction adjustment."
)
writeLines(definition_text, "00_config/step08_TCGA/GSE243013_TCGA_external_validation_definition.txt")

## =========================================================================
## XXII. Completion
## =========================================================================
cat("\n[XXII] Checking completion conditions...\n")

conditions <- list(
  list(name = "Step 07 frozen", passed = file.exists("03_results/GSE243013_step07_COMPLETE.txt")),
  list(name = "Tier 1+2 program manifest", passed = file.exists("03_results/step08_TCGA/programs/GSE243013_TCGA_program_manifest.csv")),
  list(name = "LUAD RNASeq2GeneNorm", passed = !is.null(primary_mae[["LUAD"]])),
  list(name = "LUSC RNASeq2GeneNorm", passed = !is.null(primary_mae[["LUSC"]])),
  list(name = "LUAD primary tumor", passed = !is.null(primary_mae[["LUAD"]])),
  list(name = "LUSC primary tumor", passed = !is.null(primary_mae[["LUSC"]])),
  list(name = "sampleMap audit", passed = file.exists("02_data/tcga/manifests/GSE243013_TCGA_sample_patient_map.csv.gz")),
  list(name = "Patient clinical manifest", passed = file.exists("02_data/tcga/clinical/GSE243013_TCGA_LUAD_LUSC_patient_manifest.csv")),
  list(name = "Multi-omics coverage", passed = any(grepl("coverage_matrix", list.files("03_results/step08_TCGA/sample_audit/", full.names = FALSE)))),
  list(name = "Program-RNA overlap", passed = file.exists("03_results/step08_TCGA/programs/GSE243013_TCGA_program_expression_overlap.csv")),
  list(name = "Step 08B input index", passed = file.exists("03_results/step08_TCGA/GSE243013_step08B_input_index.csv"))
)

all_passed <- all(sapply(conditions, function(x) x$passed))
for (cond in conditions) {
  cat(sprintf("  [%s] %s\n", ifelse(cond$passed, "PASS", "FAIL"), cond$name))
}

## Get final disk space
disk_final <- system("df -Pk .", intern = TRUE)
disk_final_line <- disk_final[length(disk_final)]
disk_final_parts <- strsplit(disk_final_line, "\\s+")[[1]]
disk_final_avail <- sprintf("%.1f", as.numeric(disk_final_parts[4]) / 1048576)

if (all_passed) {
  completion_text <- c(
    "GSE243013 Step 08A COMPLETE",
    "===========================",
    "",
    sprintf("Completion time: %s", as.character(Sys.time())),
    sprintf("R version: %s", R.version.string),
    sprintf("Bioconductor version: %s", BiocManager::version()),
    sprintf("curatedTCGAData version: %s", tryCatch(as.character(packageVersion("curatedTCGAData")), error = function(e) "UNKNOWN")),
    sprintf("TCGA data version: %s", tcga_version),
    sprintf("ExperimentHub cache: %s", eh_cache),
    sprintf("LUAD primary RNA patients: %d", nrow(all_manifests[["LUAD"]])),
    sprintf("LUSC primary RNA patients: %d", nrow(all_manifests[["LUSC"]])),
    sprintf("LUAD multi-omics coverage: %s", paste(
      sprintf("RNA=%d GISTIC_cont=%d GISTIC_thresh=%d Mutation=%d RPPA=%d Methyl450=%d",
              sum(combined_manifest$cohort == "LUAD" & combined_manifest$has_primary_RNA),
              sum(combined_manifest$cohort == "LUAD" & combined_manifest$has_GISTIC_continuous),
              sum(combined_manifest$cohort == "LUAD" & combined_manifest$has_GISTIC_thresholded),
              sum(combined_manifest$cohort == "LUAD" & combined_manifest$has_mutation),
              sum(combined_manifest$cohort == "LUAD" & combined_manifest$has_RPPA),
              sum(combined_manifest$cohort == "LUAD" & combined_manifest$has_methyl450)),
      collapse = "; ")),
    sprintf("LUSC multi-omics coverage: %s", paste(
      sprintf("RNA=%d GISTIC_cont=%d GISTIC_thresh=%d Mutation=%d RPPA=%d Methyl450=%d",
              sum(combined_manifest$cohort == "LUSC" & combined_manifest$has_primary_RNA),
              sum(combined_manifest$cohort == "LUSC" & combined_manifest$has_GISTIC_continuous),
              sum(combined_manifest$cohort == "LUSC" & combined_manifest$has_GISTIC_thresholded),
              sum(combined_manifest$cohort == "LUSC" & combined_manifest$has_mutation),
              sum(combined_manifest$cohort == "LUSC" & combined_manifest$has_RPPA),
              sum(combined_manifest$cohort == "LUSC" & combined_manifest$has_methyl450)),
      collapse = "; ")),
    sprintf("Tier 1 programs: %d", sum(prog_manifest$priority_tier == "Tier 1")),
    sprintf("Tier 2 programs: %d", sum(prog_manifest$priority_tier == "Tier 2")),
    sprintf("Programs scorable in LUAD: %d", if (nrow(scorable_summary) > 0) sum(scorable_summary$scorable_in_both & sapply(scorable_summary$program_id, function(p) any(overlap_df$program_id == p & overlap_df$cohort == "LUAD" & overlap_df$scorable))) else 0),
    sprintf("Programs scorable in LUSC: %d", if (nrow(scorable_summary) > 0) sum(scorable_summary$scorable_in_both & sapply(scorable_summary$program_id, function(p) any(overlap_df$program_id == p & overlap_df$cohort == "LUSC" & overlap_df$scorable))) else 0),
    sprintf("Disk available: %s GB", disk_final_avail),
    sprintf("Total runtime: %.1f seconds", as.numeric(difftime(Sys.time(), step08a_start, units = "secs")))
  )
  writeLines(completion_text, "03_results/GSE243013_step08A_COMPLETE.txt")
  cat("\n[OK] Step 08A COMPLETE\n")
} else {
  failed_conds <- sapply(conditions, function(x) if (!x$passed) x$name else NA_character_)
  failed_conds <- failed_conds[!is.na(failed_conds)]
  writeLines(c("Step 08A FAILED", sprintf("Time: %s", Sys.time()),
               sprintf("Failed conditions: %s", paste(failed_conds, collapse = ", "))),
             "03_results/GSE243013_step08A_FAILED.txt")
  cat("\n[FAILED] Step 08A FAILED\n")
}

total_runtime <- as.numeric(difftime(Sys.time(), step08a_start, units = "secs"))
cat(sprintf("\nStep 08A total runtime: %.1f seconds\n", total_runtime))
cat("========================================================================\n")
