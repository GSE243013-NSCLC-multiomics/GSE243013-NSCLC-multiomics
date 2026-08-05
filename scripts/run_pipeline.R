#!/usr/bin/env Rscript
# run_pipeline.R — GSE243013 NSCLC Multi-Omics Analysis Pipeline Runner
# Usage:
#   Rscript scripts/run_pipeline.R
#   Rscript scripts/run_pipeline.R --from-step 05 --to-step 07
#   Rscript scripts/run_pipeline.R --step 06
#   Rscript scripts/run_pipeline.R --dry-run
#   Rscript scripts/run_pipeline.R --config config/paths.yml

suppressPackageStartupMethods(library(optparse))

option_list <- list(
  make_option("--from-step", type="character", default="00",
              help="Start from this step [default: %default]"),
  make_option("--to-step", type="character", default="99",
              help="Stop after this step [default: %default]"),
  make_option("--step", type="character", default=NULL,
              help="Run only this single step"),
  make_option("--config", type="character", default="config/paths.yml",
              help="Path to config file [default: %default]"),
  make_option("--dry-run", action="store_true", default=FALSE,
              help="Show what would run without executing")
)

parser <- OptionParser(
  usage = "Rscript scripts/run_pipeline.R [options]",
  option_list = option_list
)
args <- parse_args(parser)

# Set project root
project_root <- Sys.getenv("GSE243013_PROJECT_ROOT",
                           unset = normalizePath(".", mustWork = FALSE))
cat("Project root:", project_root, "\n")

# Check config
config_path <- file.path(project_root, args$config)
if (!file.exists(config_path)) {
  cat("WARNING: Config file not found:", config_path, "\n")
  cat("Copying from paths.example.yml\n")
  example <- file.path(project_root, "config/paths.example.yml")
  if (file.exists(example)) file.copy(example, config_path)
}

# Script registry
scripts <- list(
  "00"  = "00_environment_and_manifest.R",
  "01"  = "01_download_and_inspect_GSE243013_metadata.R",
  "02"  = "02_build_GSE243013_patient_manifest.R",
  "03"  = "03_repair_GSE243013_cohort_definition.R",
  "04"  = "04_download_and_import_GSE243013_counts.R",
  "04A" = "04A_install_BPCells_direct_binary.R",
  "05"  = "05_build_GSE243013_patient_celltype_pseudobulk.R",
  "06"  = "06_edgeR_patient_level_differential_expression.R",
  "07"  = "07_pathway_TF_program_integration.R",
  "07A" = "07A_install_extracted_msigdb_and_resume.R",
  "08A" = "08A_download_audit_TCGA_multiomics.R",
  "08B1"= "08B1_QC2_reconcile_and_rebuild_TCGA_clinical_models.R",
  "08B2"= "08B2_TCGA_multiomics_integration.R",
  "09"  = "09_finalize_project_and_build_evidence_report.R",
  "09A" = "09A_close_final_project_gaps.R"
)

step_order <- c("00","01","02","03","04","04A","05","06","07","07A",
                "08A","08B1","08B2","09","09A")

# Determine steps to run
if (!is.null(args$step)) {
  steps_to_run <- args$step
} else {
  from_idx <- match(args$from.step, step_order)
  to_idx <- match(args$to.step, step_order)
  if (is.na(from_idx)) from_idx <- 1
  if (is.na(to_idx)) to_idx <- length(step_order)
  steps_to_run <- step_order[seq(from_idx, to_idx)]
}

cat("\n========================================\n")
cat("  GSE243013 NSCLC Multi-Omics Pipeline\n")
cat("========================================\n")
cat("  Steps:", paste(steps_to_run, collapse=", "), "\n")
cat("  Dry run:", args$dry.run, "\n")
cat("========================================\n\n")

run_count <- 0
fail_count <- 0

for (step in steps_to_run) {
  script_name <- scripts[[step]]
  if (is.null(script_name)) {
    cat("[SKIP] Step", step, ": Unknown step\n")
    next
  }

  script_path <- file.path(project_root, "scripts/analysis", script_name)
  if (!file.exists(script_path)) {
    cat("[SKIP] Step", step, ": Script not found\n")
    next
  }

  if (args$dry.run) {
    cat("[DRY-RUN] Step", step, ":", script_name, "\n")
    next
  }

  cat("[RUN] Step", step, ":", script_name, "\n")
  start_time <- Sys.time()

  result <- tryCatch({
    source(script_path, local = new.env(parent = globalenv()))
    TRUE
  }, error = function(e) {
    cat("[ERROR]", conditionMessage(e), "\n")
    FALSE
  })

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
  if (result) {
    cat("[DONE] Step", step, "completed in", round(elapsed), "s\n")
    run_count <- run_count + 1
  } else {
    cat("[FAIL] Step", step, "failed after", round(elapsed), "s\n")
    fail_count <- fail_count + 1
  }
  cat("\n")
}

cat("========================================\n")
cat("  Pipeline complete\n")
cat("  Steps run:", run_count, "\n")
cat("  Steps failed:", fail_count, "\n")
cat("========================================\n")
