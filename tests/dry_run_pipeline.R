#!/usr/bin/env Rscript
# Test: Dry-run pipeline simulation
cat("=== Dry Run Pipeline Test ===\n")

repo <- Sys.getenv("GSE243013_PROJECT_ROOT", unset = normalizePath("..", mustWork = FALSE))

steps <- list(
  c("00", "00_environment_and_manifest.R"),
  c("01", "01_download_and_inspect_GSE243013_metadata.R"),
  c("02", "02_build_GSE243013_patient_manifest.R"),
  c("03", "03_repair_GSE243013_cohort_definition.R"),
  c("04", "04_download_and_import_GSE243013_counts.R"),
  c("04A", "04A_install_BPCells_direct_binary.R"),
  c("05", "05_build_GSE243013_patient_celltype_pseudobulk.R"),
  c("06", "06_edgeR_patient_level_differential_expression.R"),
  c("07", "07_pathway_TF_program_integration.R"),
  c("07A", "07A_install_extracted_msigdb_and_resume.R"),
  c("08A", "08A_download_audit_TCGA_multiomics.R"),
  c("08B1", "08B1_QC2_reconcile_and_rebuild_TCGA_clinical_models.R"),
  c("08B2", "08B2_TCGA_multiomics_integration.R"),
  c("09", "09_finalize_project_and_build_evidence_report.R"),
  c("09A", "09A_close_final_project_gaps.R")
)

cat("Step-by-step dry run:\n")
all_ok <- TRUE
for (s in steps) {
  sp <- file.path(repo, "scripts/analysis", s[2])
  exists <- file.exists(sp)
  status <- ifelse(exists, "OK", "MISSING")
  if (!exists) all_ok <- FALSE
  cat(sprintf("  Step %-5s %-55s [%s]\n", s[1], s[2], status))
}

if (all_ok) {
  cat("\nPASS: All", length(steps), "scripts found\n")
  quit(status=0)
} else {
  cat("\nFAIL: Some scripts missing\n")
  quit(status=1)
}
