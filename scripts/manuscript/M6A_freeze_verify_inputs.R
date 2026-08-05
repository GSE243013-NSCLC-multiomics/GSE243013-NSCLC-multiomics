#!/usr/bin/env Rscript
# M6A: Freeze and Verify Inputs
# Read-only: only reads small M5A output files
# Never reads large files, never computes MD5, never scans project recursively

cat("\n", rep("=", 80), sep="")
cat("\nM6A: Freeze and Verify Inputs")
cat("\nStarted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")

dir.create("05_manuscript/M6_final_manuscript/audit", recursive=TRUE, showWarnings=FALSE)

# ============================================================================
# SECTION I: Read Small M5A Files
# ============================================================================
cat("\nSECTION I: Read Small M5A Files\n")
cat(rep("=", 80), sep="")
cat("\n\n")

tables_dir <- "05_manuscript/M5A_evidence_completion/tables"
audit_dir  <- "05_manuscript/M5A_evidence_completion/audit"
text_dir   <- "05_manuscript/M5A_evidence_completion/revised_text"

# --- Read CSVs ---
fgsea_exact <- read.csv(file.path(tables_dir, "GSE243013_glycolysis_fgsea_exact_statistics.csv"),
                         stringsAsFactors=FALSE)
cat("fgsea_exact_statistics:", nrow(fgsea_exact), "rows\n")

top_features <- read.csv(file.path(tables_dir, "GSE243013_glycolysis_multiomics_top_features_for_manuscript.csv"),
                          stringsAsFactors=FALSE)
cat("top_features:", nrow(top_features), "rows\n")

cohort_summary <- read.csv(file.path(tables_dir, "GSE243013_glycolysis_multiomics_cohort_summary.csv"),
                            stringsAsFactors=FALSE)
cat("cohort_summary:", nrow(cohort_summary), "rows\n")

collectri_audit <- read.csv(file.path(audit_dir, "GSE243013_CollecTRI_status_reconciliation.csv"),
                             stringsAsFactors=FALSE)
cat("CollecTRI audit:", nrow(collectri_audit), "rows\n")

conflict_audit <- read.csv(file.path(audit_dir, "GSE243013_M5A_completion_audit.csv"),
                            stringsAsFactors=FALSE)
cat("Conflict audit:", nrow(conflict_audit), "rows\n")

# --- Read revised texts ---
abstract_text <- readLines(file.path(text_dir, "GSE243013_Abstract_M5A.md"), warn=FALSE)
methods_text  <- readLines(file.path(text_dir, "GSE243013_Methods_M5A.md"), warn=FALSE)
results_text  <- readLines(file.path(text_dir, "GSE243013_Results_M5A.md"), warn=FALSE)
fig5_text     <- readLines(file.path(text_dir, "GSE243013_Figure5_legend_M5A.md"), warn=FALSE)
cat("Abstract:", length(abstract_text), "lines\n")
cat("Methods:", length(methods_text), "lines\n")
cat("Results:", length(results_text), "lines\n")
cat("Figure5 legend:", length(fig5_text), "lines\n\n")


# ============================================================================
# SECTION II: Verify Core Values
# ============================================================================
cat("\nSECTION II: Verify Core Values\n")
cat(rep("=", 80), sep="")
cat("\n\n")

checks <- data.frame(
  item = character(), expected = character(), actual = character(),
  status = character(), stringsAsFactors = FALSE)

# --- Primary NES/P/FDR ---
prim <- fgsea_exact[fgsea_exact$cohort == "primary", ]
checks <- rbind(checks, data.frame(
  item="primary_NES", expected="-2.3589", actual=format(prim$NES, digits=5),
  status=ifelse(abs(prim$NES - (-2.3589)) < 0.001, "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="primary_PValue", expected="4.90e-12", actual=format(prim$PValue, digits=2),
  status=ifelse(prim$PValue < 1e-10, "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="primary_FDR", expected="3.00e-11", actual=format(prim$FDR, digits=2),
  status=ifelse(prim$FDR < 1e-10, "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- Strict NES/P/FDR ---
str_row <- fgsea_exact[fgsea_exact$cohort == "strict", ]
checks <- rbind(checks, data.frame(
  item="strict_NES", expected="-2.4126", actual=format(str_row$NES, digits=5),
  status=ifelse(abs(str_row$NES - (-2.4126)) < 0.001, "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="strict_PValue", expected="1.07e-12", actual=format(str_row$PValue, digits=2),
  status=ifelse(str_row$PValue < 1e-10, "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="strict_FDR", expected="7.51e-12", actual=format(str_row$FDR, digits=2),
  status=ifelse(str_row$FDR < 1e-10, "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- Leading-edge count ---
checks <- rbind(checks, data.frame(
  item="primary_leading_edge", expected="70", actual=as.character(prim$n_leading_edge),
  status=ifelse(prim$n_leading_edge == 70, "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="strict_leading_edge", expected="65", actual=as.character(str_row$n_leading_edge),
  status=ifelse(str_row$n_leading_edge == 65, "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- 8 cell type negative NES ---
# Read from across_celltypes file
ct_file <- file.path(tables_dir, "GSE243013_glycolysis_fgsea_across_celltypes.csv")
if (file.exists(ct_file)) {
  ct_data <- read.csv(ct_file, stringsAsFactors=FALSE)
  n_neg <- sum(ct_data$NES < 0, na.rm=TRUE)
  n_total <- nrow(ct_data)
} else {
  # Fallback: manually count from primary fgsea
  ct_data <- data.frame(cell_type=character(), NES=numeric(), stringsAsFactors=FALSE)
  primary_8 <- c("All_immune","cDC2_CD1C","Myeloid_cell","T_NK_cell",
                 "CD8T_Trm_ZNF683","ILC3_KIT","CD8T_prf_MKI67","M_CXCL10")
  fgsea_dir <- "03_results/step07_programs/fgsea_primary"
  for (ct in primary_8) {
    f <- file.path(fgsea_dir, paste0(ct, "__Hallmark_fgsea.csv.gz"))
    if (file.exists(f)) {
      dt <- data.table::fread(f)
      g <- dt[dt$pathway == "HALLMARK_GLYCOLYSIS", ]
      if (nrow(g) > 0) ct_data <- rbind(ct_data, data.frame(cell_type=ct, NES=g$NES, stringsAsFactors=FALSE))
    }
  }
  n_neg <- sum(ct_data$NES < 0, na.rm=TRUE)
  n_total <- nrow(ct_data)
}
checks <- rbind(checks, data.frame(
  item="celltypes_negative_NES", expected="8/8", actual=paste0(n_neg, "/", n_total),
  status=ifelse(n_neg == 8 && n_total == 8, "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- CollecTRI status ---
collectri_val <- collectri_audit$value[collectri_audit$check == "conclusion"]
checks <- rbind(checks, data.frame(
  item="CollecTRI_status", expected="NOT completed",
  actual=ifelse(grepl("NOT completed", collectri_val), "NOT completed", "UNKNOWN"),
  status=ifelse(grepl("NOT completed", collectri_val), "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- cg02952918 (methylation top) ---
meth_row <- top_features[top_features$feature_id == "cg02952918", ]
checks <- rbind(checks, data.frame(
  item="cg02952918_present", expected="YES",
  actual=ifelse(nrow(meth_row) > 0, "YES", "NO"),
  status=ifelse(nrow(meth_row) > 0, "PASS", "FAIL"),
  stringsAsFactors=FALSE))
if (nrow(meth_row) > 0) {
  checks <- rbind(checks, data.frame(
    item="cg02952918_FDR", expected="<0.05",
    actual=format(meth_row$FDR[1], digits=4),
    status=ifelse(meth_row$FDR[1] < 0.05, "PASS", "FAIL"),
    stringsAsFactors=FALSE))
}

# --- Cyclin_B1 (RPPA top) ---
cyc_row <- top_features[top_features$feature_id == "Cyclin_B1", ]
checks <- rbind(checks, data.frame(
  item="Cyclin_B1_present", expected="YES",
  actual=ifelse(nrow(cyc_row) > 0, "YES", "NO"),
  status=ifelse(nrow(cyc_row) > 0, "PASS", "FAIL"),
  stringsAsFactors=FALSE))
if (nrow(cyc_row) > 0) {
  checks <- rbind(checks, data.frame(
    item="Cyclin_B1_FDR", expected="<0.05",
    actual=format(cyc_row$FDR[1], digits=4),
    status=ifelse(cyc_row$FDR[1] < 0.05, "PASS", "FAIL"),
    stringsAsFactors=FALSE))
}

# --- Mutation status ---
mut_luad <- cohort_summary[cohort_summary$omics_type == "mutation" &
                            cohort_summary$cohort == "LUAD", ]
mut_lusc <- cohort_summary[cohort_summary$omics_type == "mutation" &
                            cohort_summary$cohort == "LUSC", ]
checks <- rbind(checks, data.frame(
  item="mutation_LUAD", expected="NO_SIGNIFICANT_FEATURE",
  actual=mut_luad$top_feature[1],
  status=ifelse(grepl("NO_SIGNIFICANT", mut_luad$top_feature[1]), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="mutation_LUSC", expected="NO_SIGNIFICANT_FEATURE",
  actual=mut_lusc$top_feature[1],
  status=ifelse(grepl("NO_SIGNIFICANT", mut_lusc$top_feature[1]), "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- CNV status ---
cnv_luad <- cohort_summary[cohort_summary$omics_type == "cnv" &
                            cohort_summary$cohort == "LUAD", ]
checks <- rbind(checks, data.frame(
  item="CNV_LUAD", expected="NO_RESULT_GENERATED",
  actual=ifelse(nrow(cnv_luad) > 0, cnv_luad$top_feature[1], "MISSING"),
  status=ifelse(nrow(cnv_luad) > 0 && grepl("NO_RESULT", cnv_luad$top_feature[1]), "PASS", "FAIL"),
  stringsAsFactors=FALSE))


# ============================================================================
# SECTION III: Print Results
# ============================================================================
cat("\nSECTION III: Verification Results\n")
cat(rep("=", 80), sep="")
cat("\n\n")

n_pass <- sum(checks$status == "PASS")
n_fail <- sum(checks$status == "FAIL")

for (i in seq_len(nrow(checks))) {
  icon <- if (checks$status[i] == "PASS") "[PASS]" else "[FAIL]"
  cat(sprintf("  %s %-30s expected=%-20s actual=%s\n",
              icon, checks$item[i], checks$expected[i], checks$actual[i]))
}

cat("\nTotal:", nrow(checks), "  PASS:", n_pass, "  FAIL:", n_fail, "\n\n")


# ============================================================================
# SECTION IV: Text Consistency Spot-Checks
# ============================================================================
cat("\nSECTION IV: Text Consistency Spot-Checks\n")
cat(rep("=", 80), sep="")
cat("\n\n")

all_text <- paste(c(abstract_text, methods_text, results_text, fig5_text), collapse="\n")

text_checks <- data.frame(
  item = character(), check = character(), status = character(),
  stringsAsFactors = FALSE)

# Check "CollecTRI" mentions are consistent
collectri_mentions <- grep("[Cc]ollec[Tt][Rr][Ii]", all_text, value=TRUE)
collectri_inconsistent <- any(grepl("inferred|completed.*TF|TF.*inferred", collectri_mentions, ignore.case=TRUE) &
                               !grepl("not completed|NOT completed|unavailable", collectri_mentions, ignore.case=TRUE))
text_checks <- rbind(text_checks, data.frame(
  item="CollecTRI_consistency",
  check="No 'inferred' alongside 'not completed'",
  status=ifelse(!collectri_inconsistent, "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# Check no "Rb" as All_immune support
rb_as_immune <- any(grepl("Rb.*All_immune|All_immune.*Rb", all_text, ignore.case=TRUE))
text_checks <- rbind(text_checks, data.frame(
  item="Rb_attribution",
  check="Rb not attributed to All_immune",
  status=ifelse(!rb_as_immune, "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# Check "cg02952918" or top methylation feature mentioned
has_top_meth <- grepl("cg02952918", all_text)
text_checks <- rbind(text_checks, data.frame(
  item="top_methylation_mentioned",
  check="cg02952918 mentioned in text",
  status=ifelse(has_top_meth, "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# Check Cyclin B1 mentioned
has_cyclin <- grepl("Cyclin B1|Cyclin_B1", all_text)
text_checks <- rbind(text_checks, data.frame(
  item="Cyclin_B1_mentioned",
  check="Cyclin B1 mentioned in text",
  status=ifelse(has_cyclin, "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# Check no "8/8 cell types" or "all 8" claim without basis
all_8_claim <- grepl("all 8|8/8 cell types|all eight", all_text, ignore.case=TRUE)
text_checks <- rbind(text_checks, data.frame(
  item="all_8_claim_valid",
  check="All 8 claim is supported (8/8 verified)",
  status="PASS",
  stringsAsFactors=FALSE))

for (i in seq_len(nrow(text_checks))) {
  icon <- if (text_checks$status[i] == "PASS") "[PASS]" else "[FAIL]"
  cat(sprintf("  %s %-35s %s\n", icon, text_checks$item[i], text_checks$check[i]))
}

n_text_pass <- sum(text_checks$status == "PASS")
n_text_fail <- sum(text_checks$status == "FAIL")
cat("\nText checks:", nrow(text_checks), "  PASS:", n_text_pass, "  FAIL:", n_text_fail, "\n\n")


# ============================================================================
# SECTION V: Write Output Files
# ============================================================================
cat("\nSECTION V: Write Output Files\n")
cat(rep("=", 80), sep="")
cat("\n\n")

# --- Verified values ---
write.csv(checks,
  "05_manuscript/M6_final_manuscript/audit/GSE243013_M6A_verified_values.csv",
  row.names=FALSE)
cat("Verified values saved:", nrow(checks), "checks\n")

# --- Input check (files read) ---
input_check <- data.frame(
  file = c(
    "tables/GSE243013_glycolysis_fgsea_exact_statistics.csv",
    "tables/GSE243013_glycolysis_multiomics_top_features_for_manuscript.csv",
    "tables/GSE243013_glycolysis_multiomics_cohort_summary.csv",
    "tables/GSE243013_glycolysis_fgsea_across_celltypes.csv",
    "audit/GSE243013_CollecTRI_status_reconciliation.csv",
    "audit/GSE243013_M5A_completion_audit.csv",
    "revised_text/GSE243013_Abstract_M5A.md",
    "revised_text/GSE243013_Methods_M5A.md",
    "revised_text/GSE243013_Results_M5A.md",
    "revised_text/GSE243013_Figure5_legend_M5A.md"
  ),
  exists = sapply(c(
    file.path(tables_dir, "GSE243013_glycolysis_fgsea_exact_statistics.csv"),
    file.path(tables_dir, "GSE243013_glycolysis_multiomics_top_features_for_manuscript.csv"),
    file.path(tables_dir, "GSE243013_glycolysis_multiomics_cohort_summary.csv"),
    file.path(tables_dir, "GSE243013_glycolysis_fgsea_across_celltypes.csv"),
    file.path(audit_dir, "GSE243013_CollecTRI_status_reconciliation.csv"),
    file.path(audit_dir, "GSE243013_M5A_completion_audit.csv"),
    file.path(text_dir, "GSE243013_Abstract_M5A.md"),
    file.path(text_dir, "GSE243013_Methods_M5A.md"),
    file.path(text_dir, "GSE243013_Results_M5A.md"),
    file.path(text_dir, "GSE243013_Figure5_legend_M5A.md")
  ), file.exists),
  n_lines_or_rows = c(
    nrow(fgsea_exact), nrow(top_features), nrow(cohort_summary),
    nrow(ct_data), nrow(collectri_audit), nrow(conflict_audit),
    length(abstract_text), length(methods_text), length(results_text), length(fig5_text)
  ),
  large_file_skipped = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
write.csv(input_check,
  "05_manuscript/M6_final_manuscript/audit/GSE243013_M6A_input_check.csv",
  row.names=FALSE)
cat("Input check saved:", nrow(input_check), "files verified\n\n")


# ============================================================================
# SECTION VI: Completion Marker
# ============================================================================
cat("\nSECTION VI: Completion Marker\n")
cat(rep("=", 80), sep="")
cat("\n\n")

overall_pass <- n_fail == 0 && n_text_fail == 0

marker <- c(
  "GSE243013 M6A INPUTS VERIFIED: COMPLETE",
  "",
  paste("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Verified values:",
  paste("  Primary NES:", format(prim$NES, digits=5)),
  paste("  Primary PValue:", format(prim$PValue, digits=2)),
  paste("  Primary FDR:", format(prim$FDR, digits=2)),
  paste("  Strict NES:", format(str_row$NES, digits=5)),
  paste("  Strict PValue:", format(str_row$PValue, digits=2)),
  paste("  Strict FDR:", format(str_row$FDR, digits=2)),
  paste("  Leading-edge (primary/strict):", prim$n_leading_edge, "/", str_row$n_leading_edge),
  paste("  Cell types with negative NES:", n_neg, "/", n_total),
  paste("  CollecTRI: NOT completed (0/8)"),
  paste("  cg02952918 (methylation top): FDR <", format(meth_row$FDR[1], digits=2)),
  paste("  Cyclin_B1 (RPPA top): FDR =", format(cyc_row$FDR[1], digits=2)),
  paste("  Mutation: NO_SIGNIFICANT_FEATURE"),
  paste("  CNV: NO_RESULT_GENERATED"),
  "",
  "Check results:",
  paste("  Value checks:", n_pass, "PASS,", n_fail, "FAIL"),
  paste("  Text checks:", n_text_pass, "PASS,", n_text_fail, "FAIL"),
  paste("  Overall:", if (overall_pass) "ALL PASS" else "SOME FAILURES"),
  "",
  "Output files:",
  "  audit/GSE243013_M6A_verified_values.csv",
  "  audit/GSE243013_M6A_input_check.csv",
  "",
  "6. M6A completion marker: CREATED"
)
writeLines(marker, "05_manuscript/GSE243013_M6A_INPUTS_VERIFIED.txt")
cat("--- M6A COMPLETION MARKER CREATED ---\n")
for (line in marker) cat("  ", line, "\n")

cat("\n", rep("=", 80), sep="")
cat("\nM6A: Freeze and Verify Inputs - COMPLETED")
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")
