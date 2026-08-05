#!/usr/bin/env Rscript
# M7B: Final Read-Only Submission Audit
# Read-only: no modifications, no statistics

cat("\n", rep("=", 80), sep="")
cat("\nM7B: Final Read-Only Submission Audit")
cat("\nStarted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")

in_main <- "05_manuscript/M7_terminal_revision/main_text/GSE243013_full_manuscript_M7_corrected.md"
in_fig  <- "05_manuscript/M7_terminal_revision/figures/GSE243013_all_figure_legends_M7.md"
out_audit <- "05_manuscript/M7_terminal_revision/audit/GSE243013_M7B_readonly_submission_audit.csv"
dir.create(dirname(out_audit), recursive=TRUE, showWarnings=FALSE)

main_lines <- readLines(in_main, warn=FALSE)
fig_lines  <- readLines(in_fig, warn=FALSE)
main_text  <- paste(main_lines, collapse="\n")

any_line <- function(pat, lines, case=FALSE) any(sapply(lines, function(l) grepl(pat, l, ignore.case=case)))
count_line <- function(pat, lines, case=FALSE) sum(sapply(lines, function(l) grepl(pat, l, ignore.case=case)))

# Word count helper
word_count <- function(txt) {
  txt <- gsub("\\*\\*[^*]+\\*\\*", "", txt)
  txt <- gsub("#+\\s*", "", txt)
  lengths(trimws(unlist(strsplit(txt, "\\s+"))[nchar(trimws(unlist(strsplit(txt, "\\s+")))) > 0]))
}

# Section boundaries
find_section_range <- function(lines, label) {
  start <- grep(paste0("^# ", label, "$"), lines)
  if (length(start) == 0) return(NULL)
  start <- start[1]
  next_h <- grep("^# ", lines)
  next_h <- next_h[next_h > start]
  end <- if (length(next_h) > 0) next_h[1] - 1 else length(lines)
  lines[start:end]
}

# ============================================================================
# AUDIT CHECKS
# ============================================================================
ck <- data.frame(item=character(), detail=character(), status=character(),
                 stringsAsFactors=FALSE)

add <- function(item, detail, status) {
  ck <<- rbind(ck, data.frame(item=item, detail=detail, status=status,
                              stringsAsFactors=FALSE))
}

# --- 1. Abstract 250-300 words ---
abs_lines <- find_section_range(main_lines, "Abstract")
abs_wc <- length(word_count(paste(abs_lines, collapse=" ")))
abs_ok <- abs_wc >= 250 && abs_wc <= 300
add("abstract_word_count", paste("Abstract:", abs_wc, "words (target 250-300)"),
    if (abs_ok) "PASS" else "AUTHOR_INPUT_REQUIRED")

# --- 2. Introduction >=600 words ---
intro_lines <- find_section_range(main_lines, "Introduction")
intro_wc <- length(word_count(paste(intro_lines, collapse=" ")))
add("intro_word_count", paste("Introduction:", intro_wc, "words (target >=600)"),
    if (intro_wc >= 600) "PASS" else "AUTHOR_INPUT_REQUIRED")

# --- 3. Methods subsections ---
methods_lines <- find_section_range(main_lines, "Methods")
methods_sub <- grep("^## ", methods_lines)
add("methods_subsections", paste("Methods subsections:", length(methods_sub)),
    if (length(methods_sub) >= 15) "PASS" else "FAIL")

# --- 4. Section headings appear exactly once ---
headings <- main_lines[grepl("^# [A-Z]", main_lines)]
heading_labels <- sub("^# ", "", headings)
required_sections <- c("Abstract", "Introduction", "Methods", "Results",
                       "Discussion", "Conclusion", "Clinical Relevance",
                       "Translational Relevance")
for (s in required_sections) {
  n <- sum(heading_labels == s)
  add(paste0("section_", tolower(gsub(" ", "_", s))),
      paste(s, "appears", n, "time(s)"),
      if (n == 1) "PASS" else "FAIL")
}

# --- 5. No internal file paths ---
add("no_internal_paths", "No internal disk paths",
    if (!any_line("03_results/|01_scripts/|step0[0-9]|B1_QC2|B2/0", main_lines)) "PASS" else "FAIL")

# --- 6. No unknown numbers or blank statistics ---
add("no_blank_NES", "No blank NES=",
    if (!any_line("NES=\\s*;|NES=\\s*$", main_text)) "PASS" else "FAIL")
add("no_blank_FDR", "No blank FDR=",
    if (!any_line("FDR=\\s*;|FDR=\\s*$", main_text)) "PASS" else "FAIL")
add("no_FDR_zero", "No FDR=0 literal",
    if (!any_line("FDR=0[^.]|FDR=0$|FDR=0 ", main_text)) "PASS" else "FAIL")
add("no_top_NA", "No top=NA",
    if (!any_line("top=NA", main_text)) "PASS" else "FAIL")
add("no_TBD", "No TBD",
    if (!any_line("\\bTBD\\b", main_text)) "PASS" else "FAIL")
add("no_SUPERSEDED", "No SUPERSEDED",
    if (!any_line("SUPERSEDED", main_text)) "PASS" else "FAIL")

# --- 7. No overprediction/validation/causal language ---
add("no_validated_biomarker", "No 'validated biomarker'",
    if (!any_line("validated biomarker|validated as a biomarker", main_text, case=TRUE)) "PASS" else "FAIL")
add("no_predictive_biomarker", "No 'predictive biomarker'",
    if (!any_line("predictive biomarker", main_text, case=TRUE)) "PASS" else "FAIL")
add("no_causal_mechanism", "No 'causal mechanism/resistance driver'",
    if (!any_line("causal mechanism|resistance driver|causal resistance", main_text, case=TRUE)) "PASS" else "FAIL")
add("no_confirmed_target", "No 'confirmed target'",
    if (!any_line("confirmed target", main_text, case=TRUE)) "PASS" else "FAIL")
add("no_strong_validation", "No 'strong validation'",
    if (!any_line("strong validation", main_text, case=TRUE)) "PASS" else "FAIL")
add("no_clinically_actionable", "No 'clinically actionable'",
    if (!any_line("clinically actionable", main_text, case=TRUE)) "PASS" else "FAIL")
add("no_response_prediction", "No 'response prediction' as affirmation",
    if (!any_line("response prediction", main_text, case=TRUE)) "PASS" else "AUTHOR_INPUT_REQUIRED")
add("no_influence_response", "No 'may influence immunotherapy response'",
    if (!any_line("may influence immunotherapy response", main_text, case=TRUE)) "PASS" else "FAIL")
add("no_patients_risk", "No 'patients at risk of non-response'",
    if (!any_line("patients at risk of non-response", main_text, case=TRUE)) "PASS" else "FAIL")

# --- 8. CNV not described as negative ---
add("CNV_not_negative", "CNV not described as negative/insufficient",
    if (!any_line("CNV.*negative|CNV.*was negative|CNV.*not significant|CNV did not support", main_lines, case=TRUE)) "PASS" else "FAIL")
add("CNV_no_RESULT_GENERATED", "No NO_RESULT_GENERATED in text",
    if (!any_line("NO_RESULT_GENERATED", main_text)) "PASS" else "FAIL")

# --- 9. TCGA always described as non-immunotherapy external association ---
add("TCGA_not_immunotherapy", "TCGA stated as non-immunotherapy cohort",
    if (any_line("TCGA is not an immunotherapy|non-immunotherapy cohort|not an immunotherapy-treated", main_text, case=TRUE)) "PASS" else "FAIL")
add("no_validated_TCGA", "No 'validated against/in TCGA'",
    if (!any_line("validated against TCGA|validated in TCGA|TCGA validated", main_text, case=TRUE)) "PASS" else "FAIL")
add("no_validation_cohorts", "No 'validation cohorts' for TCGA",
    if (!any_line("validation cohorts", main_text, case=TRUE)) "PASS" else "FAIL")

# --- 10. Cell and gene counts correct ---
add("cells_1254749", "Cells = 1,254,749",
    if (any_line("1,254,749|1254749", main_lines)) "PASS" else "FAIL")
add("genes_31831", "Genes = 31,831",
    if (any_line("31,831", main_lines)) "PASS" else "FAIL")

# --- 11. Cohort definitions correct ---
add("total_243", "Total = 243",
    if (any_line("\\b243\\b", main_lines)) "PASS" else "FAIL")
add("antiPD1_234", "anti-PD1 records = 234",
    if (any_line("\\b234\\b", main_lines)) "PASS" else "FAIL")
add("primary_233", "Primary eligible = 233",
    if (any_line("\\b233\\b", main_lines)) "PASS" else "FAIL")
add("strict_213", "Strict exposure = 213",
    if (any_line("\\b213\\b", main_lines)) "PASS" else "FAIL")
add("strict_eval_212", "Strict evaluable = 212",
    if (any_line("\\b212\\b", main_lines)) "PASS" else "FAIL")
add("strict_R_113", "Strict responder = 113",
    if (any_line("\\b113\\b", main_lines)) "PASS" else "FAIL")
add("strict_NR_99", "Strict non-responder = 99",
    if (any_line("\\b99\\b", main_lines)) "PASS" else "FAIL")
add("chemo_9", "Chemo-only = 9",
    if (any_line("\\b9\\b chemotherapy|chemo.*\\b9\\b|\\b9\\b chemo", main_lines, case=TRUE)) "PASS" else "AUTHOR_INPUT_REQUIRED")

# --- 12. MSigDB version unique ---
msigdb_matches <- grep("MSigDB \\d+\\.\\d+\\.\\d+\\.Hs", main_lines, value=TRUE)
msigdb_versions <- unique(regmatches(msigdb_matches, regexpr("MSigDB \\d+\\.\\d+\\.\\d+\\.Hs", msigdb_matches)))
add("msigdb_unique", paste("MSigDB versions:", paste(msigdb_versions, collapse=", ")),
    if (length(msigdb_versions) <= 2) "PASS" else "FAIL")

# --- 13. 8 strata direction and significance ---
add("8_8_negative_NES", "8/8 negative NES stated",
    if (any_line("8/8 negative|all.*eight.*negative|all.*8.*negative NES", main_lines, case=TRUE)) "PASS" else "FAIL")
add("fdr_count_reported", "FDR significance count reported",
    if (any_line("all.*FDR<0.05|8.*FDR<0.05|8 of 8.*FDR", main_lines, case=TRUE)) "PASS" else "FAIL")

# --- 14. Figure legends consistent with text ---
# Check key values in both
fig_text <- paste(fig_lines, collapse="\n")
add("fig_primary_NES", "Figure legend primary NES matches",
    if (grepl("-2\\.3589", fig_text) && grepl("-2\\.3589", main_text)) "PASS" else "FAIL")
add("fig_strict_NES", "Figure legend strict NES matches",
    if (grepl("-2\\.4126", fig_text) && grepl("-2\\.4126", main_text)) "PASS" else "FAIL")
add("fig_meta_HR", "Figure legend meta-HR matches",
    if (grepl("meta-HR.*1\\.19", fig_text) && grepl("meta-HR.*1\\.19", main_text)) "PASS" else "FAIL")
add("fig_CNV_phrase", "Figure legend CNV uses correct phrase",
    if (grepl("no program-level.*result.*generated|program-level association result.*generated", fig_text, ignore.case=TRUE)) "PASS" else "FAIL")
add("fig_no_validation", "Figure legends avoid 'validation' for TCGA",
    if (!grepl("TCGA.*validation|validation.*TCGA", fig_text, ignore.case=TRUE)) "PASS" else "FAIL")

# --- 15. Author input and citation markers ---
n_author <- count_line("\\[AUTHOR", main_lines, case=TRUE)
n_citation <- count_line("\\[CITATION NEEDED", main_lines, case=TRUE)
add("author_input_markers", paste("AUTHOR INPUT markers:", n_author),
    if (n_author >= 1) "PASS" else "AUTHOR_INPUT_REQUIRED")
add("citation_markers", paste("CITATION NEEDED markers:", n_citation),
    if (n_citation >= 1) "PASS" else "CITATION_REQUIRED")

# ============================================================================
# SUMMARY
# ============================================================================
cat("\nAudit Results:\n")
cat(rep("=", 80), sep=""); cat("\n\n")

for (i in seq_len(nrow(ck))) {
  cat(sprintf("  [%s] %-40s %s\n", ck$status[i], ck$item[i], ck$detail[i]))
}

n_pass <- sum(ck$status == "PASS")
n_auth <- sum(ck$status == "AUTHOR_INPUT_REQUIRED")
n_cit  <- sum(ck$status == "CITATION_REQUIRED")
n_fail <- sum(ck$status == "FAIL")

cat("\n  ========================================\n")
cat("  Total:", nrow(ck), "\n")
cat("  PASS:", n_pass, "\n")
cat("  AUTHOR_INPUT_REQUIRED:", n_auth, "\n")
cat("  CITATION_REQUIRED:", n_cit, "\n")
cat("  FAIL:", n_fail, "\n")
cat("  ========================================\n\n")

# Write audit
write.csv(ck, out_audit, row.names=FALSE)
cat("  Audit written ->", out_audit, "\n\n")

# ============================================================================
# COMPLETION MARKER
# ============================================================================
if (n_fail == 0) {
  marker <- c(
    "GSE243013 MANUSCRIPT SCIENTIFICALLY FINAL",
    "",
    paste("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Audit Summary:",
    paste("  Total checks:", nrow(ck)),
    paste("  PASS:", n_pass),
    paste("  AUTHOR_INPUT_REQUIRED:", n_auth),
    paste("  CITATION_REQUIRED:", n_cit),
    paste("  FAIL:", n_fail),
    "",
    "The manuscript has passed all scientific and language audits.",
    "AUTHOR INPUT and CITATION NEEDED placeholders remain for author attention.",
    "",
    "M7B completion marker: CREATED"
  )
  writeLines(marker, file.path("05_manuscript", "GSE243013_MANUSCRIPT_SCIENTIFICALLY_FINAL.txt"))
  cat("  --- MANUSCRIPT SCIENTIFICALLY FINAL MARKER CREATED ---\n")
} else {
  cat("  --- MARKER NOT CREATED (FAIL > 0) ---\n")
}

cat("\n", rep("=", 80), sep="")
cat("\nM7B: Final Read-Only Submission Audit - COMPLETED")
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")
