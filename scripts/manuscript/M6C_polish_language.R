#!/usr/bin/env Rscript
# M6C: English Language Polish
# Read-only: reads M6B integrated manuscript + figure legends
# Only modifies English grammar, tense, redundancy, transitions, terminology, punctuation, abbreviations
# No statistical analysis, no number changes, no new biological claims

cat("\n", rep("=", 80), sep="")
cat("\nM6C: English Language Polish")
cat("\nStarted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")

in_main  <- "05_manuscript/M6_final_manuscript/main_text/GSE243013_manuscript_M6B_integrated.md"
in_fig   <- "05_manuscript/M6_final_manuscript/figures/GSE243013_figure_legends_M6B.md"
out_main <- "05_manuscript/M6_final_manuscript/main_text/GSE243013_full_manuscript_M6_polished.md"
out_fig  <- "05_manuscript/M6_final_manuscript/figures/GSE243013_all_figure_legends_M6_polished.md"
out_log  <- "05_manuscript/M6_final_manuscript/audit/GSE243013_M6C_language_change_log.csv"

# ============================================================================
# SECTION I: Read Input Files
# ============================================================================
cat("\nSECTION I: Read Input Files\n")
cat(rep("=", 80), sep=""); cat("\n\n")

main_lines <- readLines(in_main, warn=FALSE)
fig_lines  <- readLines(in_fig, warn=FALSE)
cat("  Main manuscript:", length(main_lines), "lines\n")
cat("  Figure legends: ", length(fig_lines), "lines\n\n")

# ============================================================================
# SECTION II: Line-by-Line Language Corrections
# ============================================================================
cat("\nSECTION II: Apply Language Corrections\n")
cat(rep("=", 80), sep=""); cat("\n\n")

changes <- data.frame(
  section = character(), line_no = integer(),
  original_text = character(), revised_text = character(),
  change_type = character(), numeric_unchanged = logical(),
  stringsAsFactors = FALSE
)

log_change <- function(section, line_no, original, revised, change_type) {
  changes <<- rbind(changes, data.frame(
    section = section, line_no = line_no,
    original_text = original, revised_text = revised,
    change_type = change_type, numeric_unchanged = TRUE,
    stringsAsFactors = FALSE
  ))
}

detect_section <- function(lines, idx) {
  for (i in idx:1) {
    if (grepl("^# ", lines[i])) return(trimws(gsub("^#+\\s*", "", lines[i])))
  }
  return("Unknown")
}

# Helper: exact string replacement on a single line (no numbers touched)
fix_line <- function(lines, line_no, old, new, change_type) {
  if (line_no > length(lines)) return(lines)
  if (!grepl(old, lines[line_no], fixed=TRUE)) return(lines)
  sec <- detect_section(lines, line_no)
  orig <- lines[line_no]
  lines[line_no] <- gsub(old, new, lines[line_no], fixed=TRUE)
  log_change(sec, line_no, orig, lines[line_no], change_type)
  lines
}

# Apply fixes to main manuscript
m <- main_lines

# --- Abstract ---
m <- fix_line(m, 14, "determinants of response remain incompletely characterized",
              "determinants of response remain poorly understood",
              "redundancy")
m <- fix_line(m, 22, "validated against clinical outcomes in",
              "assessed for association with clinical outcomes in",
              "terminology")
m <- fix_line(m, 31, "meta-HR=1.19, meta-FDR=0.0499), driven primarily by",
              "meta-HR = 1.19, meta-FDR = 0.0499), driven primarily by",
              "formatting")
m <- fix_line(m, 32, "I2=90.4%, P=0.0012).",
              "I2 = 90.4%, P = 0.0012).",
              "formatting")

# --- Introduction ---
m <- fix_line(m, 44, "Analyses that treat individual cells as independent biological replicates can",
              "Analytical approaches that treat individual cells as independent biological replicates can",
              "grammar")
m <- fix_line(m, 47, "pathological-response-associated in non-small cell lung cancer (NSCLC) patients",
              "associated with pathological response in non-small cell lung cancer (NSCLC) patients",
              "grammar")

# --- Methods ---
m <- fix_line(m, 57, "sign(log2FC) * sqrt(F)",
              "sign(log2FC) \u00d7 sqrt(F)",
              "formatting")
m <- fix_line(m, 62, "HALLMARK_GLYCOLYSIS was tested in All_immune as the core program.",
              "HALLMARK_GLYCOLYSIS was tested in the All_immune aggregate as the primary program of interest.",
              "precision")
m <- fix_line(m, 65, "520 patients, 515 with RNA",
              "520 patients, 515 with RNA-seq data",
              "clarity")
m <- fix_line(m, 66, "504 patients, 501 with RNA",
              "504 patients, 501 with RNA-seq data",
              "clarity")
m <- fix_line(m, 67, "ties=efron).",
              "ties = efron).",
              "formatting")
m <- fix_line(m, 74, "somatic mutation burden (Spearman correlation and Cohen's d for driver mutations)",
              "somatic mutation burden (Spearman correlation and Cohen\u2019s d for individual driver mutations)",
              "clarity")
m <- fix_line(m, 80, "TF activity was therefore not incorporated into the final program interpretation.",
              "TF activity was therefore excluded from the final program interpretation.",
              "redundancy")

# --- Results ---
m <- fix_line(m, 87, "P=4.90 x 10^-12, FDR=3.00 x 10^-11, ES=-0.4859,",
              "P = 4.90 \u00d7 10^-12, FDR = 3.00 \u00d7 10^-11, ES = -0.4859,",
              "formatting")
m <- fix_line(m, 88, "FDR=7.51 x 10^-12, leading-edge=65 genes.",
              "FDR = 7.51 \u00d7 10^-12, leading-edge = 65 genes.",
              "formatting")
m <- fix_line(m, 91, "Directionality was consistent across all 8 analyzed strata (8/8 negative NES).",
              "Directionality was consistent across all 8 analyzed strata (8/8 negative NES).",
              "no change")
m <- fix_line(m, 93, "All 8 strata showed negative NES, indicating higher glycolysis in non-responders.",
              "All 8 strata exhibited negative NES, indicating higher glycolysis-associated transcription in non-responders.",
              "precision")
m <- fix_line(m, 94, "We note that negative NES indicates enrichment in the non-responder direction;",
              "Negative NES indicates enrichment in the non-responder direction,",
              "redundancy")
m <- fix_line(m, 95, "this does not imply statistical significance for all individual strata (e.g., ILC3_KIT FDR=0.034).",
              "which does not imply statistical significance for all individual strata (e.g., ILC3_KIT FDR = 0.034).",
              "grammar")
m <- fix_line(m, 104, "The absence of association in LUSC and I2 above 90% argue against",
              "The absence of an association in LUSC and I2 exceeding 90% argue against",
              "grammar")

# --- Discussion ---
m <- fix_line(m, 130, "direction in both primary (NES=-2.3589, FDR=3.00 x 10^-11) and strict (NES=-2.4126, FDR=7.51 x 10^-12) cohorts.",
              "direction in both primary (NES = -2.3589, FDR = 3.00 \u00d7 10^-11) and strict (NES = -2.4126, FDR = 7.51 \u00d7 10^-12) cohorts.",
              "formatting")
m <- fix_line(m, 133, "with similar effect sizes (ES approximately -0.49).",
              "with similar effect sizes (ES approximately \u22120.49).",
              "formatting")
m <- fix_line(m, 140, "but the cell-type origin of this signal remains uncertain.",
              "but the cell-type origin of this signal cannot be resolved from the All_immune aggregate.",
              "precision")
m <- fix_line(m, 144, "the cell-type origin of this signal remains uncertain.",
              "the cell-type origin of this signal cannot be resolved from the All_immune aggregate.",
              "precision")
m <- fix_line(m, 148, "Previous studies have reported associations between tumor glycolysis and immunotherapy resistance.",
              "Prior studies have reported associations between tumor glycolysis and immunotherapy resistance.",
              "terminology")
m <- fix_line(m, 150, "but causality cannot be established from this observational design.",
              "but causality cannot be inferred from this observational study.",
              "terminology")
m <- fix_line(m, 157, "The I2 statistic was approximately 90%, indicating substantial heterogeneity.",
              "The I2 statistic was approximately 90%, indicating substantial between-cohort heterogeneity.",
              "precision")
m <- fix_line(m, 158, "The fixed-effect summary therefore has descriptive rather than inferential significance:",
              "The fixed-effect summary therefore has descriptive rather than inferential value:",
              "precision")
m <- fix_line(m, 159, "the pooled estimate assumes a common true effect across cohorts, which is unlikely given the histology-specific pattern.",
              "the pooled estimate assumes a common true effect across histologies, which is unsupported given the LUAD-specific pattern.",
              "precision")
m <- fix_line(m, 162, "TCGA is not an immunotherapy-treated cohort; these associations reflect general cancer biology",
              "TCGA is not an immunotherapy-treated cohort; these associations reflect general cancer biology",
              "no change")
m <- fix_line(m, 166, "Exploratory feature-level multi-omics integration revealed histology-specific patterns.",
              "Exploratory feature-level multi-omics integration revealed histology-specific association patterns.",
              "precision")
m <- fix_line(m, 171, "Cyclin B1 is a proliferation marker and may reflect tumor cell proliferation, purity, or compositional confounding rather than a direct glycolysis-proteome relationship.",
              "Cyclin B1 is a proliferation marker; the correlation may reflect tumor cell proliferation, purity, or compositional confounding rather than a direct glycolysis\u2013proteome relationship.",
              "grammar")
m <- fix_line(m, 173, "The cg02952918 CpG requires independent annotation and functional validation before biological interpretation.",
              "The cg02952918 CpG requires independent annotation and functional follow-up before biological interpretation.",
              "terminology")
m <- fix_line(m, 177, "this reflects absence of analysis output, not a statistically negative finding.",
              "this reflects the absence of analysis output rather than a statistically negative finding.",
              "grammar")
m <- fix_line(m, 183, "TF activity data are therefore unavailable and are not reported.",
              "TF activity data are therefore unavailable and are not presented.",
              "terminology")
m <- fix_line(m, 188, "the use of preranked fgsea (avoiding ssGSEA normalization artifacts),",
              "the use of preranked fgsea (avoiding normalization artifacts inherent to ssGSEA),",
              "clarity")
m <- fix_line(m, 197, "the I2 of approximately 90% limits the interpretability of the meta-analytic summary",
              "the I2 of approximately 90% limits the interpretability of the fixed-effect meta-analytic summary",
              "precision")
m <- fix_line(m, 206, "the TCGA-LUAD prognostic association provide convergent evidence, though substantial",
              "the TCGA-LUAD prognostic association provide converging evidence, though substantial",
              "terminology")
m <- fix_line(m, 210, "These findings require prospective validation with cell-type-resolved profiling",
              "These findings require prospective validation using cell-type-resolved profiling",
              "grammar")
m <- fix_line(m, 211, "and pre-treatment samples to establish clinical utility.",
              "and pre-treatment samples to establish potential clinical utility.",
              "precision")
m <- fix_line(m, 215, "The immune-compartment glycolysis program may serve as a candidate biomarker for immunotherapy response prediction",
              "The immune-compartment glycolysis program represents a candidate correlate of immunotherapy response",
              "terminology")
m <- fix_line(m, 216, "though current evidence is limited by the post-treatment sampling design and histology-specific effects.",
              "though current evidence is limited by the post-treatment sampling design, histology-specific effects, and the absence of pre-treatment profiling.",
              "completeness")
m <- fix_line(m, 218, "Clinical application would require pre-treatment validation and cell-type-resolved measurement.",
              "Clinical application would require pre-treatment validation and cell-type-resolved quantification.",
              "terminology")
m <- fix_line(m, 219, "The TCGA prognostic association (meta-FDR=0.0499) is hypothesis-generating but should not be interpreted as clinical validation",
              "The TCGA prognostic association (meta-FDR = 0.0499) is hypothesis-generating but should not be interpreted as clinical validation",
              "formatting")
m <- fix_line(m, 225, "These findings suggest that immune-cell glycolytic metabolism may influence immunotherapy response in NSCLC.",
              "These findings suggest that immune-cell glycolytic metabolism may be associated with immunotherapy response in NSCLC.",
              "precision")
m <- fix_line(m, 227, "metabolic profiling of pre-treatment biopsies may identify patients at risk of non-response",
              "metabolic profiling of pre-treatment biopsies may identify patients at higher risk of non-response",
              "precision")

# --- Figure 5 ---
m <- fix_line(m, 245, "P=4.90 x 10^-12, FDR=3.00 x 10^-11,",
              "P = 4.90 \u00d7 10^-12, FDR = 3.00 \u00d7 10^-11,",
              "formatting")
m <- fix_line(m, 246, "P=1.07 x 10^-12, FDR=7.51 x 10^-12,",
              "P = 1.07 \u00d7 10^-12, FDR = 7.51 \u00d7 10^-12,",
              "formatting")
m <- fix_line(m, 255, "TCGA-LUAD (HR=1.47, 95% CI 1.24-1.73) and TCGA-LUSC (HR=1.02, 95% CI 0.89-1.18).",
              "TCGA-LUAD (HR = 1.47, 95% CI 1.24\u20131.73) and TCGA-LUSC (HR = 1.02, 95% CI 0.89\u20131.18).",
              "formatting")
m <- fix_line(m, 256, "meta-HR=1.19, meta-FDR=0.0499.",
              "meta-HR = 1.19, meta-FDR = 0.0499.",
              "formatting")
m <- fix_line(m, 259, "FDR<2.6 x 10^-7);",
              "FDR < 2.6 \u00d7 10^-7);",
              "formatting")
m <- fix_line(m, 261, "rho=0.565, FDR below machine-reportable precision);",
              "rho = 0.565, FDR below machine-reportable precision);",
              "formatting")

# --- Figure Legends ---
m <- fix_line(m, 282, "EdgeR differential expression summary across 47 cell types.",
              "EdgeR differential expression summary across 47 cell types.",
              "no change")
m <- fix_line(m, 309, "I2=90.4%).",
              "I2 = 90.4%).",
              "formatting")
m <- fix_line(m, 315, "core glycolysis-related program",
              "core glycolysis-associated program",
              "terminology")
m <- fix_line(m, 317, "top: Cyclin B1, rho=0.565)",
              "top: Cyclin B1, rho = 0.565)",
              "formatting")

# --- Figure Legends (separate file) ---
f <- fig_lines

f <- fix_line(f, 13, "P=4.90 x 10^-12, FDR=3.00 x 10^-11,",
              "P = 4.90 \u00d7 10^-12, FDR = 3.00 \u00d7 10^-11,",
              "formatting")
f <- fix_line(f, 14, "P=1.07 x 10^-12, FDR=7.51 x 10^-12,",
              "P = 1.07 \u00d7 10^-12, FDR = 7.51 \u00d7 10^-12,",
              "formatting")
f <- fix_line(f, 23, "TCGA-LUAD (HR=1.47, 95% CI 1.24-1.73) and TCGA-LUSC (HR=1.02, 95% CI 0.89-1.18).",
              "TCGA-LUAD (HR = 1.47, 95% CI 1.24\u20131.73) and TCGA-LUSC (HR = 1.02, 95% CI 0.89\u20131.18).",
              "formatting")
f <- fix_line(f, 24, "meta-HR=1.19, meta-FDR=0.0499.",
              "meta-HR = 1.19, meta-FDR = 0.0499.",
              "formatting")
f <- fix_line(f, 27, "FDR<2.6 x 10^-7);",
              "FDR < 2.6 \u00d7 10^-7);",
              "formatting")
f <- fix_line(f, 29, "rho=0.565, FDR below machine-reportable precision);",
              "rho = 0.565, FDR below machine-reportable precision);",
              "formatting")
f <- fix_line(f, 50, "EdgeR differential expression summary across 47 cell types.",
              "EdgeR differential expression summary across 47 cell types.",
              "no change")
f <- fix_line(f, 77, "I2=90.4%).",
              "I2 = 90.4%).",
              "formatting")
f <- fix_line(f, 83, "core glycolysis-related program",
              "core glycolysis-associated program",
              "terminology")
f <- fix_line(f, 85, "top: Cyclin B1, rho=0.565)",
              "top: Cyclin B1, rho = 0.565)",
              "formatting")

# Remove the "no change" entries
changes <- changes[changes$change_type != "no change", ]

cat("  Total changes applied:", nrow(changes), "\n\n")

# ============================================================================
# SECTION III: Write Output Files
# ============================================================================
cat("\nSECTION III: Write Output Files\n")
cat(rep("=", 80), sep=""); cat("\n\n")

dir.create(dirname(out_main), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(out_fig), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(out_log), recursive=TRUE, showWarnings=FALSE)

writeLines(m, out_main)
writeLines(f, out_fig)
write.csv(changes, out_log, row.names=FALSE)

cat("  -> main_text/GSE243013_full_manuscript_M6_polished.md  (", length(m), " lines)\n", sep="")
cat("  -> figures/GSE243013_all_figure_legends_M6_polished.md  (", length(f), " lines)\n", sep="")
cat("  -> audit/GSE243013_M6C_language_change_log.csv  (", nrow(changes), " changes)\n", sep="")

# ============================================================================
# SECTION IV: Verification
# ============================================================================
cat("\nSECTION IV: Verification\n")
cat(rep("=", 80), sep=""); cat("\n\n")

# Extract all numeric strings from original and polished
# Normalize Unicode separators first for comparison
normalize_nums <- function(txt) {
  txt <- gsub("\u00d7", "x", txt, fixed=TRUE)  # × -> x
  txt <- gsub("\u2212", "-", txt, fixed=TRUE)  # − -> -
  txt <- gsub("\u2013", "-", txt, fixed=TRUE)  # – -> -
  txt <- gsub("\u2014", "-", txt, fixed=TRUE)  # — -> -
  txt
}
extract_nums <- function(lines) {
  txt <- normalize_nums(paste(lines, collapse=" "))
  m <- gregexpr("-?\\d+\\.?\\d*(?:e[+-]?\\d+)?|\\d+\\.?\\d*\\s*[xX]\\s*10\\^\\s*-?\\d+", txt)
  unlist(regmatches(txt, m))
}

orig_nums <- extract_nums(main_lines)
pol_nums  <- extract_nums(m)

# Check numeric preservation
nums_match <- all(sort(orig_nums) == sort(pol_nums))

# Check line counts
lines_ok <- length(m) == length(main_lines)

# Check section counts
orig_sections <- sum(grepl("^# ", main_lines))
pol_sections  <- sum(grepl("^# ", m))
sections_ok   <- orig_sections == pol_sections

# Check AUTHOR INPUT preserved
has_author <- any(grepl("AUTHOR INPUT|\\[AUTHOR", m, ignore.case=TRUE))

# Check no new statistics
no_new_stats <- !any(grepl("P\\s*=\\s*0\\.\\d{4,}|FDR\\s*=\\s*0\\.\\d{4,}", m) &
                     !grepl("P\\s*=\\s*0\\.\\d{4,}|FDR\\s*=\\s*0\\.\\d{4,}", main_lines))

# Check no forbidden terms
forbidden <- c("validated biomarker", "predictive biomarker", "causal mechanism",
               "resistance driver", "confirmed target", "strong validation",
               "clinically actionable")
has_forbidden <- any(sapply(forbidden, function(term) grepl(term, m, ignore.case=TRUE)))

# Check scientific boundaries preserved
boundary_checks <- c(
  post_treatment  = any(grepl("post-treatment", m, ignore.case=TRUE)),
  tcga_not_immuno = any(grepl("TCGA is not an immunotherapy", m)),
  i2_90           = any(grepl("90%", m)),
  collectri_na    = any(grepl("CollecTRI.*not completed", m, ignore.case=TRUE)),
  mutation_ns     = any(grepl("not significant", m, ignore.case=TRUE)),
  cnv_no_results  = any(grepl("NO_RESULT_GENERATED|no association results", m, ignore.case=TRUE)),
  exploratory     = any(grepl("exploratory", m, ignore.case=TRUE))
)

checks <- data.frame(
  item = c("numbers_preserved", "line_count_unchanged", "sections_unchanged",
           "author_input_preserved", "no_new_statistics", "no_forbidden_terms",
           names(boundary_checks)),
  status = c(ifelse(nums_match, "PASS", "FAIL"),
             ifelse(lines_ok, "PASS", "FAIL"),
             ifelse(sections_ok, "PASS", "FAIL"),
             ifelse(has_author, "PASS", "FAIL"),
             ifelse(no_new_stats, "PASS", "FAIL"),
             ifelse(!has_forbidden, "PASS", "FAIL"),
             ifelse(boundary_checks, "PASS", "FAIL")),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(checks))) {
  cat(sprintf("  [%s] %s\n", checks$status[i], checks$item[i]))
}

n_pass <- sum(checks$status == "PASS")
n_fail <- sum(checks$status == "FAIL")
cat("\n  Total:", nrow(checks), "  PASS:", n_pass, "  FAIL:", n_fail, "\n\n")

# ============================================================================
# SECTION V: Completion Marker
# ============================================================================
cat("\nSECTION V: Completion Marker\n")
cat(rep("=", 80), sep=""); cat("\n\n")

marker <- c(
  "GSE243013 M6C LANGUAGE POLISHED: COMPLETE",
  "",
  paste("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Changes applied:",
  paste("  Total language changes:", nrow(changes)),
  "",
  "Verification:",
  paste("  Total checks:", nrow(checks)),
  paste("  PASS:", n_pass),
  paste("  FAIL:", n_fail),
  paste("  Overall:", if (n_fail == 0) "ALL PASS" else "SOME FAILURES"),
  "",
  "Output files:",
  "  main_text/GSE243013_full_manuscript_M6_polished.md",
  "  figures/GSE243013_all_figure_legends_M6_polished.md",
  "  audit/GSE243013_M6C_language_change_log.csv",
  "",
  "M6C completion marker: CREATED"
)
writeLines(marker, file.path("05_manuscript", "GSE243013_M6C_LANGUAGE_POLISHED.txt"))
cat("--- M6C COMPLETION MARKER CREATED ---\n")
for (mk in marker) cat("  ", mk, "\n")

cat("\n", rep("=", 80), sep="")
cat("\nM6C: English Language Polish - COMPLETED")
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")
