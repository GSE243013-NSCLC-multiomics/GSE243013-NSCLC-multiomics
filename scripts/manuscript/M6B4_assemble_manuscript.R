#!/usr/bin/env Rscript
# M6B4: Assemble Final Integrated Manuscript
# Pure concatenation of M6B1/M6B2/M6B3 outputs — no analysis, no re-polishing

cat("\n", rep("=", 80), sep="")
cat("\nM6B4: Assemble Final Integrated Manuscript")
cat("\nStarted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")

dirs <- list(
  b1  = "05_manuscript/M6_final_manuscript/work/M6B1_sections",
  b2  = "05_manuscript/M6_final_manuscript/work/M6B2_core_sections",
  b3  = "05_manuscript/M6_final_manuscript/work/M6B3_secondary_sections",
  out_main = "05_manuscript/M6_final_manuscript/main_text",
  out_fig  = "05_manuscript/M6_final_manuscript/figures",
  out_audit = "05_manuscript/M6_final_manuscript/audit"
)
invisible(lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE))

# ============================================================================
# SECTION I: Read All Section Files
# ============================================================================
cat("\nSECTION I: Read All Section Files\n")
cat(rep("=", 80), sep=""); cat("\n\n")

read_sec <- function(path, label) {
  if (!file.exists(path)) stop("MISSING: ", path)
  txt <- readLines(path, warn=FALSE)
  cat(sprintf("  %-40s %d lines\n", label, length(txt)))
  txt
}

sec <- list(
  title       = c("# Immune-Compartment Glycolysis Program Is Associated with",
                   "  Non-Response to Neoadjuvant Anti-PD-1 Therapy in NSCLC", "",
                   "**Authors:** [AUTHOR LIST]", "**Affiliations:** [AFFILIATIONS]",
                   "**Corresponding author:** [CORRESPONDING AUTHOR]", "",
                   "**Keywords:** NSCLC, neoadjuvant immunotherapy, glycolysis,",
                   "  single-cell RNA sequencing, TCGA", ""),
  abstract    = read_sec(file.path(dirs$b2, "Abstract_M6B2.md"),                "Abstract (M6B2)"),
  intro       = read_sec(file.path(dirs$b1, "02_Introduction_M5.md"),           "Introduction (M5/M6B1)"),
  methods     = read_sec(file.path(dirs$b2, "Methods_M6B2.md"),                 "Methods (M6B2)"),
  results     = read_sec(file.path(dirs$b2, "Results_M6B2.md"),                 "Results (M6B2)"),
  discussion  = read_sec(file.path(dirs$b3, "Discussion_M6B3.md"),              "Discussion (M6B3)"),
  conclusion  = read_sec(file.path(dirs$b3, "Conclusion_M6B3.md"),              "Conclusion (M6B3)"),
  clinical    = read_sec(file.path(dirs$b3, "Clinical_Relevance_M6B3.md"),      "Clinical Relevance (M6B3)"),
  transl      = read_sec(file.path(dirs$b3, "Translational_Relevance_M6B3.md"), "Translational Relevance (M6B3)"),
  fig5        = read_sec(file.path(dirs$b2, "Figure5_M6B2.md"),                 "Figure 5 (M6B2)"),
  other_figs  = read_sec(file.path(dirs$b3, "Other_Figure_Legends_M6B3.md"),    "Other Figure Legends (M6B3)")
)

# ============================================================================
# SECTION II: Concatenate
# ============================================================================
cat("\nSECTION II: Concatenate Sections\n")
cat(rep("=", 80), sep=""); cat("\n\n")

# Fixed order: Title, Abstract, Introduction, Methods, Results,
#              Discussion, Conclusion, Clinical Relevance,
#              Translational Relevance, Figure Legends
concat_order <- c("title", "abstract", "intro", "methods", "results",
                  "discussion", "conclusion", "clinical", "transl",
                  "fig5", "other_figs")

manuscript <- character()
for (key in concat_order) {
  manuscript <- c(manuscript, sec[[key]], "")
}

writeLines(manuscript, file.path(dirs$out_main, "GSE243013_manuscript_M6B_integrated.md"))
cat("  -> main_text/GSE243013_manuscript_M6B_integrated.md  (", length(manuscript), " lines)\n", sep="")

# Figure legends (Figure 5 + others)
fig_leg <- c(sec$fig5, "", sec$other_figs)
writeLines(fig_leg, file.path(dirs$out_fig, "GSE243013_figure_legends_M6B.md"))
cat("  -> figures/GSE243013_figure_legends_M6B.md  (", length(fig_leg), " lines)\n", sep="")

# ============================================================================
# SECTION III: Verification
# ============================================================================
cat("\nSECTION III: Verification\n")
cat(rep("=", 80), sep=""); cat("\n\n")

lines <- manuscript
txt   <- paste(lines, collapse="\n")

any_line <- function(pat, case=FALSE) any(sapply(lines, function(l) grepl(pat, l, ignore.case=case)))

ck <- function(item, label, pass) {
  data.frame(item=item, check=label, status=ifelse(pass, "PASS", "FAIL"),
             stringsAsFactors=FALSE)
}

checks <- rbind(
  ck("title",                    "Title section present",                  any_line("^# Immune-Compartment")),
  ck("section_abstract",         "Abstract appears exactly once",          sum(sapply(lines, function(l) grepl("^# Abstract$", l))) == 1),
  ck("section_intro",            "Introduction appears exactly once",      sum(sapply(lines, function(l) grepl("^# Introduction$", l))) == 1),
  ck("section_methods",          "Methods appears exactly once",           sum(sapply(lines, function(l) grepl("^# Methods$", l))) == 1),
  ck("section_results",          "Results appears exactly once",           sum(sapply(lines, function(l) grepl("^# Results$", l))) == 1),
  ck("section_discussion",       "Discussion appears exactly once",        sum(sapply(lines, function(l) grepl("^# Discussion$", l))) == 1),
  ck("section_conclusion",       "Conclusion appears exactly once",        sum(sapply(lines, function(l) grepl("^# Conclusion$", l))) == 1),
  ck("section_clinical",         "Clinical Relevance appears exactly once",sum(sapply(lines, function(l) grepl("^# Clinical Relevance$", l))) == 1),
  ck("section_transl",           "Translational Relevance appears once",   sum(sapply(lines, function(l) grepl("^# Translational Relevance$", l))) == 1),
  ck("section_fig5",             "Figure 5 heading present",              any_line("^# Figure 5")),
  ck("section_fig_legends",      "Figure Legends heading present",         any_line("^# Figure Legends")),
  # Content checks
  ck("no_empty_NES",             "No blank NES= values",                  !any_line("^NES=\\s*$|\\sNES=\\s*;")),
  ck("no_empty_FDR",             "No blank FDR= values",                  !any_line("^FDR=\\s*$|\\sFDR=\\s*;")),
  ck("no_FDR_zero",              "No FDR=0 literal",                      !any_line("FDR=0[^.]|FDR=0$|FDR=0 ")),
  ck("no_top_NA",                "No top=NA",                             !any_line("top=NA|top: NA")),
  ck("no_TBD",                   "No TBD",                                !any_line("\\bTBD\\b")),
  ck("no_TF_support",            "No TF support claims",                  !any_line("[Ss]upporting TF.*(show|demonstrate|indicate|suggest|reveal|were|is available)", case=TRUE)),
  ck("no_Rb_All_immune",         "Rb not attributed to All_immune",       !any_line("Rb.*All_immune|All_immune.*Rb")),
  ck("mutation_accurate",        "Mutation stated as not significant",    any_line("not significant|NO_SIGNIFICANT_FEATURE", case=TRUE)),
  ck("cnv_accurate",             "CNV stated as no results (not negative)",any_line("no association results|NO_RESULT_GENERATED|not generated", case=TRUE) & !any_line("CNV.*\\bnegative\\b|CNV.*statistically negative", case=TRUE)),
  ck("no_old_fig5",              "Uses M6B2 Figure 5 (not old)",          any_line("Integrated Evidence|M6B2|\\(A\\).*Preranked fgsea NES")),
  ck("pathway_preranked",        "GSE243013 uses preranked fgsea",        any_line("preranked fgsea|fgseaMultilevel") & !any_line("ssGSEA.*GSE243013|GSE243013.*ssGSEA", case=TRUE)),
  ck("no_internal_paths",        "No internal result paths in text",      !any_line("03_results/|01_scripts/|step0[0-9]|B1_QC2|B2/"))
)

n_pass <- sum(checks$status == "PASS")
n_fail <- sum(checks$status == "FAIL")

for (i in seq_len(nrow(checks))) {
  cat(sprintf("  [%s] %-40s %s\n", checks$status[i], checks$item[i], checks$check[i]))
}
cat("\n  Total:", nrow(checks), "  PASS:", n_pass, "  FAIL:", n_fail, "\n\n")

write.csv(checks, file.path(dirs$out_audit, "GSE243013_M6B_integration_check.csv"), row.names=FALSE)

# ============================================================================
# SECTION IV: Completion Marker
# ============================================================================
cat("\nSECTION IV: Completion Marker\n")
cat(rep("=", 80), sep=""); cat("\n\n")

marker <- c(
  "GSE243013 M6B TEXT INTEGRATED: COMPLETE",
  "",
  paste("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Assembly:",
  paste("  Total lines:", length(manuscript)),
  paste("  Sections: Title, Abstract, Introduction, Methods, Results,",
        "Discussion, Conclusion, Clinical Relevance, Translational Relevance, Figure Legends"),
  "",
  "Verification:",
  paste("  Total checks:", nrow(checks)),
  paste("  PASS:", n_pass),
  paste("  FAIL:", n_fail),
  paste("  Overall:", if (n_fail == 0) "ALL PASS" else "SOME FAILURES"),
  "",
  "Output files:",
  "  main_text/GSE243013_manuscript_M6B_integrated.md",
  "  figures/GSE243013_figure_legends_M6B.md",
  "  audit/GSE243013_M6B_integration_check.csv",
  "",
  "M6B4 completion marker: CREATED"
)
writeLines(marker, file.path("05_manuscript", "GSE243013_M6B_TEXT_INTEGRATED.txt"))
cat("--- M6B4 COMPLETION MARKER CREATED ---\n")
for (m in marker) cat("  ", m, "\n")

cat("\n", rep("=", 80), sep="")
cat("\nM6B4: Assemble Final Integrated Manuscript - COMPLETED")
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")
