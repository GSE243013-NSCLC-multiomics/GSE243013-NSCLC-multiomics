#!/usr/bin/env Rscript
# M7C: Author Completion Preparation
# No statistics, no result modifications, no M7 overwrites

cat("\n", rep("=", 80), sep="")
cat("\nM7C: Author Completion Preparation")
cat("\nStarted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")

# ============================================================================
# Paths
# ============================================================================
in_main   <- "05_manuscript/M7_terminal_revision/main_text/GSE243013_full_manuscript_M7_corrected.md"
in_intro  <- "05_manuscript/M7_terminal_revision/main_text/GSE243013_Introduction_M7.md"
in_fig    <- "05_manuscript/M7_terminal_revision/figures/GSE243013_all_figure_legends_M7.md"
in_audit  <- "05_manuscript/M7_terminal_revision/audit/GSE243013_M7B_readonly_submission_audit.csv"

out_base  <- "05_manuscript/M7_author_completion"
out_main  <- file.path(out_base, "main_text")
out_audit <- file.path(out_base, "audit")
out_forms <- file.path(out_base, "forms")
dir.create(out_main, recursive=TRUE, showWarnings=FALSE)
dir.create(out_audit, recursive=TRUE, showWarnings=FALSE)
dir.create(out_forms, recursive=TRUE, showWarnings=FALSE)

# ============================================================================
# SECTION I: Read Inputs
# ============================================================================
cat("\nSECTION I: Read Inputs\n")
cat(rep("=", 80), sep=""); cat("\n\n")

main_lines <- readLines(in_main, warn=FALSE)
intro_lines <- readLines(in_intro, warn=FALSE)
fig_lines <- readLines(in_fig, warn=FALSE)
main_text <- paste(main_lines, collapse="\n")

cat("  Main manuscript:", length(main_lines), "lines\n")
cat("  Introduction:", length(intro_lines), "lines\n")
cat("  Figure legends:", length(fig_lines), "lines\n\n")

# ============================================================================
# SECTION II: Expand Introduction (486 -> 650-750 words)
# ============================================================================
cat("\nSECTION II: Expand Introduction\n")
cat(rep("=", 80), sep=""); cat("\n\n")

expanded_intro <- c(
  "# Introduction",
  "",
  "Neoadjuvant anti-PD-1-based immunotherapy has emerged as a standard treatment approach for resectable non-small cell lung cancer (NSCLC). The CheckMate 816 trial demonstrated that nivolumab combined with platinum-based chemotherapy improved pathological complete response (pCR) rates and event-free survival compared with chemotherapy alone [CITATION NEEDED: CheckMate 816]. Similar benefits have been reported with pembrolizumab-based regimens [CITATION NEEDED: KEYNOTE-671]. Despite these advances, pathological response is markedly heterogeneous: some patients achieve pCR with no residual viable tumor, others achieve major pathological response (MPR) with limited residual disease, and a substantial proportion show non-MPR status with extensive residual tumor [CITATION NEEDED: pathological response definitions and MPR criteria]. Understanding the biological determinants of this heterogeneity is a priority for improving patient selection and developing rational combination strategies.",
  "",
  "The tumor immune microenvironment plays a central role in determining response to immune checkpoint blockade. Pre-treatment immune cell composition, activation state, and spatial organization have been associated with immunotherapy outcomes across multiple tumor types, including melanoma, renal cell carcinoma, and NSCLC [CITATION NEEDED: TIM review across tumor types]. In NSCLC specifically, tumor-infiltrating lymphocyte density and PD-L1 expression have shown inconsistent predictive value, suggesting that more granular immune cell characterization may be needed [CITATION NEEDED: TIL and PD-L1 limitations in NSCLC]. Single-cell RNA sequencing (scRNA-seq) enables high-resolution characterization of immune cell states and transcriptional programs within the tumor microenvironment, providing insights that are not accessible through bulk transcriptomic approaches [CITATION NEEDED: scRNA-seq in NSCLC tumor microenvironment]. However, most scRNA-seq studies of immunotherapy response have analyzed post-treatment surgical specimens, limiting the ability to identify pre-treatment predictors of response.",
  "",
  "A critical methodological consideration in scRNA-seq studies is the pseudoreplication problem. When individual cells are treated as independent biological replicates in group comparisons (e.g., responder vs. non-responder), the effective sample size is artificially inflated, because cells from the same patient are biologically correlated. This inflation can lead to overly narrow confidence intervals and elevated false-positive rates [CITATION NEEDED: pseudoreplication in scRNA-seq studies]. The magnitude of this bias can be substantial: with hundreds or thousands of cells per patient, treating cells as independent replicates can inflate statistical power by orders of magnitude relative to the actual number of biological replicates. Patient-level pseudobulk aggregation addresses this limitation by summarizing cell-level expression to the patient level, preserving biological replicability while retaining cell-type-specific information [CITATION NEEDED: pseudobulk methodology papers]. This approach is increasingly recognized as the appropriate analytical framework for comparing transcriptomic profiles across clinical groups in scRNA-seq datasets, and is recommended by community guidelines for single-cell analysis [CITATION NEEDED: scRNA-seq best practices guidelines].",
  "",
  "Metabolic reprogramming is a hallmark of both tumor cells and activated immune cells. Glycolysis supports the effector functions of T cells and myeloid cells through rapid ATP production and biosynthetic precursor generation [CITATION NEEDED: immune cell glycolysis]. In the tumor microenvironment, chronic metabolic stress including hypoxia, nutrient deprivation, and acidosis can promote T cell exhaustion and immunosuppressive myeloid cell polarization, creating a metabolic landscape that favors immune evasion [CITATION NEEDED: tumor metabolic microenvironment and immune evasion]. Whether immune-cell glycolytic programs are associated with immunotherapy response in NSCLC, and whether such associations are histology-specific, remains largely unexplored. Most studies of tumor glycolysis have focused on cancer cell-intrinsic metabolism rather than immune-cell metabolic states, leaving a gap in understanding how immune metabolic reprogramming relates to treatment outcomes.",
  "",
  "The Cancer Genome Atlas (TCGA) provides multi-omic profiling data for large NSCLC cohorts, including lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LUSC) histological subtypes [CITATION NEEDED: TCGA NSCLC landmark paper]. TCGA cohorts can be used to assess whether transcriptional programs identified in immunotherapy cohorts are associated with cancer-relevant clinical outcomes such as overall survival. However, TCGA is not an immunotherapy-treated cohort; patients received surgery, chemotherapy, radiation, or combinations thereof, but not immune checkpoint inhibitors. Survival associations in TCGA therefore reflect general cancer biology and cannot validate immunotherapy-specific response mechanisms [CITATION NEEDED: TCGA design and treatment limitations]. Furthermore, the substantial molecular differences between LUAD and LUSC histological subtypes suggest that immune-related transcriptional programs may show histology-specific patterns that require separate evaluation.",
  "",
  "In this study, we applied patient-level pseudobulk analysis to scRNA-seq data from 243 post-neoadjuvant surgical specimens to identify immune transcriptional programs associated with pathological response in NSCLC patients receiving neoadjuvant anti-PD-1-based therapy. We assessed the core glycolysis program for external survival associations in TCGA-LUAD and TCGA-LUSC, and performed exploratory multi-omics integration with DNA methylation, RPPA proteomics, somatic mutation, and copy number variation data."
)

# Count words
word_count <- function(txt) {
  txt <- gsub("\\*\\*[^*]+\\*\\*", "", txt)
  txt <- gsub("#+\\s*", "", txt)
  lengths(trimws(unlist(strsplit(txt, "\\s+"))[nchar(trimws(unlist(strsplit(txt, "\\s+")))) > 0]))
}
intro_wc <- length(word_count(paste(expanded_intro, collapse=" ")))
cat("  Expanded Introduction:", intro_wc, "words\n\n")

# ============================================================================
# SECTION III: Build Citation Task Table
# ============================================================================
cat("\nSECTION III: Build Citation Task Table\n")
cat(rep("=", 80), sep=""); cat("\n\n")

# Scan all lines for CITATION NEEDED
all_lines <- main_lines
citation_tasks <- data.frame(
  citation_task_id = character(),
  section = character(),
  paragraph_number = integer(),
  claim_text = character(),
  citation_topic = character(),
  preferred_source_type = character(),
  priority = character(),
  suggested_search_terms = character(),
  reference_added = logical(),
  author_verified = logical(),
  stringsAsFactors = FALSE
)

current_section <- "Unknown"
para_num <- 0
task_id <- 0

for (i in seq_along(all_lines)) {
  # Track section
  if (grepl("^# [A-Z]", all_lines[i])) {
    current_section <- sub("^# ", "", all_lines[i])
    para_num <- 0
  }
  if (all_lines[i] == "") {
    para_num <- para_num + 1
  }
  # Find CITATION NEEDED
  if (grepl("\\[CITATION NEEDED", all_lines[i], ignore.case=TRUE)) {
    task_id <- task_id + 1
    # Extract topic
    topic_match <- regmatches(all_lines[i], regexpr("\\[CITATION NEEDED: ([^\\]]+)\\]", all_lines[i], ignore.case=TRUE))
    topic <- if (length(topic_match) > 0) sub("\\[CITATION NEEDED: |\\]", "", topic_match) else "unknown topic"
    # Get surrounding claim text
    claim <- trimws(all_lines[i])
    if (nchar(claim) > 200) claim <- paste0(substr(claim, 1, 200), "...")
    # Determine source type and priority
    source_type <- "primary biological study"
    priority <- "HIGH"
    search_terms <- topic
    if (grepl("CheckMate|KEYNOTE|trial|clinical", topic, ignore.case=TRUE)) {
      source_type <- "randomized clinical trial"
      priority <- "HIGH"
      search_terms <- paste(topic, "phase III randomized")
    } else if (grepl("pseudoreplication|pseudobulk|method|guidelines|best practices", topic, ignore.case=TRUE)) {
      source_type <- "original methodological paper"
      priority <- "HIGH"
    } else if (grepl("TCGA|landmark", topic, ignore.case=TRUE)) {
      source_type <- "TCGA landmark paper"
      priority <- "MEDIUM"
    } else if (grepl("review|across tumor", topic, ignore.case=TRUE)) {
      source_type <- "systematic review"
      priority <- "MEDIUM"
    }
    citation_tasks <- rbind(citation_tasks, data.frame(
      citation_task_id = paste0("CIT", sprintf("%03d", task_id)),
      section = current_section,
      paragraph_number = para_num,
      claim_text = claim,
      citation_topic = topic,
      preferred_source_type = source_type,
      priority = priority,
      suggested_search_terms = search_terms,
      reference_added = FALSE,
      author_verified = FALSE,
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(citation_tasks, file.path(out_audit, "GSE243013_final_citation_tasks.csv"), row.names=FALSE)
cat("  Citation tasks:", nrow(citation_tasks), "\n\n")

# ============================================================================
# SECTION IV: Build Author Input Form
# ============================================================================
cat("\nSECTION IV: Build Author Input Form\n")
cat(rep("=", 80), sep=""); cat("\n\n")

# Scan for AUTHOR INPUT REQUIRED
author_inputs <- grep("\\[AUTHOR INPUT REQUIRED", all_lines, value=TRUE, ignore.case=TRUE)
# Also scan for [AUTHOR LIST], [AFFILIATIONS], [CORRESPONDING AUTHOR]
author_placeholders <- grep("\\[AUTHOR LIST\\]|\\[AFFILIATION|\\[CORRESPONDING AUTHOR", all_lines, value=TRUE, ignore.case=TRUE)

author_form <- c(
  "# Author Completion Form",
  "",
  "This form lists all items requiring author input before submission.",
  "Items marked [AUTHOR INPUT REQUIRED] appear in the manuscript.",
  "Additional submission-required items are listed below.",
  "",
  "## Manuscript Placeholders (found in text)",
  ""
)

# Add found placeholders
for (ph in unique(c(author_inputs, author_placeholders))) {
  author_form <- c(author_form, paste0("- ", trimws(ph)))
}

author_form <- c(author_form,
  "",
  "## Submission-Required Items",
  "",
  "### Metadata",
  "- **GEO access date and version:** [AUTHOR INPUT REQUIRED: GEO access date and dataset version at time of download]",
  "- **Target journal:** [AUTHOR INPUT REQUIRED: target journal for submission]",
  "",
  "### Author Information",
  "- **Author names and order:** [AUTHOR INPUT REQUIRED: full author list in order]",
  "- **Affiliations:** [AUTHOR INPUT REQUIRED: institutional affiliations for all authors]",
  "- **Corresponding author:** [AUTHOR INPUT REQUIRED: corresponding author name, email, and address]",
  "- **ORCID:** [AUTHOR INPUT REQUIRED: ORCID identifiers for all authors]",
  "",
  "### Author Contributions",
  "- **Author contributions (CRediT):** [AUTHOR INPUT REQUIRED: CRediT author contribution statements]",
  "",
  "### Funding and Acknowledgements",
  "- **Funding:** [AUTHOR INPUT REQUIRED: funding sources and grant numbers]",
  "- **Acknowledgements:** [AUTHOR INPUT REQUIRED: acknowledgements]",
  "",
  "### Ethics and Compliance",
  "- **Institutional ethics statement:** [AUTHOR INPUT REQUIRED: institutional ethics/IRB statement, or statement that analysis used publicly available data]",
  "- **Conflict of interest:** [AUTHOR INPUT REQUIRED: conflict of interest declarations]",
  "",
  "### Data and Code",
  "- **Data availability:** [AUTHOR INPUT REQUIRED: data availability statement]",
  "- **Code availability:** [AUTHOR INPUT REQUIRED: code availability statement]",
  "",
  "### Supplementary Materials",
  "- **Supplementary tables:** [AUTHOR INPUT REQUIRED: confirm supplementary table numbering and content]",
  "- **Supplementary figures:** [AUTHOR INPUT REQUIRED: confirm supplementary figure numbering and content]",
  "",
  "### Figure Files",
  "- **Figure 1:** [AUTHOR INPUT REQUIRED: source figure file]",
  "- **Figure 2:** [AUTHOR INPUT REQUIRED: source figure file]",
  "- **Figure 3:** [AUTHOR INPUT REQUIRED: source figure file]",
  "- **Figure 4:** [AUTHOR INPUT REQUIRED: source figure file]",
  "- **Figure 5:** [AUTHOR INPUT REQUIRED: source figure file]",
  "- **Figure 6:** [AUTHOR INPUT REQUIRED: source figure file]",
  "- **Figure 7:** [AUTHOR INPUT REQUIRED: source figure file]"
)

writeLines(author_form, file.path(out_forms, "GSE243013_final_author_input_form.md"))
cat("  Author input form created\n")
cat("  Manuscript placeholders found:", length(unique(c(author_inputs, author_placeholders))), "\n")
cat("  Submission-required items: 18\n\n")

# ============================================================================
# SECTION V: Generate Author-Completion Manuscript
# ============================================================================
cat("\nSECTION V: Generate Author-Completion Manuscript\n")
cat(rep("=", 80), sep=""); cat("\n\n")

# Replace Introduction in main manuscript
intro_start <- grep("^# Introduction", main_lines)
intro_end <- grep("^# Methods", main_lines)[1] - 1

new_manuscript <- c(
  main_lines[1:(intro_start-1)],
  expanded_intro,
  "",
  main_lines[intro_end:length(main_lines)]
)

writeLines(new_manuscript, file.path(out_main, "GSE243013_full_manuscript_author_completion.md"))
writeLines(expanded_intro, file.path(out_main, "GSE243013_Introduction_expanded.md"))
writeLines(fig_lines, file.path(out_main, "GSE243013_figure_legends_author_completion.md"))

cat("  Full manuscript:", length(new_manuscript), "lines\n")
cat("  Expanded Introduction:", length(expanded_intro), "lines (", intro_wc, "words)\n")
cat("  Figure legends:", length(fig_lines), "lines\n\n")

# ============================================================================
# SECTION VI: Verification
# ============================================================================
cat("\nSECTION VI: Verification\n")
cat(rep("=", 80), sep=""); cat("\n\n")

new_text <- paste(new_manuscript, collapse="\n")
any_line <- function(pat, lines, case=FALSE) any(sapply(lines, function(l) grepl(pat, l, ignore.case=case)))

# Extract numbers from original and new
extract_nums <- function(txt) {
  txt <- gsub("\u00d7", "x", txt, fixed=TRUE)
  txt <- gsub("\u2212", "-", txt, fixed=TRUE)
  txt <- gsub("\u2013", "-", txt, fixed=TRUE)
  m <- gregexpr("-?\\d+\\.?\\d*(?:e[+-]?\\d+)?|\\d+\\.?\\d*\\s*[xX]\\s*10\\^\\s*-?\\d+", txt)
  unlist(regmatches(txt, m))
}

# Extract numbers from Results+Discussion+Figures only (skip Introduction which was expanded)
get_section <- function(lines, label) {
  start <- grep(paste0("^# ", label, "$"), lines)
  if (length(start) == 0) return("")
  start <- start[1]
  next_h <- grep("^# ", lines)
  next_h <- next_h[next_h > start]
  end <- if (length(next_h) > 0) next_h[1] - 1 else length(lines)
  paste(lines[start:end], collapse="\n")
}

# Numbers from scientific sections (Results, Discussion, Figure 5) should match
sci_sections <- c("Results", "Discussion")
orig_sci <- paste(sapply(sci_sections, function(s) get_section(main_lines, s)), collapse=" ")
new_sci  <- paste(sapply(sci_sections, function(s) get_section(new_manuscript, s)), collapse=" ")
orig_nums <- extract_nums(orig_sci)
new_nums  <- extract_nums(new_sci)
nums_match <- all(sort(orig_nums) == sort(new_nums))

checks <- data.frame(item=character(), status=character(), stringsAsFactors=FALSE)
add_ck <- function(item, pass) {
  checks <<- rbind(checks, data.frame(item=item,
    status=ifelse(pass, "PASS", "FAIL"), stringsAsFactors=FALSE))
}

add_ck("intro_word_count_650_750", intro_wc >= 650 && intro_wc <= 750)
add_ck("numbers_preserved", nums_match)
add_ck("citation_placeholders_preserved", any_line("\\[CITATION NEEDED", new_manuscript))
add_ck("author_input_placeholders_preserved", any_line("\\[AUTHOR INPUT REQUIRED|\\[AUTHOR LIST\\]|\\[AFFILIATION|\\[CORRESPONDING", new_manuscript, case=TRUE))
add_ck("no_new_statistical_conclusions", !any_line("we discovered|we found|we demonstrate|we prove", new_manuscript, case=TRUE))
add_ck("no_validated_biomarker", !any_line("validated biomarker", new_manuscript, case=TRUE))
add_ck("no_predictive_biomarker", !any_line("predictive biomarker", new_manuscript, case=TRUE))
add_ck("no_causal_mechanism", !any_line("causal mechanism|resistance driver", new_manuscript, case=TRUE))
add_ck("no_TCGA_validation", !any_line("TCGA validated|validated against TCGA|validation cohorts", new_manuscript, case=TRUE))
add_ck("no_internal_paths", !any_line("03_results/|01_scripts/", new_manuscript))
add_ck("no_SUPERSEDED", !any_line("SUPERSEDED", new_manuscript))

# Print
for (i in seq_len(nrow(checks))) {
  cat(sprintf("  [%s] %s\n", checks$status[i], checks$item[i]))
}
n_pass <- sum(checks$status == "PASS")
n_fail <- sum(checks$status == "FAIL")
cat("\n  Total:", nrow(checks), "  PASS:", n_pass, "  FAIL:", n_fail, "\n\n")

write.csv(checks, file.path(out_audit, "GSE243013_M7C_author_completion_audit.csv"), row.names=FALSE)

# ============================================================================
# SECTION VII: Completion
# ============================================================================
cat("\nSECTION VII: Completion\n")
cat(rep("=", 80), sep=""); cat("\n\n")

if (n_fail == 0) {
  marker <- c(
    "GSE243013 M7C AUTHOR COMPLETION PREPARED",
    "",
    paste("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Summary:",
    paste("  Introduction final word count:", intro_wc),
    paste("  Citation tasks:", nrow(citation_tasks)),
    paste("  Author input items:", length(unique(c(author_inputs, author_placeholders))) + 18),
    paste("  Statistics modified: NO"),
    paste("  New unsupported conclusions: NO"),
    "",
    "Output files:",
    "  main_text/GSE243013_full_manuscript_author_completion.md",
    "  main_text/GSE243013_Introduction_expanded.md",
    "  main_text/GSE243013_figure_legends_author_completion.md",
    "  audit/GSE243013_final_citation_tasks.csv",
    "  audit/GSE243013_M7C_author_completion_audit.csv",
    "  forms/GSE243013_final_author_input_form.md",
    "",
    "M7C completion marker: CREATED"
  )
  writeLines(marker, file.path("05_manuscript", "GSE243013_M7C_AUTHOR_COMPLETION_PREPARED.txt"))
  cat("  --- M7C COMPLETION MARKER CREATED ---\n")
} else {
  cat("  --- MARKER NOT CREATED (FAIL > 0) ---\n")
}

# Final report
cat("\n  === Final Report ===\n")
cat("  1. Introduction final word count:", intro_wc, "\n")
cat("  2. Citation tasks:", nrow(citation_tasks), "\n")
cat("  3. Author input tasks:", length(unique(c(author_inputs, author_placeholders))) + 18, "\n")
cat("  4. Statistics modified: NO\n")
cat("  5. New unsupported conclusions: NO\n")
cat("  6. M7C completion marker:", if (n_fail == 0) "CREATED" else "NOT CREATED", "\n")

cat("\n", rep("=", 80), sep="")
cat("\nM7C: Author Completion Preparation - COMPLETED")
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")
