#!/usr/bin/env Rscript
# M6B1: Prepare and Verify Sections for Integration
# Read-only: reads M5 manuscript and M5A revised texts only
# No statistical analysis, no recursive scanning, no large files

cat("\n", rep("=", 80), sep="")
cat("\nM6B1: Prepare and Verify Sections")
cat("\nStarted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")

work_dir <- "05_manuscript/M6_final_manuscript/work/M6B1_sections"
dir.create(work_dir, recursive=TRUE, showWarnings=FALSE)

m5_file  <- "05_manuscript/M5_scientific_revision/main_text/GSE243013_manuscript_M5_clean.md"
m5a_abs  <- "05_manuscript/M5A_evidence_completion/revised_text/GSE243013_Abstract_M5A.md"
m5a_meth <- "05_manuscript/M5A_evidence_completion/revised_text/GSE243013_Methods_M5A.md"
m5a_res  <- "05_manuscript/M5A_evidence_completion/revised_text/GSE243013_Results_M5A.md"
m5a_fig5 <- "05_manuscript/M5A_evidence_completion/revised_text/GSE243013_Figure5_legend_M5A.md"
m6a_vals <- "05_manuscript/M6_final_manuscript/audit/GSE243013_M6A_verified_values.csv"


# ============================================================================
# SECTION I: Check Input Files
# ============================================================================
cat("\nSECTION I: Check Input Files\n")
cat(rep("=", 80), sep="")
cat("\n\n")

input_files <- data.frame(
  name = c("M5_clean", "M5A_Abstract", "M5A_Methods", "M5A_Results", "M5A_Figure5", "M6A_values"),
  path = c(m5_file, m5a_abs, m5a_meth, m5a_res, m5a_fig5, m6a_vals),
  exists = logical(6), nonempty = logical(6), utf8 = logical(6),
  stringsAsFactors = FALSE)

for (i in seq_len(nrow(input_files))) {
  p <- input_files$path[i]
  input_files$exists[i] <- file.exists(p)
  if (input_files$exists[i]) {
    txt <- tryCatch(readLines(p, warn=FALSE, encoding="UTF-8"), error=function(e) NULL)
    input_files$nonempty[i] <- !is.null(txt) && length(txt) > 0
    input_files$utf8[i] <- !is.null(txt)
  }
}

for (i in seq_len(nrow(input_files))) {
  icon <- if (input_files$exists[i] && input_files$nonempty[i] && input_files$utf8[i]) "[OK]" else "[FAIL]"
  cat(sprintf("  %s %-20s exists=%s  nonempty=%s  utf8=%s\n",
              icon, input_files$name[i],
              input_files$exists[i], input_files$nonempty[i], input_files$utf8[i]))
}

all_ok <- all(input_files$exists & input_files$nonempty & input_files$utf8)
cat("\nAll inputs OK:", all_ok, "\n\n")

if (!all_ok) {
  stop("Input file check failed. Stopping.")
}


# ============================================================================
# SECTION II: Split M5 Manuscript into Sections
# ============================================================================
cat("\nSECTION II: Split M5 Manuscript into Sections\n")
cat(rep("=", 80), sep="")
cat("\n\n")

m5_lines <- readLines(m5_file, warn=FALSE, encoding="UTF-8")
cat("M5 manuscript:", length(m5_lines), "lines\n")

# Find heading lines (lines starting with # or ##)
heading_pattern <- "^#{1,2}\\s+"
is_heading <- grepl(heading_pattern, m5_lines)
heading_lines <- which(is_heading)

cat("Headings found:", length(heading_lines), "\n")
for (hl in heading_lines) {
  cat(sprintf("  L%04d: %s\n", hl, substr(trimws(m5_lines[hl]), 1, 70)))
}

# Define section extraction function
extract_section <- function(lines, heading_line, next_heading_line) {
  if (is.na(next_heading_line)) {
    lines[heading_line:length(lines)]
  } else {
    lines[heading_line:(next_heading_line - 1)]
  }
}

# Find specific section boundaries
find_heading <- function(lines, pattern, heading_lines) {
  matches <- which(grepl(pattern, lines[heading_lines], ignore.case=TRUE))
  if (length(matches) > 0) as.integer(heading_lines[matches[1]]) else NA_integer_
}

# Map sections
sec_intro    <- find_heading(m5_lines, "^#+\\s+Introduction", heading_lines)
sec_methods  <- find_heading(m5_lines, "^#+\\s+Methods", heading_lines)
sec_results  <- find_heading(m5_lines, "^#+\\s+Results", heading_lines)
sec_discuss  <- find_heading(m5_lines, "^#+\\s+Discussion", heading_lines)
sec_conclude <- find_heading(m5_lines, "^#+\\s+Conclusion", heading_lines)
sec_clinical <- find_heading(m5_lines, "^#+\\s+Clinical\\s+Relevance", heading_lines)
sec_transl   <- find_heading(m5_lines, "^#+\\s+Translational\\s+Relevance", heading_lines)

# Figure legends: find "## Figures" or similar
sec_figures  <- find_heading(m5_lines, "^#+\\s+Figures|^#+\\s+Figure\\s+Legends|^#+\\s+Figure\\s+1", heading_lines)

cat("\nSection boundaries (line numbers):\n")
cat("  Introduction:", sec_intro, "\n")
cat("  Methods:", sec_methods, "\n")
cat("  Results:", sec_results, "\n")
cat("  Discussion:", sec_discuss, "\n")
cat("  Conclusion:", sec_conclude, "\n")
cat("  Clinical Relevance:", sec_clinical, "\n")
cat("  Translational Relevance:", sec_transl, "\n")
cat("  Figures:", sec_figures, "\n")


# ============================================================================
# SECTION III: Extract Each Section
# ============================================================================
cat("\nSECTION III: Extract Sections\n")
cat(rep("=", 80), sep="")
cat("\n\n")

# Build ordered list of all section headings
all_sections <- data.frame(
  heading = character(), line = integer(), stringsAsFactors = FALSE)

section_defs <- list(
  list(heading="Introduction",    line=sec_intro),
  list(heading="Methods",         line=sec_methods),
  list(heading="Results",         line=sec_results),
  list(heading="Discussion",      line=sec_discuss),
  list(heading="Conclusion",      line=sec_conclude),
  list(heading="Clinical Relevance", line=sec_clinical),
  list(heading="Translational Relevance", line=sec_transl),
  list(heading="Figures",         line=sec_figures)
)

for (sd in section_defs) {
  if (!is.na(sd$line) && sd$line > 0) {
    all_sections <- rbind(all_sections, data.frame(
      heading=sd$heading, line=sd$line, stringsAsFactors=FALSE))
  }
}
all_sections <- all_sections[order(all_sections$line), ]

cat("Ordered sections:\n")
for (i in seq_len(nrow(all_sections))) {
  cat(sprintf("  L%04d: %s\n", all_sections$line[i], all_sections$heading[i]))
}

# Extract Introduction
extract_and_save <- function(lines, start, end, filename, section_name) {
  if (is.na(start) || start == 0) {
    cat(sprintf("  WARNING: %s not found\n", section_name))
    return(data.frame(section_name=section_name, source_file="M5_clean.md",
                      start_line=NA_integer_, end_line=NA_integer_,
                      line_count=0L, character_count=0L,
                      status="NOT_FOUND", stringsAsFactors=FALSE))
  }
  end_line <- if (is.na(end)) length(lines) else end - 1
  section_lines <- lines[start:end_line]
  # Remove leading blank lines
  while (length(section_lines) > 0 && trimws(section_lines[1]) == "") {
    section_lines <- section_lines[-1]
  }
  writeLines(section_lines, file.path(work_dir, filename))
  cat(sprintf("  %s: lines %d-%d (%d lines, %d chars) -> %s\n",
              section_name, start, end_line, length(section_lines),
              sum(nchar(section_lines)), filename))
  return(data.frame(section_name=section_name, source_file="M5_clean.md",
                    start_line=start, end_line=end_line,
                    line_count=length(section_lines),
                    character_count=sum(nchar(section_lines)),
                    status="OK", stringsAsFactors=FALSE))
}

# Get next heading line after a given heading
get_next_heading <- function(current_line, section_lines) {
  future <- section_lines[section_lines > current_line]
  if (length(future) > 0) future[1] else NA_integer_
}

# Extract sections in order
manifest_rows <- list()

# Introduction
next_after_intro <- get_next_heading(sec_intro, all_sections$line)
manifest_rows[[length(manifest_rows)+1]] <- extract_and_save(
  m5_lines, sec_intro, next_after_intro, "02_Introduction_M5.md", "Introduction")

# Methods -> replaced by M5A
# Results -> replaced by M5A
# Discussion
next_after_discuss <- get_next_heading(sec_discuss, all_sections$line)
manifest_rows[[length(manifest_rows)+1]] <- extract_and_save(
  m5_lines, sec_discuss, next_after_discuss, "05_Discussion_M5.md", "Discussion")

# Conclusion
next_after_conclude <- get_next_heading(sec_conclude, all_sections$line)
manifest_rows[[length(manifest_rows)+1]] <- extract_and_save(
  m5_lines, sec_conclude, next_after_conclude, "06_Conclusion_M5.md", "Conclusion")

# Clinical Relevance
next_after_clinical <- get_next_heading(sec_clinical, all_sections$line)
manifest_rows[[length(manifest_rows)+1]] <- extract_and_save(
  m5_lines, sec_clinical, next_after_clinical, "07_Clinical_Relevance_M5.md", "Clinical Relevance")

# Translational Relevance
next_after_transl <- get_next_heading(sec_transl, all_sections$line)
manifest_rows[[length(manifest_rows)+1]] <- extract_and_save(
  m5_lines, sec_transl, next_after_transl, "08_Translational_Relevance_M5.md", "Translational Relevance")

# Figure legends section
next_after_figs <- get_next_heading(sec_figures, all_sections$line)
fig_section_lines <- extract_section(m5_lines, sec_figures, next_after_figs)

# Split figures: Figure 5 vs others
fig5_lines <- grep("^.*[Ff]igure\\s*5[:.\\s]|^\\*\\*Figure 5", fig_section_lines, value=TRUE)
other_fig_lines <- fig_section_lines[!grepl("^.*[Ff]igure\\s*5[:.\\s]|^\\*\\*Figure 5", fig_section_lines)]

# Find Figure 5 block within figures section
fig5_start_in_section <- grep("[Ff]igure\\s*5[:.\\s]|\\*\\*Figure 5", fig_section_lines)
if (length(fig5_start_in_section) > 0) {
  # Find end of Figure 5 block (next figure heading or end of section)
  fig5_block_start <- fig5_start_in_section[1]
  next_fig <- grep("^\\*\\*Figure [0-9]|^## Figure|^# Figure", fig_section_lines[-(1:fig5_block_start)])
  fig5_block_end <- if (length(next_fig) > 0) fig5_block_start + next_fig[1] - 1 else length(fig_section_lines)
  fig5_extracted <- fig_section_lines[fig5_block_start:fig5_block_end]
  other_figs_extracted <- fig_section_lines[c(1:(fig5_block_start-1),
                                               (fig5_block_end+1):length(fig_section_lines))]
  other_figs_extracted <- other_figs_extracted[other_figs_extracted != ""]
} else {
  fig5_extracted <- "# Figure 5 (M5A version will replace this)"
  other_figs_extracted <- fig_section_lines
}

# Save Figure 5 section placeholder
writeLines(fig5_extracted, file.path(work_dir, "09_Figure5_M5.md"))
cat(sprintf("  Figure 5 placeholder: %d lines\n", length(fig5_extracted)))

# Save other figure legends
writeLines(other_figs_extracted, file.path(work_dir, "10_Other_Figure_Legends_M5.md"))
cat(sprintf("  Other figure legends: %d lines\n", length(other_figs_extracted)))

manifest_rows[[length(manifest_rows)+1]] <- data.frame(
  section_name="Figure 5", source_file="M5_clean.md",
  start_line=sec_figures + fig5_block_start - 1,
  end_line=sec_figures + fig5_block_end - 1,
  line_count=length(fig5_extracted),
  character_count=sum(nchar(fig5_extracted)),
  status="PLACEHOLDER_M5A_WILL_REPLACE", stringsAsFactors=FALSE)

manifest_rows[[length(manifest_rows)+1]] <- data.frame(
  section_name="Other Figure Legends", source_file="M5_clean.md",
  start_line=sec_figures, end_line=sec_figures + length(fig_section_lines) - 1,
  line_count=length(other_figs_extracted),
  character_count=sum(nchar(other_figs_extracted)),
  status="OK", stringsAsFactors=FALSE)


# ============================================================================
# SECTION IV: Copy M5A Files as Standard Sections
# ============================================================================
cat("\nSECTION IV: Copy M5A Files\n")
cat(rep("=", 80), sep="")
cat("\n\n")

# Abstract M5A
abs_lines <- readLines(m5a_abs, warn=FALSE, encoding="UTF-8")
writeLines(abs_lines, file.path(work_dir, "01_Abstract_M5A.md"))
cat(sprintf("  Abstract M5A: %d lines -> 01_Abstract_M5A.md\n", length(abs_lines)))
manifest_rows[[length(manifest_rows)+1]] <- data.frame(
  section_name="Abstract", source_file="GSE243013_Abstract_M5A.md",
  start_line=1L, end_line=length(abs_lines),
  line_count=length(abs_lines), character_count=sum(nchar(abs_lines)),
  status="M5A_DIRECT_COPY", stringsAsFactors=FALSE)

# Methods M5A
meth_lines <- readLines(m5a_meth, warn=FALSE, encoding="UTF-8")
writeLines(meth_lines, file.path(work_dir, "03_Methods_M5A.md"))
cat(sprintf("  Methods M5A: %d lines -> 03_Methods_M5A.md\n", length(meth_lines)))
manifest_rows[[length(manifest_rows)+1]] <- data.frame(
  section_name="Methods", source_file="GSE243013_Methods_M5A.md",
  start_line=1L, end_line=length(meth_lines),
  line_count=length(meth_lines), character_count=sum(nchar(meth_lines)),
  status="M5A_DIRECT_COPY", stringsAsFactors=FALSE)

# Results M5A
res_lines <- readLines(m5a_res, warn=FALSE, encoding="UTF-8")
writeLines(res_lines, file.path(work_dir, "04_Results_M5A.md"))
cat(sprintf("  Results M5A: %d lines -> 04_Results_M5A.md\n", length(res_lines)))
manifest_rows[[length(manifest_rows)+1]] <- data.frame(
  section_name="Results", source_file="GSE243013_Results_M5A.md",
  start_line=1L, end_line=length(res_lines),
  line_count=length(res_lines), character_count=sum(nchar(res_lines)),
  status="M5A_DIRECT_COPY", stringsAsFactors=FALSE)

# Figure 5 M5A
fig5a_lines <- readLines(m5a_fig5, warn=FALSE, encoding="UTF-8")
writeLines(fig5a_lines, file.path(work_dir, "09_Figure5_M5A.md"))
cat(sprintf("  Figure 5 M5A: %d lines -> 09_Figure5_M5A.md\n", length(fig5a_lines)))
manifest_rows[[length(manifest_rows)+1]] <- data.frame(
  section_name="Figure 5 (M5A)", source_file="GSE243013_Figure5_legend_M5A.md",
  start_line=1L, end_line=length(fig5a_lines),
  line_count=length(fig5a_lines), character_count=sum(nchar(fig5a_lines)),
  status="M5A_DIRECT_COPY", stringsAsFactors=FALSE)


# ============================================================================
# SECTION V: Build and Save Manifest
# ============================================================================
cat("\nSECTION V: Build Manifest\n")
cat(rep("=", 80), sep="")
cat("\n\n")

manifest <- do.call(rbind, manifest_rows)

# Add section_id for ordering
manifest$section_id <- seq_len(nrow(manifest))

# Reorder columns
manifest <- manifest[, c("section_id", "section_name", "source_file",
                          "start_line", "end_line", "line_count",
                          "character_count", "status")]

write.csv(manifest,
  file.path(work_dir, "00_section_manifest.csv"),
  row.names=FALSE)

cat("Manifest saved:", nrow(manifest), "sections\n\n")

# Print summary
for (i in seq_len(nrow(manifest))) {
  cat(sprintf("  [%02d] %-30s  lines=%-5d  chars=%-6d  %s\n",
              manifest$section_id[i], manifest$section_name[i],
              manifest$line_count[i], manifest$character_count[i],
              manifest$status[i]))
}


# ============================================================================
# SECTION VI: Verify Outputs
# ============================================================================
cat("\nSECTION VI: Verify Outputs\n")
cat(rep("=", 80), sep="")
cat("\n\n")

expected_files <- c(
  "00_section_manifest.csv",
  "01_Abstract_M5A.md",
  "02_Introduction_M5.md",
  "03_Methods_M5A.md",
  "04_Results_M5A.md",
  "05_Discussion_M5.md",
  "06_Conclusion_M5.md",
  "07_Clinical_Relevance_M5.md",
  "08_Translational_Relevance_M5.md",
  "09_Figure5_M5A.md",
  "10_Other_Figure_Legends_M5.md"
)

verify <- data.frame(
  file = expected_files,
  exists = file.exists(file.path(work_dir, expected_files)),
  nonempty = sapply(expected_files, function(f) {
    fp <- file.path(work_dir, f)
    if (file.exists(fp)) file.info(fp)$size > 0 else FALSE
  }),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(verify))) {
  icon <- if (verify$exists[i] && verify$nonempty[i]) "[OK]" else "[FAIL]"
  cat(sprintf("  %s %s\n", icon, verify$file[i]))
}

n_ok <- sum(verify$exists & verify$nonempty)
n_total <- nrow(verify)
cat(sprintf("\nFiles: %d/%d OK\n", n_ok, n_total))

# Check no duplicates
dup_check <- any(duplicated(manifest$section_name))
cat("Duplicate sections:", dup_check, "\n")

# Check M5A sections present
m5a_present <- all(c("01_Abstract_M5A.md", "03_Methods_M5A.md",
                      "04_Results_M5A.md", "09_Figure5_M5A.md") %in%
                    verify$file[verify$exists])
cat("M5A Abstract/Methods/Results/Figure5:", m5a_present, "\n")

all_pass <- n_ok == n_total && !dup_check && m5a_present
cat("Overall:", if (all_pass) "ALL PASS" else "SOME FAILURES", "\n\n")


# ============================================================================
# SECTION VII: Completion Marker
# ============================================================================
cat("\nSECTION VII: Completion Marker\n")
cat(rep("=", 80), sep="")
cat("\n\n")

marker <- c(
  "GSE243013 M6B1 SECTIONS PREPARED: COMPLETE",
  "",
  paste("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Input files checked: 6/6",
  paste("M5 manuscript lines:", length(m5_lines)),
  paste("M5 headings found:", length(heading_lines)),
  "",
  "Section manifest:",
  paste("  Total sections:", nrow(manifest)),
  paste("  M5 extracted:", sum(manifest$status == "OK"), "(Introduction, Discussion, Conclusion, Clinical, Translational, Other Figs)"),
  paste("  M5A direct copy:", sum(manifest$status == "M5A_DIRECT_COPY"), "(Abstract, Methods, Results, Figure 5 M5A)"),
  paste("  M5 placeholder:", sum(grepl("PLACEHOLDER", manifest$status)), "(Figure 5 M5 placeholder)"),
  "",
  "Output files:",
  paste("  Directory:", work_dir),
  paste("  Files:", n_ok, "/", n_total, "OK"),
  paste("  Duplicates:", if (dup_check) "YES (problem)" else "NONE"),
  paste("  M5A sections:", if (m5a_present) "ALL PRESENT" else "MISSING"),
  "",
  "8. M6B1 completion marker: CREATED"
)
writeLines(marker, "05_manuscript/GSE243013_M6B1_SECTIONS_PREPARED.txt")
cat("--- M6B1 COMPLETION MARKER CREATED ---\n")
for (line in marker) cat("  ", line, "\n")

cat("\n", rep("=", 80), sep="")
cat("\nM6B1: Prepare and Verify Sections - COMPLETED")
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")
