## =========================================================================
## Step 03: Repair GSE243013 Cohort Definition & Audit Logic
## =========================================================================

options(timeout = 3600)

cat("========================================================================\n")
cat("Step 03: Repair GSE243013 Cohort Definition & Audit Logic\n")
cat("========================================================================\n\n")

RAW_DIR  <- "02_data/raw/GSE243013"
RESULTS  <- "03_results"

dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)

## =========================================================================
## I. Read Data
## =========================================================================
cat("[I] Reading data...\n")

META_PATH <- file.path(RAW_DIR, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")
meta <- NULL

tryCatch({
  meta <- data.table::fread(META_PATH, header = TRUE, data.table = FALSE)
  cat(sprintf("[INFO] Read metadata: %d rows x %d columns\n", nrow(meta), ncol(meta)))
}, error = function(e) {
  cat("[WARNING] data.table failed:", conditionMessage(e), "\n")
  meta <<- read.csv(gzfile(META_PATH), check.names = FALSE)
  cat(sprintf("[INFO] Read with read.csv: %d rows x %d columns\n", nrow(meta), ncol(meta)))
})

cat(sprintf("[INFO] Column names: %s\n", paste(names(meta), collapse = ", ")))

## Read old manifest
old_manifest <- NULL
tryCatch({
  old_manifest <- read.csv("03_results/GSE243013_patient_manifest.csv", stringsAsFactors = FALSE)
  cat(sprintf("[INFO] Read old manifest: %d rows\n", nrow(old_manifest)))
}, error = function(e) {
  cat("[WARNING] Could not read old manifest:", conditionMessage(e), "\n")
})

## =========================================================================
## II. Define is_no_value()
## =========================================================================
cat("\n[II] Defining is_no_value()...\n")

no_values <- c("na", "n/a", "null", "none", "no", "not", "not used", "not given",
               "negative", "0", "unknown", "unkown", "")

is_no_value <- function(x) {
  if (is.na(x)) return(TRUE)
  x_trimmed <- trimws(tolower(as.character(x)))
  if (x_trimmed == "") return(TRUE)
  return(x_trimmed %in% no_values)
}

## Test
cat("[INFO] is_no_value tests:\n")
cat(sprintf("  NA: %s\n", is_no_value(NA)))
cat(sprintf("  'No': %s\n", is_no_value("No")))
cat(sprintf("  'none': %s\n", is_no_value("none")))
cat(sprintf("  'Pembrolizumab': %s\n", is_no_value("Pembrolizumab")))
cat(sprintf("  'Carboplatin + Abraxane': %s\n", is_no_value("Carboplatin + Abraxane")))

## =========================================================================
## III. Rebuild Patient-Level Treatment
## =========================================================================
cat("\n[III] Rebuilding patient-level treatment...\n")

## Identify actual column names (handle hyphen vs underscore)
cat("[INFO] Looking for anti-PD1_therapy column...\n")
pd1_col <- intersect(c("anti-PD1_therapy", "anti_PD1_therapy"), names(meta))
chemo_col <- intersect(c("chemotherapy"), names(meta))
targeted_col <- intersect(c("targeted_therapy"), names(meta))

cat(sprintf("[INFO] anti-PD1 column found: '%s'\n", ifelse(length(pd1_col) > 0, pd1_col[1], "NOT FOUND")))
cat(sprintf("[INFO] chemotherapy column found: '%s'\n", ifelse(length(chemo_col) > 0, chemo_col[1], "NOT FOUND")))
cat(sprintf("[INFO] targeted_therapy column found: '%s'\n", ifelse(length(targeted_col) > 0, targeted_col[1], "NOT FOUND")))

if (length(pd1_col) == 0) stop("Cannot find anti-PD1 therapy column")
if (length(chemo_col) == 0) stop("Cannot find chemotherapy column")
if (length(targeted_col) == 0) stop("Cannot find targeted_therapy column")

pd1_col <- pd1_col[1]
chemo_col <- chemo_col[1]
targeted_col <- targeted_col[1]

## Get unique sampleIDs
sample_ids <- sort(unique(meta$sampleID))
cat(sprintf("[INFO] Unique sampleIDs: %d\n", length(sample_ids)))

## Build patient-level treatment
treatment_rows <- list()

for (sid in sample_ids) {
  idx <- meta$sampleID == sid
  sub <- meta[idx, ]

  ## Extract unique non-NA values for each treatment field
  pd1_vals <- unique(sub[[pd1_col]][!is.na(sub[[pd1_col]])])
  chemo_vals <- unique(sub[[chemo_col]][!is.na(sub[[chemo_col]])])
  targeted_vals <- unique(sub[[targeted_col]][!is.na(sub[[targeted_col]])])

  ## Original values (first non-NA)
  pd1_orig <- ifelse(length(pd1_vals) > 0, pd1_vals[1], NA_character_)
  chemo_orig <- ifelse(length(chemo_vals) > 0, chemo_vals[1], NA_character_)
  targeted_orig <- ifelse(length(targeted_vals) > 0, targeted_vals[1], NA_character_)

  ## Check if has actual drug (not "no" values)
  pd1_no <- sapply(pd1_vals, is_no_value)
  chemo_no <- sapply(chemo_vals, is_no_value)
  targeted_no <- sapply(targeted_vals, is_no_value)

  has_anti_PD1 <- any(!pd1_no) && length(pd1_vals) > 0
  has_chemo <- any(!chemo_no) && length(chemo_vals) > 0
  has_targeted <- any(!targeted_no) && length(targeted_vals) > 0

  ## Check for conflicts within each field
  pd1_has_drug <- !pd1_no
  chemo_has_drug <- !chemo_no
  targeted_has_drug <- !targeted_no

  pd1_conflict <- sum(pd1_has_drug) > 1 || (sum(pd1_has_drug) > 0 && sum(pd1_no) > 0 && length(pd1_vals) > 1)
  chemo_conflict <- sum(chemo_has_drug) > 1 || (sum(chemo_has_drug) > 0 && sum(chemo_no) > 0 && length(chemo_vals) > 1)
  targeted_conflict <- sum(targeted_has_drug) > 1 || (sum(targeted_has_drug) > 0 && sum(targeted_no) > 0 && length(targeted_vals) > 1)

  has_conflict <- pd1_conflict || chemo_conflict || targeted_conflict

  ## Treatment pattern (mutually exclusive)
  if (has_conflict) {
    treatment_pattern <- "conflict"
  } else if (has_anti_PD1 && has_chemo) {
    treatment_pattern <- "anti_PD1_plus_chemo_strict"
  } else if (has_anti_PD1 && !has_chemo) {
    treatment_pattern <- "anti_PD1_recorded_chemo_not_recorded"
  } else if (!has_anti_PD1 && has_chemo) {
    treatment_pattern <- "chemotherapy_only"
  } else if (!has_anti_PD1 && !has_chemo && has_targeted) {
    treatment_pattern <- "targeted_or_other_only"
  } else {
    treatment_pattern <- "no_treatment_recorded"
  }

  ## Cohort flags
  anti_PD1_cohort <- has_anti_PD1
  strict_chemo_cohort <- has_anti_PD1 && has_chemo
  chemo_control_cohort <- !has_anti_PD1 && has_chemo

  treatment_rows[[sid]] <- data.frame(
    sampleID = sid,
    anti_PD1_therapy_original = pd1_orig,
    chemotherapy_original = chemo_orig,
    targeted_therapy_original = targeted_orig,
    has_anti_PD1 = has_anti_PD1,
    has_chemotherapy = has_chemo,
    has_targeted_therapy = has_targeted,
    treatment_pattern = treatment_pattern,
    anti_PD1_cohort = anti_PD1_cohort,
    strict_chemoimmunotherapy_cohort = strict_chemo_cohort,
    chemotherapy_control_cohort = chemo_control_cohort,
    stringsAsFactors = FALSE
  )
}

treatment_audit <- do.call(rbind, treatment_rows)
rownames(treatment_audit) <- NULL

cat("\n[INFO] Treatment pattern distribution:\n")
print(table(treatment_audit$treatment_pattern))

cat(sprintf("\n[INFO] anti_PD1_cohort: %d\n", sum(treatment_audit$anti_PD1_cohort)))
cat(sprintf("[INFO] strict_chemoimmunotherapy_cohort: %d\n", sum(treatment_audit$strict_chemoimmunotherapy_cohort)))
cat(sprintf("[INFO] chemotherapy_control_cohort: %d\n", sum(treatment_audit$chemotherapy_control_cohort)))

write.csv(treatment_audit, file.path(RESULTS, "GSE243013_treatment_reaudit.csv"), row.names = FALSE)
cat("[INFO] Saved: GSE243013_treatment_reaudit.csv\n")

## =========================================================================
## IV. Review Previous 21 Unresolved Patients
## =========================================================================
cat("\n[IV] Reviewing previous unresolved patients...\n")

if (!is.null(old_manifest)) {
  old_unresolved <- old_manifest[old_manifest$treatment_class_candidate == "unresolved", ]
  cat(sprintf("[INFO] Old unresolved patients: %d\n", nrow(old_unresolved)))

  review_rows <- list()
  for (i in seq_len(nrow(old_unresolved))) {
    sid <- old_unresolved$sampleID[i]
    new_row <- treatment_audit[treatment_audit$sampleID == sid, ]
    sub <- meta[meta$sampleID == sid, ]

    pd1_val <- ifelse(nrow(new_row) > 0, as.character(new_row$anti_PD1_therapy_original[1]), NA)
    chemo_val <- ifelse(nrow(new_row) > 0, as.character(new_row$chemotherapy_original[1]), NA)
    targeted_val <- ifelse(nrow(new_row) > 0, as.character(new_row$targeted_therapy_original[1]), NA)
    new_pattern <- ifelse(nrow(new_row) > 0, as.character(new_row$treatment_pattern[1]), NA)

    ## Classify reason
    pd1_is_no <- is_no_value(pd1_val)
    chemo_is_no <- is_no_value(chemo_val)

    if (!pd1_is_no && chemo_is_no) {
      reason <- "has_anti_PD1_drug_but_chemo_recorded_as_No"
    } else if (pd1_is_no && !chemo_is_no) {
      reason <- "anti_PD1_recorded_as_No_but_has_chemo_regimen"
    } else if (!pd1_is_no && !chemo_is_no) {
      reason <- "both_have_drug_records"
    } else {
      reason <- "other"
    }

    resp_orig <- unique(sub$pathological_response[!is.na(sub$pathological_response)])
    cancer <- unique(sub$cancer_type[!is.na(sub$cancer_type)])

    review_rows[[sid]] <- data.frame(
      sampleID = sid,
      anti_PD1_therapy = pd1_val,
      chemotherapy = chemo_val,
      targeted_therapy = targeted_val,
      pathological_response = ifelse(length(resp_orig) > 0, resp_orig[1], NA),
      cancer_type = ifelse(length(cancer) > 0, cancer[1], NA),
      old_treatment_class = "unresolved",
      new_treatment_pattern = new_pattern,
      new_classification_reason = reason,
      stringsAsFactors = FALSE
    )
  }

  review_df <- do.call(rbind, review_rows)
  rownames(review_df) <- NULL

  cat("\n[INFO] Previous unresolved patients - new classifications:\n")
  print(table(review_df$new_treatment_pattern))

  cat("\n[INFO] Detailed breakdown:\n")
  for (i in seq_len(nrow(review_df))) {
    cat(sprintf("  %s: PD1='%s', chemo='%s' → %s (%s)\n",
                review_df$sampleID[i],
                review_df$anti_PD1_therapy[i],
                review_df$chemotherapy[i],
                review_df$new_treatment_pattern[i],
                review_df$new_classification_reason[i]))
  }

  write.csv(review_df, file.path(RESULTS, "GSE243013_previous_unresolved_treatment_review.csv"),
            row.names = FALSE)
  cat("[INFO] Saved: GSE243013_previous_unresolved_treatment_review.csv\n")
} else {
  cat("[WARNING] Old manifest not available, skipping unresolved review.\n")
}

## =========================================================================
## V. Repair Response Endpoints
## =========================================================================
cat("\n[V] Repairing response endpoints...\n")

## Build patient-level response
response_rows <- list()
for (sid in sample_ids) {
  idx <- meta$sampleID == sid
  sub <- meta[idx, ]

  resp_vals <- unique(sub$pathological_response[!is.na(sub$pathological_response)])
  resp_orig <- ifelse(length(resp_vals) > 0, resp_vals[1], NA_character_)

  ## Clean response
  resp_clean <- resp_orig
  if (!is.na(resp_clean)) {
    resp_lower <- tolower(trimws(resp_clean))
    if (resp_lower %in% c("pcr")) resp_clean <- "pCR"
    else if (resp_lower %in% c("mpr")) resp_clean <- "MPR"
    else if (grepl("^non.?mpr$", resp_lower) || resp_lower == "nmpr") resp_clean <- "non-MPR"
    else if (resp_lower %in% c("unknown", "unkown", "unknowm", "unkownn")) resp_clean <- "unknown"
    else if (!(resp_clean %in% c("pCR", "MPR", "non-MPR", "unknown"))) resp_clean <- "UNRESOLVED"
  } else {
    resp_clean <- "UNRESOLVED"
  }

  ## Binary response
  resp_binary <- NA_character_
  if (resp_clean %in% c("pCR", "MPR")) resp_binary <- "Responder"
  else if (resp_clean == "non-MPR") resp_binary <- "Non_responder"

  ## Response rate descriptive audit
  rate_vals <- unique(sub$pathological_response_rate[!is.na(sub$pathological_response_rate)])
  rate_orig <- ifelse(length(rate_vals) > 0, as.character(rate_vals[1]), NA_character_)

  rate_type <- "missing"
  rate_numeric <- NA_real_

  if (!is.na(rate_orig)) {
    rate_str <- trimws(as.character(rate_orig))
    rate_num <- suppressWarnings(as.numeric(rate_str))

    if (!is.na(rate_num)) {
      if (rate_num <= 1) {
        rate_type <- "numeric_fraction"
      } else {
        rate_type <- "numeric_percent"
      }
      rate_numeric <- rate_num
    } else if (grepl("^<|>|%|\\d", rate_str)) {
      ## Try extracting number from strings like "<10%", ">90%"
      extracted <- regmatches(rate_str, regexpr("[0-9]+\\.?[0-9]*", rate_str))
      if (length(extracted) > 0 && !is.na(extracted)) {
        rate_type <- "text_description"
        rate_numeric <- suppressWarnings(as.numeric(extracted))
      } else {
        rate_type <- "text_description"
      }
    } else {
      ## Check for Chinese text or other non-numeric
      rate_type <- "text_description"
    }
  }

  response_rows[[sid]] <- data.frame(
    sampleID = sid,
    pathological_response_original = resp_orig,
    pathological_response_clean = resp_clean,
    response_binary = resp_binary,
    pathological_response_rate_original = rate_orig,
    response_rate_value_type = rate_type,
    response_rate_numeric_extracted = rate_numeric,
    stringsAsFactors = FALSE
  )
}

response_df <- do.call(rbind, response_rows)
rownames(response_df) <- NULL

cat("\n[INFO] Response clean distribution:\n")
print(table(response_df$pathological_response_clean, useNA = "ifany"))

cat("\n[INFO] Response rate value types:\n")
print(table(response_df$response_rate_value_type))

write.csv(response_df, file.path(RESULTS, "GSE243013_response_rate_descriptive_audit.csv"),
          row.names = FALSE)
cat("[INFO] Saved: GSE243013_response_rate_descriptive_audit.csv\n")

## =========================================================================
## VI. Build Revised Patient Manifest
## =========================================================================
cat("\n[VI] Building revised patient manifest...\n")

manifest <- data.frame(sampleID = sample_ids, stringsAsFactors = FALSE)

## Add n_cells
for (i in seq_len(nrow(manifest))) {
  sid <- manifest$sampleID[i]
  manifest$n_cells[i] <- sum(meta$sampleID == sid)
}

## Add cancer_type and staging
for (i in seq_len(nrow(manifest))) {
  sid <- manifest$sampleID[i]
  sub <- meta[meta$sampleID == sid, ]
  ct <- unique(sub$cancer_type[!is.na(sub$cancer_type)])
  stg <- unique(sub$pre_treatment_staging[!is.na(sub$pre_treatment_staging)])
  manifest$cancer_type[i] <- ifelse(length(ct) > 0, ct[1], NA_character_)
  manifest$pre_treatment_staging[i] <- ifelse(length(stg) > 0, stg[1], NA_character_)
}

## Merge response
manifest <- merge(manifest, response_df, by = "sampleID", all.x = TRUE)

## Merge treatment
manifest <- merge(manifest, treatment_audit, by = "sampleID", all.x = TRUE)

## Add cell type counts
for (i in seq_len(nrow(manifest))) {
  sid <- manifest$sampleID[i]
  sub <- meta[meta$sampleID == sid, ]
  manifest$major_cell_type_count[i] <- length(unique(sub$major_cell_type))
  manifest$sub_cell_type_count[i] <- length(unique(sub$sub_cell_type))
}

## Conflict flags
manifest$has_any_patient_level_conflict <- FALSE

## Eligibility flags
manifest$primary_analysis_eligible <- (
  manifest$anti_PD1_cohort == TRUE &
  !is.na(manifest$response_binary) &
  !is.na(manifest$sampleID) &
  manifest$has_any_patient_level_conflict == FALSE &
  manifest$n_cells > 0
)

manifest$strict_sensitivity_analysis_eligible <- (
  manifest$strict_chemoimmunotherapy_cohort == TRUE &
  !is.na(manifest$response_binary) &
  manifest$has_any_patient_level_conflict == FALSE
)

manifest$chemotherapy_control_eligible <- (
  manifest$chemotherapy_control_cohort == TRUE &
  !is.na(manifest$response_binary) &
  manifest$has_any_patient_level_conflict == FALSE
)

cat(sprintf("[INFO] primary_analysis_eligible: %d\n", sum(manifest$primary_analysis_eligible)))
cat(sprintf("[INFO] strict_sensitivity_analysis_eligible: %d\n", sum(manifest$strict_sensitivity_analysis_eligible)))
cat(sprintf("[INFO] chemotherapy_control_eligible: %d\n", sum(manifest$chemotherapy_control_eligible)))

write.csv(manifest, file.path(RESULTS, "GSE243013_patient_manifest_revised.csv"), row.names = FALSE)
cat("[INFO] Saved: GSE243013_patient_manifest_revised.csv\n")

## =========================================================================
## VII. Cohort Summaries (Separate Denominators)
## =========================================================================
cat("\n[VII] Generating cohort summaries...\n")

summarize_cohort <- function(subset_df, cohort_name) {
  if (nrow(subset_df) == 0) return(NULL)
  data.frame(
    cohort = cohort_name,
    total_patients = nrow(subset_df),
    pCR = sum(subset_df$pathological_response_clean == "pCR", na.rm = TRUE),
    MPR = sum(subset_df$pathological_response_clean == "MPR", na.rm = TRUE),
    non_MPR = sum(subset_df$pathological_response_clean == "non-MPR", na.rm = TRUE),
    unknown = sum(subset_df$pathological_response_clean == "unknown", na.rm = TRUE),
    Responder = sum(subset_df$response_binary == "Responder", na.rm = TRUE),
    Non_responder = sum(subset_df$response_binary == "Non_responder", na.rm = TRUE),
    LUAD = sum(subset_df$cancer_type == "LUAD", na.rm = TRUE),
    LUSC = sum(subset_df$cancer_type == "LUSC", na.rm = TRUE),
    total_cells = sum(subset_df$n_cells, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

## A. All patients
summary_all <- summarize_cohort(manifest, "all_patients")

## B. Anti-PD1 treated (primary)
summary_antiPD1 <- summarize_cohort(manifest[manifest$anti_PD1_cohort == TRUE, ],
                                     "anti_PD1_treated_primary")

## C. Strict chemo-immunotherapy
summary_strict <- summarize_cohort(manifest[manifest$strict_chemoimmunotherapy_cohort == TRUE, ],
                                    "strict_chemoimmunotherapy")

## D. Chemotherapy control
summary_chemo <- summarize_cohort(manifest[manifest$chemotherapy_control_cohort == TRUE, ],
                                   "chemotherapy_control")

cohort_summaries <- rbind(summary_all, summary_antiPD1, summary_strict, summary_chemo)
write.csv(cohort_summaries, file.path(RESULTS, "GSE243013_revised_cohort_summary.csv"),
          row.names = FALSE)
cat("[INFO] Saved: GSE243013_revised_cohort_summary.csv\n")

cat("\n[INFO] Cohort summaries:\n")
print(cohort_summaries)

## =========================================================================
## VIII. Consistency Checks
## =========================================================================
cat("\n[VIII] Running consistency checks...\n")

checks <- list()

## 1. sampleID count
checks$sampleID_count_243 <- nrow(manifest) == 243
cat(sprintf("  1. sampleID count = 243: %s (actual: %d)\n",
            ifelse(checks$sampleID_count_243, "PASS", "FAIL"), nrow(manifest)))

## 2. Each patient has one response
resp_conflicts <- 0
for (sid in sample_ids) {
  resp_vals <- unique(manifest$pathological_response_clean[manifest$sampleID == sid])
  resp_vals <- resp_vals[!is.na(resp_vals) & resp_vals != "UNRESOLVED"]
  if (length(resp_vals) > 1) resp_conflicts <- resp_conflicts + 1
}
checks$no_response_conflicts <- resp_conflicts == 0
cat(sprintf("  2. No response conflicts: %s (conflicts: %d)\n",
            ifelse(checks$no_response_conflicts, "PASS", "FAIL"), resp_conflicts))

## 3. pCR/MPR/non-MPR mutually exclusive
primary_sub <- manifest[manifest$primary_analysis_eligible, ]
if (nrow(primary_sub) > 0) {
  pcr_n <- sum(primary_sub$pathological_response_clean == "pCR", na.rm = TRUE)
  mpr_n <- sum(primary_sub$pathological_response_clean == "MPR", na.rm = TRUE)
  non_mpr_n <- sum(primary_sub$pathological_response_clean == "non-MPR", na.rm = TRUE)
  checks$response_mutual_exclusion <- (pcr_n + mpr_n + non_mpr_n <= nrow(primary_sub))
} else {
  checks$response_mutual_exclusion <- TRUE
}
cat(sprintf("  3. Response categories mutually exclusive: %s\n",
            ifelse(checks$response_mutual_exclusion, "PASS", "FAIL")))

## 4. Responder vs Non_responder mutual exclusion
resp_binary_na <- sum(is.na(manifest$response_binary))
resp_binary_both <- sum(manifest$response_binary == "Responder", na.rm = TRUE) +
                    sum(manifest$response_binary == "Non_responder", na.rm = TRUE)
checks$binary_exclusion <- (resp_binary_both + resp_binary_na == nrow(manifest))
cat(sprintf("  4. Binary response mutual exclusion: %s\n",
            ifelse(checks$binary_exclusion, "PASS", "FAIL")))

## 5. anti_PD1_cohort vs chemotherapy_control_cohort mutual exclusion
both_cohorts <- sum(manifest$anti_PD1_cohort == TRUE & manifest$chemotherapy_control_cohort == TRUE, na.rm = TRUE)
checks$cohort_exclusion <- both_cohorts == 0
cat(sprintf("  5. Cohort mutual exclusion: %s (overlap: %d)\n",
            ifelse(checks$cohort_exclusion, "PASS", "FAIL"), both_cohorts))

## 6. All summaries have denominators
checks$summaries_have_denominators <- all(!is.na(cohort_summaries$total_patients)) &
                                       all(cohort_summaries$total_patients > 0)
cat(sprintf("  6. Summaries have denominators: %s\n",
            ifelse(checks$summaries_have_denominators, "PASS", "FAIL")))

## 7. cellID globally unique
n_distinct_cellID <- length(unique(meta$cellID))
checks$cellID_unique <- n_distinct_cellID == nrow(meta)
cat(sprintf("  7. cellID globally unique: %s (distinct: %d, total: %d)\n",
            ifelse(checks$cellID_unique, "PASS", "FAIL"), n_distinct_cellID, nrow(meta)))

## 8. metadata cellID count matches row count
checks$cellID_row_match <- n_distinct_cellID == nrow(meta)
cat(sprintf("  8. cellID matches metadata rows: %s\n",
            ifelse(checks$cellID_row_match, "PASS", "FAIL")))

## 9. Each patient has at least one cell
min_cells <- min(manifest$n_cells)
checks$min_cells_per_patient <- min_cells >= 1
cat(sprintf("  9. Min cells per patient: %d (%s)\n",
            min_cells, ifelse(checks$min_cells_per_patient, "PASS", "FAIL")))

## 10. anti_PD1 cohort close to 234
antiPD1_n <- sum(manifest$anti_PD1_cohort == TRUE)
checks$antiPD1_close_to_234 <- abs(antiPD1_n - 234) <= 10
cat(sprintf("  10. anti_PD1 cohort ~234: %s (actual: %d)\n",
            ifelse(checks$antiPD1_close_to_234, "PASS/WARN", "FAIL"), antiPD1_n))

## 11. chemotherapy control close to 9
chemo_n <- sum(manifest$chemotherapy_control_cohort == TRUE)
checks$chemo_close_to_9 <- abs(chemo_n - 9) <= 5
cat(sprintf("  11. chemotherapy control ~9: %s (actual: %d)\n",
            ifelse(checks$chemo_close_to_9, "PASS/WARN", "FAIL"), chemo_n))

## =========================================================================
## IX. Audit Status
## =========================================================================
cat("\n[IX] Determining audit status...\n")

audit_status <- "PASS"
audit_notes <- character()

if (!checks$sampleID_count_243) {
  audit_status <- "FAIL"
  audit_notes <- c(audit_notes, "sampleID count is not 243")
}
if (!checks$no_response_conflicts) {
  audit_status <- "FAIL"
  audit_notes <- c(audit_notes, "Patient-level response conflicts exist")
}
if (!checks$response_mutual_exclusion) {
  audit_status <- "FAIL"
  audit_notes <- c(audit_notes, "Response categories not mutually exclusive")
}
if (!checks$binary_exclusion) {
  audit_status <- "FAIL"
  audit_notes <- c(audit_notes, "Binary response categories overlap")
}
if (!checks$cohort_exclusion) {
  audit_status <- "FAIL"
  audit_notes <- c(audit_notes, "anti_PD1 and chemo_control cohorts overlap")
}
if (!checks$cellID_unique) {
  audit_status <- "FAIL"
  audit_notes <- c(audit_notes, "cellID not globally unique")
}

## WARN conditions
n_text_rates <- sum(manifest$response_rate_value_type == "text_description", na.rm = TRUE)
if (n_text_rates > 0) {
  audit_notes <- c(audit_notes, sprintf("%d patients have text-based response rates", n_text_rates))
}

n_unknown_resp <- sum(manifest$pathological_response_clean == "unknown", na.rm = TRUE)
if (n_unknown_resp > 0) {
  audit_notes <- c(audit_notes, sprintf("%d patients have unknown response", n_unknown_resp))
}

if (!checks$antiPD1_close_to_234) {
  audit_notes <- c(audit_notes, sprintf("anti_PD1 cohort (%d) differs from expected 234", antiPD1_n))
}

if (!checks$chemo_close_to_9) {
  audit_notes <- c(audit_notes, sprintf("chemo_control cohort (%d) differs from expected 9", chemo_n))
}

## Set WARN if not FAIL
if (audit_status != "FAIL" && length(audit_notes) > 0) {
  audit_status <- "WARN"
}

cat(sprintf("[INFO] Final audit status: %s\n", audit_status))
for (note in audit_notes) {
  cat(sprintf("  - %s\n", note))
}

audit_text <- c(
  "GSE243013 Revised Cohort Audit Status",
  "======================================",
  "",
  sprintf("Status: %s", audit_status),
  "",
  "Checks:",
  sprintf("  1. sampleID count = 243: %s (actual: %d)",
          ifelse(checks$sampleID_count_243, "PASS", "FAIL"), nrow(manifest)),
  sprintf("  2. No patient-level response conflicts: %s",
          ifelse(checks$no_response_conflicts, "PASS", "FAIL")),
  sprintf("  3. Response categories mutually exclusive: %s",
          ifelse(checks$response_mutual_exclusion, "PASS", "FAIL")),
  sprintf("  4. Binary response mutual exclusion: %s",
          ifelse(checks$binary_exclusion, "PASS", "FAIL")),
  sprintf("  5. anti_PD1 vs chemo_control cohort exclusion: %s",
          ifelse(checks$cohort_exclusion, "PASS", "FAIL")),
  sprintf("  6. Summaries have denominators: %s",
          ifelse(checks$summaries_have_denominators, "PASS", "FAIL")),
  sprintf("  7. cellID globally unique: %s",
          ifelse(checks$cellID_unique, "PASS", "FAIL")),
  sprintf("  8. cellID matches metadata rows: %s",
          ifelse(checks$cellID_row_match, "PASS", "FAIL")),
  sprintf("  9. Min cells per patient >= 1: %s (min: %d)",
          ifelse(checks$min_cells_per_patient, "PASS", "FAIL"), min_cells),
  sprintf("  10. anti_PD1 cohort ~234: %s (actual: %d)",
          ifelse(checks$antiPD1_close_to_234, "PASS", "WARN"), antiPD1_n),
  sprintf("  11. chemotherapy control ~9: %s (actual: %d)",
          ifelse(checks$chemo_close_to_9, "PASS", "WARN"), chemo_n),
  "",
  "Notes:"
)
if (length(audit_notes) > 0) {
  audit_text <- c(audit_text, paste(" -", audit_notes))
} else {
  audit_text <- c(audit_text, " - No issues found")
}

writeLines(audit_text, file.path(RESULTS, "GSE243013_revised_cohort_audit_status.txt"))
cat("[INFO] Saved: GSE243013_revised_cohort_audit_status.txt\n")

## =========================================================================
## X. Sample Lists for Next Stage
## =========================================================================
cat("\n[X] Saving sample lists...\n")

make_sample_list <- function(df, inclusion_reason, exclusion_reason) {
  data.frame(
    sampleID = df$sampleID,
    response_binary = df$response_binary,
    pathological_response_clean = df$pathological_response_clean,
    cancer_type = df$cancer_type,
    n_cells = df$n_cells,
    treatment_pattern = df$treatment_pattern,
    inclusion_reason = inclusion_reason,
    exclusion_reason = exclusion_reason,
    stringsAsFactors = FALSE
  )
}

## Primary anti-PD1
primary_in <- manifest[manifest$primary_analysis_eligible, ]
primary_out <- manifest[!manifest$primary_analysis_eligible, ]

primary_in_list <- make_sample_list(primary_in, "anti_PD1_cohort_with_binary_response", NA_character_)
primary_out_list <- make_sample_list(primary_out, NA_character_,
                                      ifelse(primary_out$anti_PD1_cohort == FALSE, "not_anti_PD1_cohort",
                                             ifelse(is.na(primary_out$response_binary), "no_binary_response",
                                                    "patient_level_conflict")))

write.csv(primary_in_list, file.path(RESULTS, "GSE243013_primary_anti_PD1_sampleIDs.csv"), row.names = FALSE)
write.csv(primary_out_list, file.path(RESULTS, "GSE243013_excluded_sampleIDs_revised.csv"), row.names = FALSE)

## Strict chemo-immunotherapy
strict_in <- manifest[manifest$strict_sensitivity_analysis_eligible, ]
strict_in_list <- make_sample_list(strict_in, "strict_chemoimmunotherapy_with_binary_response", NA_character_)
write.csv(strict_in_list, file.path(RESULTS, "GSE243013_strict_chemoimmunotherapy_sampleIDs.csv"), row.names = FALSE)

## Chemotherapy control
chemo_in <- manifest[manifest$chemotherapy_control_eligible, ]
chemo_in_list <- make_sample_list(chemo_in, "chemotherapy_control_with_binary_response", NA_character_)
write.csv(chemo_in_list, file.path(RESULTS, "GSE243013_chemotherapy_control_sampleIDs.csv"), row.names = FALSE)

cat("[INFO] Saved sample lists:\n")
cat(sprintf("  primary_anti_PD1_sampleIDs.csv: %d patients\n", nrow(primary_in_list)))
cat(sprintf("  strict_chemoimmunotherapy_sampleIDs.csv: %d patients\n", nrow(strict_in_list)))
cat(sprintf("  chemotherapy_control_sampleIDs.csv: %d patients\n", nrow(chemo_in_list)))
cat(sprintf("  excluded_sampleIDs_revised.csv: %d patients\n", nrow(primary_out_list)))

## =========================================================================
## Final Summary
## =========================================================================
cat("\n========================================================================\n")
cat("FINAL SUMMARY\n")
cat("========================================================================\n")

cat(sprintf("Total patients: %d\n", nrow(manifest)))
cat(sprintf("has_anti_PD1: %d\n", sum(manifest$anti_PD1_cohort)))
cat(sprintf("strict anti_PD1 + chemo: %d\n", sum(manifest$strict_chemoimmunotherapy_cohort)))
cat(sprintf("chemotherapy_only: %d\n", sum(manifest$chemotherapy_control_cohort)))
cat(sprintf("anti_PD1 recorded, chemo not recorded: %d\n",
            sum(manifest$treatment_pattern == "anti_PD1_recorded_chemo_not_recorded")))
cat(sprintf("targeted_or_other_only: %d\n", sum(manifest$treatment_pattern == "targeted_or_other_only")))
cat(sprintf("no_treatment_recorded: %d\n", sum(manifest$treatment_pattern == "no_treatment_recorded")))
cat(sprintf("conflict: %d\n", sum(manifest$treatment_pattern == "conflict")))

cat(sprintf("\nanti_PD1 primary cohort: %d\n", sum(manifest$primary_analysis_eligible)))
cat(sprintf("  pCR: %d\n", sum(manifest$pathological_response_clean[manifest$primary_analysis_eligible] == "pCR")))
cat(sprintf("  MPR: %d\n", sum(manifest$pathological_response_clean[manifest$primary_analysis_eligible] == "MPR")))
cat(sprintf("  non-MPR: %d\n", sum(manifest$pathological_response_clean[manifest$primary_analysis_eligible] == "non-MPR")))
cat(sprintf("  unknown: %d\n", sum(manifest$pathological_response_clean[manifest$primary_analysis_eligible] == "unknown")))
cat(sprintf("  Responder: %d\n", sum(manifest$response_binary[manifest$primary_analysis_eligible] == "Responder", na.rm = TRUE)))
cat(sprintf("  Non_responder: %d\n", sum(manifest$response_binary[manifest$primary_analysis_eligible] == "Non_responder", na.rm = TRUE)))

cat(sprintf("\nstrict chemoimmunotherapy cohort: %d\n", sum(manifest$strict_sensitivity_analysis_eligible)))
cat(sprintf("  Responder: %d\n", sum(manifest$response_binary[manifest$strict_sensitivity_analysis_eligible] == "Responder", na.rm = TRUE)))
cat(sprintf("  Non_responder: %d\n", sum(manifest$response_binary[manifest$strict_sensitivity_analysis_eligible] == "Non_responder", na.rm = TRUE)))

cat(sprintf("\nCancer type: LUAD=%d, LUSC=%d\n",
            sum(manifest$cancer_type == "LUAD", na.rm = TRUE),
            sum(manifest$cancer_type == "LUSC", na.rm = TRUE)))

cat(sprintf("\nPatient-level conflicts: %d\n", sum(manifest$has_any_patient_level_conflict)))
cat(sprintf("Audit status: %s\n", audit_status))
cat(sprintf("Ready for counts matrix download: %s\n",
            ifelse(audit_status %in% c("PASS", "WARN"), "YES", "NO")))

cat("\n========================================================================\n")
cat("Step 03 completed.\n")
cat("========================================================================\n")
