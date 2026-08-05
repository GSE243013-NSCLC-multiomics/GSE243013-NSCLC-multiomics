## =========================================================================
## Step 02: Patient-Level Cohort Audit & Analysis Cohort Construction
## =========================================================================

options(timeout = 3600)

cat("========================================================================\n")
cat("Step 02: Patient-Level Cohort Audit & Analysis Cohort Construction\n")
cat("========================================================================\n\n")

RAW_DIR  <- "02_data/raw/GSE243013"
RESULTS  <- "03_results"
CONFIG   <- "00_config"

dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG, recursive = TRUE, showWarnings = FALSE)

## =========================================================================
## I. Read Data
## =========================================================================
cat("[I] Reading metadata...\n")

META_PATH <- file.path(RAW_DIR, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")
meta <- NULL

tryCatch({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cat("[INFO] Installing data.table from CRAN...\n")
    install.packages("data.table", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  meta <- data.table::fread(META_PATH, header = TRUE, data.table = FALSE)
  cat(sprintf("[INFO] Read with data.table::fread: %d rows x %d columns\n",
              nrow(meta), ncol(meta)))
}, error = function(e) {
  cat("[WARNING] data.table::fread failed:", conditionMessage(e), "\n")
  meta <<- read.csv(gzfile(META_PATH), check.names = FALSE)
  cat(sprintf("[INFO] Read with read.csv: %d rows x %d columns\n",
              nrow(meta), ncol(meta)))
})

cat(sprintf("[INFO] Columns: %s\n", paste(names(meta), collapse = ", ")))

## Verify required fields
required_fields <- c("sampleID", "cellID", "pathological_response",
                      "pathological_response_rate", "anti-PD1_therapy",
                      "chemotherapy", "targeted_therapy", "cancer_type",
                      "pre_treatment_staging", "major_cell_type", "sub_cell_type")

cat("\n[INFO] Checking required fields:\n")
for (f in required_fields) {
  found <- f %in% names(meta)
  cat(sprintf("  %s: %s\n", f, ifelse(found, "FOUND", "MISSING")))
}

## =========================================================================
## II. Clean Strings
## =========================================================================
cat("\n[II] Cleaning string fields...\n")

## Trim whitespace and convert empty strings to NA for all character columns
for (cn in names(meta)) {
  if (is.character(meta[[cn]])) {
    meta[[cn]] <- trimws(meta[[cn]])
    meta[[cn]][meta[[cn]] == ""] <- NA_character_
  }
}

## Create lowercase version for comparison
meta$pathological_response_lower <- tolower(meta$pathological_response)

## Create cleaned pathological_response
meta$pathological_response_clean <- meta$pathological_response

## Normalization rules (case-insensitive via lower)
pr <- meta$pathological_response_clean
pr_lower <- tolower(pr)

## pCR
idx <- grepl("^pcr$", pr_lower)
pr[idx] <- "pCR"

## MPR (but not non-MPR)
idx <- grepl("^mpr$", pr_lower) & !grepl("non", pr_lower)
pr[idx] <- "MPR"

## non-MPR variants
idx <- grepl("non.?mpr", pr_lower) | grepl("^nmpR$", pr_lower) | grepl("^non mpr$", pr_lower)
pr[idx] <- "non-MPR"

## unknown variants
idx <- grepl("^unk[n]?own$", pr_lower) | grepl("^unkown$", pr_lower) |
       grepl("^unknowm$", pr_lower) | grepl("^unkownn$", pr_lower) |
       grepl("^unknown$", pr_lower)
pr[idx] <- "unknown"

## Mark UNRESOLVED
valid_values <- c("pCR", "MPR", "non-MPR", "unknown")
idx_unresolved <- !(pr %in% valid_values) & !is.na(pr)
pr[idx_unresolved] <- "UNRESOLVED"

meta$pathological_response_clean <- pr

cat("[INFO] pathological_response_clean value counts:\n")
print(table(meta$pathological_response_clean, useNA = "ifany"))

## =========================================================================
## III. Sample-Level Consistency Audit
## =========================================================================
cat("\n[III] Building sample-level consistency audit...\n")

audit_fields <- c("pathological_response_clean", "pathological_response_rate",
                   "anti-PD1_therapy", "chemotherapy", "targeted_therapy",
                   "cancer_type", "pre_treatment_staging")

sample_ids <- sort(unique(meta$sampleID))
cat(sprintf("[INFO] Unique sampleIDs: %d\n", length(sample_ids)))

audit_rows <- list()
conflict_rows <- list()

for (sid in sample_ids) {
  idx <- meta$sampleID == sid
  sub <- meta[idx, ]

  row_data <- data.frame(
    sampleID = sid,
    n_cells = sum(idx),
    n_cellID_unique = length(unique(sub$cellID)),
    n_major_cell_types = length(unique(sub$major_cell_type)),
    n_sub_cell_types = length(unique(sub$sub_cell_type)),
    stringsAsFactors = FALSE
  )

  for (af in audit_fields) {
    vals <- sub[[af]]
    vals_no_na <- vals[!is.na(vals)]
    n_unique <- length(unique(vals_no_na))
    all_vals <- paste(sort(unique(vals_no_na)), collapse = " | ")
    has_conflict <- n_unique > 1

    row_data[[paste0(af, "_n_unique")]] <- n_unique
    row_data[[paste0(af, "_all_values")]] <- all_vals
    row_data[[paste0(af, "_conflict")]] <- has_conflict

    if (has_conflict) {
      conflict_rows[[length(conflict_rows) + 1]] <- data.frame(
        sampleID = sid,
        field = af,
        n_unique_values = n_unique,
        all_values = all_vals,
        stringsAsFactors = FALSE
      )
    }
  }

  audit_rows[[sid]] <- row_data
}

audit_df <- do.call(rbind, audit_rows)
rownames(audit_df) <- NULL

write.csv(audit_df, file.path(RESULTS, "GSE243013_sample_level_consistency_audit.csv"),
          row.names = FALSE)
cat(sprintf("[INFO] Saved: GSE243013_sample_level_consistency_audit.csv (%d rows)\n",
            nrow(audit_df)))

if (length(conflict_rows) > 0) {
  conflicts_df <- do.call(rbind, conflict_rows)
  rownames(conflicts_df) <- NULL
} else {
  conflicts_df <- data.frame(sampleID = character(), field = character(),
                              n_unique_values = integer(), all_values = character(),
                              stringsAsFactors = FALSE)
}
write.csv(conflicts_df, file.path(RESULTS, "GSE243013_sample_level_conflicts.csv"),
          row.names = FALSE)
cat(sprintf("[INFO] Saved: GSE243013_sample_level_conflicts.csv (%d conflicts)\n",
            nrow(conflicts_df)))

## =========================================================================
## IV. Response Cross Table
## =========================================================================
cat("\n[IV] Building response cross table...\n")

## Use unique sampleID + response_clean combinations
sample_response <- unique(meta[, c("sampleID", "pathological_response_clean")])
sample_response <- sample_response[!is.na(sample_response$pathological_response_clean), ]

response_categories <- c("pCR", "MPR", "non-MPR", "unknown", "UNRESOLVED")
cross_table <- data.frame(sampleID = sample_ids, stringsAsFactors = FALSE)
for (rc in response_categories) {
  cross_table[[rc]] <- FALSE
}

for (i in seq_len(nrow(sample_response))) {
  sid <- sample_response$sampleID[i]
  resp <- sample_response$pathological_response_clean[i]
  if (resp %in% response_categories) {
    cross_table[cross_table$sampleID == sid, resp] <- TRUE
  }
}

write.csv(cross_table, file.path(RESULTS, "GSE243013_sample_response_cross_table.csv"),
          row.names = FALSE)
cat(sprintf("[INFO] Saved: GSE243013_sample_response_cross_table.csv\n"))

## Response mapping conflicts
cat("[INFO] Checking for response mapping conflicts...\n")
conflict_response <- data.frame(sampleID = character(), issue = character(),
                                 details = character(), stringsAsFactors = FALSE)

for (sid in sample_ids) {
  row <- cross_table[cross_table$sampleID == sid, ]
  true_cats <- response_categories[sapply(response_categories, function(rc) row[[rc]])]

  if (length(true_cats) > 1) {
    conflict_response <- rbind(conflict_response, data.frame(
      sampleID = sid,
      issue = "multiple_response_categories",
      details = paste(true_cats, collapse = " | "),
      stringsAsFactors = FALSE
    ))
  } else if (length(true_cats) == 0) {
    conflict_response <- rbind(conflict_response, data.frame(
      sampleID = sid,
      issue = "missing_response",
      details = "",
      stringsAsFactors = FALSE
    ))
  }

  ## Check for UNRESOLVED in original data
  sub <- meta[meta$sampleID == sid, ]
  unresolved_vals <- unique(sub$pathological_response_clean[sub$pathological_response_clean == "UNRESOLVED"])
  if (length(unresolved_vals) > 0) {
    orig_vals <- unique(sub$pathological_response[!is.na(sub$pathological_response)])
    conflict_response <- rbind(conflict_response, data.frame(
      sampleID = sid,
      issue = "unresolved_response_value",
      details = paste(orig_vals, collapse = " | "),
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(conflict_response, file.path(RESULTS, "GSE243013_response_mapping_conflicts.csv"),
          row.names = FALSE)
cat(sprintf("[INFO] Saved: GSE243013_response_mapping_conflicts.csv (%d rows)\n",
            nrow(conflict_response)))

## =========================================================================
## V. Response Rate Cross-Validation
## =========================================================================
cat("\n[V] Cross-validating response rate...\n")

rate_summary <- list(
  min = min(meta$pathological_response_rate, na.rm = TRUE),
  max = max(meta$pathological_response_rate, na.rm = TRUE),
  median = median(meta$pathological_response_rate, na.rm = TRUE),
  n_unique = length(unique(meta$pathological_response_rate[!is.na(meta$pathological_response_rate)])),
  n_missing = sum(is.na(meta$pathological_response_rate))
)

cat(sprintf("[INFO] Rate summary: min=%s, max=%s, median=%s, unique=%d, missing=%d\n",
            rate_summary$min, rate_summary$max, rate_summary$median,
            rate_summary$n_unique, rate_summary$n_missing))

## Determine scale
rate_vals <- as.numeric(meta$pathological_response_rate[!is.na(meta$pathological_response_rate)])
rate_vals <- rate_vals[!is.nan(rate_vals)]
max_rate <- max(rate_vals, na.rm = TRUE)
if (max_rate <= 1.0) {
  rate_scale <- "0_to_1"
  cat("[INFO] Rate scale: 0 to 1 (proportion)\n")
} else {
  rate_scale <- "0_to_100"
  cat("[INFO] Rate scale: 0 to 100 (percentage)\n")
}

## Build cross-check at sample level
sample_rate <- unique(meta[, c("sampleID", "pathological_response_clean", "pathological_response_rate")])
sample_rate <- sample_rate[!is.na(sample_rate$pathological_response_clean), ]

sample_rate$rate_scale <- rate_scale
sample_rate$expected_category_from_rate <- NA_character_
sample_rate$category_rate_consistent <- NA_character_
sample_rate$issue <- ""

for (i in seq_len(nrow(sample_rate))) {
  resp <- sample_rate$pathological_response_clean[i]
  rate <- suppressWarnings(as.numeric(sample_rate$pathological_response_rate[i]))

  if (is.na(rate)) {
    sample_rate$issue[i] <- "rate_missing"
    sample_rate$category_rate_consistent[i] <- "unknown"
    next
  }

  if (rate_scale == "0_to_1") {
    if (abs(rate) < 1e-8) {
      sample_rate$expected_category_from_rate[i] <- "pCR"
    } else if (rate > 0 && rate <= 0.10 + 1e-8) {
      sample_rate$expected_category_from_rate[i] <- "MPR"
    } else if (rate > 0.10 + 1e-8) {
      sample_rate$expected_category_from_rate[i] <- "non-MPR"
    } else {
      sample_rate$expected_category_from_rate[i] <- "uninterpretable"
    }
  } else {
    if (abs(rate) < 1e-8) {
      sample_rate$expected_category_from_rate[i] <- "pCR"
    } else if (rate > 0 && rate <= 10 + 1e-8) {
      sample_rate$expected_category_from_rate[i] <- "MPR"
    } else if (rate > 10 + 1e-8) {
      sample_rate$expected_category_from_rate[i] <- "non-MPR"
    } else {
      sample_rate$expected_category_from_rate[i] <- "uninterpretable"
    }
  }

  exp_cat <- sample_rate$expected_category_from_rate[i]
  actual <- resp
  if (actual %in% c("pCR", "MPR", "non-MPR")) {
    consistent <- (exp_cat == actual)
    sample_rate$category_rate_consistent[i] <- ifelse(consistent, "consistent", "INCONSISTENT")
    if (!consistent) {
      sample_rate$issue[i] <- sprintf("expected_%s_but_got_%s", exp_cat, actual)
    } else {
      sample_rate$issue[i] <- "ok"
    }
  } else {
    sample_rate$category_rate_consistent[i] <- "cannot_check"
    sample_rate$issue[i] <- "non_standard_category"
  }
}

write.csv(sample_rate, file.path(RESULTS, "GSE243013_response_rate_crosscheck.csv"),
          row.names = FALSE)
cat(sprintf("[INFO] Saved: GSE243013_response_rate_crosscheck.csv\n"))

## =========================================================================
## VI. Treatment Classification
## =========================================================================
cat("\n[VI] Classifying treatment...\n")

## Stats for anti-PD1_therapy
cat("\n[INFO] anti-PD1_therapy value counts (by sampleID):\n")
pd1_vals <- unique(meta[, c("sampleID", "anti-PD1_therapy")])
pd1_counts <- table(pd1_vals[["anti-PD1_therapy"]], useNA = "ifany")
print(sort(pd1_counts, decreasing = TRUE))

cat("\n[INFO] chemotherapy value counts (by sampleID):\n")
chemo_vals <- unique(meta[, c("sampleID", "chemotherapy")])
chemo_counts <- table(chemo_vals$chemotherapy, useNA = "ifany")
print(sort(chemo_counts, decreasing = TRUE))

## Classify at sample level
sample_treatment <- unique(meta[, c("sampleID", "anti-PD1_therapy", "chemotherapy")])

sample_treatment$treatment_class_candidate <- "unresolved"

for (i in seq_len(nrow(sample_treatment))) {
  pd1 <- sample_treatment[["anti-PD1_therapy"]][i]
  chemo <- sample_treatment[["chemotherapy"]][i]

  ## Check if anti-PD1 has a drug name
  pd1_has_drug <- !is.na(pd1) && !tolower(pd1) %in% c("no", "none", "na")
  chemo_has_drug <- !is.na(chemo) && !tolower(chemo) %in% c("no", "none", "na")

  if (pd1_has_drug && chemo_has_drug) {
    sample_treatment$treatment_class_candidate[i] <- "chemo_immunotherapy"
  } else if (!pd1_has_drug && chemo_has_drug) {
    sample_treatment$treatment_class_candidate[i] <- "chemotherapy_only"
  } else {
    sample_treatment$treatment_class_candidate[i] <- "unresolved"
  }
}

cat("\n[INFO] Treatment classification results:\n")
print(table(sample_treatment$treatment_class_candidate))

write.csv(sample_treatment, file.path(RESULTS, "GSE243013_treatment_classification.csv"),
          row.names = FALSE)
cat(sprintf("[INFO] Saved: GSE243013_treatment_classification.csv\n"))

## Treatment value counts
treatment_value_counts <- data.frame(
  field = character(), value = character(), n_sampleIDs = integer(),
  stringsAsFactors = FALSE
)
for (field in c("anti-PD1_therapy", "chemotherapy", "targeted_therapy")) {
  vals <- unique(meta[, c("sampleID", field)])
  tbl <- table(vals[[field]], useNA = "ifany")
  for (v in names(tbl)) {
    treatment_value_counts <- rbind(treatment_value_counts, data.frame(
      field = field, value = ifelse(is.na(v), "<NA>", v),
      n_sampleIDs = as.integer(tbl[v]),
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(treatment_value_counts, file.path(RESULTS, "GSE243013_treatment_value_counts.csv"),
          row.names = FALSE)
cat(sprintf("[INFO] Saved: GSE243013_treatment_value_counts.csv\n"))

## =========================================================================
## VII. Patient Manifest
## =========================================================================
cat("\n[VII] Creating patient manifest...\n")

## Merge audit, treatment, and response
patient_manifest <- data.frame(
  sampleID = sample_ids,
  stringsAsFactors = FALSE
)

## Add cell counts from audit
for (i in seq_len(nrow(patient_manifest))) {
  sid <- patient_manifest$sampleID[i]
  a_row <- audit_df[audit_df$sampleID == sid, ]
  if (nrow(a_row) > 0) {
    patient_manifest$n_cells[i] <- a_row$n_cells
    patient_manifest$n_major_cell_types[i] <- a_row$n_major_cell_types
    patient_manifest$n_sub_cell_types[i] <- a_row$n_sub_cell_types
  }
}

## Add response fields
patient_manifest$pathological_response_original <- NA_character_
patient_manifest$pathological_response_clean <- NA_character_
patient_manifest$pathological_response_rate <- NA_character_

for (i in seq_len(nrow(patient_manifest))) {
  sid <- patient_manifest$sampleID[i]
  sub <- meta[meta$sampleID == sid, ]
  ## Take first non-NA value for original (should be consistent within sample)
  orig_vals <- unique(sub$pathological_response[!is.na(sub$pathological_response)])
  patient_manifest$pathological_response_original[i] <- ifelse(length(orig_vals) > 0, orig_vals[1], NA_character_)

  clean_vals <- unique(sub$pathological_response_clean[!is.na(sub$pathological_response_clean)])
  patient_manifest$pathological_response_clean[i] <- ifelse(length(clean_vals) > 0, clean_vals[1], NA_character_)

  rate_vals <- unique(sub$pathological_response_rate[!is.na(sub$pathological_response_rate)])
  patient_manifest$pathological_response_rate[i] <- ifelse(length(rate_vals) > 0, as.character(rate_vals[1]), NA_character_)
}

## Add treatment fields
for (i in seq_len(nrow(patient_manifest))) {
  sid <- patient_manifest$sampleID[i]
  t_row <- sample_treatment[sample_treatment$sampleID == sid, ]
  if (nrow(t_row) > 0) {
    patient_manifest$anti_PD1_therapy[i] <- t_row[["anti-PD1_therapy"]][1]
    patient_manifest$chemotherapy[i] <- t_row$chemotherapy[1]
    patient_manifest$treatment_class_candidate[i] <- t_row$treatment_class_candidate[1]
  }
}

## Add targeted_therapy, cancer_type, pre_treatment_staging
for (i in seq_len(nrow(patient_manifest))) {
  sid <- patient_manifest$sampleID[i]
  sub <- meta[meta$sampleID == sid, ]

  tt_vals <- unique(sub$targeted_therapy[!is.na(sub$targeted_therapy)])
  patient_manifest$targeted_therapy[i] <- ifelse(length(tt_vals) > 0, tt_vals[1], NA_character_)

  ct_vals <- unique(sub$cancer_type[!is.na(sub$cancer_type)])
  patient_manifest$cancer_type[i] <- ifelse(length(ct_vals) > 0, ct_vals[1], NA_character_)

  stg_vals <- unique(sub$pre_treatment_staging[!is.na(sub$pre_treatment_staging)])
  patient_manifest$pre_treatment_staging[i] <- ifelse(length(stg_vals) > 0, stg_vals[1], NA_character_)
}

## Add conflict flags
patient_manifest$has_response_conflict <- FALSE
patient_manifest$has_treatment_conflict <- FALSE
patient_manifest$has_cancer_type_conflict <- FALSE
patient_manifest$has_stage_conflict <- FALSE

for (i in seq_len(nrow(patient_manifest))) {
  sid <- patient_manifest$sampleID[i]
  a_row <- audit_df[audit_df$sampleID == sid, ]
  if (nrow(a_row) > 0) {
    patient_manifest$has_response_conflict[i] <- a_row$pathological_response_clean_conflict
    patient_manifest$has_treatment_conflict[i] <- a_row$`anti-PD1_therapy_conflict` |
                                                   a_row$chemotherapy_conflict |
                                                   a_row$targeted_therapy_conflict
    patient_manifest$has_cancer_type_conflict[i] <- a_row$cancer_type_conflict
    patient_manifest$has_stage_conflict[i] <- a_row$pre_treatment_staging_conflict
  }
}

patient_manifest$has_any_patient_level_conflict <- patient_manifest$has_response_conflict |
  patient_manifest$has_treatment_conflict | patient_manifest$has_cancer_type_conflict |
  patient_manifest$has_stage_conflict

write.csv(patient_manifest, file.path(RESULTS, "GSE243013_patient_manifest.csv"),
          row.names = FALSE)
cat(sprintf("[INFO] Saved: GSE243013_patient_manifest.csv (%d rows)\n",
            nrow(patient_manifest)))

## =========================================================================
## VIII. Analysis Endpoints
## =========================================================================
cat("\n[VIII] Creating analysis endpoints...\n")

## response_three_level
patient_manifest$response_three_level <- patient_manifest$pathological_response_clean
invalid <- !(patient_manifest$response_three_level %in% c("pCR", "MPR", "non-MPR", "unknown"))
patient_manifest$response_three_level[invalid] <- "UNRESOLVED"
patient_manifest$response_three_level[is.na(patient_manifest$response_three_level)] <- "UNRESOLVED"

## response_binary
patient_manifest$response_binary <- NA_character_
patient_manifest$response_binary[patient_manifest$response_three_level %in% c("pCR", "MPR")] <- "Responder"
patient_manifest$response_binary[patient_manifest$response_three_level == "non-MPR"] <- "Non_responder"

## primary_cohort_eligible
patient_manifest$primary_cohort_eligible <- (
  patient_manifest$treatment_class_candidate == "chemo_immunotherapy" &
  !is.na(patient_manifest$response_binary) &
  !patient_manifest$has_response_conflict &
  !patient_manifest$has_treatment_conflict &
  !patient_manifest$has_cancer_type_conflict &
  !is.na(patient_manifest$sampleID)
)

cat(sprintf("[INFO] primary_cohort_eligible TRUE count: %d\n",
            sum(patient_manifest$primary_cohort_eligible, na.rm = TRUE)))

## Save updated manifest
write.csv(patient_manifest, file.path(RESULTS, "GSE243013_patient_manifest.csv"),
          row.names = FALSE)

## =========================================================================
## IX. Cohort Statistics
## =========================================================================
cat("\n[IX] Generating cohort statistics...\n")

## A. Overall stats
cat("\n--- A. Overall 243 sampleID stats ---\n")
cat(sprintf("Total sampleIDs: %d\n", nrow(patient_manifest)))
cat(sprintf("Total cells: %d\n", sum(patient_manifest$n_cells)))
cat(sprintf("chemo_immunotherapy: %d\n",
            sum(patient_manifest$treatment_class_candidate == "chemo_immunotherapy")))
cat(sprintf("chemotherapy_only: %d\n",
            sum(patient_manifest$treatment_class_candidate == "chemotherapy_only")))
cat(sprintf("unresolved treatment: %d\n",
            sum(patient_manifest$treatment_class_candidate == "unresolved")))

overall_stats <- data.frame(
  metric = c("total_sampleIDs", "total_cells", "chemo_immunotherapy", "chemotherapy_only", "unresolved_treatment"),
  value = c(nrow(patient_manifest), sum(patient_manifest$n_cells),
            sum(patient_manifest$treatment_class_candidate == "chemo_immunotherapy"),
            sum(patient_manifest$treatment_class_candidate == "chemotherapy_only"),
            sum(patient_manifest$treatment_class_candidate == "unresolved")),
  stringsAsFactors = FALSE
)

## B. Primary immunotherapy cohort
cat("\n--- B. Primary immunotherapy cohort ---\n")
primary <- patient_manifest[patient_manifest$primary_cohort_eligible, ]
cat(sprintf("Primary cohort size: %d\n", nrow(primary)))
cat(sprintf("pCR: %d\n", sum(primary$response_three_level == "pCR")))
cat(sprintf("MPR: %d\n", sum(primary$response_three_level == "MPR")))
cat(sprintf("non-MPR: %d\n", sum(primary$response_three_level == "non-MPR")))
cat(sprintf("unknown: %d\n", sum(primary$response_three_level == "unknown")))
cat(sprintf("Responder (pCR+MPR): %d\n", sum(primary$response_binary == "Responder", na.rm = TRUE)))
cat(sprintf("Non_responder: %d\n", sum(primary$response_binary == "Non_responder", na.rm = TRUE)))

excluded <- patient_manifest[!patient_manifest$primary_cohort_eligible, ]
cat(sprintf("Excluded: %d\n", nrow(excluded)))

exclusion_reasons <- data.frame(
  sampleID = character(), reason = character(), stringsAsFactors = FALSE
)
for (i in seq_len(nrow(excluded))) {
  sid <- excluded$sampleID[i]
  reasons <- character()
  if (excluded$treatment_class_candidate[i] != "chemo_immunotherapy") {
    reasons <- c(reasons, excluded$treatment_class_candidate[i])
  }
  if (is.na(excluded$response_binary[i])) {
    reasons <- c(reasons, "no_binary_response")
  }
  if (excluded$has_response_conflict[i]) reasons <- c(reasons, "response_conflict")
  if (excluded$has_treatment_conflict[i]) reasons <- c(reasons, "treatment_conflict")
  if (excluded$has_cancer_type_conflict[i]) reasons <- c(reasons, "cancer_type_conflict")
  if (is.na(excluded$sampleID[i])) reasons <- c(reasons, "missing_sampleID")
  exclusion_reasons <- rbind(exclusion_reasons, data.frame(
    sampleID = sid, reason = paste(reasons, collapse = "; "),
    stringsAsFactors = FALSE
  ))
}

## Save cohort summary
cohort_summary <- data.frame(
  category = c("total_sampleIDs", "total_cells", "chemo_immunotherapy", "chemotherapy_only",
               "unresolved_treatment", "primary_cohort_size", "pCR", "MPR", "non-MPR",
               "unknown_primary", "Responder", "Non_responder", "excluded"),
  count = c(nrow(patient_manifest), sum(patient_manifest$n_cells),
            sum(patient_manifest$treatment_class_candidate == "chemo_immunotherapy"),
            sum(patient_manifest$treatment_class_candidate == "chemotherapy_only"),
            sum(patient_manifest$treatment_class_candidate == "unresolved"),
            nrow(primary),
            sum(primary$response_three_level == "pCR"),
            sum(primary$response_three_level == "MPR"),
            sum(primary$response_three_level == "non-MPR"),
            sum(primary$response_three_level == "unknown"),
            sum(primary$response_binary == "Responder", na.rm = TRUE),
            sum(primary$response_binary == "Non_responder", na.rm = TRUE),
            nrow(excluded)),
  stringsAsFactors = FALSE
)
write.csv(cohort_summary, file.path(RESULTS, "GSE243013_analysis_cohort_summary.csv"),
          row.names = FALSE)
write.csv(exclusion_reasons, file.path(RESULTS, "GSE243013_exclusion_reasons.csv"),
          row.names = FALSE)

## C. Response by cancer_type
cat("\n--- C. Response by cancer_type ---\n")
ct_response <- table(patient_manifest$cancer_type, patient_manifest$response_three_level,
                     useNA = "ifany")
print(ct_response)
write.csv(as.data.frame.matrix(ct_response),
          file.path(RESULTS, "GSE243013_response_by_cancer_type.csv"))

## D. Response by PD1 drug
cat("\n--- D. Response by PD1 drug ---\n")
pd1_response <- table(patient_manifest$anti_PD1_therapy, patient_manifest$response_binary,
                      useNA = "ifany")
print(pd1_response)
write.csv(as.data.frame.matrix(pd1_response),
          file.path(RESULTS, "GSE243013_response_by_PD1_drug.csv"))

## E. Response by stage
cat("\n--- E. Response by stage ---\n")
stg_response <- table(patient_manifest$pre_treatment_staging, patient_manifest$response_binary,
                      useNA = "ifany")
print(stg_response)
write.csv(as.data.frame.matrix(stg_response),
          file.path(RESULTS, "GSE243013_response_by_stage.csv"))

## =========================================================================
## X. Cohort Definition Document
## =========================================================================
cat("\n[X] Writing cohort definition...\n")

cohort_def <- c(
  "GSE243013 Cohort Definition",
  "==========================",
  "",
  "Dataset: GSE243013 - A single-cell atlas of immune heterogeneity in anti-PD1-treated NSCLC",
  sprintf("Total scRNA samples: %d (one per patient)", nrow(patient_manifest)),
  "Total cells: 1,254,749",
  "",
  "Primary study population: Patients receiving neoadjuvant chemo-immunotherapy",
  "Exploratory control: Patients receiving chemotherapy only",
  "",
  "Response classification (three-level):",
  "  - pCR: pathological complete response",
  "  - MPR: major pathological response",
  "  - non-MPR: non-major pathological response",
  "  - unknown: response category unknown",
  "  - UNRESOLVED: conflicting or uninterpretable",
  "",
  "Response classification (binary):",
  "  - Responder: pCR or MPR",
  "  - Non_responder: non-MPR",
  "  - NA: unknown or UNRESOLVED",
  "",
  "Primary cohort eligibility criteria:",
  "  1. treatment_class_candidate == chemo_immunotherapy",
  "  2. response_binary is not NA",
  "  3. No patient-level conflicts in response, treatment, cancer_type",
  "  4. sampleID is not missing",
  "",
  "Statistical unit: sampleID (biological replicate, one per patient)",
  "Cell count is NOT used as independent sample size",
  "",
  "Data scope: Immune cells only (T/NK, B, Myeloid)",
  "  - Cannot be used for malignant cell state analysis",
  "  - Post-treatment surgical samples only",
  "",
  "Conclusions should be framed as:",
  "  - Response-associated mechanisms",
  "  - Residual tumor-associated mechanisms",
  ""
)
writeLines(cohort_def, file.path(CONFIG, "GSE243013_cohort_definition.txt"))
cat("[INFO] Saved: 00_config/GSE243013_cohort_definition.txt\n")

## =========================================================================
## XI. Audit Status
## =========================================================================
cat("\n[XI] Generating audit status...\n")

audit_status <- "PASS"
audit_notes <- character()

## Check 1: 243 sampleIDs
if (nrow(patient_manifest) != 243) {
  audit_status <- "FAIL"
  audit_notes <- c(audit_notes, sprintf("sampleID count is %d, expected 243", nrow(patient_manifest)))
}

## Check 2: No sampleID with multiple response categories
n_response_conflicts <- sum(patient_manifest$has_response_conflict)
if (n_response_conflicts > 0) {
  audit_status <- "FAIL"
  audit_notes <- c(audit_notes, sprintf("%d sampleIDs have multiple response categories", n_response_conflicts))
}

## Check 3: No sampleID with multiple treatment categories
n_treatment_conflicts <- sum(patient_manifest$has_treatment_conflict)
if (n_treatment_conflicts > 0) {
  audit_status <- "FAIL"
  audit_notes <- c(audit_notes, sprintf("%d sampleIDs have treatment conflicts", n_treatment_conflicts))
}

## Check 4: Primary cohort pCR/MPR/non-MPR mutual exclusion
if (nrow(primary) > 0) {
  pcr_n <- sum(primary$response_three_level == "pCR")
  mpr_n <- sum(primary$response_three_level == "MPR")
  non_mpr_n <- sum(primary$response_three_level == "non-MPR")
  if (pcr_n + mpr_n + non_mpr_n > nrow(primary)) {
    audit_status <- "FAIL"
    audit_notes <- c(audit_notes, "Primary cohort response categories overlap")
  }
}

## Check 5: Rate consistency
n_inconsistent <- sum(sample_rate$category_rate_consistent == "INCONSISTENT", na.rm = TRUE)
if (n_inconsistent > 10) {
  audit_status <- "FAIL"
  audit_notes <- c(audit_notes, sprintf("%d samples have rate-category inconsistencies", n_inconsistent))
} else if (n_inconsistent > 0) {
  if (audit_status != "FAIL") audit_status <- "WARN"
  audit_notes <- c(audit_notes, sprintf("%d samples have rate-category inconsistencies", n_inconsistent))
}

## Check 6: Treatment classification
n_unresolved_treatment <- sum(patient_manifest$treatment_class_candidate == "unresolved")
if (n_unresolved_treatment > 0) {
  if (audit_status != "FAIL") audit_status <- "WARN"
  audit_notes <- c(audit_notes, sprintf("%d patients with unresolved treatment", n_unresolved_treatment))
}

## Check 7: Unknown response
n_unknown <- sum(patient_manifest$response_three_level == "unknown", na.rm = TRUE)
if (n_unknown > 0) {
  if (audit_status != "FAIL") audit_status <- "WARN"
  audit_notes <- c(audit_notes, sprintf("%d patients with unknown response", n_unknown))
}

## Check 8: Cohort size vs expected 234/9
n_chemo_immuno <- sum(patient_manifest$treatment_class_candidate == "chemo_immunotherapy")
n_chemo_only <- sum(patient_manifest$treatment_class_candidate == "chemotherapy_only")
if (abs(n_chemo_immuno - 234) > 5 || abs(n_chemo_only - 9) > 5) {
  if (audit_status != "FAIL") audit_status <- "WARN"
  audit_notes <- c(audit_notes,
                   sprintf("chemo_immunotherapy=%d (expected ~234), chemotherapy_only=%d (expected ~9)",
                           n_chemo_immuno, n_chemo_only))
}

cat(sprintf("[INFO] Audit status: %s\n", audit_status))
for (note in audit_notes) {
  cat(sprintf("  - %s\n", note))
}

audit_text <- c(
  "GSE243013 Cohort Audit Status",
  "=============================",
  "",
  sprintf("Status: %s", audit_status),
  "",
  "Checks performed:",
  sprintf("  1. sampleID count = 243: %s (actual: %d)",
          ifelse(nrow(patient_manifest) == 243, "PASS", "FAIL"), nrow(patient_manifest)),
  sprintf("  2. No response conflicts: %s (conflicts: %d)",
          ifelse(n_response_conflicts == 0, "PASS", "FAIL"), n_response_conflicts),
  sprintf("  3. No treatment conflicts: %s (conflicts: %d)",
          ifelse(n_treatment_conflicts == 0, "PASS", "FAIL"), n_treatment_conflicts),
  sprintf("  4. Response categories mutually exclusive in primary: %s",
          ifelse(pcr_n + mpr_n + non_mpr_n <= nrow(primary) || nrow(primary) == 0, "PASS", "FAIL")),
  sprintf("  5. Rate-category consistency: %s (inconsistent: %d)",
          ifelse(n_inconsistent <= 10, "PASS/WARN", "FAIL"), n_inconsistent),
  sprintf("  6. Treatment classification: %s (unresolved: %d)",
          ifelse(n_unresolved_treatment == 0, "PASS", "WARN"), n_unresolved_treatment),
  sprintf("  7. Unknown response: %s (count: %d)",
          ifelse(n_unknown == 0, "PASS", "WARN"), n_unknown),
  sprintf("  8. Cohort size vs 234/9: %s (actual: %d/%d)",
          ifelse(abs(n_chemo_immuno - 234) <= 5 && abs(n_chemo_only - 9) <= 5, "PASS", "WARN"),
          n_chemo_immuno, n_chemo_only),
  "",
  "Notes:"
)
if (length(audit_notes) > 0) {
  audit_text <- c(audit_text, paste(" -", audit_notes))
} else {
  audit_text <- c(audit_text, " - No issues found")
}

writeLines(audit_text, file.path(RESULTS, "GSE243013_cohort_audit_status.txt"))
cat("[INFO] Saved: 03_results/GSE243013_cohort_audit_status.txt\n")

## =========================================================================
## Final Summary
## =========================================================================
cat("\n========================================================================\n")
cat("FINAL SUMMARY\n")
cat("========================================================================\n")

cat(sprintf("sampleID total: %d\n", nrow(patient_manifest)))
cat(sprintf("Total cells: %d\n", sum(patient_manifest$n_cells)))
cat(sprintf("chemo_immunotherapy: %d\n", n_chemo_immuno))
cat(sprintf("chemotherapy_only: %d\n", n_chemo_only))
cat(sprintf("unresolved treatment: %d\n", n_unresolved_treatment))
cat(sprintf("pCR: %d\n", pcr_n))
cat(sprintf("MPR: %d\n", mpr_n))
cat(sprintf("non-MPR: %d\n", non_mpr_n))
cat(sprintf("unknown: %d\n", n_unknown))
cat(sprintf("Responder: %d\n", sum(patient_manifest$response_binary == "Responder", na.rm = TRUE)))
cat(sprintf("Non_responder: %d\n", sum(patient_manifest$response_binary == "Non_responder", na.rm = TRUE)))
cat(sprintf("Primary cohort: %d\n", nrow(primary)))
cat(sprintf("Excluded: %d\n", nrow(excluded)))
cat(sprintf("Response conflicts: %d\n", n_response_conflicts))
cat(sprintf("Treatment conflicts: %d\n", n_treatment_conflicts))
cat(sprintf("Rate inconsistencies: %d\n", n_inconsistent))
cat(sprintf("Audit status: %s\n", audit_status))

## Cancer type distribution
cat("\nCancer type distribution:\n")
ct_tab <- table(patient_manifest$cancer_type, useNA = "ifany")
print(sort(ct_tab, decreasing = TRUE))

cat("\n========================================================================\n")
cat("Step 02 completed.\n")
cat("========================================================================\n")
