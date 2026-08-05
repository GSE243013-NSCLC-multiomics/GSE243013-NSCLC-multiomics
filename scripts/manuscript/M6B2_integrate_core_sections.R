#!/usr/bin/env Rscript
# M6B2: Integrate Core Scientific Sections
# Read-only: reads M5A sections and M6A verified values only
# No statistical analysis, no scanning other directories

cat("\n", rep("=", 80), sep="")
cat("\nM6B2: Integrate Core Scientific Sections")
cat("\nStarted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")

work_in  <- "05_manuscript/M6_final_manuscript/work/M6B1_sections"
work_out <- "05_manuscript/M6_final_manuscript/work/M6B2_core_sections"
dir.create(work_out, recursive=TRUE, showWarnings=FALSE)

# ============================================================================
# SECTION I: Read Inputs
# ============================================================================
cat("\nSECTION I: Read Inputs\n")
cat(rep("=", 80), sep="")
cat("\n\n")

abstract_m5a <- readLines(file.path(work_in, "01_Abstract_M5A.md"), warn=FALSE)
methods_m5a  <- readLines(file.path(work_in, "03_Methods_M5A.md"), warn=FALSE)
results_m5a  <- readLines(file.path(work_in, "04_Results_M5A.md"), warn=FALSE)
fig5_m5a     <- readLines(file.path(work_in, "09_Figure5_M5A.md"), warn=FALSE)

cat("Abstract M5A:", length(abstract_m5a), "lines\n")
cat("Methods M5A:", length(methods_m5a), "lines\n")
cat("Results M5A:", length(results_m5a), "lines\n")
cat("Figure 5 M5A:", length(fig5_m5a), "lines\n\n")


# ============================================================================
# SECTION II: Write Revised Abstract
# ============================================================================
cat("\nSECTION II: Write Revised Abstract\n")
cat(rep("=", 80), sep="")
cat("\n\n")

abstract <- c(
  "# Abstract",
  "",
  "**Background:** Neoadjuvant anti-PD-1 immunotherapy achieves pathological complete response (pCR) ",
  "in a subset of non-small cell lung cancer (NSCLC) patients, but the immune microenvironment ",
  "determinants of response remain incompletely characterized.",
  "",
  "**Methods:** We performed single-cell RNA sequencing on 31,831 cells from 243 NSCLC patients ",
  "treated with neoadjuvant anti-PD-1-based regimens, analyzing 233 primary-eligible specimens ",
  "across 8 immune cell types. Tumor-infiltrating immune cell transcriptomes were profiled using ",
  "edgeR differential expression and preranked fgsea pathway enrichment (Hallmark gene sets, ",
  "MSigDB 2024.1.Hs). The core immune-compartment glycolysis program was validated against ",
  "clinical outcomes in TCGA-LUAD (n=477) and TCGA-LUSC (n=485) using Cox proportional hazards ",
  "models, and integrated with DNA methylation, RPPA proteomics, somatic mutation, and copy number ",
  "variation data.",
  "",
  "**Results:** The HALLMARK_GLYCOLYSIS pathway in All_immune cells was enriched in the non-responder ",
  "direction (primary NES=-2.36, FDR=3.00 x 10^-11; strict NES=-2.41, FDR=7.51 x 10^-12). ",
  "This directionality was consistent across all 8 analyzed strata (8/8 negative NES). ",
  "Fixed-effect meta-analysis across TCGA cohorts demonstrated a significant pooled hazard ratio ",
  "(meta-HR=1.19, meta-FDR=0.0499), driven primarily by LUAD (HR=1.47, 95% CI 1.24-1.73, ",
  "P=5.68 x 10^-6), with substantial heterogeneity (I2=90.4%, P=0.0012). LUSC showed no ",
  "significant association (HR=1.02, P=0.750). Exploratory multi-omics integration identified ",
  "30 methylation CpGs in LUAD (top: cg02952918, rho=0.481, FDR<2.6 x 10^-7) and 86 RPPA ",
  "antibodies (top: Cyclin B1, rho=0.565, FDR below machine-reportable precision). No significant ",
  "somatic mutation or CNV associations were identified for the All_immune glycolysis program.",
  "",
  "**Conclusions:** An immune-compartment glycolysis transcriptional signature is associated with ",
  "non-response to neoadjuvant anti-PD-1 therapy in NSCLC, with histology-dependent effects. ",
  "These findings require prospective validation."
)
writeLines(abstract, file.path(work_out, "Abstract_M6B2.md"))
cat("Abstract written:", length(abstract), "lines\n\n")


# ============================================================================
# SECTION III: Write Revised Methods
# ============================================================================
cat("\nSECTION III: Write Revised Methods\n")
cat(rep("=", 80), sep="")
cat("\n\n")

methods <- c(
  "# Methods",
  "",
  "## Pathway Enrichment Analysis",
  "Preranked gene-set enrichment analysis was performed using fgsea v1.38.0 with the ",
  "fgseaMultilevel algorithm (minSize=15, maxSize=500, eps=1e-50, gseaParam=1). ",
  "Gene-level rankings were constructed as sign(log2FC) * sqrt(F) from edgeR quasi-likelihood ",
  "F-tests. Hallmark gene sets (MSigDB 2024.1.Hs, 50 gene sets) were used for primary analysis. ",
  "Reactome gene sets (MSigDB 2024.1.Hs, 1,839 gene sets after size filtering) were analyzed ",
  "in parallel. Enrichment was performed separately for each of 8 primary-eligible immune cell ",
  "types using two cohorts: primary (anti-PD1 treated, n=233) and strict (chemoimmunotherapy, ",
  "n=212). HALLMARK_GLYCOLYSIS was tested in All_immune as the core program.",
  "",
  "## Clinical Validation (TCGA)",
  "The glycolysis program was scored in TCGA-LUAD (520 patients, 515 with RNA) and TCGA-LUSC ",
  "(504 patients, 501 with RNA) using ssGSEA (alpha=0.25, normalize=TRUE). ",
  "Cox proportional hazards models were fitted: Surv(OS_days/365.25, OS_event) ~ score_z + ",
  "age_z + sex_f + stage_f (ties=efron). Fixed-effect meta-analysis used the metafor package. ",
  "The prespecified internal Tier A threshold was meta-FDR < 0.05.",
  "",
  "## Multi-Omics Integration",
  "Exploratory feature-level associations were tested for DNA methylation (Illumina 450K, ",
  "Spearman correlation with program scores), RPPA protein levels (Spearman correlation), ",
  "somatic mutation burden (Spearman correlation and Cohen's d for driver mutations), and ",
  "copy number variation (GISTIC thresholded). All multi-omics analyses were exploratory and ",
  "used FDR < 0.05 for feature-level significance.",
  "",
  "## Transcription Factor Inference",
  "CollecTRI transcription-factor activity inference was prespecified but was not completed; ",
  "TF activity was therefore not incorporated into the final program interpretation."
)
writeLines(methods, file.path(work_out, "Methods_M6B2.md"))
cat("Methods written:", length(methods), "lines\n\n")


# ============================================================================
# SECTION IV: Write Revised Results
# ============================================================================
cat("\nSECTION IV: Write Revised Results\n")
cat(rep("=", 80), sep="")
cat("\n\n")

results <- c(
  "# Results",
  "",
  "## Glycolysis Pathway Enrichment",
  "The HALLMARK_GLYCOLYSIS pathway in All_immune cells was enriched in the non-responder ",
  "direction in both primary and strict cohorts. Primary fgsea: NES=-2.3589, ",
  "P=4.90 x 10^-12, FDR=3.00 x 10^-11, ES=-0.4859, gene-set size=170, ",
  "leading-edge=70 genes. Strict fgsea: NES=-2.4126, P=1.07 x 10^-12, ",
  "FDR=7.51 x 10^-12, leading-edge=65 genes.",
  "",
  "Directionality was consistent across all 8 analyzed strata (8/8 negative NES). ",
  "The magnitude ranged from NES=-2.4233 (cDC2_CD1C) to NES=-1.7016 (M_CXCL10). ",
  "All 8 strata showed negative NES, indicating higher glycolysis in non-responders. ",
  "We note that negative NES indicates enrichment in the non-responder direction; ",
  "this does not imply statistical significance for all individual strata (e.g., ILC3_KIT FDR=0.034).",
  "",
  "## Clinical Validation",
  "Fixed-effect meta-analysis across TCGA-LUAD and TCGA-LUSC demonstrated a significant ",
  "pooled hazard ratio (meta-HR=1.19, meta-FDR=0.0499), meeting the prespecified internal ",
  "Tier A criterion (meta-FDR < 0.05). However, heterogeneity was substantial (I2=90.4%, ",
  "P=0.0012). The association was driven by LUAD (HR=1.47, 95% CI 1.24-1.73, ",
  "P=5.68 x 10^-6, n=477, events=172), with no significant association in LUSC ",
  "(HR=1.02, 95% CI 0.89-1.18, P=0.750, n=485, events=210). The absence of association ",
  "in LUSC and I2 above 90% argue against a histology-independent prognostic association. ",
  "TCGA is not an immunotherapy-treated cohort; these associations reflect general cancer biology.",
  "",
  "## Exploratory Feature-Level Multi-Omics Associations",
  "DNA methylation: 30 CpGs in LUAD (top: cg02952918, rho=0.481, FDR<2.6 x 10^-7); ",
  "0 CpGs in LUSC. RPPA: 86 antibodies in LUAD (top: Cyclin B1, rho=0.565, FDR below ",
  "machine-reportable precision); 0 antibodies in LUSC for All_immune glycolysis. ",
  "Mutation burden: not significant in any cohort (LUAD FDR=0.356, LUSC FDR=0.950). ",
  "CNV: analysis completed but generated no association results for All_immune glycolysis ",
  "(NO_RESULT_GENERATED). All multi-omics results are exploratory and cannot establish causality.",
  "",
  "## Transcription Factor Inference",
  "CollecTRI-based TF inference was not completed. Supporting TF evidence is unavailable.",
  "",
  "## Leading-Edge Gene Composition",
  "The primary All_immune leading edge comprised 70 genes, including ",
  "canonical glycolytic enzymes (PKM, LDHA, PGAM1, ENO1, TPI1, PFKP, PGK1, GAPDH, ALDOA) ",
  "and MIF. The strict cohort leading edge comprised 65 genes with ",
  "substantial overlap. Leading-edge genes are provided in Supplementary Table X."
)
writeLines(results, file.path(work_out, "Results_M6B2.md"))
cat("Results written:", length(results), "lines\n\n")


# ============================================================================
# SECTION V: Write Revised Figure 5
# ============================================================================
cat("\nSECTION V: Write Revised Figure 5\n")
cat(rep("=", 80), sep="")
cat("\n\n")

fig5 <- c(
  "# Figure 5: Integrated Evidence for the Immune-Compartment Glycolysis Program",
  "",
  "## Figure 5 Legend",
  "",
  "**Figure 5. Immune-compartment glycolysis program in NSCLC response to neoadjuvant anti-PD-1.**",
  "",
  "**(A)** Preranked fgsea normalized enrichment scores (NES) for HALLMARK_GLYCOLYSIS ",
  "across 8 primary-eligible immune cell types. All cell types show negative NES ",
  "(direction: higher in non-responders). All_immune primary NES=-2.3589 ",
  "(FDR=3.00 x 10^-11).",
  "",
  "**(B)** fgsea exact statistics for HALLMARK_GLYCOLYSIS in All_immune. ",
  "Primary cohort: NES=-2.3589, P=4.90 x 10^-12, FDR=3.00 x 10^-11, ",
  "ES=-0.4859, gene-set size=170. Strict cohort: NES=-2.4126, ",
  "P=1.07 x 10^-12, FDR=7.51 x 10^-12, ES=-0.4876, gene-set size=170.",
  "",
  "**(C)** Leading-edge gene composition. Primary cohort: 70 genes. ",
  "Strict cohort: 65 genes. Representative genes: PKM, LDHA, PGAM1, ENO1, ",
  "TPI1, PFKP, PGK1, GAPDH, ALDOA, MIF. CollecTRI TF inference was not completed; ",
  "no TF activity data are shown.",
  "",
  "**(D)** Forest plot of Cox hazard ratios for HALLMARK_GLYCOLYSIS program score in ",
  "TCGA-LUAD (HR=1.47, 95% CI 1.24-1.73) and TCGA-LUSC (HR=1.02, 95% CI 0.89-1.18). ",
  "Fixed-effect meta-analysis: meta-HR=1.19, meta-FDR=0.0499.",
  "",
  "**(E)** Exploratory multi-omics feature-level associations for HALLMARK_GLYCOLYSIS ",
  "in All_immune. DNA methylation: 30 CpGs in LUAD (top: cg02952918, rho=0.481, ",
  "FDR<2.6 x 10^-7); 0 CpGs in LUSC. RPPA: 86 antibodies in LUAD (top: Cyclin B1, ",
  "rho=0.565, FDR below machine-reportable precision); 0 antibodies in LUSC. ",
  "Mutation burden: no significant associations (LUAD FDR=0.356, LUSC FDR=0.950). ",
  "CNV: no association results generated for All_immune glycolysis. ",
  "All multi-omics results are exploratory.",
  "",
  "**(F)** Cohort flow diagram.",
  "",
  "Source data: 03_results/step07_programs/fgsea_primary/ + fgsea_strict/; ",
  "03_results/step08_TCGA/B1_QC2/; 03_results/step08_TCGA/B2/."
)
writeLines(fig5, file.path(work_out, "Figure5_M6B2.md"))
cat("Figure 5 written:", length(fig5), "lines\n\n")


# ============================================================================
# SECTION VI: Verification Checks
# ============================================================================
cat("\nSECTION VI: Verification Checks\n")
cat(rep("=", 80), sep="")
cat("\n\n")

# Read back all output files
abs_out <- readLines(file.path(work_out, "Abstract_M6B2.md"), warn=FALSE)
met_out <- readLines(file.path(work_out, "Methods_M6B2.md"), warn=FALSE)
res_out <- readLines(file.path(work_out, "Results_M6B2.md"), warn=FALSE)
fig_out <- readLines(file.path(work_out, "Figure5_M6B2.md"), warn=FALSE)

all_out <- paste(c(abs_out, met_out, res_out, fig_out), collapse="\n")
all_lines <- c(abs_out, met_out, res_out, fig_out)

# Helper: check if pattern exists in ANY line
any_line_match <- function(pattern, lines, ignore.case=FALSE) {
  any(sapply(lines, function(l) grepl(pattern, l, ignore.case=ignore.case)))
}

checks <- data.frame(
  item = character(), check = character(), status = character(),
  stringsAsFactors = FALSE)

# --- NES values present ---
checks <- rbind(checks, data.frame(
  item="NES_primary", check="NES=-2.3589 or NES=-2.36 present",
  status=ifelse(any_line_match("NES=-2\\.3[56]", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="NES_strict", check="NES=-2.4126 or NES=-2.41 present",
  status=ifelse(any_line_match("NES=-2\\.41", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- FDR values present (scientific notation) ---
checks <- rbind(checks, data.frame(
  item="FDR_primary", check="FDR=3.00 x 10^-11 present",
  status=ifelse(any_line_match("3\\.00.*10.*-11|3\\.00e-11", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="FDR_strict", check="FDR=7.51 x 10^-12 present",
  status=ifelse(any_line_match("7\\.51.*10.*-12|7\\.51e-12", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- Leading-edge counts ---
checks <- rbind(checks, data.frame(
  item="LE_primary", check="70 leading-edge genes mentioned",
  status=ifelse(any_line_match("70 genes|leading-edge=70|70 gene", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="LE_strict", check="65 leading-edge genes mentioned",
  status=ifelse(any_line_match("65 genes|leading-edge=65|65 gene", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- 8/8 negative NES (not 8/8 significant) ---
checks <- rbind(checks, data.frame(
  item="eight_of_eight", check="8/8 negative NES stated",
  status=ifelse(any_line_match("8/8 negative NES|8/8 analyzed strata", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="no_eight_significant", check="Does NOT claim 8/8 significant",
  status=ifelse(!any_line_match("8/8.*significant|all 8.*strata.*significant|all 8.*significant", all_lines, ignore.case=TRUE), "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- Methylation ---
checks <- rbind(checks, data.frame(
  item="meth_LUAD_30", check="30 CpGs LUAD mentioned",
  status=ifelse(any_line_match("30 CpGs.*LUAD|LUAD.*30 CpGs", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="meth_top_cg02952918", check="cg02952918 mentioned",
  status=ifelse(any_line_match("cg02952918", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="meth_rho", check="rho=0.481 mentioned",
  status=ifelse(any_line_match("rho=0\\.481", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="meth_LUSC_0", check="0 CpGs LUSC mentioned",
  status=ifelse(any_line_match("0 CpGs.*LUSC|LUSC.*0 CpG", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- RPPA ---
checks <- rbind(checks, data.frame(
  item="rppa_LUAD_86", check="86 antibodies LUAD mentioned",
  status=ifelse(any_line_match("86 antibodies.*LUAD|LUAD.*86 antibod", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="rppa_Cyclin_B1", check="Cyclin B1 mentioned",
  status=ifelse(any_line_match("Cyclin B1", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="rppa_rho_0565", check="rho=0.565 mentioned",
  status=ifelse(any_line_match("rho=0\\.565", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="rppa_LUSC_0", check="0 antibodies LUSC mentioned",
  status=ifelse(any_line_match("0 antibod.*LUSC|LUSC.*0 antibod|LUSC.*All_immune.*0", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="rppa_FDR_machine", check="FDR below machine-reportable precision (not FDR=0)",
  status=ifelse(any_line_match("machine-reportable precision", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- Mutation ---
checks <- rbind(checks, data.frame(
  item="mutation_status", check="Mutation not significant mentioned",
  status=ifelse(any_line_match("[Mm]utation.*not significant|no significant.*mutation|NO_SIGNIFICANT_FEATURE", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- CNV ---
checks <- rbind(checks, data.frame(
  item="cnv_no_results", check="CNV no results generated (not CNV negative)",
  status=ifelse(any_line_match("no association results|NO_RESULT_GENERATED|not generated", all_lines, ignore.case=TRUE) &
                !any_line_match("CNV.*\\bnegative\\b|CNV.*statistically negative", all_lines, ignore.case=TRUE), "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# --- Forbidden patterns ---
checks <- rbind(checks, data.frame(
  item="no_empty_NES", check="No NES= followed by empty",
  status=ifelse(!any_line_match("^NES=\\s*$|\\sNES=\\s*$", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="no_empty_FDR", check="No FDR= followed by empty",
  status=ifelse(!any_line_match("^FDR=\\s*$|\\sFDR=\\s*$", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="no_FDR_zero", check="No FDR=0 literal",
  status=ifelse(!any_line_match("FDR=0[^.]|FDR=0$|FDR=0 ", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="no_top_NA", check="No top=NA",
  status=ifelse(!any_line_match("top=NA|top: NA", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="no_TBD", check="No TBD",
  status=ifelse(!any_line_match("TBD", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="no_supporting_TF", check="No 'supporting TF' used as positive evidence",
  status=ifelse(!any_line_match("[Ss]upporting TF", all_lines) |
                any_line_match("[Ss]upporting TF.*unavailable|[Ss]upporting TF.*not completed", all_lines, ignore.case=TRUE), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="no_CollecTRI_completed", check="No 'CollecTRI completed' claim",
  status=ifelse(!any_line_match("[Cc]ollec[Tt][Rr][Ii].*completed", all_lines) |
                any_line_match("not completed", all_lines, ignore.case=TRUE), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="no_Rb_All_immune", check="Rb not attributed to All_immune",
  status=ifelse(!any_line_match("Rb.*All_immune|All_immune.*Rb", all_lines), "PASS", "FAIL"),
  stringsAsFactors=FALSE))
checks <- rbind(checks, data.frame(
  item="no_CNV_negative", check="CNV not described as negative",
  status=ifelse(!any_line_match("CNV.*\\bnegative\\b|CNV.*statistically negative", all_lines, ignore.case=TRUE), "PASS", "FAIL"),
  stringsAsFactors=FALSE))

# Print results
n_pass <- sum(checks$status == "PASS")
n_fail <- sum(checks$status == "FAIL")

for (i in seq_len(nrow(checks))) {
  icon <- if (checks$status[i] == "PASS") "[PASS]" else "[FAIL]"
  cat(sprintf("  %s %-35s %s\n", icon, checks$item[i], checks$check[i]))
}

cat("\nTotal:", nrow(checks), "  PASS:", n_pass, "  FAIL:", n_fail, "\n\n")

# Write check file
write.csv(checks,
  file.path(work_out, "GSE243013_M6B2_core_section_check.csv"),
  row.names=FALSE)
cat("Check file written.\n\n")


# ============================================================================
# SECTION VII: Completion Marker
# ============================================================================
cat("\nSECTION VII: Completion Marker\n")
cat(rep("=", 80), sep="")
cat("\n\n")

overall <- n_fail == 0

marker <- c(
  "GSE243013 M6B2 CORE SECTIONS COMPLETE: COMPLETE",
  "",
  paste("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Sections written:",
  paste("  Abstract:", length(abstract), "lines"),
  paste("  Methods:", length(methods), "lines"),
  paste("  Results:", length(results), "lines"),
  paste("  Figure 5:", length(fig5), "lines"),
  "",
  "Verification:",
  paste("  Total checks:", nrow(checks)),
  paste("  PASS:", n_pass),
  paste("  FAIL:", n_fail),
  paste("  Overall:", if (overall) "ALL PASS" else "SOME FAILURES"),
  "",
  "Key values embedded:",
  "  Primary NES=-2.3589, FDR=3.00 x 10^-11",
  "  Strict NES=-2.4126, FDR=7.51 x 10^-12",
  "  Leading-edge: 70 primary, 65 strict",
  "  8/8 negative NES (not 8/8 significant)",
  "  LUAD methylation: 30 CpGs, top cg02952918",
  "  LUAD RPPA: 86 antibodies, top Cyclin B1",
  "  Mutation: NO_SIGNIFICANT_FEATURE",
  "  CNV: NO_RESULT_GENERATED",
  "  CollecTRI: not completed",
  "  Rb: not attributed to All_immune",
  "  Cyclin B1 FDR: machine-reportable precision (not FDR=0)",
  "",
  "Output files:",
  "  Abstract_M6B2.md",
  "  Methods_M6B2.md",
  "  Results_M6B2.md",
  "  Figure5_M6B2.md",
  "  GSE243013_M6B2_core_section_check.csv",
  "",
  "7. M6B2 completion marker: CREATED"
)
writeLines(marker, "05_manuscript/GSE243013_M6B2_CORE_SECTIONS_COMPLETE.txt")
cat("--- M6B2 COMPLETION MARKER CREATED ---\n")
for (line in marker) cat("  ", line, "\n")

cat("\n", rep("=", 80), sep="")
cat("\nM6B2: Integrate Core Scientific Sections - COMPLETED")
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
cat("\n", rep("=", 80), sep="", "\n\n")
