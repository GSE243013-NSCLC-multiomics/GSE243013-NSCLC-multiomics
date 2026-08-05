# An Immune-Compartment Glycolysis-Related Program Associated with Pathological Nonresponse after Neoadjuvant Anti-PD-1-Based Treatment in NSCLC

**Authors:** [AUTHOR LIST]
**Affiliations:** [AFFILIATIONS]
**Corresponding author:** [CORRESPONDING AUTHOR]

**Keywords:** NSCLC, neoadjuvant immunotherapy, glycolysis,
  single-cell RNA sequencing, patient-level pseudobulk, TCGA


# Abstract

**Background:** Neoadjuvant anti-PD-1-based immunotherapy achieves pathological complete response (pCR) in a subset of non-small cell lung cancer (NSCLC) patients, but the immune microenvironment determinants of response remain poorly understood.

**Methods:** We performed single-cell RNA sequencing on 1,254,749 cells from 243 post-neoadjuvant surgical specimens, including 234 with anti-PD-1 treatment records and 9 chemotherapy-only controls. The primary response analysis included 233 anti-PD-1-treated patients with evaluable binary pathological response, analyzed across 8 immune cell types using patient-level pseudobulk aggregation. Tumor-infiltrating immune cell transcriptomes were profiled using edgeR differential expression and preranked fgsea pathway enrichment (Hallmark gene sets, MSigDB 2024.1.Hs). The core immune-compartment glycolysis program was assessed for association with clinical outcomes in TCGA-LUAD (n=477) and TCGA-LUSC (n=485) using Cox proportional hazards models, and integrated with DNA methylation, RPPA proteomics, somatic mutation, and copy number variation data.

**Results:** The HALLMARK_GLYCOLYSIS pathway in All_immune cells was enriched in the non-responder direction (primary NES=-2.3589, FDR=3.00e-11; strict NES=-2.4126, FDR=7.51e-12). Negative NES was observed in all eight analyzed strata, and all eight met FDR<0.05. Fixed-effect meta-analysis across TCGA cohorts demonstrated a pooled hazard ratio (meta-HR=1.19, meta-FDR=0.0499), driven primarily by LUAD (HR=1.47, 95% CI 1.24-1.73, P=5.68 x 10^-6), with substantial heterogeneity (I2=90.4%, P=0.0012). LUSC showed no significant association (HR=1.02, P=0.750). Exploratory multi-omics integration identified 30 methylation CpGs in LUAD (top: cg02952918, rho=0.481, FDR<2.6 x 10^-7) and 86 RPPA antibodies (top: Cyclin B1, rho=0.565, FDR below machine-reportable precision). No significant somatic mutation associations were identified. CNV processing was completed but no program-level association result was generated for the All_immune glycolysis program.

**Conclusions:** An immune-compartment glycolysis transcriptional program is associated with the non-MPR pathological-response group in post-treatment surgical specimens, with histology-dependent effects. These findings require prospective evaluation in pretreatment immunotherapy cohorts.

# Introduction

Neoadjuvant anti-PD-1-based immunotherapy has emerged as a standard treatment approach for resectable non-small cell lung cancer (NSCLC). The CheckMate 816 trial demonstrated that nivolumab combined with platinum-based chemotherapy improved pathological complete response (pCR) rates and event-free survival compared with chemotherapy alone [CITATION NEEDED: CheckMate 816]. Similar benefits have been reported with pembrolizumab-based regimens [CITATION NEEDED: KEYNOTE-671]. Despite these advances, pathological response is markedly heterogeneous: some patients achieve pCR with no residual viable tumor, others achieve major pathological response (MPR) with limited residual disease, and a substantial proportion show non-MPR status with extensive residual tumor [CITATION NEEDED: pathological response definitions and MPR criteria]. Understanding the biological determinants of this heterogeneity is a priority for improving patient selection and developing rational combination strategies.

The tumor immune microenvironment plays a central role in determining response to immune checkpoint blockade. Pre-treatment immune cell composition, activation state, and spatial organization have been associated with immunotherapy outcomes across multiple tumor types, including melanoma, renal cell carcinoma, and NSCLC [CITATION NEEDED: TIM review across tumor types]. In NSCLC specifically, tumor-infiltrating lymphocyte density and PD-L1 expression have shown inconsistent predictive value, suggesting that more granular immune cell characterization may be needed [CITATION NEEDED: TIL and PD-L1 limitations in NSCLC]. Single-cell RNA sequencing (scRNA-seq) enables high-resolution characterization of immune cell states and transcriptional programs within the tumor microenvironment, providing insights that are not accessible through bulk transcriptomic approaches [CITATION NEEDED: scRNA-seq in NSCLC tumor microenvironment]. However, most scRNA-seq studies of immunotherapy response have analyzed post-treatment surgical specimens, limiting the ability to identify pre-treatment predictors of response.

A critical methodological consideration in scRNA-seq studies is the pseudoreplication problem. When individual cells are treated as independent biological replicates in group comparisons (e.g., responder vs. non-responder), the effective sample size is artificially inflated, because cells from the same patient are biologically correlated. This inflation can lead to overly narrow confidence intervals and elevated false-positive rates [CITATION NEEDED: pseudoreplication in scRNA-seq studies]. The magnitude of this bias can be substantial: with hundreds or thousands of cells per patient, treating cells as independent replicates can inflate statistical power by orders of magnitude relative to the actual number of biological replicates. Patient-level pseudobulk aggregation addresses this limitation by summarizing cell-level expression to the patient level, preserving biological replicability while retaining cell-type-specific information [CITATION NEEDED: pseudobulk methodology papers]. This approach is increasingly recognized as the appropriate analytical framework for comparing transcriptomic profiles across clinical groups in scRNA-seq datasets, and is recommended by community guidelines for single-cell analysis [CITATION NEEDED: scRNA-seq best practices guidelines].

Metabolic reprogramming is a hallmark of both tumor cells and activated immune cells. Glycolysis supports the effector functions of T cells and myeloid cells through rapid ATP production and biosynthetic precursor generation [CITATION NEEDED: immune cell glycolysis]. In the tumor microenvironment, chronic metabolic stress including hypoxia, nutrient deprivation, and acidosis can promote T cell exhaustion and immunosuppressive myeloid cell polarization, creating a metabolic landscape that favors immune evasion [CITATION NEEDED: tumor metabolic microenvironment and immune evasion]. Whether immune-cell glycolytic programs are associated with immunotherapy response in NSCLC, and whether such associations are histology-specific, remains largely unexplored. Most studies of tumor glycolysis have focused on cancer cell-intrinsic metabolism rather than immune-cell metabolic states, leaving a gap in understanding how immune metabolic reprogramming relates to treatment outcomes.

The Cancer Genome Atlas (TCGA) provides multi-omic profiling data for large NSCLC cohorts, including lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LUSC) histological subtypes [CITATION NEEDED: TCGA NSCLC landmark paper]. TCGA cohorts can be used to assess whether transcriptional programs identified in immunotherapy cohorts are associated with cancer-relevant clinical outcomes such as overall survival. However, TCGA is not an immunotherapy-treated cohort; patients received surgery, chemotherapy, radiation, or combinations thereof, but not immune checkpoint inhibitors. Survival associations in TCGA therefore reflect general cancer biology and cannot validate immunotherapy-specific response mechanisms [CITATION NEEDED: TCGA design and treatment limitations]. Furthermore, the substantial molecular differences between LUAD and LUSC histological subtypes suggest that immune-related transcriptional programs may show histology-specific patterns that require separate evaluation.

In this study, we applied patient-level pseudobulk analysis to scRNA-seq data from 243 post-neoadjuvant surgical specimens to identify immune transcriptional programs associated with pathological response in NSCLC patients receiving neoadjuvant anti-PD-1-based therapy. We assessed the core glycolysis program for external survival associations in TCGA-LUAD and TCGA-LUSC, and performed exploratory multi-omics integration with DNA methylation, RPPA proteomics, somatic mutation, and copy number variation data.


# Methods

## Study Design and Public Datasets
This study analyzed publicly available single-cell RNA sequencing data from the GEO accession GSE243013. The dataset comprised post-neoadjuvant surgical specimens from 243 NSCLC patients treated with neoadjuvant anti-PD-1-based regimens or chemotherapy alone.

## GEO Supplementary-File Acquisition
Raw count matrices and cell metadata were downloaded from GEO supplementary files. Cell-level annotations, patient identifiers, and clinical metadata were extracted from the provided metadata. No FASTQ files were downloaded; all analyses were performed on the provided count matrices.

## Cohort Definitions
The full dataset included 243 patients. Of these, 234 had documented anti-PD-1 treatment records and 9 received chemotherapy without anti-PD-1 therapy (chemotherapy-only controls). The primary response analysis included 233 anti-PD-1-treated patients with evaluable binary pathological response (responder vs. non-responder). A strict chemoimmunotherapy exposure subgroup included 213 patients with documented concurrent anti-PD-1 and chemotherapy; of these, 212 had evaluable binary response (113 responders, 99 non-responders).

## Pathological-Response Definitions
Pathological response was classified as binary: responder (pCR or MPR) versus non-responder (non-MPR), based on the provided clinical metadata. Pathological response rates were extracted as numeric values when available.

## Count-Matrix Validation
The count matrix was validated for dimension consistency (1,254,749 cells x 31,831 genes), barcodes, gene names, non-negativity, integer values, and non-zero content.

## BPCells Processing
Count matrices were converted to BPCells format for memory-efficient storage and access. Data were organized in column-major orientation (cells as columns, genes as rows) for efficient column slicing by patient.

## Patient-Level Pseudobulk Aggregation
For each patient and cell type, cell-level counts were aggregated to patient-level pseudobulk profiles by summing expression across all cells of a given type within each patient. This preserves biological replicability at the patient level while retaining cell-type-specific information.

## Cell-Type Eligibility
Cell types were eligible for primary analysis if they had sufficient cells across both responder and non-responder groups to support patient-level differential expression testing. A total of 8 immune cell types met eligibility criteria for the primary analysis.

## edgeR Filtering, Normalization, and Model
Patient-level pseudobulk counts were analyzed using edgeR v4.10.1. Genes with low expression were filtered using the filterByExpr function with default settings. Normalization factors were calculated using the TMM method (normLibSizes). A quasi-likelihood negative binomial generalized linear model was fitted with the design formula: ~ cancer_type + response_binary, with response_binaryResponder as the target coefficient. The contrast was set using makeContrasts for the response comparison.

## QL and glmTreat Testing
Quasi-likelihood F-tests were performed using glmQLFTest. A treat-like filtering was applied using glmTreat with a log2-fold-change threshold of log2(1.2), testing for biologically meaningful differences rather than merely statistical significance.

## Multiple-Testing Correction
False discovery rates were estimated using the Benjamini-Hochberg method within edgeR. For pathway-level analysis, fgsea reported nominal P values and FDR-adjusted P values for each gene set.

## Preranked fgsea
Gene-set enrichment analysis was performed using fgsea v1.38.0 with the fgseaMultilevel algorithm (minSize=15, maxSize=500, eps=1e-50, gseaParam=1). Gene-level rankings were constructed as sign(log2FC) x sqrt(F) from edgeR quasi-likelihood F-tests, combining effect direction and statistical evidence.

## Hallmark and Reactome Resources
Hallmark gene sets (MSigDB 2024.1.Hs, 50 gene sets) were used for primary analysis. Reactome gene sets (MSigDB 2024.1.Hs, 1,839 gene sets after size filtering) were analyzed in parallel. Enrichment was performed separately for each of 8 primary-eligible immune cell types using two cohorts: primary (anti-PD1 treated, n=233) and strict (chemoimmunotherapy, n=212).

## Leading-Edge Extraction
Leading-edge gene sets were extracted from fgsea results for significant gene sets. The leading edge comprises the subset of genes contributing most strongly to the enrichment score.

## CollecTRI Transcription Factor Inference
CollecTRI-based transcription factor activity inference was prespecified but was not completed due to software dependency issues. TF activity data are therefore unavailable and are not reported.

## TCGA Data Acquisition
TCGA-LUAD and TCGA-LUSC RNA-seq, clinical, methylation, RPPA, mutation, and copy number data were obtained from the Genomic Data Commons. TCGA-LUAD included 520 patients (515 with RNA-seq); TCGA-LUSC included 504 patients (501 with RNA-seq).

## ssGSEA and GSVA Scoring
The glycolysis program was scored in TCGA cohorts using single-sample Gene Set Enrichment Analysis (ssGSEA, alpha=0.25, normalize=TRUE) as the primary method. Gene Set Variation Analysis (GSVA, Gaussian kernel, tau=1) was used as a sensitivity analysis. Scores were z-scored within each cohort-program combination.

## Cox Models and Covariates
Cox proportional hazards models were fitted: Surv(OS_days/365.25, OS_event) ~ score_z + age_z + sex_f + stage_f, with ties handled using the Efron method. Models were fitted separately for TCGA-LUAD and TCGA-LUSC.

## Fixed-Effect Meta-Analysis and Heterogeneity
Fixed-effect meta-analysis was performed using the metafor package to pool hazard ratios across TCGA-LUAD and TCGA-LUSC. The I2 statistic and Cochran Q test were used to assess between-study heterogeneity. With I2 approximately 90%, the fixed-effect pooled estimate has descriptive rather than inferential significance, as it assumes a common true effect across histologies.

## Exploratory Multi-Omics Procedures
Exploratory feature-level associations were tested for DNA methylation (Illumina 450K, Spearman correlation with program scores), RPPA protein levels (Spearman correlation), somatic mutation burden (Spearman correlation and Cohen's d for individual driver mutations), and copy number variation (GISTIC thresholded). All multi-omics analyses were exploratory and used FDR<0.05 for feature-level significance.

## Internal Integrated-Evidence Criterion
A prespecified internal Tier A criterion was defined as meta-FDR<0.05 across TCGA cohorts. This criterion was met for the HALLMARK_GLYCOLYSIS program (meta-FDR=0.0499).

## Software and Reproducibility
All analyses were performed in R 4.6.1 on macOS Sonoma (ARM64). Key packages: edgeR v4.10.1, fgsea v1.38.0, metafor, BPCells. All scripts and parameters are provided in the supplementary code repository.

## Ethics Statement
[AUTHOR INPUT REQUIRED: institutional ethics statement]

# Results

## Cohort and Matrix Overview
The dataset included post-neoadjuvant surgical specimens from 243 NSCLC patients. Of these, 234 had anti-PD-1 treatment records and 9 received chemotherapy alone. The count matrix comprised 1,254,749 cells and 31,831 genes.

## Patient-Level Pseudobulk and Model Eligibility
The primary response analysis included 233 anti-PD-1-treated patients with evaluable binary pathological response. A strict chemoimmunotherapy exposure subgroup included 213 patients, of whom 212 had evaluable response (113 responders, 99 non-responders). Patient-level pseudobulk profiles were generated for 8 primary-eligible immune cell types.

## Cell-Type Differential Expression
edgeR quasi-likelihood tests identified differentially expressed genes between responders and non-responders for each cell type. Gene-level rankings were constructed as sign(log2FC) x sqrt(F) for downstream pathway enrichment.

## Glycolysis Enrichment in Primary and Strict Analyses
The HALLMARK_GLYCOLYSIS pathway in All_immune cells was enriched in the non-responder direction in both primary (NES=-2.3589, P=4.90e-12, FDR=3.00e-11, ES=-0.4859, gene-set size=170, leading-edge=70 genes) and strict (NES=-2.4126, P=1.07e-12, FDR=7.51e-12, ES=-0.4876, leading-edge=65 genes) cohorts.

## Eight-Strata Direction and Significance
Negative NES was observed in all 8 analyzed strata (8/8 negative NES). All 8 of 8 strata met FDR<0.05.

## Leading-Edge Genes
The primary All_immune leading edge comprised 70 genes, including canonical glycolytic enzymes (PKM, LDHA, PGAM1, ENO1, TPI1, PFKP, PGK1, GAPDH, ALDOA) and MIF. The strict cohort leading edge comprised 65 genes with substantial overlap.

## LUAD and LUSC Survival Associations
Fixed-effect meta-analysis across TCGA-LUAD and TCGA-LUSC demonstrated a pooled hazard ratio (meta-HR=1.19, meta-FDR=0.0499), meeting the prespecified Tier A criterion. The association was driven by LUAD (HR=1.47, 95% CI 1.24-1.73, P=5.68 x 10^-6, n=477, events=172), with no significant association in LUSC (HR=1.02, 95% CI 0.89-1.18, P=0.750, n=485, events=210).

## Fixed-Effect Summary and High Heterogeneity
Heterogeneity between TCGA-LUAD and TCGA-LUSC was substantial (I2=90.4%, P=0.0012). The fixed-effect pooled estimate therefore has descriptive rather than inferential significance, as it assumes a common true effect across histologies, which is unsupported by the histology-specific pattern. TCGA is not an immunotherapy-treated cohort; these associations reflect general cancer biology.

## Exploratory Methylation and RPPA
DNA methylation: 30 CpGs in LUAD (top: cg02952918, rho=0.481, FDR<2.6 x 10^-7); 0 CpGs in LUSC. RPPA: 86 antibodies in LUAD (top: Cyclin B1, rho=0.565, FDR below machine-reportable precision); 0 antibodies in LUSC for All_immune glycolysis.

## Mutation: No Significant Feature
Somatic mutation burden was not associated with glycolysis program scores in any cohort (LUAD FDR=0.356, LUSC FDR=0.950).

## CNV: No Program-Level Result Generated
CNV processing was completed, but no program-level CNV association result was generated for the All_immune glycolysis program.

## CollecTRI Not Completed
CollecTRI-based TF inference was not completed. TF activity data are therefore unavailable.

# Discussion

## Principal Finding
This study identifies an immune-compartment glycolysis-related transcriptional program associated with the non-MPR pathological-response group in post-treatment surgical specimens from NSCLC patients receiving neoadjuvant anti-PD-1-based therapy. The HALLMARK_GLYCOLYSIS pathway in All_immune cells was enriched in the non-responder direction in both primary (NES=-2.3589, FDR=3.00e-11) and strict (NES=-2.4126, FDR=7.51e-12) cohorts. Negative NES was observed in all 8 analyzed strata, and all 8 met FDR<0.05.

## Biological Interpretation
The enrichment of glycolysis-related genes in immune cells of non-responders is consistent with the metabolic reprogramming observed in activated T cells and myeloid cells. However, the All_immune aggregation precludes cell-subtype attribution: we cannot determine whether the signal reflects tumor-infiltrating lymphocyte metabolism, myeloid cell glycolysis, or compositional shifts in the immune microenvironment. Transcriptomic enrichment does not establish metabolic flux.

## Relationship to Existing Literature
Previous studies have reported associations between tumor glycolysis and immunotherapy resistance. Our findings extend this to the immune compartment specifically, though the All_immune aggregation limits mechanistic interpretation. The directionality (higher glycolysis in non-responders) is consistent with an immunosuppressive metabolic microenvironment, but causality cannot be inferred from this observational study.

## TCGA External Survival Associations
Fixed-effect meta-analysis across TCGA-LUAD and TCGA-LUSC demonstrated a pooled hazard ratio (meta-HR=1.19, meta-FDR=0.0499). However, the I2 statistic was approximately 90%, indicating substantial between-histology heterogeneity. The fixed-effect summary therefore has descriptive rather than inferential value: the pooled estimate assumes a common true effect across histologies, which is unsupported given the LUAD-specific pattern. The association was driven by LUAD (HR=1.47, 95% CI 1.24-1.73, P=5.68 x 10^-6), with no significant association in LUSC (HR=1.02, P=0.750). TCGA is not an immunotherapy-treated cohort; these associations reflect general cancer biology and cannot establish the neoadjuvant immunotherapy response association.

## Multi-Omics Integration
Exploratory feature-level multi-omics integration revealed histology-specific association patterns. In LUAD, 30 methylation CpGs were correlated with glycolysis program scores (top: cg02952918, rho=0.481, FDR<2.6 x 10^-7), and 86 RPPA antibodies showed correlation (top: Cyclin B1, rho=0.565, FDR below machine-reportable precision). Neither finding was recapitulated in LUSC (0 methylation CpGs, 0 RPPA antibodies for All_immune glycolysis at FDR<0.05). Cyclin B1 is a proliferation marker; the correlation may reflect tumor cell proliferation, purity, or compositional confounding rather than a direct glycolysis-proteome relationship. The cg02952918 CpG requires independent annotation and functional follow-up before biological interpretation. Somatic mutation burden was not associated with glycolysis program scores in any cohort (LUAD FDR=0.356, LUSC FDR=0.950). CNV processing was completed but no program-level association result was generated for the All_immune glycolysis program.

## Transcription Factor Inference
CollecTRI-based transcription factor inference was prespecified but was not completed. TF activity data are therefore unavailable and are not presented. We cannot propose TF-mediated mechanisms based on these data.

## Strengths
Key strengths include the large single-cell cohort (243 patients, 1,254,749 cells), the prespecified patient-level pseudobulk framework, separate assessment in TCGA-LUAD and TCGA-LUSC, explicit reporting of between-histology heterogeneity, and exploratory multi-omics integration. The consistent directionality across 8 cell types and 2 cohorts supports the robustness of the direction finding.

## Limitations
Important limitations include: (1) the All_immune aggregation prevents cell-subtype attribution; (2) compositional shifts in the immune microenvironment may confound the glycolysis signal; (3) CollecTRI TF inference was not completed; (4) the LUAD-specific multi-omics findings were not recapitulated in LUSC; (5) TCGA is not an immunotherapy cohort; (6) the I2 of approximately 90% limits the interpretability of the meta-analytic summary; (7) post-treatment samples cannot establish whether glycolysis preceded or resulted from treatment; and (8) the cg02952918 CpG requires functional validation.

# Conclusion

An immune-compartment glycolysis transcriptional program is associated with the non-MPR pathological-response group in post-treatment surgical specimens from NSCLC patients receiving neoadjuvant anti-PD-1-based therapy, with consistent directionality across immune cell strata and cohorts. The association with pathological response in the primary cohort and the TCGA-LUAD prognostic association provide converging evidence, though substantial heterogeneity (I2 approximately 90%) and histology-specific effects limit generalizability. CollecTRI transcription factor inference was not completed. Exploratory multi-omics integration identified methylation and RPPA associations in LUAD that were not recapitulated in LUSC. These findings require prospective evaluation in pretreatment immunotherapy cohorts with cell-type-resolved profiling.

# Clinical Relevance

This post-treatment, response-associated transcriptional program provides a hypothesis for evaluation in independent pretreatment immunotherapy cohorts. The current evidence does not support clinical prediction or biomarker use. Clinical application would require pretreatment validation and cell-type-resolved quantification.

# Translational Relevance

The findings nominate an aggregate immune-compartment glycolysis-related transcriptional state for cell-type-resolved, spatial, metabolic, and functional investigation. They do not establish a treatment target, predictive assay, or causal metabolic mechanism. The histology-specific pattern (LUAD associations not recapitulated in LUSC) suggests that any future translational application would need to account for histological heterogeneity.

# Figure 5: Integrated Evidence for the Immune-Compartment Glycolysis Program

## Figure 5 Legend

**Figure 5. Immune-compartment glycolysis program in NSCLC response to neoadjuvant anti-PD-1.**

**(A)** Preranked fgsea normalized enrichment scores (NES) for HALLMARK_GLYCOLYSIS across 8 primary-eligible immune cell types. All cell types show negative NES (direction: higher in non-responders). All_immune primary NES=-2.3589 (FDR=3.00e-11).

**(B)** fgsea exact statistics for HALLMARK_GLYCOLYSIS in All_immune. Primary cohort: NES=-2.3589, P=4.90e-12, FDR=3.00e-11, ES=-0.4859, gene-set size=170. Strict cohort: NES=-2.4126, P=1.07e-12, FDR=7.51e-12, ES=-0.4876, gene-set size=170.

**(C)** Leading-edge gene composition. Primary cohort: 70 genes. Strict cohort: 65 genes. Representative genes: PKM, LDHA, PGAM1, ENO1, TPI1, PFKP, PGK1, GAPDH, ALDOA, MIF. CollecTRI TF inference was not completed; no TF activity data are shown.

**(D)** Forest plot of Cox hazard ratios for HALLMARK_GLYCOLYSIS program score in TCGA-LUAD (HR=1.47, 95% CI 1.24-1.73) and TCGA-LUSC (HR=1.02, 95% CI 0.89-1.18). Fixed-effect meta-analysis: meta-HR=1.19, meta-FDR=0.0499, I2=90.4%, heterogeneity P=0.0012. The I2 of approximately 90% indicates substantial between-histology heterogeneity; the fixed-effect summary has descriptive rather than inferential significance.

**(E)** Exploratory multi-omics feature-level associations for HALLMARK_GLYCOLYSIS in All_immune. DNA methylation: 30 CpGs in LUAD (top: cg02952918, rho=0.481, FDR<2.6 x 10^-7); 0 CpGs in LUSC. RPPA: 86 antibodies in LUAD (top: Cyclin B1, rho=0.565, FDR below machine-reportable precision); 0 antibodies in LUSC. Mutation burden: no significant associations (LUAD FDR=0.356, LUSC FDR=0.950). CNV processing was completed but no program-level association result was generated for All_immune glycolysis. All multi-omics results are exploratory.

**(F)** Cohort flow diagram.

# Figure Legends

## Figure 1: Cohort Overview and Single-Cell Transcriptomic Profiling
**Figure 1.** Cohort overview and single-cell transcriptomic profiling of post-neoadjuvant NSCLC surgical specimens. **(A)** Study design and sample collection timeline. **(B)** Patient cohort composition: 243 total, 234 with anti-PD-1 treatment records, 233 primary-eligible for response analysis. **(C)** Cell-type annotation and distribution across samples. **(D)** Quality control metrics. Pathway enrichment for Figures 1-3 used preranked fgsea (fgseaMultilevel, MSigDB 2024.1.Hs Hallmark gene sets), not ssGSEA.

## Figure 2: Cell-Type-Specific Differential Expression and Pathway Enrichment
**Figure 2.** Cell-type-specific differential expression and preranked fgsea pathway enrichment. 
**(A)** EdgeR differential expression summary across 47 cell types. 
**(B)** Hallmark pathway enrichment heatmap. 
**(C)** Reactome pathway enrichment heatmap. 
Enrichment was performed using preranked fgsea (fgseaMultilevel) with gene-level rankings constructed as sign(log2FC) x sqrt(F) from edgeR quasi-likelihood F-tests.

## Figure 3: Glycolysis Program Across Cell Types and Cohorts
**Figure 3.** Glycolysis program enrichment across cell types and cohorts. 
**(A)** Primary cohort fgsea NES for HALLMARK_GLYCOLYSIS across 8 cell types. 
**(B)** Strict cohort fgsea NES. 
**(C)** Direction concordance between primary and strict cohorts. 
All cell types show negative NES (higher in non-responders). 
Preranked fgsea was used for all pathway enrichment analyses.

## Figure 4: Program Prioritization and External Clinical Association
**Figure 4.** Program prioritization and external clinical association. 
**(A)** Preranked fgsea NES for HALLMARK_GLYCOLYSIS in All_immune cells. 
**(B)** Leading-edge gene overlap between primary (70 genes) and strict (65 genes) cohorts. 
**(C)** Cohort direction concordance across primary and strict cohorts. 
**(D)** Program prioritization: Tier 2 pathway evidence, external TCGA clinical association, 
and multi-omics integration status. CollecTRI TF inference was not completed.

## Figure 6: TCGA Histology-Specific Survival Associations
**Figure 6.** TCGA histology-specific survival associations for the glycolysis program. **(A)** KM survival curves by program score median split (visualization only; no statistical test). **(B)** Cox hazard ratios for TCGA-LUAD and TCGA-LUSC. **(C)** Forest plot of fixed-effect meta-analysis (meta-HR=1.19, meta-FDR=0.0499, I2=90.4%, heterogeneity P=0.0012). The I2 of approximately 90% indicates substantial between-histology heterogeneity; the fixed-effect summary has descriptive rather than inferential significance.

## Figure 7: Exploratory Multi-Omics Feature-Level Associations
**Figure 7.** Exploratory multi-omics feature-level associations for the 
core glycolysis-associated program (HALLMARK_GLYCOLYSIS) in All_immune cells. 
**(A)** DNA methylation: 30 CpGs in LUAD (top: cg02952918), 0 in LUSC. 
**(B)** RPPA: 86 antibodies in LUAD (top: Cyclin B1, rho=0.565), 0 in LUSC. 
**(C)** Mutation burden: not significant in any cohort. 
**(D)** CNV processing was completed but no program-level association result was generated 
for All_immune glycolysis. All multi-omics results are exploratory and cannot establish causality.

