# Data Dictionary

## GSE243013 scRNA-seq Data

### Cell Metadata Columns

| Column | Description | Values |
|--------|-------------|--------|
| sampleID | Patient identifier | e.g., "P001" |
| barcode | Cell barcode | e.g., "AAACGAAAGTGGTATG-1" |
| cancer_type | Histological subtype | LUAD, LUSC |
| response | Pathological response (3-level) | pCR, MPR, non-MPR, unknown |
| response_binary | Binary response | Responder, Non_responder |
| treatment | Treatment regimen | Anti-PD1, Chemo-immunotherapy, Chemo-only |
| treatment_class | Treatment classification | anti_PD1, chemo_immunotherapy, chemotherapy |
| cell_type_major | Major cell type | T_NK_cell, B_cell, Myeloid_cell |
| cell_type_sub | Sub cell type | e.g., CD8T_Trm_ZNF683, cDC2_CD1C |

### Pseudobulk Matrix Format

- **Dimensions:** genes × patients
- **Values:** Raw UMI counts (summed per patient)
- **Orientation:** Transposed from original cells × genes BPCells format
- **File format:** CSV or RDS

## TCGA Data

### Clinical Variables

| Variable | Description | Type |
|----------|-------------|------|
| OS_days | Overall survival (days) | Numeric |
| OS_event | Death event (0=censored, 1=dead) | Binary |
| age | Age at diagnosis | Numeric |
| sex | Patient sex | Male, Female |
| pathologic_stage | Pathological stage | Stage I, II, III, IV |

### Program Scoring

| Method | Parameters | Use |
|--------|-----------|-----|
| ssGSEA | alpha=0.25, normalized=TRUE | Primary scoring |
| GSVA | tau=1, Gaussian kernel, maxDiff=TRUE | Sensitivity |
| Z-score | Within cohort × program | Normalization |

### Clinical Model Formula

```
Surv(OS_days/365.25, OS_event) ~ score_z + age_z + sex_f + stage_f
```

- **ties:** efron
- **Score:** Continuous z-scored program score
- **HR interpretation:** Risk per 1 SD increase in score

## Evidence Tiers

| Tier | Definition | Criteria |
|------|-----------|----------|
| Level A | Core program | meta-FDR < 0.05 AND direction consistent |
| Level B | Supportive program | FDR < 0.05 in ≥1 histology, meta-FDR ≥ 0.05 |
| Level C | Exploratory | No FDR < 0.05 in either histology |

## Key Result Tables

### Table 1: Cohort Characteristics
- Patient count by cancer type, response, treatment
- Cell count summary

### Table 3: Primary edgeR Results
- Gene-level DE statistics per cell type
- Includes logFC, F-statistic, P-value, FDR

### Table 5: TCGA Clinical Validation
- Cox model results per program per histology
- Includes HR, CI, P-value, FDR

### Supplementary Table S6: Meta-Analysis
- Fixed-effect meta-analysis across LUAD and LUSC
- Includes meta-HR, meta-FDR, heterogeneity statistics
