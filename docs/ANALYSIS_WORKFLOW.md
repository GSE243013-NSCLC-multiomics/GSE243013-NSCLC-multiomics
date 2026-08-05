# Analysis Workflow

## Overview

This document describes the computational workflow for the GSE243013 NSCLC multi-omics analysis. The pipeline processes single-cell RNA-seq data from 243 NSCLC patients through 16 analysis steps.

## Pipeline Architecture

```
Step 00: Environment Check
    ↓
Step 01: Metadata Download & Inspection
    ↓
Step 02: Patient Manifest Construction
    ↓
Step 03: Cohort Definition Repair
    ↓
Step 04: Count Matrix Import (BPCells)
    ↓
Step 05: Patient-Level Pseudobulk
    ↓
Step 06: edgeR Differential Expression
    ↓
Step 07: Pathway Enrichment (fgsea)
    ↓
Step 08A: TCGA Data Download
    ↓
Step 08B1: TCGA Clinical Modeling (Cox)
    ↓
Step 08B2: TCGA Multi-Omics Integration
    ↓
Step 09: Final Audit & Evidence Report
    ↓
Step 09A: Gap Closure & Revision
```

## Step Details

### Step 00: Environment Check
- Verifies R version, packages, and disk space
- Creates project directory structure
- Logs sessionInfo()

### Step 01: Metadata Download
- Downloads GSE243013 metadata from GEO
- Profiles metadata columns
- Identifies patient-level and cell-level annotations

### Step 02: Patient Manifest
- Constructs patient-level manifest (243 patients)
- Cross-validates response, treatment, and cancer type
- Identifies conflicts and unresolved records

### Step 03: Cohort Repair
- Fixes treatment classification logic
- Defines primary (anti-PD1), strict chemoimmunotherapy, and chemo-only cohorts
- Excludes patients with conflicting annotations

### Step 04: BPCells Import
- Downloads count matrix from GEO
- Converts to BPCells column-major format
- Validates matrix dimensions (1,254,749 × 31,831)

### Step 05: Pseudobulk
- Aggregates cell-level counts to patient-level pseudobulk
- Separates by cell type (All_immune, major, sub)
- Applies minimum cell threshold (20 cells/patient/celltype)

### Step 06: edgeR DE
- Fits patient-level GLM: ~ cancer_type + response_binary
- Tests for response-associated expression (TREAT, lfc > log2(1.2))
- Generates volcano plots and MDS plots

### Step 07: Pathway Enrichment
- Ranks genes by sign(logFC) × sqrt(F)
- Runs fgseaMultilevel on Hallmark and Reactome collections
- Identifies leading-edge genes for significant pathways

### Step 08A: TCGA Download
- Downloads TCGA-LUAD and TCGA-LUSC via curatedTCGAData
- Audits data completeness and sample types
- Prepares for clinical modeling

### Step 08B1: TCGA Clinical Models
- Scores TCGA RNA with ssGSEA (primary) and GSVA (sensitivity)
- Fits Cox PH models: Surv(OS) ~ score + age + sex + stage
- Meta-analyzes across LUAD and LUSC (fixed-effect)

### Step 08B2: Multi-Omics Integration
- Tests program associations with mutation, CNV, methylation, RPPA
- Identifies top features per program per omics layer
- Generates multi-omics summary tables

### Step 09: Final Audit
- Validates all result tables
- Checks evidence tier assignments
- Generates completion markers

### Step 09A: Gap Closure
- Addresses remaining gaps from M7B audit
- Updates manuscript with final corrections
- Generates author completion package

## Output Structure

```
results/
├── step05_pseudobulk/     # Patient-level pseudobulk matrices
├── step06_edgeR/          # DE results and plots
├── step07_programs/       # Pathway enrichment results
├── step08_TCGA/           # TCGA analysis results
└── final/                 # Final audit tables
```

## Key Design Decisions

1. **Patient as biological replicate:** All statistics use patient sampleID, not cell counts
2. **Post-treatment interpretation:** Results describe response-associated mechanisms, not pre-treatment prediction
3. **TCGA as external assessment:** TCGA provides expression-clinical correlation, not immunotherapy validation
4. **Separate histology analysis:** LUAD and LUSC are never merged
