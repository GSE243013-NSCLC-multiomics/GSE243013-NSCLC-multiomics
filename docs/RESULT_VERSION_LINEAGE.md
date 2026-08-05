# Result Version Lineage

## Overview

This document tracks the version history of key results and identifies which versions are current vs superseded.

## Analysis Pipeline Versions

### Step 06: edgeR Differential Expression
- **Version:** Original (no known issues)
- **Status:** CURRENT
- **Notes:** Patient-level pseudobulk; TREAT with lfc > log2(1.2)

### Step 07: Pathway Enrichment
- **Version:** Original with extracted MSigDB v2024.1.Hs
- **Status:** CURRENT
- **Notes:** fgseaMultilevel; Hallmark (50) + Reactome (1839) gene sets

### Step 08B1: TCGA Clinical Modeling

#### Version History
1. **Original (08B1):** SUPERSEDED
   - **Issue:** logHR extraction bug at line 772
   - **Bug:** `hr_row[2]` stored `exp(-coef)` instead of true `coef`
   - **Impact:** All logHR values incorrect

2. **QC Audit (08B1-QC):** SUPERSEDED
   - **Result:** FAIL_FDR_CALCULATION
   - **Issue:** 0% Cox comparison pass rate
   - **Impact:** Could not validate recalculated results

3. **QC2 (08B1-QC2):** CURRENT
   - **Result:** PASS
   - **Fix:** Corrected extraction using `summary(fit)$coefficients["score_z", "coef"]`
   - **Status:** 290 canonical Cox models; Level A=0→1 (after 09A), Level B=5, Level C=140

### Step 08B2: TCGA Multi-Omics Integration
- **Version:** Original
- **Status:** CURRENT (EXPLORATORY)
- **Notes:** 25 programs across 4 omics; CNV no results for approved programs

## Evidence Tier Changes

### Before 09A
- Level A: 0 programs
- Level B: 24 programs
- Level C: 1 program

### After 09A
- Level A: 1 program (HALLMARK_GLYCOLYSIS)
- Level B: 24 programs
- Level C: 1 program

### Change Reason
- 09A corrected the meta-FDR calculation for HALLMARK_GLYCOLYSIS
- Meta-FDR = 0.0499205324269691 < 0.05 threshold
- Upgraded from Level B to Level A

## Manuscript Versions

### M1-M4: Initial Drafts
- **Status:** Superseded by M5-M7

### M5: Scientific Revision
- **Status:** Superseded by M5A-M7

### M5A: Evidence Completion
- **Status:** Superseded by M6-M7

### M6: Text Assembly
- **Status:** Superseded by M7

### M7: Terminal Correction
- **Status:** CURRENT (Scientifically Final)
- **Changes:** Fixed cell counts (1,254,749), cohort numbers, TCGA language

### M7C: Author Completion
- **Status:** CURRENT (Author-Ready)
- **Changes:** Expanded Introduction (727 words), citation tasks, author input form

## Key Verified Values

| Metric | Value | Source |
|--------|-------|--------|
| Total cells | 1,254,749 | BPCells validation CSV |
| Total genes | 31,831 | BPCells validation CSV |
| Primary NES (glycolysis) | -2.3589 | fgsea exact statistics |
| Primary FDR | 3.00e-11 | fgsea exact statistics |
| Meta-FDR | 0.0499 | Canonical meta-analysis |
| LUAD HR | 1.465 | Cox model |
| LUSC HR | 1.023 | Cox model |
| I² heterogeneity | 90.44% | Meta-analysis |
