# Reproducibility Guide

## Overview

This document describes how to reproduce the analysis results from the GSE243013 NSCLC multi-omics study.

## Prerequisites

1. **R ≥ 4.6.1** with the following packages:
   - edgeR (≥ 4.10.1)
   - fgsea
   - BPCells
   - metafor
   - survminer
   - digest

2. **Raw data:**
   - GSE243013 count matrix and metadata from GEO
   - TCGA data (downloaded automatically by Step 08A)

3. **Disk space:** ~110 GB (including raw data and intermediate results)

4. **RAM:** 16 GB minimum (32 GB recommended)

## Reproduction Steps

### 1. Clone and Configure

```bash
git clone https://github.com/GSE243013-NSCLC-multiomics/GSE243013-NSCLC-multiomics.git
cd GSE243013-NSCLC-multiomics
export GSE243013_PROJECT_ROOT="$(pwd)"
cp config/paths.example.yml config/paths.yml
# Edit config/paths.yml for your environment
```

### 2. Download Data

```bash
# Download GSE243013 from GEO
# Place files in data/raw/GSE243013/

# TCGA data is downloaded automatically by Step 08A
```

### 3. Run Full Pipeline

```bash
bash reproduce.sh
```

### 4. Run Individual Steps

```bash
# Run specific steps
bash reproduce.sh --from-step 05 --to-step 07

# Run single step
bash reproduce.sh --step 06

# Dry run
bash reproduce.sh --dry-run
```

## Key Results Verification

### Glycolysis Program (Level A)

Expected values:
- Primary NES: -2.3589
- Primary FDR: 3.00e-11
- Strict NES: -2.4126
- Strict FDR: 7.51e-12
- Meta-FDR: 0.0499

### TCGA Survival Association

Expected values:
- LUAD HR: 1.465 (CI: 1.242-1.727, P = 5.68e-06)
- LUSC HR: 1.023 (CI: 0.888-1.179, P = 0.750)
- Meta HR: 1.192 (CI: 1.059-1.342)

## Known Limitations

1. **renv.lock not available:** Binary package dependencies prevent full renv lock file generation
2. **BPCells binary:** Requires compilation or binary installation for ARM64 macOS
3. **CollecTRI:** TF inference not completed
4. **PROGENy:** Pathway activity inference not completed
5. **MSigDB version:** Local extraction used for v2024.1.Hs; may differ from latest MSigDB

## Troubleshooting

### BPCells Installation

```bash
# For ARM64 macOS
Rscript scripts/analysis/04A_install_BPCells_direct_binary.R
```

### MSigDB Extraction

```bash
# Download MSigDB v2024.1.Hs from MSigDB website
# Extract to ~/Downloads/msigdbr_manual/msigdb/
# Or set MSIGDB_EXTRACTED_DIR environment variable
export MSIGDB_EXTRACTED_DIR="/path/to/extracted/msigdb"
```

### TCGA Download Issues

```bash
# If curatedTCGAData fails, check ExperimentHub cache
# Set custom cache directory
export RUSER_CACHE_DIR="/path/to/cache"
```
