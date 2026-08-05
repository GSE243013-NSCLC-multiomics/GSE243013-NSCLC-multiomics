# GSE243013 Data Download Manifest

## Overview

This repository does **not** redistribute raw patient-level data.
Raw data must be downloaded separately from GEO and TCGA before running the pipeline.

The patient-level manifest is reconstructed locally from the downloaded GEO metadata and is not redistributed in this repository.

## GSE243013 (GEO)

**Access date:** 20 July 2026
**Source:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013

### Required Files

| File | Purpose | Destination | Checksum | Redistributed | Download Source |
|------|---------|-------------|----------|---------------|-----------------|
| `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz` | Sparse count matrix (cells x genes) | `data/raw/GSE243013/` | See GEO | No | GEO supplementary |
| `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz` | Cell-level annotations with patient IDs | `data/raw/GSE243013/` | See GEO | No | GEO supplementary |
| `GSE243013_genes.csv.gz` | Gene names and IDs (for matrix assembly) | `data/raw/GSE243013/` | See GEO | No | GEO supplementary |
| `GSE243013_barcodes.csv.gz` | Cell barcode identifiers (for matrix assembly) | `data/raw/GSE243013/` | See GEO | No | GEO supplementary |

### Download Steps

1. Visit https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013
2. Download all 4 supplementary files listed above
3. Place them in `data/raw/GSE243013/`
4. Run `scripts/analysis/04_download_and_import_GSE243013_counts.R` to build BPCells format

## TCGA Multi-Omics

| Field | Value |
|-------|-------|
| Source | TCGA via curatedTCGAData (R package) |
| Data version | 2.1.1 |
| curatedTCGAData version | 1.34.0 |
| Cohorts | TCGA-LUAD, TCGA-LUSC (analyzed separately) |
| Access method | Automated download via `curatedTCGAData::curatedTCGAData()` in Step 08A |
| Expected destination | `data/raw/TCGA/` (via ExperimentHub cache) |
| Redistributed in repository | No |

### Download Steps

1. Run `scripts/analysis/08A_download_audit_TCGA_multiomics.R`
2. Data will be cached in the ExperimentHub cache directory
3. Set `RUSER_CACHE_DIR` environment variable if needed

## Expected Data Directory Structure

```
data/
├── raw/
│   ├── GSE243013/
│   │   ├── GSE243013_NSCLC_immune_scRNA_counts.mtx.gz
│   │   ├── GSE243013_NSCLC_immune_scRNA_metadata.csv.gz
│   │   ├── GSE243013_genes.csv.gz
│   │   ├── GSE243013_barcodes.csv.gz
│   │   └── bpcells/
│   │       └── GSE243013_counts_colmajor/  (created by Step 04)
│   └── TCGA/
│       ├── LUAD/
│       └── LUSC/
├── manifests/
│   ├── GSE243013_download_manifest.tsv
│   ├── TCGA_resource_manifest.tsv
│   └── GSE243013_patient_manifest_schema.tsv
└── processed/
    └── (created by pipeline)
```

## Patient Manifest

The patient-level manifest is not redistributed in this repository. It is reconstructed locally from the downloaded GEO metadata by running Steps 01-03. See `data/manifests/GSE243013_patient_manifest_schema.tsv` for the column schema.
