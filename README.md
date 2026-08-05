# GSE243013 NSCLC Multi-Omics Analysis

## 1. Study Overview

This repository contains the analysis code for a single-cell RNA sequencing (scRNA-seq) study of immune heterogeneity in non-small cell lung cancer (NSCLC) patients treated with neoadjuvant anti-PD-1 immunotherapy. The study uses GSE243013 data (243 patients, 1,254,749 cells) to identify immune-compartment expression programs associated with pathological nonresponse.

**Key finding:** The glycolysis-related program met the prespecified internal integrated-evidence criterion, with negative NES across all 8 scored cell types (FDR < 0.05) and external survival association in TCGA-LUAD (HR = 1.46, P = 5.68e-06).

## 2. Repository Scope

This repository includes:
- All analysis scripts (Step 00-09A) for the computational pipeline
- Manuscript generation and audit scripts (M1-M7C)
- Configuration templates and parameter definitions
- Small final result tables supporting the manuscript
- Documentation of the analysis workflow and data dictionary

This repository does **not** include:
- Raw patient-level scRNA-seq data (must be downloaded from GEO)
- TCGA multi-omics data (downloaded via curatedTCGAData)
- Large intermediate results (BPCells matrices, edgeR objects, scoring matrices)
- Patient-level manifest (reconstructed locally from GEO metadata)
- Complete renv.lock (see [Dependency Limitations](docs/DEPENDENCY_LIMITATIONS.md))

## 3. Main Findings

1. **Glycolysis program:** Immune-compartment glycolysis-related expression program met the prespecified internal integrated-evidence criterion (meta-FDR < 0.05)
2. **Eight cell types:** Negative NES across all 8 scored cell types (FDR < 0.05)
3. **TCGA-LUAD external association:** HR = 1.46 (CI: 1.24-1.73, P = 5.68e-06)
4. **TCGA-LUSC:** No significant association (HR = 1.02, P = 0.750)
5. **Heterogeneity:** I-squared approximately 90%; fixed-effect pooled estimate is descriptive
6. **Multi-omics:** LUAD methylation (30 CpGs, top cg02952918), RPPA (86 proteins, top Cyclin B1); no program-level CNV association result was generated
7. **CollecTRI TF inference:** Not completed for any cell type

## 4. Directory Structure

```
GSE243013-NSCLC-multiomics/
├── config/                    # Configuration templates
│   ├── paths.example.yml      # Project path template
│   ├── analysis_parameters.yml # All analysis parameters
│   └── *.txt                  # Step-specific analysis definitions
├── scripts/
│   ├── analysis/              # Step 00-09A analysis scripts
│   ├── manuscript/            # M1-M7C manuscript scripts
│   └── run_pipeline.R         # Pipeline runner
├── data/
│   ├── README.md              # Data download instructions
│   └── manifests/             # Download manifests and schema
├── resources/
│   └── resource_versions.tsv  # External resource versions
├── results/
│   └── final_tables/          # Small final result tables
├── manuscript/                # Manuscript output files
├── docs/                      # Documentation
├── tests/                     # Pipeline validation tests
├── reproduce.sh               # Main reproduction script
├── README.md                  # This file
├── LICENSE                    # MIT License (code)
├── CITATION.cff               # Citation metadata
└── .gitignore                 # Git ignore rules
```

## 5. Data Availability

- **GSE243013:** Available from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013). Four supplementary files are required (see [data/README.md](data/README.md)).
- **TCGA:** Available via curatedTCGAData R package (version 1.34.0). LUAD and LUSC are analyzed separately.
- **MSigDB:** Hallmark and Reactome gene sets from MSigDB v2024.1.Hs (local extraction).
- **Patient manifest:** Reconstructed locally from downloaded GEO metadata; not redistributed in this repository.

## 6. Hardware and Storage Requirements

- **RAM:** 16 GB minimum (32 GB recommended)
- **Disk:** ~110 GB free space (including raw data)
- **CPU:** Multi-core recommended for edgeR and fgsea
- **OS:** macOS (ARM64) or Linux (x86_64/ARM64)

See [docs/HARDWARE_REQUIREMENTS.md](docs/HARDWARE_REQUIREMENTS.md) for details.

## 7. Software Requirements

Tested with R 4.6.1. Key packages:

| Package | Purpose | Source |
|---------|---------|--------|
| edgeR | Patient-level pseudobulk DE | Bioconductor |
| fgsea | Preranked pathway enrichment | Bioconductor |
| BPCells | Column-major sparse matrix | GitHub |
| metafor | Fixed-effect meta-analysis | CRAN |
| survival | Cox PH models | CRAN |
| curatedTCGAData | TCGA data access | Bioconductor |

See [docs/DEPENDENCY_LIMITATIONS.md](docs/DEPENDENCY_LIMITATIONS.md) for full list and renv status.

## 8. Environment Restoration

```bash
git clone https://github.com/GSE243013-NSCLC-multiomics/GSE243013-NSCLC-multiomics.git
cd GSE243013-NSCLC-multiomics
export GSE243013_PROJECT_ROOT="$(pwd)"
cp config/paths.example.yml config/paths.yml
# Edit config/paths.yml for your environment
```

## 9. Data Download

Before running the pipeline, download 4 files from GEO:
1. `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz`
2. `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz`
3. `GSE243013_genes.csv.gz`
4. `GSE243013_barcodes.csv.gz`

See [data/README.md](data/README.md) for detailed instructions.

TCGA data is downloaded automatically by Step 08A via `curatedTCGAData`.

## 10. Full Execution Order

**Optional prerequisites** (run only if needed):
```bash
bash reproduce.sh --step 04A    # BPCells binary install
bash reproduce.sh --step 07A   # MSigDB extraction install
# Or include both:
bash reproduce.sh --include-install-steps
```

**Main pipeline:**
```bash
bash reproduce.sh              # Run all steps (00-09A)
bash reproduce.sh --from-step 05 --to-step 07  # Partial run
bash reproduce.sh --dry-run    # Show what would run
```

| Step | Description | Runtime |
|------|-------------|---------|
| 00 | R environment check | Short |
| 01 | Metadata download and inspection | Short |
| 02 | Patient manifest construction | Long |
| 03 | Cohort definition repair | Long |
| 04 | Count matrix import (BPCells) | Long |
| 05 | Patient-level pseudobulk | Long |
| 06 | edgeR differential expression | Long |
| 07 | Pathway enrichment (fgsea) | Long |
| 08A | TCGA data download | Long |
| 08B1 | TCGA clinical modeling (Cox) | Long |
| 08B2 | TCGA multi-omics integration | Long |
| 09 | Final audit and evidence report | Long |
| 09A | Gap closure and revision | Long |

## 11. Expected Outputs

Key outputs in `results/final_tables/`:
- `Table_1_cohort_characteristics.csv` - Patient cohort summary (aggregate counts only)
- `Table_3_primary_edgeR_results.csv` - Primary DE results
- `Table_4_pathway_TF_programs.csv` - Pathway and TF programs
- `Table_5_TCGA_clinical_validation.csv` - TCGA clinical models (external association)
- `GSE243013_final_evidence_tiers_revised.csv` - Evidence tier classification
- `GSE243013_cohort_counts_public.csv` - Aggregate cohort counts

## 12. Reproducing Manuscript Tables and Figures

```bash
# Run manuscript generation scripts (M1-M7C)
bash reproduce.sh --from-step M1 --to-step M7C
```

## 13. Known Limitations

- **CollecTRI:** TF inference not completed for any cell type
- **PROGENy:** Pathway activity inference not completed
- **renv:** Full renv.lock not generated due to binary package dependencies
- **BPCells:** Binary installation required for ARM64 macOS
- **TCGA:** Not an immunotherapy cohort; provides expression-clinical correlation data only
- **CNV:** No program-level association result was generated
- **Heterogeneity:** I-squared approximately 90%; fixed-effect pooled estimate is descriptive

## 14. Citation

```bibtex
@software{Xie2026GSE243013,
  title = {GSE243013 NSCLC Multi-Omics Analysis},
  author = {Xie, Xiaowei and Guo, Ying and Zhang, Yang},
  year = {2026},
  url = {https://github.com/GSE243013-NSCLC-multiomics/GSE243013-NSCLC-multiomics}
}
```

See [CITATION.cff](CITATION.cff) for full citation metadata.

## 15. License

- **Code:** MIT License (see [LICENSE](LICENSE))
- **Data:** Original data providers retain all rights. See GEO and TCGA data use agreements.
- **Third-party resources:** MSigDB, curatedTCGAData, and other packages retain their original licenses.

## 16. Contact

- **Corresponding author:** Xiaowei Xie
- **Email:** xiexiaowei@fjmu.edu.cn
- **ORCID:** 0009-0009-4221-2777
