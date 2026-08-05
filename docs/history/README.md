# History Directory

This directory contains **superseded scripts** that are kept for version audit purposes only.

## Purpose

- Document the evolution of the analysis pipeline
- Provide audit trail for methodological decisions
- Not part of the default reproduction workflow

## Scripts

| Script | Reason for Supersession |
|--------|------------------------|
| `08B1_TCGA_program_scoring_SUPERSeded.R` | logHR extraction bug (hr_row[2]=exp(-coef)) |
| `08B1_QC_audit_FAIL_FDR.R` | FAIL_FDR_CALCULATION - 0% pass rate |
| `07_pathway_TF_program_resume_SUPERSeded.R` | Resume script; use main 07 instead |

## Usage

These scripts are **NOT** called by `reproduce.sh` or `run_pipeline.R`.
They should **NOT** be used as sources for statistical results.
They are retained solely for reproducibility audit purposes.
