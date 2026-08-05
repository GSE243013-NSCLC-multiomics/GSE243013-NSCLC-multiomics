#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# reproduce.sh — GSE243013 NSCLC Multi-Omics Analysis Pipeline
# =============================================================================
# Usage:
#   bash reproduce.sh                              # Run main pipeline (00-09A)
#   bash reproduce.sh --from-step 05               # Start from Step 05
#   bash reproduce.sh --to-step 07                 # Stop after Step 07
#   bash reproduce.sh --dry-run                    # Show what would run
#   bash reproduce.sh --include-install-steps      # Include 04A and 07A
#   bash reproduce.sh --config config/paths.yml    # Use custom config
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="${SCRIPT_DIR}/scripts/analysis"
LOG_DIR="${SCRIPT_DIR}/logs"
CONFIG_FILE="${SCRIPT_DIR}/config/paths.yml"

FROM_STEP="00"
TO_STEP="99"
DRY_RUN=false
INCLUDE_INSTALL=false
SINGLE_STEP=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --from-step) FROM_STEP="$2"; shift 2 ;;
    --to-step) TO_STEP="$2"; shift 2 ;;
    --step) SINGLE_STEP="$2"; shift 2 ;;
    --config) CONFIG_FILE="$2"; shift 2 ;;
    --dry-run) DRY_RUN=true; shift ;;
    --include-install-steps) INCLUDE_INSTALL=true; shift ;;
    -h|--help)
      echo "Usage: bash reproduce.sh [options]"
      echo "  --from-step N         Start from step N"
      echo "  --to-step N           Stop after step N"
      echo "  --step N              Run only step N"
      echo "  --config PATH         Use config file PATH"
      echo "  --dry-run             Show what would run"
      echo "  --include-install-steps  Include install steps (04A, 07A)"
      exit 0
      ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

export GSE243013_PROJECT_ROOT="${SCRIPT_DIR}"

if ! command -v Rscript &> /dev/null; then
  echo "ERROR: Rscript not found in PATH"; exit 1
fi

R_VER=$(Rscript -e 'cat(R.version$major, R.version$minor, sep=".")')
echo "R version: ${R_VER}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "WARNING: Config not found, copying from paths.example.yml"
  cp "${SCRIPT_DIR}/config/paths.example.yml" "${CONFIG_FILE}"
fi

mkdir -p "${LOG_DIR}"

# Main pipeline steps (install steps excluded by default)
MAIN_STEPS=(00 01 02 03 04 05 06 07 08A 08B1 08B2 09 09A)
INSTALL_STEPS=(04A 07A)

declare -A SCRIPT_MAP
SCRIPT_MAP[00]="00_environment_and_manifest.R"
SCRIPT_MAP[01]="01_download_and_inspect_GSE243013_metadata.R"
SCRIPT_MAP[02]="02_build_GSE243013_patient_manifest.R"
SCRIPT_MAP[03]="03_repair_GSE243013_cohort_definition.R"
SCRIPT_MAP[04]="04_download_and_import_GSE243013_counts.R"
SCRIPT_MAP[04A]="04A_install_BPCells_direct_binary.R"
SCRIPT_MAP[05]="05_build_GSE243013_patient_celltype_pseudobulk.R"
SCRIPT_MAP[06]="06_edgeR_patient_level_differential_expression.R"
SCRIPT_MAP[07]="07_pathway_TF_program_integration.R"
SCRIPT_MAP[07A]="07A_install_extracted_msigdb_and_resume.R"
SCRIPT_MAP[08A]="08A_download_audit_TCGA_multiomics.R"
SCRIPT_MAP[08B1]="08B1_QC2_reconcile_and_rebuild_TCGA_clinical_models.R"
SCRIPT_MAP[08B2]="08B2_TCGA_multiomics_integration.R"
SCRIPT_MAP[09]="09_finalize_project_and_build_evidence_report.R"
SCRIPT_MAP[09A]="09A_close_final_project_gaps.R"

if [[ -n "${SINGLE_STEP}" ]]; then
  STEP_ORDER=("${SINGLE_STEP}")
elif [[ "${INCLUDE_INSTALL}" == true ]]; then
  STEP_ORDER=("${MAIN_STEPS[@]}" "${INSTALL_STEPS[@]}")
else
  STEP_ORDER=("${MAIN_STEPS[@]}")
fi

echo ""
echo "=========================================="
echo "  GSE243013 NSCLC Multi-Omics Pipeline"
echo "=========================================="
echo "  Project root: ${SCRIPT_DIR}"
echo "  Steps: ${STEP_ORDER[*]}"
echo "  Dry run: ${DRY_RUN}"
echo "=========================================="
echo ""

RUN_COUNT=0; FAIL_COUNT=0

for STEP in "${STEP_ORDER[@]}"; do
  if [[ "${STEP}" < "${FROM_STEP}" ]]; then continue; fi
  if [[ "${STEP}" > "${TO_STEP}" ]]; then continue; fi

  SCRIPT="${ANALYSIS_DIR}/${SCRIPT_MAP[$STEP]}"
  LOGFILE="${LOG_DIR}/step${STEP}_$(date +%Y%m%d_%H%M%S).log"

  if [[ ! -f "${SCRIPT}" ]]; then
    echo "[SKIP] Step ${STEP}: Script not found"; continue
  fi

  if [[ "${DRY_RUN}" == true ]]; then
    echo "[DRY-RUN] Step ${STEP}: ${SCRIPT_MAP[$STEP]}"; continue
  fi

  echo "[RUN] Step ${STEP}: ${SCRIPT_MAP[$STEP]}"
  START_TIME=$(date +%s)

  if Rscript "${SCRIPT}" 2>&1 | tee "${LOGFILE}"; then
    ELAPSED=$(( $(date +%s) - START_TIME ))
    echo "[DONE] Step ${STEP} completed in ${ELAPSED}s"
    RUN_COUNT=$((RUN_COUNT + 1))
  else
    ELAPSED=$(( $(date +%s) - START_TIME ))
    echo "[FAIL] Step ${STEP} failed after ${ELAPSED}s"
    FAIL_COUNT=$((FAIL_COUNT + 1))
  fi
  echo ""
done

echo "=========================================="
echo "  Pipeline complete: ${RUN_COUNT} run, ${FAIL_COUNT} failed"
echo "=========================================="
