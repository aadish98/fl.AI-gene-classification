#!/bin/bash
#SBATCH --job-name=backfill_ref_metadata
#SBATCH --mail-user=aadishms@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --account=rallada0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=24:00:00
#SBATCH --output=/nfs/turbo/umms-rallada/UM\ Lab\ Users/Aadish/fl-ai_gene_classification/Logs/%x-%A.log

set -euo pipefail

PROJECT_DIR="/nfs/turbo/umms-rallada/UM Lab Users/Aadish/fl-ai_gene_classification"
SCRIPT_PATH="${PROJECT_DIR}/HelperScripts/backfill_reference_metadata.py"
VENV_DIR="${PROJECT_DIR}/.venv-py310"
INPUT_DIR="${1:-${PROJECT_DIR}/Data/Tx-Omics_Shiju-Corr_Sets}"
MODE="${2:---dry-run}"

if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "Error: input directory not found: ${INPUT_DIR}"
  exit 1
fi

if [[ "${MODE}" != "--dry-run" && "${MODE}" != "--apply" ]]; then
  echo "Error: mode must be --dry-run or --apply"
  exit 1
fi

if ! type module >/dev/null 2>&1; then
  source /etc/profile.d/modules.sh
fi
module purge
module load python/3.10.4

source "${VENV_DIR}/bin/activate"

REPORT_PATH="${PROJECT_DIR}/Logs/backfill_reference_metadata_${SLURM_JOB_ID:-manual}.json"
echo "[$(date)] Running backfill script on ${INPUT_DIR} (${MODE})"
python "${SCRIPT_PATH}" "${INPUT_DIR}" ${MODE} --report-path "${REPORT_PATH}"
echo "[$(date)] Done. Report: ${REPORT_PATH}"
