#!/bin/bash
#SBATCH --job-name=sanitize_pubmed_cache
#SBATCH --mail-user=aadishms@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --account=rallada0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/turbo/umms-rallada/UM\ Lab\ Users/Aadish/fl-ai_gene_classification/Logs/%x-%A.log

set -euo pipefail

PROJECT_DIR="/nfs/turbo/umms-rallada/UM Lab Users/Aadish/fl-ai_gene_classification"
SCRIPT_PATH="${PROJECT_DIR}/HelperScripts/sanitize_pubmed_caches.py"
VENV_DIR="${PROJECT_DIR}/.venv-py310"
MODE="${1:---dry-run}"

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

REPORT_PATH="${PROJECT_DIR}/Logs/sanitize_pubmed_caches_${SLURM_JOB_ID:-manual}.json"
echo "[$(date)] Running cache sanitation (${MODE})"
python "${SCRIPT_PATH}" ${MODE} --report-path "${REPORT_PATH}"
echo "[$(date)] Done. Report: ${REPORT_PATH}"
