#!/bin/bash
#SBATCH --job-name=flai_txomics_corr_sets
#SBATCH --mail-user=aadishms@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --account=rallada0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=24:00:00
#SBATCH --array=0-1
#SBATCH --output=/nfs/turbo/umms-rallada/UM\ Lab\ Users/Aadish/fl-ai_gene_classification/Logs/%x-%A_%a.log

set -euo pipefail

PROJECT_DIR="/nfs/turbo/umms-rallada/UM Lab Users/Aadish/fl-ai_gene_classification"
SCRIPT_PATH="${PROJECT_DIR}/flai-gene-classification.py"
VENV_DIR="${PROJECT_DIR}/.venv-py310"
PARENT_DIR="${PROJECT_DIR}/Data/Tx-Omics_Shiju-Corr_Sets"

if [[ ! -d "${PARENT_DIR}" ]]; then
  echo "Error: parent directory not found: ${PARENT_DIR}"
  exit 1
fi

# Discover immediate subfolders under Tx-Omics_Shiju-Corr_Sets.
mapfile -t INPUT_DIRS < <(
  python3 - "${PARENT_DIR}" <<'PY'
import os
import sys

base = sys.argv[1]
subdirs = [
    os.path.join(base, name)
    for name in sorted(os.listdir(base))
    if os.path.isdir(os.path.join(base, name))
]
for path in subdirs:
    print(path)
PY
)

if [[ "${#INPUT_DIRS[@]}" -eq 0 ]]; then
  echo "Error: no subfolders found in ${PARENT_DIR}"
  exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
if (( TASK_ID < 0 || TASK_ID >= ${#INPUT_DIRS[@]} )); then
  echo "Error: SLURM_ARRAY_TASK_ID=${TASK_ID} is out of range."
  echo "Valid range: 0-$(( ${#INPUT_DIRS[@]} - 1 ))"
  echo "Discovered subfolders:"
  printf '  - %s\n' "${INPUT_DIRS[@]}"
  exit 1
fi

INPUT_DIR="${INPUT_DIRS[$TASK_ID]}"

echo "[$(date)] Starting array task ${TASK_ID} on ${INPUT_DIR}"
echo "[$(date)] Total discovered subfolders: ${#INPUT_DIRS[@]}"

# Load Python 3.10.4 from module system.
if ! type module >/dev/null 2>&1; then
  source /etc/profile.d/modules.sh
fi
module purge
module load python/3.10.4

echo "[$(date)] Python from module: $(python --version 2>&1)"

source "${VENV_DIR}/bin/activate"

python "${SCRIPT_PATH}" \
  "${INPUT_DIR}" \
  --keywords "sleep,circadian" \
  --reference-limit 500

echo "[$(date)] Finished array task ${TASK_ID}"

