#!/bin/bash
#SBATCH --job-name=flai_CSW_FC0.5
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

# Inputs for each array task index
INPUT_DIRS=(
  "${PROJECT_DIR}/Data/Test3"
  "${PROJECT_DIR}/Data/Test4"
)

INPUT_DIR="${INPUT_DIRS[$SLURM_ARRAY_TASK_ID]}"

echo "[$(date)] Starting array task ${SLURM_ARRAY_TASK_ID} on ${INPUT_DIR}"

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

echo "[$(date)] Finished array task ${SLURM_ARRAY_TASK_ID}"