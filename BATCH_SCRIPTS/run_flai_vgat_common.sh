#!/bin/bash
#SBATCH --job-name=flai_vgat_common
#SBATCH --mail-user=aadishms@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --account=rallada0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=20:00:00
#SBATCH --output=/nfs/turbo/umms-rallada/UM\ Lab\ Users/Aadish/fl-ai_gene_classification/Logs/%x-%j.log

set -euo pipefail

PROJECT_DIR="/nfs/turbo/umms-rallada/UM Lab Users/Aadish/fl-ai_gene_classification"
SCRIPT_PATH="${PROJECT_DIR}/flai-gene-classification.py"
VENV_DIR="${PROJECT_DIR}/.venv-py310"
INPUT_DIR="${PROJECT_DIR}/Data/Test-vGAT-extra-small"

echo "[$(date)] Starting job on ${INPUT_DIR}"

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

echo "[$(date)] Finished job"
