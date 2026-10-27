#!/bin/bash
#SBATCH --job-name=prepcmip6
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/homedata/%u/lagranto_run/slurm_logs/%x-%j.out
# ---------------------------------------------------------------------------
# CMIP6 P- and S-files (prepare_cmip6.py) for the run in config.sh.
# Normally submitted by submit_cmip6.sh; PREP_ARGS=--test writes one time step.
# ---------------------------------------------------------------------------
set -euo pipefail
HERE="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "${HERE}/config.sh"          # activates the conda env, sets PYTHON
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

"${PYTHON}" "${HERE}/prepare_cmip6.py" ${PREP_ARGS:-}
