#!/bin/bash
#SBATCH --job-name=trajpp6h
#SBATCH --partition=zen16           # [CHECK] spiritx partition (sinfo -s)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
# ---------------------------------------------------------------------------
# Stage 05: trajectory process decomposition (advection/adiabatic/diabatic/
# seasonal) for the ERA5T 2026 event, from the results_6h/ trajectories.
# Submit:  sbatch 05_traj_postprocess.slurm
# ---------------------------------------------------------------------------
set -euo pipefail
# Slurm runs a spool copy of this script, so BASH_SOURCE is useless here:
# use PIPELINE_DIR (exported by submit_all.sh), else the sbatch submit dir.
HERE="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "${HERE}/config.sh"          # activates the conda env, sets PYTHON

NJOBS="${SLURM_CPUS_PER_TASK:-32}"
"${PYTHON}" "${HERE}/05_traj_postprocess_6h.py" "${NJOBS}"
