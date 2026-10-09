#!/bin/bash
#SBATCH --job-name=origin6h
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/homedata/%u/lagranto_run/slurm_logs/%x-%j.out
# ---------------------------------------------------------------------------
# Stage 06: origin contributions (TA + per-process genesis + exact arrival)
# for the ERA5T 2026 event, from results_6h/traj_netcdfs/lagrangian_terms_*.nc.
# Writes $NC_RESULTS_DIR/origin_contributions.nc (on /scratchx) and copies it
# to $RESULTS_DIR (on /homedata).
# Submit:  sbatch 06_origin_contributions.slurm
# ---------------------------------------------------------------------------
set -euo pipefail
# Slurm runs a spool copy of this script, so BASH_SOURCE is useless here:
# use PIPELINE_DIR (exported by submit_all.sh), else the sbatch submit dir.
HERE="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "${HERE}/config.sh"          # activates the conda env, sets PYTHON

NJOBS="${SLURM_CPUS_PER_TASK:-32}"
"${PYTHON}" "${HERE}/06_origin_contributions_6h.py" "${NJOBS}"
