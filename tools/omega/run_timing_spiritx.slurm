#!/bin/bash
#SBATCH --job-name=timing_chk
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --array=0-1
#SBATCH --output=/homedata/%u/omega_check/logs/%x_%A_%a.log
# ---------------------------------------------------------------------------
# Time stamps of 6hrLev: 6hrLev ua vs instantaneous 6hrPlevPt ua and 3hr uas of the same
# model at offsets -12..+12 h (timing_check.py), June 2000. MRI-ESM2-0 is the control.
# Submit from the repo root:  sbatch tools/omega/run_timing_spiritx.slurm
# ---------------------------------------------------------------------------
set -u
SRC="${SLURM_SUBMIT_DIR:-$(pwd)}/tools/omega"
PYTHON="${CONDA_ENV:-/homedata/${USER}/conda/thesis}/bin/python"
export CMIP6_POOL="${CMIP6_POOL:-/bdd/CMIP6}"
MODELS=(EC-Earth3 MRI-ESM2-0)
MODEL=${MODELS[${SLURM_ARRAY_TASK_ID:-0}]}
echo "=== ${MODEL} on $(hostname), $(date)"
"${PYTHON}" "${SRC}/timing_check.py" "${MODEL}" --start 2000-06-01 --end 2000-06-30
echo "=== ${MODEL} finished $(date)"
