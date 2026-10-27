#!/bin/bash
#SBATCH --job-name=window_scan
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=0-1
#SBATCH --output=/homedata/%u/omega_check/logs/%x_%A_%a.log
# ---------------------------------------------------------------------------
# Which averaging window do the models' daily-mean omega (day wap) cover? omega_check.py --daily
# --window-scan, JJA 2000: computed daily means over 00-24 UTC shifted by -12..+12 h vs day wap.
# MRI-ESM2-0 is the control (expected best offset 0 h).
# Submit from the repo root:  sbatch tools/omega/run_windowscan_spiritx.slurm
# ---------------------------------------------------------------------------
set -u
SRC="${SLURM_SUBMIT_DIR:-$(pwd)}/tools/omega"
PYTHON="${CONDA_ENV:-/homedata/${USER}/conda/thesis}/bin/python"
export CMIP6_POOL="${CMIP6_POOL:-/bdd/CMIP6}" MPLBACKEND=Agg
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
MODELS=(EC-Earth3 MRI-ESM2-0)
MODEL=${MODELS[${SLURM_ARRAY_TASK_ID:-0}]}
OUT=/homedata/${USER}/omega_check/results/${MODEL}/window_scan
mkdir -p "${OUT}"
echo "=== ${MODEL} on $(hostname), $(date)"
"${PYTHON}" "${SRC}/omega_check.py" "${MODEL}" --daily --window-scan --no-figure \
    --start 2000-06-01 --end 2000-08-31 --outdir "${OUT}"
echo "=== ${MODEL} finished $(date)"
