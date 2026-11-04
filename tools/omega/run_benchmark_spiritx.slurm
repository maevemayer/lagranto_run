#!/bin/bash
#SBATCH --job-name=omega_bench
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-12%4
#SBATCH --output=/homedata/%u/omega_check/logs/%x_%A_%a.log
# ---------------------------------------------------------------------------
# Omega benchmark on /bdd/CMIP6, historical JJA 2000, every model with 6hrLev: first the instantaneous check
# against CFsubhr wap where the model has it, then daily means against its day/Eday wap (the only comparison
# MRI-ESM2-0, the daily reference on spiritx, has). CNRM family with the reduced-grid reconstruction.
# A model without the needed data fails here and stays NOT_BENCHMARKED.
# Submit from the repo root, then build the table:
#   sbatch tools/omega/run_benchmark_spiritx.slurm
#   python tools/omega/omega_benchmark.py --site spiritx --results /homedata/$USER/omega_check/results
# ---------------------------------------------------------------------------
set -u
SRC="${SLURM_SUBMIT_DIR:-$(pwd)}/tools/omega"
PYTHON="${CONDA_ENV:-/homedata/${USER}/conda/thesis}/bin/python"
export CMIP6_POOL="${CMIP6_POOL:-/bdd/CMIP6}" CMIP6_EXTRA_POOL="" MPLBACKEND=Agg
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
MODELS=(MRI-ESM2-0 MIROC6 CNRM-ESM2-1 CNRM-CM6-1 CNRM-CM6-1-HR EC-Earth3 NorESM2-MM CESM2
        GISS-E2-1-G GISS-E3-G IPSL-CM6A-LR IPSL-CM6A-LR-INCA MPI-ESM1-2-HR)
MODEL=${MODELS[${SLURM_ARRAY_TASK_ID:-0}]}
OUT=/homedata/${USER}/omega_check/results/${MODEL}/benchmark
mkdir -p "${OUT}"
echo "=== ${MODEL} on $(hostname), $(date)"
RANGE=(--daily --no-figure --start 2000-06-01 --end 2000-08-31 --outdir "${OUT}")
case "${MODEL}" in
  CNRM-*)
    JSON=/homedata/${USER}/omega_check/results/${MODEL}/regrid/regrid_${MODEL}.json
    [ -f "${JSON}" ] || "${PYTHON}" "${SRC}/omega_check_cnrm.py" "${MODEL}" --regrid-check --outdir "$(dirname "${JSON}")"
    "${PYTHON}" "${SRC}/omega_check_sites.py" "${MODEL}" --unregrid "${JSON}" --outdir "${OUT}"
    "${PYTHON}" "${SRC}/omega_check_cnrm.py" "${MODEL}" --unregrid "${JSON}" "${RANGE[@]}" ;;
  *)
    "${PYTHON}" "${SRC}/omega_check_sites.py" "${MODEL}" --outdir "${OUT}"
    "${PYTHON}" "${SRC}/omega_check.py" "${MODEL}" "${RANGE[@]}" ;;
esac
echo "=== ${MODEL} finished with status $? $(date)"
