#!/bin/bash
#SBATCH --job-name=omega_chk
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-2
#SBATCH --output=/homedata/%u/omega_check/logs/%x_%A_%a.log
# ---------------------------------------------------------------------------
# Validation of the diagnosed vertical velocity (omega_check.py 1.7.8 /
# omega_check_cnrm.py cnrm-1.11) on the ESPRI mirror /bdd/CMIP6, historical, JJA 2000:
# computed omega (continuity equation from 6hrLev ua, va, ps) vs the model's own omega.
#   0 MRI-ESM2-0   native Gaussian grid (gn)            -> standard
#   1 CNRM-CM6-1   regridded from ARPEGE's reduced grid -> regrid check + unregrid
#   2 EC-Earth3    regridded IFS output (gr)            -> regrid check; unregrid if the
#                                                          rows come from fewer native points
# Submit from the repo root:
#   mkdir -p /homedata/$USER/omega_check/logs && sbatch tools/omega/run_omega_spiritx.slurm
#   sbatch --array=2 tools/omega/run_omega_spiritx.slurm        # one model again
# Results: /homedata/$USER/omega_check/results/<model>/ (summary .txt, score .csv, figures)
# ---------------------------------------------------------------------------
set -u
SRC="${SLURM_SUBMIT_DIR:-$(pwd)}/tools/omega"
PYTHON="${CONDA_ENV:-/homedata/${USER}/conda/thesis}/bin/python"
export CMIP6_POOL="${CMIP6_POOL:-/bdd/CMIP6}"      # same directory layout as the DKRZ pool
export CMIP6_EXTRA_POOL="${CMIP6_EXTRA_POOL:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MPLBACKEND=Agg

# model           start       end         treatment
MODELS=(
  "MRI-ESM2-0     2000-06-01  2000-08-31  std"
  "CNRM-CM6-1     2000-06-01  2000-08-31  cnrm"
  "EC-Earth3      2000-06-01  2000-08-31  auto"
)
read -r MODEL START END KIND <<< "${MODELS[${SLURM_ARRAY_TASK_ID:-0}]}"
EXP=historical
OUT=/homedata/${USER}/omega_check/results/${MODEL}
mkdir -p "${OUT}"
echo "=== ${MODEL}  ${EXP}  ${START} .. ${END}  [${KIND}]  on $(hostname), $(date)"
"${PYTHON}" -c "import sys; sys.path.insert(0, '${SRC}'); import omega_check, omega_check_cnrm; \
print('omega_check', omega_check.__version__, '| omega_check_cnrm', omega_check_cnrm.__version__, '| pool', omega_check.POOL)" \
  || { echo "=== cannot import the omega scripts with ${PYTHON}"; exit 1; }

echo "--- native omega (wap) on the mirror for ${MODEL} ${EXP}:"
ls -d ${CMIP6_POOL}/*/*/${MODEL}/${EXP}/*/*/wap 2>/dev/null | sed "s|${CMIP6_POOL}/||" || true
echo "--- 6hrLev members with ua:"
ls -d ${CMIP6_POOL}/*/*/${MODEL}/${EXP}/*/6hrLev/ua 2>/dev/null | sed "s|${CMIP6_POOL}/||" || true

step () {   # step <label> <command ...>: run and report, never stop the job
  local label=$1; shift
  echo "--- ${MODEL}: ${label}   ($(date +%H:%M))"
  if "$@"; then echo "=== ${MODEL}: ${label} done"; return 0
  else echo "=== ${MODEL}: ${label} FAILED or not available"; return 1; fi
}

run_std () {
  step "pressure levels (E3hrPt / 6hrPlev / day)" \
       "${PYTHON}" "${SRC}/omega_check.py" "${MODEL}" --exp ${EXP} --start ${START} --end ${END} --outdir "${OUT}" \
  || step "daily means vs day/Eday (no sub-daily reference)" \
       "${PYTHON}" "${SRC}/omega_check.py" "${MODEL}" --daily --exp ${EXP} --start ${START} --end ${END} --outdir "${OUT}"
  step "model levels vs CFday" \
       "${PYTHON}" "${SRC}/omega_check.py" "${MODEL}" --cfday --exp ${EXP} --start ${START} --end ${END} --outdir "${OUT}"
}

run_unregrid () {   # needs ${JSON}
  step "model levels vs CFday, regridding undone" \
       "${PYTHON}" "${SRC}/omega_check_cnrm.py" "${MODEL}" --cfday --unregrid "${JSON}" --exp ${EXP} --start ${START} --end ${END} --outdir "${OUT}"
  step "pressure levels (6hrPlev, day), regridding undone" \
       "${PYTHON}" "${SRC}/omega_check_cnrm.py" "${MODEL}" --unregrid "${JSON}" --exp ${EXP} --start ${START} --end ${END} --outdir "${OUT}" \
  || step "daily means vs day/Eday, regridding undone" \
       "${PYTHON}" "${SRC}/omega_check_cnrm.py" "${MODEL}" --daily --unregrid "${JSON}" --exp ${EXP} --start ${START} --end ${END} --outdir "${OUT}"
}

if [ "${KIND}" = std ]; then
  run_std
else
  JSON=${OUT}/regrid/regrid_${MODEL}.json
  [ -f "${JSON}" ] || step "regrid check (native points per latitude row)" \
       "${PYTHON}" "${SRC}/omega_check_cnrm.py" "${MODEL}" --regrid-check --exp ${EXP} --outdir "${OUT}/regrid"
  NRED=0
  if [ -f "${JSON}" ]; then
    NRED=$("${PYTHON}" -c "import json; d=json.load(open('${JSON}')); print(sum(1 for r in d['rows'] if r['scheme'] != 'none' and r['n'] < d['nlon']))")
  fi
  echo "=== ${MODEL}: ${NRED} latitude rows interpolated from fewer native points"
  if [ "${KIND}" = cnrm ] || [ "${NRED}" -gt 0 ]; then
    if [ -f "${JSON}" ]; then run_unregrid; else echo "=== ${MODEL}: no regrid table -> standard run instead"; run_std; fi
  else
    run_std
  fi
fi
echo "=== ${MODEL} finished $(date)"
