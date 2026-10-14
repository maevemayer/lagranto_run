#!/bin/bash
#SBATCH --job-name=cmip6chk
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=160G
#SBATCH --time=06:00:00
#SBATCH --signal=B:USR1@600
#SBATCH --output=/homedata/%u/cmip6_check/%x-%j.out
# ---------------------------------------------------------------------------
# Read-only integrity check of /bdd/CMIP6 6hrLev hus/ta/ua/va.
# Submit:  mkdir -p /homedata/$USER/cmip6_check && sbatch tools/check_cmip6_6hrlev.slurm
# Default = all activities/experiments. Extra args go to the script, e.g.
#   --quick                              filename inventory only (fast)
#   --best                               longest clean member per model/experiment
#   --activity ScenarioMIP --exp 'ssp*'  one slice per job
#   --fresh                              ignore the checkpoint, check everything again
#
# Resuming: every finished check is saved to checkpoint_<tag>.jsonl in the
# output dir. Resubmitting the same command skips what is already done.
# 10 min before the time limit the job resubmits itself (at most
# CMIP6CHK_MAXCHAIN=5 times), so a timeout continues on its own. A crash or
# scancel does not resubmit; sbatch the same command again to resume.
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(pwd)}"
PYTHON="${CONDA_ENV:-/homedata/${USER}/conda/thesis}/bin/python"
CHAIN="${CMIP6CHK_CHAIN:-0}"
MAXCHAIN="${CMIP6CHK_MAXCHAIN:-5}"
ARGS=("$@")

resubmit() {
    if (( CHAIN < MAXCHAIN )); then
        sbatch --dependency=afterany:"${SLURM_JOB_ID}" --chdir="${HERE}" \
               --export=ALL,CMIP6CHK_CHAIN=$((CHAIN + 1)) \
               "${HERE}/tools/check_cmip6_6hrlev.slurm" ${ARGS[@]+"${ARGS[@]}"}
        echo "time limit near: continuation $((CHAIN + 1))/${MAXCHAIN} queued, it resumes from the checkpoint"
    else
        echo "time limit near: ${MAXCHAIN} continuations used, sbatch the same command again to resume"
    fi
}
trap resubmit USR1

# python runs in the background so bash can act on USR1 while it works
"${PYTHON}" "${HERE}/tools/check_cmip6_6hrlev.py" --njobs "${SLURM_CPUS_PER_TASK:-64}" ${ARGS[@]+"${ARGS[@]}"} &
PID=$!
wait "${PID}"; rc=$?
while kill -0 "${PID}" 2>/dev/null; do        # wait was interrupted by the signal
    wait "${PID}"; rc=$?
done
exit "${rc}"
