#!/bin/bash
# ===========================================================================
# Read-only sanity check of the spiritx environment for the LAGRANTO pipeline.
# Resolves every [CHECK] item in case_studies/*/config.sh. Writes nothing.
#
#   bash tools/check_spiritx_env.sh [YYYYMMDD]      (default 20260525)
# ===========================================================================
DAY="${1:-20260525}"
CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../case_studies/June_2026" && pwd)"
# config.sh would mkdir the output tree; read the few values we need instead
eval "$(grep -E '^export (ERA5_DIR|LAGRANTO_ROOT|CONDA_ENV)=' "${CASE_DIR}/config.sh")"

hr() { printf '\n==== %s ====\n' "$*"; }

hr "Slurm partitions (name, cpus/node, mem/node MB, time limit, nodes)"
sinfo -o "%P %c %m %l %D" 2>&1 | column -t

hr "conda env '${CONDA_ENV}' + python packages"
if command -v conda >/dev/null; then
    source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate "${CONDA_ENV}"
fi
which python3
python3 - <<'EOF'
for m in ("cdo", "xarray", "netCDF4", "pandas", "joblib", "numpy"):
    try:
        mod = __import__(m); print(f"  ok   {m:8s} {getattr(mod, '__version__', '')}")
    except Exception as e:
        print(f"  MISS {m:8s} {e}")
EOF

hr "cdo binary (Levante used 2.2.2)"
which cdo && cdo --version 2>&1 | head -2
command -v module >/dev/null && module -t avail cdo 2>&1 | head

hr "storage: WORK_ROOT is /homedata/$USER/lagranto_run"
for d in /homedata/$USER; do
    [[ -d $d ]] && df -h "$d" | tail -1 | awk -v d="$d" '{printf "  %-24s free %s of %s\n", d, $4, $2}'
done

hr "LAGRANTO datelist interval (must be 6)"
grep -n "interval" "${LAGRANTO_ROOT}/goodies/datelist.sh" | head -5

F="${ERA5_DIR}/${DAY:0:4}/ERA5EN${DAY}.grb"
hr "ERA5 input: $F"
ls -lh "$F" || exit 1

hr "variables / levels / grid (cdo sinfon)"
cdo -s sinfon "$F" 2>&1 | head -60

hr "time steps in the file"
cdo -s showtime "$F" 2>&1 | tr ' ' '\n' | grep -v '^$' | sort -u | tr '\n' ' '; echo

hr "GRIB summary: shortName paramId typeOfLevel #levels edition"
if command -v grib_ls >/dev/null; then
    grib_ls -p shortName,paramId,typeOfLevel,edition "$F" 2>/dev/null \
        | awk 'NF==4 && $1!="shortName"' | sort | uniq -c
    grib_ls -p gridType,Ni,Nj,latitudeOfFirstGridPointInDegrees,longitudeOfFirstGridPointInDegrees,latitudeOfLastGridPointInDegrees,longitudeOfLastGridPointInDegrees,iDirectionIncrementInDegrees "$F" 2>/dev/null \
        | sed -n '2,3p'
else
    echo "  grib_ls not on PATH (conda install -c conda-forge eccodes)"
fi

hr "years available for the 9-year climatology"
ls -d "${ERA5_DIR}"/[12][0-9][0-9][0-9] 2>/dev/null | xargs -n1 basename | tr '\n' ' '; echo

cat <<'EOF'

==== what to look for ====
 * w (paramId 135, omega) on hybrid levels. If only etadot (77) is there,
   stop: the P-file OMEGA has to be derived differently (see README).
 * sp (134) or lnsp (152): either works.
 * 137 hybrid levels (or note the count), hours incl. 00/06/12/18.
 * Grid: global, or at least the whole NH. A regional domain cuts off
   15-day back-trajectories.
EOF
