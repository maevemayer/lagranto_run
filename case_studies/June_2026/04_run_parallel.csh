#!/bin/csh
# ===========================================================================
# Stage 04 - single-datetime LAGRANTO backward run (15 days).
# Adapted from 05_run_parallel-7.sh. Math/trajectory setup UNCHANGED.
# All paths (data, results, LAGRANTO, tracevars) come from config.sh.
#
# REMINDER: set  set interval = "6"  in
#   $LAGRANTO/goodies/datelist.sh   (was "3")
# so LAGRANTO steps on the 6-hourly P/S files.
#
# Usage: csh 04_run_parallel.csh YYYYMMDD_HH
# ===========================================================================
if ($#argv != 1) then
    echo "Usage: $0 YYYYMMDD_HH"
    exit 1
endif

set datetime = $argv[1]
echo "Processing $datetime ..."

set source_dir = ${DATA_DIR}      # from config.sh (via 04_run_lagranto.slurm)
set results_dir = ${RESULTS_DIR}
set tmpdir = ${source_dir}/tmp/tmp_$datetime
mkdir -p $tmpdir

# tracevars + LAGRANTO path come from config.sh (via 04_run_lagranto.slurm);
# create_startf.param is written by create_startf.sh itself on every call
\cp ${LAGRANTO_PARAM_DIR}/tracevars.shellscript $tmpdir/
cd $tmpdir

if ("${LAGRANTO_LD_LIBRARY_PATH}" != "") then
    if (! $?LD_LIBRARY_PATH) setenv LD_LIBRARY_PATH ""
    setenv LD_LIBRARY_PATH ${LAGRANTO_LD_LIBRARY_PATH}:$LD_LIBRARY_PATH
endif

setenv MODEL era5
setenv NETCDF_FORMAT CF
setenv LAGRANTO ${LAGRANTO_ROOT}

\ln -sf tracevars.shellscript tracevars

set backward_start = `${PYTHON} -c "from datetime import datetime, timedelta; dt=datetime.strptime('${datetime}','%Y%m%d_%H'); print((dt - timedelta(days=15)).strftime('%Y%m%d_%H'))"`

foreach date (`${LAGRANTO}/goodies/datelist.sh stdout -create $backward_start $datetime`)
    ln -sf ${source_dir}/P${date} ${tmpdir}/P${date}
    ln -sf ${source_dir}/S${date} ${tmpdir}/S${date}
end

set start_file = "${datetime}_startf.2"
set lonlat_file = "start_lonlat_${datetime}.txt"

# arrival box + levels from the EVENT block of config.sh
${PYTHON} ${PIPELINE_DIR}/make_start_points.py "$lonlat_file"
if ($status != 0) exit 1

${LAGRANTO}/create_startf/create_startf.sh $datetime $start_file \
"file(${lonlat_file}) @profile(${ARRIVAL_P_BOTTOM},${ARRIVAL_P_TOP},${ARRIVAL_NLEV}) @hPa" -changet
if ($status != 0) exit 1

set traj_file = "trajectory_${datetime}"
${LAGRANTO}/caltra/caltra.sh $datetime $backward_start $start_file $traj_file -j
if ($status != 0) exit 1
${LAGRANTO}/trace/trace.sh $traj_file $traj_file
if ($status != 0) exit 1

mkdir -p $results_dir
mv -f $traj_file $results_dir/
if ($status != 0) exit 1

# CLEANUP_INTERMEDIATES=1: drop this run folder (symlinks + start/param files)
# now that the trajectory is safely in $results_dir; failed runs keep theirs
if ($?CLEANUP_INTERMEDIATES) then
    if ("$CLEANUP_INTERMEDIATES" == "1") then
        cd ${source_dir}
        \rm -rf $tmpdir
    endif
endif
exit 0
