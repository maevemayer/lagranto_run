#!/usr/bin/env python3
"""
Stage 01 - prepare ERA5 P-files (CORRECTED 6-hourly snapshots).

Same math as 01_prepare_ERA5.ipynb EXCEPT:
  * -selhour,0,6,12,18 is applied BEFORE the spectral transform,
  * the 3-hourly mean (cdo timselmean(3)) is REMOVED, and
  * ERA5 is read from the spiritx all-variable daily GRIB files (era5_input.py).
The 6-hourly instantaneous snapshots are used directly, so the resulting
P-files carry the exact valid times 00/06/12/18 (no +1 h offset -> no rename).

Usage:
    python 01_prepare_era5_6h.py START_DATE END_DATE [N_JOBS]
    e.g. python 01_prepare_era5_6h.py 2026-05-15 2026-06-01 30

Paths, HOURS, LEVELS, grid come from the environment (source config.sh first).
"""
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr
from joblib import Parallel, delayed

# ------------------------------------------------------------------ config ---
from era5_input import regrid_daily, FILE_TAG   # spiritx reader (all fields in one daily GRIB)
from lagranto_io import write_lagranto_netcdf

SCRATCH   = os.path.join(os.environ["SCRATCH_BASE"], "prepare")
FINAL     = os.environ["DATA_DIR"]
os.makedirs(SCRATCH, exist_ok=True)
os.makedirs(FINAL, exist_ok=True)

era5_vars = ["sp", "t", "u", "v", "w"]
CLEANUP = os.environ.get("CLEANUP_INTERMEDIATES", "0") == "1"
SLOT_HOURS = [int(h) for h in os.environ.get("HOURS_CSV", "0,6,12,18").split(",")]


# ------------------------------------------------------- daily regridding ----
def process_era5_daily(year, month, day, var):
    """Regrid one day of `var` (REMAP_GRID), hours 00/06/12/18 only (see era5_input)."""
    regrid_file = regrid_daily(year, month, day, var, SCRATCH, FILE_TAG)
    if regrid_file:
        print(f"Processed {var}: {regrid_file}", flush=True)
    return regrid_file


# --------------------------------------------------- pressure & dp (as-is) ---
def compute_pressure_and_dp(sp_file, t_file, output_dir, date_str):
    """Same math as the notebook (hybrid p and dp from T-file coefficients).

    Changes: the coefficient dimensions are looked up instead of being
    hard-coded as nhym/nhyi (cdo 2.2.2 on Levante), so any cdo version works;
    and the levels actually present in the T file are used. With LEVELS
    (e.g. 18/137) cdo keeps all 137/138 coefficients, so model level k
    (1-based, the lev value) takes hyam[k], and its thickness is the
    difference of interfaces k and k-1. With all levels this is exactly
    the original p = hyam + hybm*sp, dp = p_half[1:] - p_half[:-1].
    """
    os.makedirs(output_dir, exist_ok=True)

    ds_t = xr.open_dataset(t_file)
    ds_sp = xr.open_dataset(sp_file)
    sp = ds_sp["sp"]

    missing = [v for v in ("hyam", "hybm", "hyai", "hybi") if v not in ds_t]
    if missing:
        raise KeyError(f"{t_file} lacks {missing} (has: {list(ds_t.variables)})")
    mid_dim, = ds_t["hyam"].dims        # nhym with cdo 2.2.2
    half_dim, = ds_t["hyai"].dims       # nhyi with cdo 2.2.2
    # positional differences below: drop any coordinate labels on these dims
    coef = ds_t[["hyam", "hybm", "hyai", "hybi"]].drop_vars([mid_dim, half_dim], errors="ignore")

    lev = ds_t["lev"].values
    k = lev.astype(int) - 1                       # 0-based model level index
    nmid = coef.sizes[mid_dim]
    if not (np.allclose(lev, k + 1) and k.min() >= 0 and k.max() < nmid
            and coef.sizes[half_dim] == nmid + 1):
        raise ValueError(f"{t_file}: lev {lev[0]}..{lev[-1]} are not model level numbers "
                         f"1..{nmid} matching hyam/hyai")

    e = {"time": sp.time, "lat": sp.lat, "lon": sp.lon}
    hyam = coef["hyam"].isel({mid_dim: k}).expand_dims(e)
    hybm = coef["hybm"].isel({mid_dim: k}).expand_dims(e)
    hyai_top = coef["hyai"].isel({half_dim: k}).expand_dims(e)       # interface above level k
    hybi_top = coef["hybi"].isel({half_dim: k}).expand_dims(e)
    hyai_bot = coef["hyai"].isel({half_dim: k + 1}).expand_dims(e)   # interface below level k
    hybi_bot = coef["hybi"].isel({half_dim: k + 1}).expand_dims(e)

    p = hyam + hybm * sp
    dp = (hyai_bot + hybi_bot * sp) - (hyai_top + hybi_top * sp)

    p = p.rename({mid_dim: "lev"}) if mid_dim != "lev" else p
    dp = dp.rename({half_dim: "lev"})
    p = p.assign_coords(lev=ds_t["lev"]); p.name = "p"; p.attrs["units"] = "Pa"
    dp = dp.assign_coords(lev=ds_t["lev"]); dp.name = "dp"; dp.attrs["units"] = "Pa"

    os.makedirs(f"{output_dir}/p", exist_ok=True)
    os.makedirs(f"{output_dir}/dp", exist_ok=True)
    p_file = os.path.join(output_dir, f"p/p_{date_str}_{FILE_TAG}.nc")
    dp_file = os.path.join(output_dir, f"dp/dp_{date_str}_{FILE_TAG}.nc")
    p.to_dataset().to_netcdf(p_file)
    dp.to_dataset().to_netcdf(dp_file)
    return (p_file, dp_file)


# ------------------------------------------------- merge & split (as-is) -----
def merge_and_split_daily_vars(date_str, var_list=None, input_dir=None, output_dir=None):
    """UNCHANGED math: merge vars, sp Pa->hPa, rename, split by time -> P-files.

    With 6-hourly snapshots the split naturally yields PYYYYMMDD_00/06/12/18
    with correct valid times (no centring offset, so no rename step needed).
    """
    if var_list is None:
        var_list = ["t", "u", "v", "w", "sp", "dp", "p"]
    if input_dir is None:
        input_dir = SCRATCH
    os.makedirs(output_dir, exist_ok=True)

    datasets = []
    for var in var_list:
        file_path = os.path.join(input_dir, var, f"{var}_{date_str}_{FILE_TAG}.nc") \
            if var in ("p", "dp") else \
            os.path.join(input_dir, f"{var}_{date_str}_{FILE_TAG}.nc")
        if os.path.exists(file_path):
            datasets.append(xr.open_dataset(file_path))
        else:
            print(f"File not found for {var}: {file_path}", flush=True)

    if not datasets:
        print("No variables found to merge!", flush=True)
        return

    merged = xr.merge(datasets)

    if "sp" in merged.data_vars:
        merged["sp"] = merged["sp"] / 100.0
        merged["sp"].attrs["units"] = "hPa"
        merged["sp"].attrs["long_name"] = "Surface Pressure"

    rename_dict = {"u": "U", "v": "V", "w": "OMEGA", "sp": "PS", "t": "T"}
    merged = merged.rename({k: v for k, v in rename_dict.items() if k in merged.data_vars})

    for t_idx in range(len(merged.time)):
        single = merged.isel(time=t_idx)
        time_value = pd.to_datetime(single.time.values)
        time_str = time_value.strftime("P%Y%m%d_%H")
        single = single.expand_dims(time=[merged.time.isel(time=t_idx).values])
        out_file = os.path.join(output_dir, time_str)
        write_lagranto_netcdf(single, out_file)      # LAGRANTO-readable layout
        single.close()
        print(f"Saved {out_file}", flush=True)

    merged.close()
    for ds in datasets:
        ds.close()


# --------------------------------------------------------------- driver ------
def process_one_day(date):
    year, month, day = date.year, date.month, date.day
    date_str = f"{year}-{month:02d}-{day:02d}"
    print(f"\n=== Processing {date_str} ===", flush=True)

    regrid_files = {}
    for var_short in era5_vars:
        f = process_era5_daily(year, month, day, var_short)
        if f:
            regrid_files[var_short] = f

    if "t" in regrid_files and "sp" in regrid_files:
        compute_pressure_and_dp(regrid_files["sp"], regrid_files["t"], SCRATCH, date_str)
    else:
        print("Skipping p/dp computation (missing t or sp).", flush=True)

    merge_and_split_daily_vars(date_str, input_dir=SCRATCH, output_dir=FINAL)

    if CLEANUP:
        cleanup_day(date, date_str)


def cleanup_day(date, date_str):
    """CLEANUP_INTERMEDIATES=1: drop the day's regridded fields and p/dp once
    all its P-files exist and were built from all 7 inputs (else keep them,
    so a rerun can resume)."""
    inputs = [os.path.join(SCRATCH, f"{v}_{date_str}_{FILE_TAG}.nc") for v in era5_vars] + \
             [os.path.join(SCRATCH, v, f"{v}_{date_str}_{FILE_TAG}.nc") for v in ("p", "dp")]
    p_files = [os.path.join(FINAL, f"P{date:%Y%m%d}_{h:02d}") for h in SLOT_HOURS]
    if not all(os.path.exists(f) for f in inputs + p_files):
        print(f"[cleanup] {date_str}: P-files or inputs incomplete - keeping intermediates", flush=True)
        return
    for f in inputs:
        os.remove(f)
    print(f"[cleanup] {date_str}: removed {len(inputs)} intermediate files", flush=True)


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: 01_prepare_era5_6h.py START_DATE END_DATE [N_JOBS]")
    start_date, end_date = sys.argv[1], sys.argv[2]
    n_jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    dates = pd.date_range(start_date, end_date, freq="D")
    Parallel(n_jobs=n_jobs, backend="loky", verbose=10)(
        delayed(process_one_day)(d) for d in dates
    )


if __name__ == "__main__":
    main()
