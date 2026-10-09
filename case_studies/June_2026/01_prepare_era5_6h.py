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
import pandas as pd
import xarray as xr
from joblib import Parallel, delayed

# ------------------------------------------------------------------ config ---
from era5_input import regrid_daily     # spiritx reader (all fields in one daily GRIB)

SCRATCH   = os.path.join(os.environ["SCRATCH_BASE"], "prepare")
FINAL     = os.environ["DATA_DIR"]
os.makedirs(SCRATCH, exist_ok=True)
os.makedirs(FINAL, exist_ok=True)

era5_vars = ["sp", "t", "u", "v", "w"]


# ------------------------------------------------------- daily regridding ----
def process_era5_daily(year, month, day, var):
    """Regrid one day of `var` to 0.5 deg, hours 00/06/12/18 only (see era5_input)."""
    regrid_file = regrid_daily(year, month, day, var, SCRATCH, "0.5deg_6h")
    if regrid_file:
        print(f"Processed {var}: {regrid_file}", flush=True)
    return regrid_file


# --------------------------------------------------- pressure & dp (as-is) ---
def compute_pressure_and_dp(sp_file, t_file, output_dir, date_str):
    """UNCHANGED from the notebook (hybrid p and dp from T-file coefficients)."""
    os.makedirs(output_dir, exist_ok=True)

    ds_t = xr.open_dataset(t_file)
    ds_sp = xr.open_dataset(sp_file)
    sp = ds_sp["sp"]

    hyam = ds_t["hyam"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})
    hybm = ds_t["hybm"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})
    hyai = ds_t["hyai"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})
    hybi = ds_t["hybi"].expand_dims({"time": sp.time, "lat": sp.lat, "lon": sp.lon})

    p = hyam + hybm * sp
    p_half = hyai + hybi * sp
    dp = p_half.isel(nhyi=slice(1, None)) - p_half.isel(nhyi=slice(None, -1))

    p = p.rename({"nhym": "lev"})
    dp = dp.rename({"nhyi": "lev"})
    p = p.assign_coords(lev=ds_t["lev"]); p.name = "p"; p.attrs["units"] = "Pa"
    dp = dp.assign_coords(lev=ds_t["lev"]); dp.name = "dp"; dp.attrs["units"] = "Pa"

    os.makedirs(f"{output_dir}/p", exist_ok=True)
    os.makedirs(f"{output_dir}/dp", exist_ok=True)
    p_file = os.path.join(output_dir, f"p/p_{date_str}_0.5deg_6h.nc")
    dp_file = os.path.join(output_dir, f"dp/dp_{date_str}_0.5deg_6h.nc")
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
        file_path = os.path.join(input_dir, var, f"{var}_{date_str}_0.5deg_6h.nc") \
            if var in ("p", "dp") else \
            os.path.join(input_dir, f"{var}_{date_str}_0.5deg_6h.nc")
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
        single.to_netcdf(out_file)
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
