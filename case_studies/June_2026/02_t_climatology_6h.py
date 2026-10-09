#!/usr/bin/env python3
"""
Stage 02 - 9-year T climatology (tbar), CORRECTED 6-hourly snapshots.

Same math as 02_S_files_T_climatology.ipynb EXCEPT:
  * -selhour,0,6,12,18 applied BEFORE the spectral transform,
  * the 3-hourly mean (timselmean) is REMOVED (instantaneous snapshots),
  * slots are 00/06/12/18 instead of 01/04/.../22,
  * ERA5 is read from the spiritx all-variable daily GRIB files (era5_input.py).
The 9-year selection, the 21-day running mean, the chunked cat/timmean, and the
final multi-year timmean are all UNCHANGED.

Usage:
    python 02_t_climatology_6h.py BASE_START BASE_END [N_JOBS]
    e.g. python 02_t_climatology_6h.py 2026-05-14 2026-06-02 28

tbar files are written to $TBAR_DIR as t_YYYY-MM-DD_6hHH_tbar.nc
NOTE: choose BASE_START/END to extend +/-1 day beyond the S-file period
(stage 03 needs the +/-6 h climatology neighbours across day boundaries).
"""
import os
import sys
import glob
import pandas as pd
import xarray as xr
from cdo import Cdo
from joblib import Parallel, delayed

from era5_input import era5_file, regrid_daily   # spiritx reader (all fields in one daily GRIB)

SCRATCH_BASE = os.path.join(os.environ["SCRATCH_BASE"], "t_clim")
dump_path = os.path.join(SCRATCH_BASE, "tbar_dump")
temp_dir = os.path.join(dump_path, "tmp")
intermediate_path = os.path.join(SCRATCH_BASE, "bar_intermediate")
final_path = os.environ["TBAR_DIR"]
for p in [dump_path, intermediate_path, temp_dir, final_path]:
    os.makedirs(p, exist_ok=True)

SLOTS = ["00", "06", "12", "18"]          # was 01/04/07/10/13/16/19/22

# base period is set by main(); referenced by year-selection helpers
base_period_start = None
base_period_end = None


def process_era5_daily(year, month, day, var_num, var, ml="ml"):
    """Regrid one daily T file, 6-hourly snapshots only (no 3h mean).

    var_num/ml are kept for call compatibility; spiritx has a single ERA5 pool
    for every year, so there is no ET/E5 choice any more.
    """
    return regrid_daily(year, month, day, var, dump_path, "0.5deg_6h")


def process_extended_period(start_date, end_date, var_num=130, var="t", n_jobs=8):
    extended_range = pd.date_range(start_date - pd.Timedelta(days=11),
                                   end_date + pd.Timedelta(days=11), freq="D")
    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(process_era5_daily)(d.year, d.month, d.day, var_num, var)
        for d in extended_range
    )
    return [r for r in results if r is not None]


def split_into_6hourly_files(regrid_6h_file, var, date_str):
    """Split the daily 6-hourly file (4 steps) into one file per slot."""
    ds = xr.open_dataset(regrid_6h_file)
    output_files = []
    for t_idx, slot in enumerate(SLOTS):
        single = ds.isel(time=t_idx)
        single = single.expand_dims(time=[single.time.values])
        out_file = os.path.join(dump_path, f"{var}_{date_str}_6h{slot}.nc")
        single.to_netcdf(out_file)
        single.close()
        output_files.append(out_file)
    ds.close()
    if os.path.exists(regrid_6h_file):
        os.remove(regrid_6h_file)
    return output_files


def compute_running_mean_split(center_date, var="t", window_days=21, chunk_size=3):
    """UNCHANGED 21-day running mean per slot (chunked cat/timmean)."""
    half_window = window_days // 2
    window = pd.date_range(center_date - pd.Timedelta(days=half_window),
                           center_date + pd.Timedelta(days=half_window), freq="D")
    outputs = []
    for slot in SLOTS:
        files = [os.path.join(dump_path, f"{var}_{d.strftime('%Y-%m-%d')}_6h{slot}.nc")
                 for d in window
                 if os.path.exists(os.path.join(dump_path, f"{var}_{d.strftime('%Y-%m-%d')}_6h{slot}.nc"))]
        if len(files) < window_days:
            continue
        out_file = os.path.join(intermediate_path,
                                f"{var}_{center_date.strftime('%Y-%m-%d')}_6h{slot}_mean.nc")
        if os.path.exists(out_file):
            continue
        tmp_chunks = []
        cdo = Cdo()
        for i in range(0, len(files), chunk_size):
            sub = files[i:i + chunk_size]
            tmp = os.path.join(temp_dir, f"chunk_{i}_{center_date.strftime('%Y-%m-%d')}_{slot}.nc")
            cdo.cat(input=" ".join(sub), output=tmp)
            tmp_m = tmp.replace(".nc", "_mean.nc")
            cdo.timmean(input=tmp, output=tmp_m)
            os.remove(tmp)
            tmp_chunks.append(tmp_m)
        tmp_cat = os.path.join(temp_dir, f"tmp_{slot}_{center_date.strftime('%Y-%m-%d')}.nc")
        cdo.cat(input=" ".join(tmp_chunks), output=tmp_cat)
        cdo.timmean(input=tmp_cat, output=out_file)
        os.remove(tmp_cat)
        for f in tmp_chunks:
            os.remove(f)
        outputs.append(out_file)
    return outputs


def era5_file_exists_for_date(d, var_num=130, ml="ml"):
    return os.path.exists(era5_file(d.year, d.month, d.day))


def parse_clim_years(spec):
    """CLIM_YEARS from config.sh: "1991-2020" (range) or "2014 2016 2018" (list)."""
    spec = spec.replace(",", " ").strip()
    if "-" in spec:
        y0, y1 = (int(x) for x in spec.split("-"))
        return list(range(y0, y1 + 1))
    return [int(y) for y in spec.split()]


def choose_climatology_years(base_start, base_end, var_num=130, ml="ml"):
    """CLIM_YEARS (config.sh) if set; else UNCHANGED auto rule: the base year
    and the 8 before it, or the previous 9 if the base year is incomplete."""
    clim_years = os.environ.get("CLIM_YEARS", "").strip()
    if clim_years:
        years = parse_clim_years(clim_years)
        print(f"Climatology years from CLIM_YEARS: {years[0]}..{years[-1]} ({len(years)} years)", flush=True)
        return years
    current_year = pd.Timestamp.today().year
    extended_range = pd.date_range(base_start - pd.Timedelta(days=11),
                                   base_end + pd.Timedelta(days=11), freq="D")
    missing = [d for d in extended_range if not era5_file_exists_for_date(d, var_num=var_num, ml=ml)]
    if base_start.year == current_year and missing:
        print(f"Current year {base_start.year} incomplete ({len(missing)} missing). "
              f"Using {base_start.year-1}..{base_start.year-9}.", flush=True)
        return [base_start.year - i for i in range(1, 10)]
    return [base_start.year - i for i in range(0, 9)]


def compute_all_years_split(base_start, base_end, var="t", var_num=130, n_jobs=8):
    base_years = choose_climatology_years(base_start, base_end, var_num=var_num)
    for year in base_years:
        print(f"\n=== Processing year {year} ===", flush=True)
        start = base_start.replace(year=year)
        end = base_end.replace(year=year)
        date_range = pd.date_range(start, end, freq="D")

        existing = [os.path.join(intermediate_path, f"{var}_{d.strftime('%Y-%m-%d')}_6h{s}_mean.nc")
                    for d in date_range for s in SLOTS]
        existing = [f for f in existing if os.path.exists(f)]
        if len(existing) == len(date_range) * len(SLOTS):
            print(f"All running-mean files for {year} already exist. Skipping.", flush=True)
            continue

        for f in glob.glob(os.path.join(dump_path, f"{var}_{year}-*_0.5deg_6h.nc")):
            os.remove(f)

        daily_files = process_extended_period(start, end, var_num, var, n_jobs=n_jobs)
        print(f"Splitting {len(daily_files)} daily files...", flush=True)
        Parallel(n_jobs=n_jobs, verbose=10)(
            delayed(split_into_6hourly_files)(
                df, var,
                date_str=os.path.basename(df).split(f"{var}_")[1].split("_0.5deg_6h.nc")[0]
            ) for df in daily_files
        )
        print(f"Computing 21-day running means for {len(date_range)} days...", flush=True)
        Parallel(n_jobs=n_jobs, verbose=10)(
            delayed(compute_running_mean_split)(d, var) for d in date_range
        )
        for f in glob.glob(os.path.join(dump_path, f"{var}_*.nc")):
            os.remove(f)
        print(f"Cleaned hourly files for {year}.", flush=True)
    print("\nAll years processed.", flush=True)
    return base_years


def process_finalclim(d, var, climatology_years):
    cdo = Cdo()
    for slot in SLOTS:
        hourly_files = []
        for year in climatology_years:
            try:
                file_date = d.replace(year=year)
            except ValueError:
                continue
            fp = os.path.join(intermediate_path,
                              f"{var}_{file_date.strftime('%Y-%m-%d')}_6h{slot}_mean.nc")
            if os.path.exists(fp):
                hourly_files.append(fp)
        if not hourly_files:
            print(f"No files for {d.date()} slot {slot}, skipping", flush=True)
            continue
        out_file = os.path.join(final_path, f"{var}_{d.strftime('%Y-%m-%d')}_6h{slot}_tbar.nc")
        if os.path.exists(out_file):
            continue
        tmp_merged = os.path.join(final_path, f"tmp_merged_{d.strftime('%Y-%m-%d')}_6h{slot}.nc")
        cdo.cat(input=" ".join(hourly_files), output=tmp_merged)
        cdo.timmean(input=tmp_merged, output=out_file)
        if os.path.exists(tmp_merged):
            os.remove(tmp_merged)
        print(f"Climatology written: {out_file}", flush=True)


def compute_final_climatology_parallel(base_start, base_end, climatology_years, var="t", n_jobs=8):
    date_range = pd.date_range(base_start, base_end, freq="D")
    Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(process_finalclim)(d, var, climatology_years) for d in date_range
    )


def main():
    global base_period_start, base_period_end
    if len(sys.argv) < 3:
        sys.exit("usage: 02_t_climatology_6h.py BASE_START BASE_END [N_JOBS]")
    base_period_start = pd.Timestamp(sys.argv[1])
    base_period_end = pd.Timestamp(sys.argv[2])
    n_jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 28

    years = compute_all_years_split(base_period_start, base_period_end,
                                    var="t", var_num=130, n_jobs=n_jobs)
    compute_final_climatology_parallel(base_period_start, base_period_end,
                                       years, var="t", n_jobs=n_jobs)


if __name__ == "__main__":
    main()
