#!/usr/bin/env python3
"""
Stage 03 - build S-files (THETA, TCLIM, DTDP, DTDT, PS), CORRECTED 6-hourly.

Same math as 03_S_files_final.ipynb EXCEPT (forced by the 6-hourly cadence):
  * slots are 00/06/12/18,
  * the climatology filenames are t_YYYY-MM-DD_6hHH_tbar.nc,
  * the DTDT central time difference uses +/-6 h (dt_seconds = 2*6*3600)
    instead of +/-3 h, because the climatology now lives on a 6 h grid.
THETA, the DTDP central difference in pressure, PS handling, dataset assembly
and attributes are all UNCHANGED.

Usage:
    python 03_build_s_files_6h.py START_DATE END_DATE [N_JOBS]
    e.g. python 03_build_s_files_6h.py 2026-05-15 2026-06-01 64
"""
import os
import sys
import pandas as pd
import xarray as xr
from joblib import Parallel, delayed

clim_dir   = os.environ["TBAR_DIR"]           # tbar files from stage 02
p_files    = os.environ["DATA_DIR"]           # P-files from stage 01
final_path = p_files                          # S-files land beside the P-files
SLOTS      = ["00", "06", "12", "18"]
DT_HOURS   = 6                                # was 3 (climatology cadence)


HYBRID_COEFF_VARS = ["hyam", "hybm"]      # folded-in 04_write_coords.sh


def compute_dp_p_theta(date_str, hour_str):
    """Read P-file, compute potential temperature, and pull hyam/hybm.

    The hyam/hybm read replaces the standalone 04_write_coords.sh
    (ncks -A -v hyam,hybm P S): the coefficients are copied straight from the
    P-file into the S dataset while it is written, so no second NCO pass and no
    `module load nco` are needed. Values are loaded into memory before the file
    closes (many workers open P-files concurrently). The theta math is unchanged.
    """
    filename = os.path.join(p_files, f"P{date_str}_{hour_str}")
    with xr.open_dataset(filename) as ds:
        p_aligned = ds.p.transpose("time", "lev", "lat", "lon")
        p0 = 100000.0
        R_cp = 287.0 / 1004.0
        theta = ds.T * (p0 / p_aligned) ** R_cp
        theta.name = "theta"
        theta.attrs["long_name"] = "Potential Temperature"
        theta.attrs["units"] = "K"
        coeffs = {}
        for name in HYBRID_COEFF_VARS:
            if name in ds:
                coeffs[name] = ds[name].load()
            else:
                print(f"[WARN] P{date_str}_{hour_str} has no {name}; S file will omit it", flush=True)
        return ds.PS.load(), p_aligned.load(), theta.load(), coeffs


def load_clim_with_overwritten_time(dt, hh):
    """UNCHANGED math; only the filename slot label is 6h instead of 3h."""
    date = pd.to_datetime(dt + hh, format="%Y%m%d%H")
    date_dashed = date.strftime("%Y-%m-%d")
    filename = f"{clim_dir}/t_{date_dashed}_6h{hh}_tbar.nc"
    with xr.open_dataset(filename) as ds:
        ds_new = ds.rename({"t": "t_clim"})
        ds_new = ds_new.assign_coords(time=("time", [date]))
        return ds_new


def _central_diff_Tp(Tprof, Pprof):
    """UNCHANGED 1D central difference along lev."""
    Tprof = Tprof.astype(float); Pprof = Pprof.astype(float)
    n = len(Tprof); dTdp = Tprof.copy()
    for i in range(1, n - 1):
        dTdp[i] = (Tprof[i + 1] - Tprof[i - 1]) / (Pprof[i + 1] - Pprof[i - 1])
    dTdp[0] = (Tprof[1] - Tprof[0]) / (Pprof[1] - Pprof[0])
    dTdp[-1] = (Tprof[-1] - Tprof[-2]) / (Pprof[-1] - Pprof[-2])
    return dTdp


def compute_dTdp(t_clim, p):
    """UNCHANGED: vertical pressure gradient of climatological temperature."""
    if t_clim.dims != p.dims:
        raise ValueError("t_clim and p must have identical dimensions")
    dTdp = xr.apply_ufunc(_central_diff_Tp, t_clim, p,
                          input_core_dims=[["lev"], ["lev"]],
                          output_core_dims=[["lev"]], vectorize=True)
    dTdp.name = "dTdp"
    dTdp.attrs["units"] = "K Pa-1"
    dTdp.attrs["long_name"] = "Pressure gradient of climatological temperature"
    return dTdp


def shift_datetime(dt, hh, offset_hours):
    base = pd.to_datetime(dt + hh, format="%Y%m%d%H")
    shifted = base + pd.Timedelta(hours=offset_hours)
    return shifted.strftime("%Y%m%d"), shifted.strftime("%H")


def load_three_times(dt, hh):
    """Neighbours at +/-DT_HOURS (=6 h) instead of +/-3 h."""
    dt_minus, hh_minus = shift_datetime(dt, hh, -DT_HOURS)
    dt_plus, hh_plus = shift_datetime(dt, hh, +DT_HOURS)
    t_minus = load_clim_with_overwritten_time(dt_minus, hh_minus).t_clim
    t_center = load_clim_with_overwritten_time(dt, hh).t_clim
    t_plus = load_clim_with_overwritten_time(dt_plus, hh_plus).t_clim
    return t_minus, t_center, t_plus


def compute_dtclim_dt_centered(t_minus, t_center, t_plus):
    """Central time difference; dt_seconds scales with the 6 h stencil."""
    dt_seconds = 2 * 3600 * DT_HOURS          # was 2*3600*3
    dTdt_values = (t_plus.drop("time") - t_minus.drop("time")) / dt_seconds
    dTdt = dTdt_values.assign_coords({"time": t_center.time}).transpose("time", "lev", "lat", "lon")
    dTdt.name = "dTdt"
    dTdt.attrs["units"] = "K s-1"
    dTdt.attrs["long_name"] = "Time derivative of T_clim"
    return dTdt


def save_compiled_dataset(dt, hh, output_dir):
    """Assembly of THETA/TCLIM/DTDP/DTDT/PS (+ hyam/hybm) -> SYYYYMMDD_hh."""
    PS, p, theta, coeffs = compute_dp_p_theta(dt, hh)
    t_minus, t_clim, t_plus = load_three_times(dt, hh)
    dTdp = compute_dTdp(t_clim, p)
    dTdt = compute_dtclim_dt_centered(t_minus, t_clim, t_plus)
    dTdp = dTdp.transpose("time", "lev", "lat", "lon")

    data_vars = {"THETA": theta, "TCLIM": t_clim, "DTDP": dTdp, "DTDT": dTdt, "PS": PS}
    data_vars.update(coeffs)  # hyam, hybm from the P file (folded-in stage 04)
    ds = xr.Dataset(data_vars)
    ds.PS.attrs["long_name"] = "Pressure"; ds.PS.attrs["units"] = "hPa"
    ds.THETA.attrs["long_name"] = "Potential Temperature"; ds.THETA.attrs["units"] = "K"
    ds.DTDP.attrs["long_name"] = "Pressure gradient of climatological temperature"; ds.DTDP.attrs["units"] = "K Pa-1"
    ds.DTDT.attrs["long_name"] = "Time derivative of T_clim"; ds.DTDT.attrs["units"] = "K s-1"

    out_file = f"{output_dir}/S{dt}_{hh}"
    ds.to_netcdf(out_file)
    print(f"Saved {out_file}", flush=True)
    return 0


def process_period_all_hours(base_start, base_end, output_dir, n_jobs=4, slots=None):
    if slots is None:
        slots = SLOTS
    date_range = pd.date_range(base_start, base_end, freq="D")
    tasks = [(d.strftime("%Y%m%d"), slot) for d in date_range for slot in slots]

    def safe_call(dt, slot):
        try:
            save_compiled_dataset(dt, slot, output_dir)
            return (dt, slot)
        except FileNotFoundError:
            print(f"File for {dt} slot {slot} not found. Skipping.", flush=True)
        except Exception as e:
            print(f"Error processing {dt} slot {slot}: {e}", flush=True)

    results = Parallel(n_jobs=n_jobs, verbose=10, backend="loky")(
        delayed(safe_call)(dt, slot) for dt, slot in tasks
    )
    return [r for r in results if r is not None]


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: 03_build_s_files_6h.py START_DATE END_DATE [N_JOBS]")
    base_start = pd.Timestamp(sys.argv[1])
    base_end = pd.Timestamp(sys.argv[2])
    n_jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    process_period_all_hours(base_start, base_end, p_files, n_jobs=n_jobs)
    print("Stage 03 done.", flush=True)


if __name__ == "__main__":
    main()
