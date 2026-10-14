#!/usr/bin/env python3
"""Stage 05 - trajectory process decomposition for the ERA5T 2026 event.

Framework port of 06_trajectories.ipynb, which is ALREADY 6-hourly
(delta_t = 6*3600; the loader asserts 61 steps = 360 h / 6 h). The
decomposition math is unchanged except for d(TCLIM)/dp: it is no longer a
precomputed S-file field but computed along the trajectory from TCLIM traced
10 hPa below/above the parcel (Hotz et al. 2024). Otherwise only the I/O is wired to the event's
results_6h/ folder (from config.sh), and every trajectory found there is
processed in parallel.

For each trajectory_<YYYYMMDD>_<HH> in $RESULTS_DIR (/homedata) it writes
lagrangian_terms_, integrated_terms_ and cumulative_terms_ NetCDFs into
$NC_RESULTS_DIR/traj_netcdfs/ (/scratchx).

Usage:
    python 05_traj_postprocess_6h.py [N_JOBS]
    (source config.sh first so RESULTS_DIR / NC_RESULTS_DIR are set)
"""
import os
import re
import sys
from pathlib import Path

import pandas as pd
import xarray as xr
from joblib import Parallel, delayed

RESULTS_DIR = os.environ["RESULTS_DIR"]
OUT_DIR = os.path.join(os.environ["NC_RESULTS_DIR"], "traj_netcdfs")
DELTA_T = 3600.0 * 6                                # 6-hourly (unchanged from nb)


# --------------------------------------------------------------------------- #
# Verbatim math from 06_trajectories.ipynb                                     #
# --------------------------------------------------------------------------- #
def vertical_gradient_tclim(df, dp_hpa=10.0):
    """d(TCLIM)/dp along the trajectory, as in Hotz et al. (2024):
    [TCLIM(p+10 hPa) - TCLIM(p-10 hPa)] / 20 hPa. Where one side is missing
    (p+10 hPa below the ground, or outside the data) the one-sided difference
    with TCLIM at the parcel over 10 hPa is used.
    Returned x1000 (K per 1000 Pa), the scaling the old traced DTDP field had,
    so DTDP_mid below is unchanged and still comes out in K Pa-1."""
    pa = 100.0 * dp_hpa
    up, mid, dn = df["TCLIM_M10"], df["TCLIM"], df["TCLIM_P10"]
    grad = (dn - up) / (2 * pa)
    grad = grad.fillna((mid - up) / pa).fillna((dn - mid) / pa)
    return 1000.0 * grad


def load_trajectory(date, hour, folder):
    hour_str = f"{int(hour):02d}"
    file = f"{folder}/trajectory_{date}_{hour_str}"
    # column order = tracevars.shellscript: TCLIM is also traced 10 hPa below
    # (TCLIM:+10HPA, p + 10 hPa) and above (TCLIM:-10HPA) the parcel
    float_cols = ["time", "lon", "lat", "p", "T", "OMEGA", "THETA", "TCLIM",
                  "TCLIM_P10", "TCLIM_M10", "DTDT"]
    df = pd.read_csv(file, sep=r"\s+", skiprows=4, names=float_cols, dtype=str)
    for c in float_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    sentinel_cols = float_cols[1:]
    df[sentinel_cols] = df[sentinel_cols].where(~df[sentinel_cols].eq(-999.0))
    df["DTDP"] = vertical_gradient_tclim(df)
    df["valid"] = df[["lon", "lat", "p", "T", "TCLIM"]].notna().all(axis=1)
    df["new_traj"] = df["time"].fillna(0) == 0
    df["traj_id"] = df["new_traj"].cumsum()
    df["traj_id"] = df["traj_id"].where(df["traj_id"] > 0, 1)
    df["step"] = -df.groupby("traj_id").cumcount()
    df["TA"] = df["T"] - df["TCLIM"]
    df = df.set_index(["traj_id", "step"])
    ds = df.drop(columns=["new_traj"]).to_xarray()
    # sign check of the +/-10HPA convention: T_clim increases with pressure
    # throughout the troposphere, so the gradient must be positive there
    trop = (df["p"] > 300) & (df["p"] < 900)
    if trop.any() and not df.loc[trop, "DTDP"].median() > 0:
        raise ValueError(f"{file}: dTCLIM/dp median {df.loc[trop, 'DTDP'].median():.3g} <= 0 in "
                         "300-900 hPa - check the sign convention of TCLIM:+/-10HPA in tracevars")
    assert ds.step.size == 61, (
        f"{ds.step.size} steps -> {360/(ds.step.size-1):.0f} h input "
        f"(expected 61/6 h): stale/3-h file {file}"
    )
    return ds


def compute_lagrangian_midpoints(ds, delta_t=3600.0 * 6):
    ds = ds.copy()
    ds["OMEGA_mid"] = 0.5 / 100.0 * (ds.OMEGA + ds.OMEGA.shift(step=1))
    ds["T_mid"] = 0.5 * (ds.T + ds.T.shift(step=1))
    ds["P_mid"] = 0.5 * (ds.p + ds.p.shift(step=1)) * 100
    ds["DTDT_mid"] = 0.5 / 1e6 * (ds.DTDT + ds.DTDT.shift(step=1))
    ds["DTDP_mid"] = 0.5 / 1e3 * (ds.DTDP + ds.DTDP.shift(step=1))
    ds["dp_dt_mid"] = (ds.p - ds.p.shift(step=1)) * 100.0 / delta_t
    ds["dTHETA_dt_mid"] = (ds.THETA - ds.THETA.shift(step=1)) / delta_t
    ds["dTCLIM_dt_mid"] = (ds.TCLIM - ds.TCLIM.shift(step=1)) / delta_t
    return ds


def compute_lagrangian_terms(ds, kappa=0.286):
    ds["advection"] = - ds["dTCLIM_dt_mid"] + ds["DTDT_mid"] + ds["dp_dt_mid"] * ds["DTDP_mid"]
    ds["adiabatic"] = ((kappa * ds["T_mid"] / ds["P_mid"]) - ds["DTDP_mid"]) * ds["dp_dt_mid"]
    ds["diabatic"] = ((ds["P_mid"] / 1e5) ** kappa) * ds["dTHETA_dt_mid"]
    ds["seasonal"] = - ds["DTDT_mid"]
    return ds


def integrate_terms_to_sign_change_vectorized(ds, delta_t=3600.0 * 6,
                                              terms=("advection", "adiabatic", "diabatic", "seasonal")):
    TA = ds["TA"]
    TA0 = TA.sel(step=0).broadcast_like(TA)
    sign_change = (TA * TA0 < 0)
    sign_change = sign_change.where(TA.step < 0)
    sorted_ds = sign_change.sortby("step", ascending=False)
    arr = sorted_ds.fillna(False).astype(int)
    cumsum = arr.cumsum("step")
    first_mask = cumsum == 1
    first_switch_step = first_mask.idxmax("step")
    step_2d = ds.step.broadcast_like(TA)
    step_end_broadcast = first_switch_step.broadcast_like(TA)
    mask = (step_2d <= 0) & (step_2d >= step_end_broadcast)
    integrated_data = {}
    for term in terms:
        integrated_data[term] = ("traj_id", ds[term].where(mask).sum(dim="step").values * delta_t)
    return xr.Dataset(integrated_data, coords={"traj_id": ds.traj_id})


def cumulative_terms_from_sign_change(ds, delta_t=6 * 3600.0,
                                      terms=("advection", "adiabatic", "diabatic", "seasonal")):
    TA = ds["TA"]
    TA0 = TA.sel(step=0).broadcast_like(TA)
    sign_change = (TA * TA0 < 0)
    sign_change = sign_change.where(TA.step < 0)
    sorted_ds = sign_change.sortby("step", ascending=False)
    arr = sorted_ds.fillna(False).astype(int)
    cumsum = arr.cumsum("step")
    first_mask = cumsum == 1
    first_switch_step = first_mask.idxmax("step")
    step_2d = ds.step.broadcast_like(TA)
    step_end_broadcast = first_switch_step.broadcast_like(TA)
    mask = (step_2d >= step_end_broadcast) & (step_2d <= 0)
    cum_integrated = {}
    for term in terms:
        vals = ds[term].where(mask) * delta_t
        vals_filled = vals.fillna(0)
        cum_vals = vals_filled.cumsum(dim="step")
        cum_vals = cum_vals.where(mask)
        cum_integrated[term] = cum_vals
    return cum_integrated


def process_and_save_trajectory(date_str, hour, delta_t, results_folder, out_folder):
    os.makedirs(out_folder, exist_ok=True)
    ds = load_trajectory(date_str, hour, folder=results_folder)
    ds_mid = compute_lagrangian_midpoints(ds, delta_t=delta_t)
    ds_terms = compute_lagrangian_terms(ds_mid)
    integrated_ds = integrate_terms_to_sign_change_vectorized(ds_terms, delta_t=delta_t)
    cum_ds = xr.Dataset(cumulative_terms_from_sign_change(ds_terms, delta_t=delta_t))

    out_str = f"{date_str}_{int(hour):02d}"
    ds_terms.to_netcdf(os.path.join(out_folder, f"lagrangian_terms_{out_str}.nc"))
    integrated_ds.to_netcdf(os.path.join(out_folder, f"integrated_terms_{out_str}.nc"))
    cum_ds.to_netcdf(os.path.join(out_folder, f"cumulative_terms_{out_str}.nc"))
    for d in (ds, ds_mid, ds_terms, integrated_ds, cum_ds):
        d.close()
    print(f"[OK] {out_str}", flush=True)
    return 0


# --------------------------------------------------------------------------- #
_TRAJ_RE = re.compile(r"^trajectory_(\d{8})_(\d{2})$")


def discover_tasks(results_dir):
    tasks = []
    for f in sorted(Path(results_dir).glob("trajectory_????????_??")):
        m = _TRAJ_RE.match(f.name)
        if m:
            tasks.append((m.group(1), int(m.group(2))))
    return tasks


def main():
    n_jobs = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    if not os.path.isdir(RESULTS_DIR):
        sys.exit(f"RESULTS_DIR not found: {RESULTS_DIR}")
    tasks = discover_tasks(RESULTS_DIR)
    print(f"Trajectories found: {len(tasks)} in {RESULTS_DIR}", flush=True)
    print(f"Decomposition output: {OUT_DIR}", flush=True)
    if not tasks:
        sys.exit(f"No trajectory files in {RESULTS_DIR}")
    Parallel(n_jobs=n_jobs, backend="loky", verbose=10)(
        delayed(process_and_save_trajectory)(date_str, hour, DELTA_T, RESULTS_DIR, OUT_DIR)
        for date_str, hour in tasks
    )
    print(f"Stage 05 done: {len(tasks)} trajectories decomposed.", flush=True)


if __name__ == "__main__":
    main()
