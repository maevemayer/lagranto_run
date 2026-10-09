#!/usr/bin/env python3
"""Stage 06 - origin contributions for the ERA5T 2026 event.

From the per-arrival lagrangian_terms_*.nc written by stage 05
($NC_RESULTS_DIR/traj_netcdfs/), build ONE origin-contribution file matching the IFS
`tanom` idea but richer. Per trajectory it stores:

  * exact ARRIVAL point : lon0, lat0, p0 (arrival_p)
  * TA                  : value at arrival + its genesis (TA_lon/lat/p_genesis,
                          TA_genesis_step) = first TA sign-change back from arrival
  * each PROCESS (advection/adiabatic/diabatic/seasonal): the cumulative
    contribution at arrival + its OWN genesis (<proc>_lon/lat/p_genesis,
    <proc>_genesis_step). A process genesis is the LAST sign-change of that
    process's cumulative contribution (closest to arrival), falling back to the
    TA genesis when it never flips -> it coincides with, or lies downstream of,
    the total-TA genesis.

Output: $NC_RESULTS_DIR/origin_contributions.nc  (dims: time, traj_id),
        copied to $RESULTS_DIR on /homedata (not purged like /scratchx)

Usage:
    python 06_origin_contributions_6h.py [N_JOBS]     # N_JOBS accepted for
    (source config.sh first so NC_RESULTS_DIR is set) # symmetry; build is serial
"""
import os
import glob
import shutil
import sys

import numpy as np
import pandas as pd
import xarray as xr

NC_RESULTS_DIR = os.environ["NC_RESULTS_DIR"]
IN_DIR = os.path.join(NC_RESULTS_DIR, "traj_netcdfs")
OUT_FILE = os.path.join(NC_RESULTS_DIR, "origin_contributions.nc")
KEEP_DIR = os.environ["RESULTS_DIR"]        # /homedata copy of the final result
DELTA_T = 6 * 3600.0
PROCS = ["advection", "adiabatic", "diabatic", "seasonal"]


# --------------------------------------------------------------------------- #
def _at(da, s):
    return da.sel(step=s).drop_vars("step", errors="ignore")


def ta_genesis_step(TA):
    """First TA sign-change going back from arrival (total-TA genesis)."""
    TA0 = TA.sel(step=0).broadcast_like(TA)
    sc = (TA * TA0 < 0).where(TA.step < 0)
    arr = sc.sortby("step", ascending=False).fillna(False).astype(int)
    return (arr.cumsum("step") == 1).idxmax("step").astype(int)


def cumulative(ds, fss):
    """Cumulative process contribution from the TA genesis to each step (K)."""
    step2 = ds.step.broadcast_like(ds["TA"])
    mask = (step2 >= fss.broadcast_like(ds["TA"])) & (step2 <= 0)
    return {p: (ds[p].where(mask) * DELTA_T).fillna(0).cumsum("step").where(mask) for p in PROCS}


def proc_genesis_step(cum_p, fss):
    """Per-process genesis: last sign-change of the cumulative contribution
    (closest to arrival); fall back to the TA genesis if it never flips."""
    sign = np.sign(cum_p)
    change = (sign.diff("step") != 0)
    last = change["step"].where(change).max("step")
    firstv = cum_p["step"].where(cum_p.notnull()).min("step")
    g = xr.where(last.isnull(), firstv, last)
    return g.where(np.isfinite(g), fss).astype(int)


def build_all(ds):
    """One arrival file -> per-trajectory Dataset: arrival + TA + per-process geneses."""
    for v in ["lon", "lat", "p"]:
        ds[v] = ds[v].where(ds[v] != -999.0)               # mask LAGRANTO sentinels
    TA = ds["TA"]
    fss = ta_genesis_step(TA)
    cum = cumulative(ds, fss)
    d = dict(
        lon0=_at(ds["lon"], 0), lat0=_at(ds["lat"], 0), p0=_at(ds["p"], 0), arrival_p=_at(ds["p"], 0),
        TA=_at(TA, 0),
        TA_lon_genesis=_at(ds["lon"], fss), TA_lat_genesis=_at(ds["lat"], fss),
        TA_p_genesis=_at(ds["p"], fss), TA_genesis_step=fss.astype("float32"),
    )
    for p in PROCS:
        gp = proc_genesis_step(cum[p], fss)
        d[p] = _at(cum[p], 0)                               # cumulative contribution at arrival
        d[f"{p}_lon_genesis"] = _at(ds["lon"], gp)
        d[f"{p}_lat_genesis"] = _at(ds["lat"], gp)
        d[f"{p}_p_genesis"] = _at(ds["p"], gp)
        d[f"{p}_genesis_step"] = gp.astype("float32")
    return xr.Dataset(d).reset_coords(drop=True)


def main():
    files = sorted(glob.glob(os.path.join(IN_DIR, "lagrangian_terms_*.nc")))
    if not files:
        sys.exit(f"No lagrangian_terms_*.nc in {IN_DIR} (run stage 05 first)")
    print(f"[06] {len(files)} arrival files from {IN_DIR}")
    parts = []
    for f in files:
        stamp = os.path.basename(f).replace("lagrangian_terms_", "").replace(".nc", "")  # YYYYMMDD_HH
        try:
            t = pd.Timestamp(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]} {stamp[9:11]}:00")
        except Exception:
            print("  skip (bad stamp):", f); continue
        ds = xr.open_dataset(f)
        parts.append(build_all(ds).expand_dims(time=[np.datetime64(t)]))
        ds.close()
        print("  done", stamp)
    pooled = xr.concat(parts, dim="time", join="outer")
    pooled.to_netcdf(OUT_FILE)
    print(f"[06] saved -> {OUT_FILE} | dims {dict(pooled.sizes)} | {len(pooled.data_vars)} vars")
    os.makedirs(KEEP_DIR, exist_ok=True)
    kept = shutil.copy2(OUT_FILE, KEEP_DIR)
    print(f"[06] copied -> {kept}")


if __name__ == "__main__":
    main()
