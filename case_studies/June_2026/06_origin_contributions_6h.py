#!/usr/bin/env python3
"""Stage 06 - origin contributions for the ERA5T 2026 event.

From the per-arrival lagrangian_terms_*.nc written by stage 05
($NC_RESULTS_DIR/traj_netcdfs/), build ONE origin-contribution file matching the IFS
`tanom` idea but richer. Per trajectory it stores:

  * exact ARRIVAL point : lon0, lat0, p0 (arrival_p)
  * TA                  : value at arrival + its genesis (TA_lon/lat/p_genesis,
                          TA_genesis_step) = last point with the arrival sign of TA
                          (whole trajectory if TA never flips)
  * each PROCESS (advection/adiabatic/diabatic/seasonal): the cumulative
    contribution at arrival + its OWN genesis (<proc>_lon/lat/p_genesis,
    <proc>_genesis_step). A process genesis is the LAST sign-change of that
    process's cumulative contribution (closest to arrival), falling back to the
    TA genesis when it never flips -> it coincides with, or lies downstream of,
    the total-TA genesis.

GENESIS_RULE (config.sh):
  hotz  (default, paper) integration from the last point with the arrival sign of TA
        (Hotz et al. 2024, as the colleague's tanom files). The anomaly gained between the
        zero crossing and that point (on average ~0.5 K) is not attributed, so the
        process terms add up to TA(arrival) - TA(genesis).
  zero  integration from the zero crossing of TA, linearly interpolated inside the 6-h
        interval where TA changes sign: each process also gets that fraction of its
        crossing-interval term, so the terms add up to TA (residual ~0.02 K). TA genesis
        position and TA_genesis_step are interpolated the same way (step may be fractional).
        Trajectories without a sign change are integrated over the whole trajectory, as in hotz.

Output: $NC_RESULTS_DIR/origin_contributions.nc (hotz) or origin_contributions_zerocross.nc
        (zero), dims (time, traj_id), copied to $RESULTS_DIR on /homedata (not purged like
        /scratchx)

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
GENESIS_RULE = os.environ.get("GENESIS_RULE", "hotz")
if GENESIS_RULE not in ("hotz", "zero"):
    sys.exit(f"GENESIS_RULE must be 'hotz' or 'zero', not '{GENESIS_RULE}'")
OUT_FILE = os.path.join(NC_RESULTS_DIR, "origin_contributions.nc" if GENESIS_RULE == "hotz"
                        else "origin_contributions_zerocross.nc")
KEEP_DIR = os.environ["RESULTS_DIR"]        # /homedata copy of the final result
DELTA_T = 6 * 3600.0
PROCS = ["advection", "adiabatic", "diabatic", "seasonal"]


# --------------------------------------------------------------------------- #
def _at(da, s):
    return da.sel(step=s).drop_vars("step", errors="ignore")


def ta_genesis_step(TA):
    """TA genesis as in Hotz et al. (2024) and the colleague's tanom files: the
    LAST point back from arrival that still has the arrival sign of TA (one step
    after the first opposite sign). If TA never changes sign, the earliest valid
    point, so the whole trajectory is integrated (her age = 360 h)."""
    TA0 = TA.sel(step=0)
    opp = (TA * TA0 < 0) & (TA.step < 0)
    first_opp = TA.step.where(opp).max("step")
    earliest = TA.step.where(TA.notnull()).min("step")
    return xr.where(first_opp.notnull(), first_opp + 1, earliest).fillna(0).astype(int)


def crossing_fraction(TA, gen):
    """GENESIS_RULE=zero: share of the interval (gen-1 -> gen) that lies after TA = 0, from
    linear interpolation of TA; 0 where TA never changes sign (and always 0 for hotz)."""
    if GENESIS_RULE == "hotz":
        return xr.zeros_like(gen, dtype=float)
    TA0 = TA.sel(step=0)
    flips = ((TA * TA0 < 0) & (TA.step < 0)).any("step")
    ta_g = _at(TA, gen)
    ta_o = _at(TA, (gen - 1).clip(min=int(TA.step.min())))
    return (ta_g / (ta_g - ta_o)).where(flips, 0.0).clip(0.0, 1.0).fillna(0.0)


def _at_frac(da, gen, frac, lon=False):
    """Value at the interpolated genesis: frac of the way from step gen back to gen-1."""
    a = _at(da, gen)
    b = _at(da, (gen - 1).clip(min=int(da.step.min())))
    d = b - a
    if lon:
        d = (d + 180.0) % 360.0 - 180.0
    return a + (frac * d).where(frac > 0, 0.0)


def cumulative(ds, gen, frac):
    """Cumulative process contribution from the TA genesis to each step (K).
    A term at step s is the interval (s-1 -> s): only steps after the genesis, plus (rule
    zero) the fraction frac of the crossing interval, i.e. the term at the genesis step."""
    step2 = ds.step.broadcast_like(ds["TA"])
    g2, f2 = gen.broadcast_like(ds["TA"]), frac.broadcast_like(ds["TA"])
    w = xr.where(step2 > g2, 1.0, 0.0) + xr.where(step2 == g2, f2, 0.0)
    mask = (w > 0) & (step2 <= 0)
    return {p: (ds[p] * w * DELTA_T).where(mask).fillna(0).cumsum("step").where(mask) for p in PROCS}


def proc_genesis_step(cum_p, gen):
    """Per-process genesis: last sign-change of the cumulative contribution
    (closest to arrival); fall back to the TA genesis if it never flips."""
    sign = np.sign(cum_p)
    change = (sign.diff("step") != 0) & sign.shift(step=1).notnull().isel(step=slice(1, None))
    last = change["step"].where(change).max("step")
    return last.fillna(gen).astype(int)


def build_all(ds):
    """One arrival file -> per-trajectory Dataset: arrival + TA + per-process geneses."""
    for v in ["lon", "lat", "p"]:
        ds[v] = ds[v].where(ds[v] != -999.0)               # mask LAGRANTO sentinels
    TA = ds["TA"]
    gen = ta_genesis_step(TA)
    frac = crossing_fraction(TA, gen)
    cum = cumulative(ds, gen, frac)
    d = dict(
        lon0=_at(ds["lon"], 0), lat0=_at(ds["lat"], 0), p0=_at(ds["p"], 0), arrival_p=_at(ds["p"], 0),
        TA=_at(TA, 0),
        TA_lon_genesis=_at_frac(ds["lon"], gen, frac, lon=True), TA_lat_genesis=_at_frac(ds["lat"], gen, frac),
        TA_p_genesis=_at_frac(ds["p"], gen, frac), TA_genesis_step=(gen - frac).astype("float32"),
    )
    for p in PROCS:
        gp = proc_genesis_step(cum[p], gen)
        # cumulative contribution at arrival; genesis AT arrival (TA flips in the
        # first 6 h back) -> nothing to integrate, 0 K
        # (rule zero: a flip in the first 6 h back still contributes its crossing fraction)
        d[p] = _at(cum[p], 0).where((gen < 0) | (frac > 0), xr.where(TA.sel(step=0).notnull(), 0.0, np.nan))
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
    pooled.attrs["genesis_rule"] = GENESIS_RULE
    print(f"[06] genesis rule: {GENESIS_RULE}")
    pooled.to_netcdf(OUT_FILE)
    print(f"[06] saved -> {OUT_FILE} | dims {dict(pooled.sizes)} | {len(pooled.data_vars)} vars")
    os.makedirs(KEEP_DIR, exist_ok=True)
    kept = shutil.copy2(OUT_FILE, KEEP_DIR)
    print(f"[06] copied -> {kept}")


if __name__ == "__main__":
    main()
