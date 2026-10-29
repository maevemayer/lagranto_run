#!/usr/bin/env python3
"""Record-breaking heatwave days of a CMIP6 model over the paper's target box.

Index as in the paper (Fig 1/5): box mean of daily-mean 2 m temperature over the target box
(-5.5..5.0E, 44.5..55.0N = the LAGRANTO arrival box), fractional cell weights x cos(lat).
Model data: day/tas of one member from the ESPRI mirror (CMIP6_POOL, default /bdd/CMIP6).

Prints
  * the warmest days (separate events, >= SEP days apart), absolute values
  * the record progression: years whose warmest day beats all earlier years (record-breaking)
  * for each candidate: is 6hrLev ua/va/ta (with ps) available for the LAGRANTO period
    (event - 5 - 15 days .. event + 6 days) and for the 9-year TCLIM window?
  * the submit command for case_studies/CMIP6 (arrivals event-5 .. event+6, as July 2022:
    13-24 July around 18 July)
and writes the daily box series to OUTDIR/<model>_<exp>_<member>_box_tas.csv.

Usage: python find_record_heatwave.py MRI-ESM2-0 [--exp historical] [--member r1i1p1f1]
                                       [--months 5,6,7,8,9] [--outdir .]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools", "omega"))
import omega_check as oc  # noqa: E402

BOX = (-5.5, 5.0, 44.5, 55.0)       # lon0, lon1, lat0, lat1 (paper target box = arrival box)
SEP = 10                            # days between separate events
BEFORE, AFTER, BACK = 5, 6, 15      # arrivals event-5 .. event+6; 15-day back-trajectories
CLIM_NYEARS, CLIM_FIRST, CLIM_LAST = 9, 1950, 2014

ap = argparse.ArgumentParser()
ap.add_argument("model")
ap.add_argument("--exp", default="historical")
ap.add_argument("--member", default="")
ap.add_argument("--months", default="5,6,7,8,9")
ap.add_argument("--top", type=int, default=10)
ap.add_argument("--outdir", default=".")
a = ap.parse_args()
months = [int(m) for m in a.months.split(",")]


def overlap(b0, b1, a0, a1):
    return np.clip(np.minimum(b1, a1) - np.maximum(b0, a0), 0, None)


def box_weights(ds):
    """Fractional overlap of each cell with BOX (lon given in -180..180, data in 0..360) x cos(lat)."""
    lat, lon = ds["lat"].values, ds["lon"].values
    lb = ds["lat_bnds"].values if "lat_bnds" in ds else np.stack([lat - np.gradient(lat) / 2, lat + np.gradient(lat) / 2], 1)
    ob = ds["lon_bnds"].values if "lon_bnds" in ds else np.stack([lon - np.gradient(lon) / 2, lon + np.gradient(lon) / 2], 1)
    lb, ob = np.sort(lb, 1), np.sort(ob, 1)
    fy = overlap(lb[:, 0], lb[:, 1], BOX[2], BOX[3]) / (lb[:, 1] - lb[:, 0])
    l0, l1 = BOX[0] % 360, BOX[1] % 360
    parts = [(l0, 360.0), (0.0, l1)] if l0 > l1 else [(l0, l1)]
    ovl = sum(overlap(ob[:, 0], ob[:, 1], p0, p1) for p0, p1 in parts)
    ovl += sum(overlap(ob[:, 0] - 360, ob[:, 1] - 360, p0, p1) for p0, p1 in parts)   # bounds beyond 360
    fx = ovl / (ob[:, 1] - ob[:, 0])
    w = fy[:, None] * fx[None, :] * np.cos(np.deg2rad(lat))[:, None]
    return xr.DataArray(w, dims=("lat", "lon"), coords={"lat": lat, "lon": lon})


# ---------------- daily box series
dirs = oc.version_dirs(a.model, "day", "tas", a.exp)
if not dirs:
    sys.exit(f"{a.model} {a.exp}: no day/tas under {oc.POOL}")
member = a.member or ("r1i1p1f1" if "r1i1p1f1" in dirs else sorted(dirs)[0])
if member not in dirs:
    sys.exit(f"{a.model} {a.exp}: member {member} has no day/tas (members: {sorted(dirs)})")
files = sorted(f for f in os.listdir(dirs[member]) if f.endswith(".nc"))
print(f"{a.model} {a.exp} {member}: day/tas, {len(files)} files in {dirs[member]}")
ds = xr.open_mfdataset([os.path.join(dirs[member], f) for f in files], combine="by_coords",
                       data_vars="minimal", coords="minimal", compat="override")
if not np.issubdtype(ds["time"].dtype, np.datetime64):
    sys.exit(f"{a.model}: non-standard calendar, not handled")
w = box_weights(ds)
sel = w > 0
w = w.where(sel, drop=True)
tas = ds["tas"].sel(lat=w["lat"], lon=w["lon"])
s = ((tas * w).sum(("lat", "lon")) / w.sum()).load() - 273.15
t = pd.DatetimeIndex(s["time"].values).normalize()
ser = pd.Series(s.values, index=t, name="tas_box_C")
os.makedirs(a.outdir, exist_ok=True)
csv = os.path.join(a.outdir, f"{a.model}_{a.exp}_{member}_box_tas.csv")
ser.to_csv(csv)
print(f"  box {BOX}: {int(sel.sum())} cells with weight; {ser.index[0].date()}..{ser.index[-1].date()}; -> {csv}")


# ---------------- candidates
def lev_ok(day):
    """6hrLev ua/va/ta of this member cover the LAGRANTO period and the TCLIM years?"""
    s0 = day - pd.Timedelta(days=BEFORE + BACK + 1)
    s1 = day + pd.Timedelta(days=AFTER + 1)
    ok = all(member in oc._covering(oc.version_dirs(a.model, "6hrLev", v, a.exp),
                                    s0.strftime("%Y-%m-%d"), s1.strftime("%Y-%m-%d")) for v in ("ua", "va", "ta"))
    y0 = min(max(day.year - CLIM_NYEARS // 2, CLIM_FIRST), CLIM_LAST - CLIM_NYEARS + 1)
    miss = [y for y in range(y0, y0 + CLIM_NYEARS)
            if member not in oc._covering(oc.version_dirs(a.model, "6hrLev", "ta", a.exp),
                                          f"{y}-{s0:%m-%d}", f"{y}-{s1:%m-%d}")]
    return ok, f"{y0}-{y0 + CLIM_NYEARS - 1}", miss


season = ser[ser.index.month.isin(months)]
top, used = [], []
for d, v in season.sort_values(ascending=False).items():
    if all(abs((d - u).days) >= SEP for u in used):
        used.append(d)
        top.append((d, v))
    if len(top) >= a.top:
        break
ann = season.groupby(season.index.year).agg(["max", "idxmax"])   # warmest day of each year
rec, best = [], -np.inf
for y, r in ann.iterrows():
    if r["max"] > best:
        best = r["max"]
        rec.append((r["idxmax"], r["max"]))
recent_rec = [r for r in rec if r[0].year >= 1900]

clim = season.groupby(season.index.dayofyear).mean().rolling(31, center=True, min_periods=1).mean()
print(f"\nWARMEST DAYS (months {months}, separate events >= {SEP} d apart), box-mean daily tas")
print("   day          tas C   anomaly vs all-years doy mean   6hrLev  TCLIM years (missing)")
for d, v in top:
    ok, cy, miss = lev_ok(d)
    print(f"  {d.date()}   {v:6.2f}      {v - clim.get(d.dayofyear, np.nan):+5.2f}                     "
          f"{'yes' if ok else 'NO '}   {cy} {miss if miss else ''}")
print("\nRECORD PROGRESSION since 1900 (warmest day of a year that beats all earlier years)")
for d, v in recent_rec[-12:]:
    print(f"  {d.date()}   {v:6.2f}")
d0, v0 = top[0]
print(f"\nAbsolute record: {d0.date()} ({v0:.2f} C).  LAGRANTO run for it:")
print(f"  cd ~/code/lagranto_run/case_studies/CMIP6 && CMIP6_MODEL={a.model} CMIP6_EXP={a.exp} "
      f"CMIP6_MEMBER={member} TRAJ_START={(d0 - pd.Timedelta(days=BEFORE)).date()} "
      f"TRAJ_END={(d0 + pd.Timedelta(days=AFTER)).date()} ./submit_cmip6.sh --test")
