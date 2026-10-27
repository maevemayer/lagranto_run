#!/usr/bin/env python3
"""Are the 6hrLev time stamps right? Compare 6hrLev ua with instantaneous winds of the same model.

  (1) 6hrLev ua, interpolated in log p to the levels of 6hrPlevPt ua (instantaneous), vs
      6hrPlevPt ua at time offsets of -12..+12 h (6-hourly steps).
  (2) 6hrLev ua of the lowest model level vs 3hr uas (10 m), offsets -12..+12 h in 3-h steps
      (every offset for which the 3hr file has a matching time stamp).

Correlations are pooled over |lat| < 70 after removing each grid point's time mean. If the
6hrLev stamps are right, the correlation peaks at offset 0 h.

Usage:  python timing_check.py MODEL [--exp historical] [--start 2000-06-01] [--end 2000-06-30]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import omega_check as oc  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("model")
ap.add_argument("--exp", default="historical")
ap.add_argument("--member", default="")
ap.add_argument("--start", default="2000-06-01")
ap.add_argument("--end", default="2000-06-30")
a = ap.parse_args()
oc.DONOR_FOR = a.model
s0, s1 = a.start, a.end


def dirs(table, var):
    return oc._covering(oc.version_dirs(a.model, table, var, a.exp), s0, s1)


def times_of(da):
    return pd.DatetimeIndex(da["time"].values)


def anomaly(x):
    return x - x.mean("time")


def corr(c, r):
    m = np.isfinite(c) & np.isfinite(r)
    c, r = c[m], r[m]
    return float(np.corrcoef(c, r)[0, 1]), float(np.sum(c * r) / np.sum(r * r)), int(m.sum())


def lag_table(comp, ref, label):
    """comp, ref: anomalies (time, ...). Correlate comp(t) with ref(t + offset)."""
    tc, tr = times_of(comp), times_of(ref)
    dt_ref = np.median(np.diff(tr.values)).astype("timedelta64[m]").astype(int)
    print(f"\n{label}\n  6hrLev times {tc[0]}, {tc[1]} ...; reference times {tr[0]}, {tr[1]} ... (every {dt_ref} min)")
    print("  offset h    r       slope   n")
    rows = []
    for off in range(-12 * 60, 12 * 60 + 1, 30):
        tgt = tc + pd.Timedelta(minutes=off)
        ok = tgt.isin(tr)
        if ok.sum() < 0.8 * len(tc):
            continue
        c = comp.isel(time=np.flatnonzero(ok)).values
        r = ref.sel(time=tgt[ok].values).values
        rr, sl, n = corr(c.ravel(), r.ravel())
        rows.append((off / 60, rr))
        print(f"  {off / 60:+6.1f}   {rr:7.4f}  {sl:7.3f}  {n}")
    if rows:
        best = max(rows, key=lambda x: x[1])
        print(f"  -> best offset {best[0]:+.1f} h (r = {best[1]:.4f})")


member_filter = (lambda m: m == a.member) if a.member else (lambda m: True)
lev = {v: dirs("6hrLev", v) for v in ("ua", "ps")}
pt = dirs("6hrPlevPt", "ua")
s3 = dirs("3hr", "uas")
common = sorted(m for m in set(lev["ua"]) & set(lev["ps"]) if member_filter(m))
if not common:
    sys.exit(f"{a.model}: no 6hrLev ua/ps for {s0}..{s1}")
m = next((x for x in common if x in pt or x in s3), common[0])
print(f"{a.model} {a.exp} {m}, {s0}..{s1}; 6hrPlevPt ua: {m in pt}; 3hr uas: {m in s3}")

dsu = oc.open_files(lev["ua"][m], s0, s1 + "T23:59")
dps = oc.open_files(lev["ps"][m], s0, s1 + "T23:59")
A, B, vdim, flip = oc.half_level_coeffs(dsu, "ua")
u = dsu["ua"]
if flip:
    u = u.isel({vdim: slice(None, None, -1)})
ps = dps["ps"].sel(time=u["time"])
Af, Bf = 0.5 * (A[1:] + A[:-1]), 0.5 * (B[1:] + B[:-1])
lat_ok = np.abs(u["lat"].values) < 70
print(f"  {u.sizes[vdim]} levels, {u.sizes['lat']}x{u.sizes['lon']}, {u.sizes['time']} time steps")

if m in pt:
    dpt = oc.open_files(pt[m], s0, s1 + "T23:59")
    upt = dpt["ua"]
    pdim = next(d for d in upt.dims if d not in ("time", "lat", "lon"))
    plevs = np.asarray(upt[pdim].values, float)
    plevs_pa = plevs * (100.0 if plevs.max() < 2000 else 1.0)
    sel = [p for p in plevs_pa if 2e4 <= p <= 9.5e4]
    print(f"  6hrPlevPt ua levels (Pa): {plevs_pa.tolist()} -> using {sel}")
    if (upt.sizes["lat"], upt.sizes["lon"]) != (u.sizes["lat"], u.sizes["lon"]):
        sys.exit("6hrPlevPt grid differs from 6hrLev grid")
    out = np.full((u.sizes["time"], len(sel), int(lat_ok.sum()), u.sizes["lon"]), np.nan, "f4")
    for i in range(u.sizes["time"]):
        uu = u.isel(time=i).transpose(vdim, "lat", "lon").values[:, lat_ok]
        pp = ps.isel(time=i).values[lat_ok]
        lnp = np.log(Af[:, None, None] + Bf[:, None, None] * pp[None])       # increasing with index
        for j, pt_ in enumerate(sel):
            k = np.clip((lnp < np.log(pt_)).sum(0), 1, len(Af) - 1)            # first level below target
            l0 = np.take_along_axis(lnp, (k - 1)[None], 0)[0]
            l1 = np.take_along_axis(lnp, k[None], 0)[0]
            w = (np.log(pt_) - l0) / (l1 - l0)
            v = (1 - w) * np.take_along_axis(uu, (k - 1)[None], 0)[0] + w * np.take_along_axis(uu, k[None], 0)[0]
            out[i, j] = np.where(pp > pt_ + 1000.0, v, np.nan)               # level above ground
    comp = xr.DataArray(out, dims=("time", "plev", "y", "lon"), coords={"time": u["time"].values})
    ref = upt.sel({pdim: [p / (100.0 if plevs.max() < 2000 else 1.0) for p in sel]})
    ref = ref.isel(lat=np.flatnonzero(lat_ok)).load().transpose("time", pdim, "lat", "lon")
    ref = xr.DataArray(ref.values, dims=("time", "plev", "y", "lon"), coords={"time": ref["time"].values})
    lag_table(anomaly(comp), anomaly(ref), "(1) 6hrLev ua interpolated to pressure levels vs 6hrPlevPt ua (instantaneous)")

if m in s3:
    d3 = oc.open_files(s3[m], s0, s1 + "T23:59")
    u3 = d3["uas"]
    cm = u3.attrs.get("cell_methods", "")
    print(f"\n  3hr uas cell_methods: '{cm}'")
    if (u3.sizes["lat"], u3.sizes["lon"]) != (u.sizes["lat"], u.sizes["lon"]):
        sys.exit("3hr grid differs from 6hrLev grid")
    low = u.isel({vdim: -1}).isel(lat=np.flatnonzero(lat_ok)).load()
    r3 = u3.isel(lat=np.flatnonzero(lat_ok)).load()
    lag_table(anomaly(low.drop_vars([c for c in low.coords if c not in ("time", "lat", "lon")])),
              anomaly(r3.drop_vars([c for c in r3.coords if c not in ("time", "lat", "lon")])),
              "(2) 6hrLev ua lowest model level vs 3hr uas (10 m)")
