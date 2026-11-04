#!/usr/bin/env python
"""Instantaneous check of the diagnosed omega against CFsubhr wap (table of ~120 model columns, 'area: point
time: point', every 30 min, on model levels, sampled on the model's native grid by nearest neighbour).

No daily mean and no vertical interpolation enter, so neither sampling nor the time definition of a day table
can affect the scores. Each 6-hourly instant is compared with the same instant in the nearest column.
Controls in the log: surface pressure of the same columns, and a time-lag scan, which must peak at 0.
The files are read by computed record, with the time stamp of each record checked: their time axes are stored
in chunks of one record, so reading a whole axis (263000 values in IPSL-CM6A-LR) takes hours on /bdd.
With a reconstruction (--unregrid, CNRM family) the computed field is on the full Gaussian grid but the site on the
native reduced grid: the columns are then up to half a native cell apart, so the scores are a lower limit.

    python omega_check_sites.py IPSL-CM6A-LR --outdir results/IPSL-CM6A-LR/benchmark
    python omega_check_sites.py CNRM-ESM2-1 --unregrid regrid_CNRM-ESM2-1.json --outdir results/CNRM-ESM2-1/benchmark
"""
import argparse
import glob
import os
import sys
import time

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("model")
ap.add_argument("--exp", default="historical")
ap.add_argument("--start", default="2000-06-01")
ap.add_argument("--end", default="2000-08-31")
ap.add_argument("--member", default=None)
ap.add_argument("--unregrid", default=None, metavar="JSON")
ap.add_argument("--outdir", default=".")
a = ap.parse_args()
if a.unregrid:
    import omega_check_cnrm as oc
    oc.UNREGRID = oc.load_unregrid(a.unregrid)
else:
    import omega_check as oc

LEVELS = (925, 850, 700, 600, 500, 250)
HALF_BAND = 40.0                                 # hPa around each level
LAGS = (-360, -60, -30, 0, 30, 60, 360)          # minutes
NAE = (-80, 40, 20, 75)
t0 = time.time()
oc.DONOR_FOR = a.model
oc.LON_DERIV = "cd"
lev, ok = oc._lev_members(a.model, a.exp, a.start, a.end)
sd = oc._covering(oc.version_dirs(a.model, "CFsubhr", "wap", a.exp), a.start, a.end)
members = sorted(ok & set(sd))
if a.member:
    members = [m for m in members if m == a.member]
if not members:
    sys.exit(f"{a.model} {a.exp}: no member with 6hrLev and CFsubhr wap for {a.start}..{a.end}")
member = members[0]
parts = oc._common_time([oc.open_files(lev[v][member], a.start, a.end) for v in ("ua", "va", "ps")])
ds = xr.merge([parts[0], parts[1][["va"]], parts[2][["ps"]]], compat="override", join="override")
A, B, vdim, flip = oc.half_level_coeffs(ds)
try:
    oc._sh_for(np.sort(ds["lat"].values), ds["lon"].values)
    deriv = "sh"
except ValueError:
    deriv = "cd"
SITE_FILES = [netCDF4.Dataset(f) for f in sorted(glob.glob(os.path.join(sd[member], "*.nc")))
              if int(f.rsplit("_", 1)[-1][:4]) <= int(a.end[:4]) and int(f.rsplit("-", 1)[-1][:4]) >= int(a.start[:4])]
site = SITE_FILES[0]
if site["wap"].dimensions[0] != "time" or site["wap"].dimensions[-1] != "site":
    sys.exit(f"unexpected dimensions of CFsubhr wap: {site['wap'].dimensions}")
print(f"{a.model} {member} {a.start}..{a.end}: {site.dimensions['site'].size} columns, {site['wap'].shape[1]} site levels, "
      f"{ds.sizes[vdim]} 6hrLev levels, deriv {deriv}, wap cell_methods '{site['wap'].cell_methods}'", flush=True)

slat, slon = site["latitude"][:].data, ((site["longitude"][:].data + 180) % 360) - 180
lats, lons = np.sort(ds["lat"].values), ds["lon"].values
ilat = np.array([np.abs(lats - x).argmin() for x in slat])
ilon = np.array([np.abs(((lons - x + 180) % 360) - 180).argmin() for x in slon])
dlat, dlon = np.abs(lats[ilat] - slat), np.abs(((lons[ilon] - slon + 180) % 360) - 180)
print(f"  nearest column: up to {dlat.max():.2f} deg latitude, {dlon.max():.2f} deg longitude away", flush=True)
C, P, S, T = [], [], [], []
for i in range(ds.sizes["time"]):
    s = ds[["ua", "va", "ps"]].isel(time=i).load()
    terms = oc.flux_terms(s, A, B, vdim, flip, deriv)
    om, p = oc.omega_from_terms(terms, A, B, vdim, "sb", deriv)
    C.append(om.sortby("lat").transpose(vdim, "lat", "lon").values[:, ilat, ilon])
    P.append(p.sortby("lat").transpose(vdim, "lat", "lon").values[:, ilat, ilon])
    S.append(terms[2].sortby("lat").transpose("lat", "lon").values[ilat, ilon])
    T.append(s["time"].values)
    if i % 40 == 0:
        print(f"    step {i + 1}/{ds.sizes['time']} ({time.time() - t0:.0f} s)", flush=True)
C, P, S, T = np.stack(C), np.stack(P), np.stack(S), pd.DatetimeIndex(np.array(T))


def at(var, minutes=0):
    """var at the instants T + minutes, NaN where a file has no record within 5 min of it."""
    shape = (len(T),) + site[var].shape[1:]
    out = np.full(shape, np.nan, dtype="float32")
    for i, t in enumerate(T + pd.Timedelta(minutes=minutes)):
        for n in SITE_FILES:
            tv = n["time"]
            first, step = float(tv[0]), float(tv[1]) - float(tv[0])
            want = netCDF4.date2num(t.to_pydatetime(), tv.units, tv.calendar)
            k = int(round((want - first) / step))
            if 0 <= k < tv.shape[0] and abs(float(tv[k]) - want) <= 5.0 / 1440:
                out[i] = np.ma.filled(n[var][k].astype("float32"), np.nan)
                break
    return out


def anom(x):                                     # time mean removed per column and level
    return x - np.nanmean(x, axis=0, keepdims=True)


def corr(c, r):
    m = np.isfinite(c) & np.isfinite(r)
    return float(np.corrcoef(c[m], r[m])[0, 1])


Rt = at("wap")
same, flipped = corr(anom(C).ravel(), anom(Rt).ravel()), corr(anom(C).ravel(), anom(Rt[:, ::-1, :]).ravel())
order = -1 if flipped > same else 1
print(f"  level order of the sites found from the data: as computed r {same:.3f}, reversed r {flipped:.3f}")
Rt = Rt[:, ::order, :]
Ps = at("ps")
print(f"  control, surface pressure of the same columns: r {corr(anom(S).ravel(), anom(Ps).ravel()):.4f}, "
      f"mean difference {np.nanmean(S - Ps):+.1f} Pa, rms {np.sqrt(np.nanmean((S - Ps) ** 2)):.1f} Pa")
print("  time lag (reference shifted), r: " + "  ".join(
    f"{m:+d} min {corr(anom(C).ravel(), anom(at('wap', m)[:, ::order, :]).ravel()):.3f}" for m in LAGS), flush=True)

Ca, Ra = anom(C), anom(Rt)
Slat, Slo = np.broadcast_to(slat[None, None, :], C.shape), np.broadcast_to(slon[None, None, :], C.shape)
AGL = S[:, None, :] - P
regions = {"global": np.abs(Slat) < 70, "midlatitudes": (np.abs(Slat) >= 30) & (np.abs(Slat) <= 60),
           "N Atlantic-Europe": (Slo >= NAE[0]) & (Slo <= NAE[1]) & (Slat >= NAE[2]) & (Slat <= NAE[3])}
rows = []


def row(reg, mask, plev, agl):
    c, r = Ca[mask], Ra[mask]
    m = np.isfinite(c) & np.isfinite(r)
    c, r = c[m], r[m]
    if m.sum() < 50:
        return
    cc = float(np.corrcoef(c, r)[0, 1])
    rows.append(dict(region=reg, grid="6-hourly instants", plev_hPa=plev, agl_bin=agl, n=int(m.sum()), r=cc, r_time=cc,
                     std_ratio=float(c.std() / r.std()), slope=float((c * r).sum() / (r * r).sum()),
                     rmse_pct=float(np.sqrt(np.mean((c - r) ** 2)) / r.std() * 100)))


for reg, mk in regions.items():
    for lv in LEVELS:
        row(reg, mk & (np.abs(P / 100.0 - lv) <= HALF_BAND), lv, "")
    row(reg, mk & (AGL / 100.0 <= 50.0), "", "0-50")
print(f"\nINSTANTANEOUS, computed vs CFsubhr wap, {a.model} {member}\n"
      "  region              hPa      n       r   slope  std ratio  RMSE % sigma")
for x in rows:
    print(f"  {x['region']:18s} {str(x['plev_hPa'] or 'sfc-50'):>6s} {x['n']:7d} {x['r']:6.3f} {x['slope']:6.2f} {x['std_ratio']:7.2f} {x['rmse_pct']:8.0f}")
os.makedirs(a.outdir, exist_ok=True)
tag = f"{a.model}_cfsubhr" + ("_unregrid" if a.unregrid else "")
out = pd.DataFrame(rows).assign(model=a.model, member=member, reference="CFsubhr", start=a.start, end=a.end, variant=f"sb/{deriv}")
out.to_csv(os.path.join(a.outdir, f"omega_check_{tag}.csv"), index=False)
print(f"-> {os.path.join(a.outdir, f'omega_check_{tag}.csv')}  ({time.time() - t0:.0f} s)")
