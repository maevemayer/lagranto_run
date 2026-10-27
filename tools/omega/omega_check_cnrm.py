#!/usr/bin/env python
"""
omega_check_cnrm.py -- EXPERIMENTAL companion of omega_check.py for models whose 6hrLev output is
regridded ('gr', e.g. CNRM-CM6-1, CNRM-CM6-1-HR, CNRM-ESM2-1). Use omega_check.py (validated on
MPI-ESM1-2-HR) for everything else; this script is for diagnosing and testing fixes only.
With default options it computes omega exactly as omega_check.py 1.7.

    python omega_check_cnrm.py CNRM-CM6-1 --regrid-check        # 1-2 min: native points per latitude
    python omega_check_cnrm.py CNRM-CM6-1 --cfday --unregrid regrid_CNRM-CM6-1.json   # undo the regridding
    python omega_check_cnrm.py CNRM-CM6-1 --cfday --scan        # compare zonal filters, no filter, cd
    python omega_check_cnrm.py CNRM-CM6-1-HR --daily --unregrid regrid_CNRM-CM6-1-HR.json   # no 6hrPlev/CFday
    python omega_check_cnrm.py CNRM-CM6-1 --cfday --zonal-filter 1.0
    python omega_check_cnrm.py CNRM-CM6-1 --cfday --trunc 106   # lower spherical-harmonic truncation

Findings so far (CNRM-CM6-1, JJA 2000, daily means vs CFday on model levels):
  * full method (as omega_check.py): free troposphere r 0.94-0.98, slope ~0.98; near the surface in
    the midlatitudes noise bands at 45-60 deg (r ~0.7, std ratio up to 1.3)
  * lower truncation (T106/T85/T68) reduces the noise but removes real signal (slope 0.68 at T85)
  * filtering v.grad(ps) is harmful (removed)

Reference, chosen automatically per model:
    1. E3hrPt  wap  (instantaneous, 3-hourly, 7 pressure levels)  -> exact instants
    2. 6hrPlev wap  (6-hour means, 4 pressure levels)              -> mean of bounding instants
The member used is the first one that has 6hrLev ua/va/ps AND the reference.

omega is diagnosed from the continuity equation in hybrid coordinates
(Simmons & Burridge 1981), on full levels, then interpolated linearly in
log(p) to the reference pressure levels (NaN below ground).

Usage on Levante
    python omega_check.py MPI-ESM1-2-HR                       # default: JJA 2000
    python omega_check.py MRI-ESM2-0 --start 2000-06-01 --end 2000-08-31
    python omega_check.py MPI-ESM1-2-HR --box -10 30 40 55    # extra region (lon0 lon1 lat0 lat1)
    python omega_check.py MPI-ESM1-2-HR --compare --end 2000-06-30   # isolate causes, one month
or in Jupyter
    %run omega_check.py MPI-ESM1-2-HR
    # or:  from omega_check import run;  scores = run("MPI-ESM1-2-HR")

Daily timescale (if the same member has native daily-mean wap in table day or Eday):
    daily means of computed omega are compared with the model's true daily mean,
    and the error is split into
      method   : computed vs native instants, both averaged to daily means
      sampling : native omega sampled 6-hourly vs native true daily mean
      total    : computed (6-hourly) vs native true daily mean
    Daily means from instants use the trapezoidal rule over 00-24 UTC.

--cfday (1.7): daily means of the computed omega on the MODEL LEVELS compared with the
model's daily-mean omega on model levels (table CFday), i.e. without any interpolation to
pressure levels. Scores at the model levels nearest to 925/850/700/600/500/250 hPa (nominal,
for ps = 1000 hPa), by region and by height above ground (actual level pressure).
    python omega_check.py CNRM-CM6-1 --cfday

1.9: CNRM-CM6-1 output is regridded ('gr') from ARPEGE's reduced Gaussian grid. Interpolation
along latitude circles leaves small kinks that derivatives amplify (noise bands at 45-60 deg).
--zonal-filter C removes zonal wavenumbers m > C*(nlon/2)*cos(lat) before the derivatives.
    python omega_check.py CNRM-CM6-1 --cfday --scan     # compare C = 1.3/1.15/1.0/0.9, none, cd
(1.8 tried a lower spherical-harmonic truncation: it removed real signal and not the noise;
default is back to the full truncation, --trunc none. The 1.8 filtering of v.grad(ps) is removed.)

Defaults since 1.6: --deriv auto (spherical harmonics on Gaussian grids, i.e. spectral models
such as MPI-ESM, AWI-ESM, CNRM, MIROC, MRI, BCC; centred differences on all other grids) and
--fullevel sb. Validated for MPI-ESM1-2-HR (June 2000): r = 0.991-1.000, slope 0.998-1.001.
Map figures at 925/850/500 hPa (PNG at 200 dpi + PDF): omega_check_<model>_maps (snapshot of
the instantaneous fields: computed, native, difference; choose the time with --snapshot) and
omega_check_<model>_skillmaps (temporal r, RMSE in % of local sigma, local regression slope);
with --box also zoomed versions (_box_maps, _box_skillmaps).

Options added in 1.5 (to isolate the causes of the amplitude underestimation)
    --fullevel mean|sb   omega on full levels as the mean of the two half levels (mean)
                         or with the exact Simmons & Burridge (1981) full-level formula (sb, default)
    --deriv auto|cd|fft|sh  horizontal derivatives (auto = sh on Gaussian grids, else cd), centred
                         differences (cd), Fourier in
                         longitude (fft), or spherical harmonics (sh; Gaussian grids only, i.e.
                         spectral models; divergence and grad(ln ps) computed as in the model:
                         div(v dp) = dp*D + dB * v.grad(ps))
    --compare            compute mean/cd, sb/cd, mean/sh, sb/sh in one pass and print a
                         side-by-side table (+ omega_check_<model>_compare.csv/.png)
    4x4 coarse scores now use only points valid in both fields before averaging

Diagnostics added in 1.4
    * all reference levels between 1000 and 100 hPa (E3hrPt: 925 ... 250 hPa)
    * regions: global (|lat|<70), tropics (<23.5), subtropics (23.5-35), midlatitudes (35-65),
      both hemispheres; NH mid-latitude land
      (30-60N, land > 50 %, mean ps >= 900 hPa), high terrain (mean ps < 850 hPa), optional box
    * height above ground: scores pooled in bins of (ps - p) = 0-50, 50-100, ... hPa
    * under/overestimation: regression slope of computed on native omega (< 1: too weak),
      mean bias, and the fraction of the native strongest ascent (lowest 5 %) and strongest
      descent (highest 5 %) reproduced by the computed omega at the same points

RMSE is put in context: as a fraction of the reference standard deviation of omega
(sigma, over all points and times), and as the vertical displacement an error of that
size would cause if it persisted over one 6-hour step.

Outputs: printed score table, omega_check_<model>.csv, omega_check_<model>.png,
         omega_check_<model>_daily.png, omega_check_<model>_summary.txt
         (plain-language interpretation + a LaTeX sentence with the numbers),
         omega_check_<model>_nearsurface_6h.png / _daily.png
"""
import argparse
import datetime
import glob
import os
import sys
import time

import numpy as np
import xarray as xr

__version__ = "cnrm-1.12"   # 1.12: 6hrLev ps read from the ua files when there is no 6hrLev/ps folder (ESPRI mirror); 1.11: daily means of 03/09/15/21 UTC instants by the midpoint rule, donor coefficients kept in data order, reference dates from time bounds, time-lag check, --ref-shift; 1.6: all dataset versions searched, ua/va/ps aligned on common times; 1.7: W Europe and N Atlantic regions, CMIP6_EXTRA_POOL; 1.8: IPSL half levels, calendar-safe window matching

POOL = os.environ.get("CMIP6_POOL", "/pool/data/CMIP6/data")
# further data roots with the same directory layout (e.g. references downloaded from ESGF), ":"-separated
EXTRA_POOLS = [x for x in os.environ.get("CMIP6_EXTRA_POOL", "").split(":") if x]
# regions evaluated level by level in addition to the zonal bands (lon0, lon1, lat0, lat1)
REGIONS = {"W Europe": (-10., 20., 40., 60.),            # target region of the heat extremes
           "N Atlantic-Europe": (-80., 40., 20., 75.)}  # domain of 10-day back-trajectories
A_EARTH = 6.371e6                       # m
TOL = datetime.timedelta(minutes=30)    # time matching tolerance
DEFAULT_LEVELS = None   # None: all reference levels between 1000 and 100 hPa
NON_VERTICAL = ("time", "lat", "lon", "bnds", "nbnd", "plev")


# ============================================================ file handling
def version_dirs(model, table, var, exp="historical"):
    """{member: newest version directory} for one model/table/variable."""
    out = {}
    for root in reversed([POOL] + EXTRA_POOLS):                  # pool last -> pool wins
        for p in sorted(glob.glob(f"{root}/*/*/{model}/{exp}/*/{table}/{var}/*/[vl]*")):   # [vl]*: version dirs v2019... or "latest" (ESPRI mirror)
            out[p.split("/")[-5]] = p                            # sorted -> newest version wins
    if not out and table == "6hrLev" and var == "ps":
        # the ESPRI mirror has no 6hrLev/ps folder; CMOR stores ps (formula term) in every 6hrLev file
        out = version_dirs(model, table, "ua", exp)
        if out and (model, exp) not in _PS_NOTE:
            _PS_NOTE.add((model, exp))
            print(f"  note: no 6hrLev/ps for {model} {exp}; ps is read from the ua files")
    return out


_PS_NOTE = set()


def files_in_range(path, start, end):
    """Files whose date range in the file name overlaps [start, end]."""
    s, e = start.replace("-", ""), end.replace("-", "")
    out = []
    for f in sorted(glob.glob(f"{path}/*.nc")):
        stamp = os.path.basename(f)[:-3].rsplit("_", 1)[-1]
        if "-" not in stamp:                       # fx-like file without time range
            out.append(f)
            continue
        t0, t1 = stamp.split("-")
        if t0[:8] <= e and t1[:8] >= s:
            out.append(f)
    if not out:
        raise FileNotFoundError(f"no files in {path} overlapping {start}..{end}")
    return out


def _covers(path, start, end):
    try:
        files_in_range(path, start, end)
        return True
    except FileNotFoundError:
        return False


def _covering(dirs, start, end):
    """Members whose files overlap [start, end]. Some groups (e.g. BCC) publish each chunk of years as
    its own dataset version, so for every member the newest version that covers the period is used."""
    out = {}
    for m, p in dirs.items():
        for v in sorted(glob.glob(os.path.join(os.path.dirname(os.path.dirname(p)), "*", "[vl]*")), reverse=True):
            if _covers(v, start, end):
                out[m] = v
                break
    return out


def _common_time(parts):
    """Restrict datasets to their common time steps (a variable may miss or repeat a time step)."""
    t = parts[0]["time"].values
    for d in parts[1:]:
        t = np.intersect1d(t, d["time"].values)
    n = [d.sizes["time"] for d in parts]
    if any(k != len(t) for k in n):
        print(f"  note: ua/va/ps have {n} time steps; using the {len(t)} common ones")
    return [d.sel(time=t) for d in parts]


def _period(path):
    """First and last date stamp of the files in a directory (for error messages)."""
    st = [os.path.basename(f)[:-3].rsplit("_", 1)[-1] for f in sorted(glob.glob(f"{path}/*.nc"))]
    st = [x for x in st if "-" in x]
    return f"{st[0].split('-')[0][:8]}..{st[-1].split('-')[1][:8]}" if st else "?"


def _lev_members(model, exp, start, end):
    """6hrLev members with ua, va, ps covering [start, end]; clear error otherwise."""
    raw = {v: version_dirs(model, "6hrLev", v, exp) for v in ("ua", "va", "ps")}
    allm = set(raw["ua"]) & set(raw["va"]) & set(raw["ps"])
    if not allm:
        miss = [v for v in raw if not raw[v]]
        raise FileNotFoundError(f"{model}: no member with ua, va, ps in 6hrLev/{exp} under {POOL}"
                                + (f" (missing: {', '.join(miss)})" if miss else ""))
    lev = {v: _covering(raw[v], start, end) for v in raw}
    ok = set(lev["ua"]) & set(lev["va"]) & set(lev["ps"])
    if not ok:
        m = sorted(allm)[0]
        raise FileNotFoundError(f"{model}: no 6hrLev member covers {start}..{end}; e.g. {m}: ua "
                                f"{_period(raw['ua'][m])}, va {_period(raw['va'][m])}, ps {_period(raw['ps'][m])} "
                                f"(newest versions) -> choose --start/--end")
    return lev, ok


def open_files(path, start, end):
    """Open and concatenate along time only the time-dependent variables
    (keeps ap_bnds, b_bnds, ... as they are)."""
    files = files_in_range(path, start, end)
    ds = xr.open_mfdataset(files, combine="by_coords", data_vars="minimal",
                           coords="minimal", compat="override")
    return ds.sel(time=slice(start, end))


# ====================================================== vertical coordinate
def _vdim(da):
    return next(d for d in da.dims if d not in NON_VERTICAL)


def _fixed(ds, name):
    """Variable without any time dimension, as numpy array."""
    da = ds[name]
    if "time" in da.dims:
        da = da.isel(time=0)
    return np.asarray(da.values, float)


def _bnds_to_half(b):
    """(nlev, 2) layer bounds -> (nlev+1,) interface values (single coefficient)."""
    b = np.asarray(b, float)
    while b.ndim > 2:
        b = b[0]
    if np.allclose(b[1:, 0], b[:-1, 1]):
        return np.append(b[:, 0], b[-1, 1])
    if np.allclose(b[1:, 1], b[:-1, 0]):
        return np.append(b[:, 1], b[-1, 0])
    raise ValueError("vertical bounds are not contiguous")


def _pair_to_half(ab, bb, scale=1.0, name="ap/b"):
    """Interface coefficients from layer bounds of A (ab, scaled by `scale`) and B (bb).
    Robust to: (2, nlev) layout, column order differing between rows, and bounds stored
    surface->top. The level order of the file is kept (it must match the data)."""
    def prep(x):
        x = np.array(x, dtype=float, copy=True)          # never modify the dataset in place
        while x.ndim > 2:
            x = x[0]
        if x.shape[0] == 2 and x.shape[1] != 2:
            x = x.T
        return x
    a, b = prep(ab) * scale, prep(bb)
    if a.shape != b.shape or a.shape[1] != 2:
        raise ValueError(f"unexpected bounds shapes {a.shape}, {b.shape}")
    p = a + b * 1e5                                   # interface pressure for ps = 1000 hPa
    swap = p[:, 0] > p[:, 1]                          # put the upper interface in column 0
    a[swap], b[swap], p[swap] = a[swap][:, ::-1], b[swap][:, ::-1], p[swap][:, ::-1]
    tol = dict(rtol=1e-4, atol=1.0)                   # 1 Pa
    if np.allclose(p[1:, 0], p[:-1, 1], **tol):       # file ordered top -> surface
        A, B = np.append(a[0, 0], a[:, 1]), np.append(b[0, 0], b[:, 1])
    elif np.allclose(p[1:, 1], p[:-1, 0], **tol):     # file ordered surface -> top
        A, B = np.append(a[0, 1], a[:, 0]), np.append(b[0, 1], b[:, 0])
    else:
        gaps = np.abs(p[1:, 0] - p[:-1, 1])
        k = int(np.argmax(gaps))
        raise ValueError(
            f"vertical bounds ({name}) are not contiguous: largest gap {gaps[k]:.1f} Pa between "
            f"layers {k} and {k + 1}.\n  first layers (p_upper, p_lower) [Pa]: {np.round(p[:3], 2).tolist()}"
            f"\n  last layers: {np.round(p[-3:], 2).tolist()}"
            f"\n  A bounds first/last: {np.round(a[:2], 4).tolist()} / {np.round(a[-2:], 4).tolist()}"
            f"\n  B bounds first/last: {np.round(b[:2], 6).tolist()} / {np.round(b[-2:], 6).tolist()}")
    return A, B


def _unique_interfaces(ab, bb, scale=1.0):
    """Rebuild interface coefficients from bounds whose entries are complete but scrambled
    (e.g. written with the wrong array ordering, as in CNRM files): every interior
    interface occurs twice, surface and top once. Returns (A, B) ordered top -> surface."""
    a = np.asarray(ab, float).ravel() * scale
    b = np.asarray(bb, float).ravel()
    p = a + b * 1e5
    order = np.argsort(p, kind="stable")
    a, b, p = a[order], b[order], p[order]
    keep = np.append(True, np.diff(p) > 1e-3 * np.maximum(1.0, np.abs(p[1:])) * 1e-3 + 1e-6)
    A, B = a[keep], b[keep]
    return A, B


# models whose files carry invalid vertical coefficients -> model with the identical vertical grid
VCOORD_DONOR = {"CNRM-CM6-1-HR": "CNRM-CM6-1"}
DONOR_FOR = None          # set to the model name in main()
_DONOR_CACHE = {}


def _donor_coeffs(donor, ds, vdim, n):
    """A, B (data order of the donor) from the first 6hrLev ua file of `donor`, if its level coordinate
    and bounds are identical to those of ds."""
    if donor in _DONOR_CACHE:
        dds = _DONOR_CACHE[donor]
    else:
        files = sorted(glob.glob(f"{POOL}/*/*/{donor}/historical/*/6hrLev/ua/*/[vl]*/*.nc"))
        if not files:
            return None
        dds = xr.open_dataset(files[0], decode_times=False)
        _DONOR_CACHE[donor] = dds
    if vdim not in ds.variables or vdim not in dds.variables or dds.sizes.get(vdim) != n:
        return None
    if not np.allclose(ds[vdim].values, dds[vdim].values, rtol=1e-6, atol=0):
        return None
    bn = ds[vdim].attrs.get("bounds")
    if bn and bn in ds.variables and bn in dds.variables and \
            not np.allclose(ds[bn].values, dds[bn].values, rtol=1e-6, atol=1e-12):
        return None
    global DONOR_FOR
    keep, DONOR_FOR = DONOR_FOR, None              # avoid recursion
    try:
        A, B, _, flip = half_level_coeffs(dds)
    finally:
        DONOR_FOR = keep
    return (A[::-1].copy(), B[::-1].copy()) if flip else (A, B)   # back to the donor's data order


def half_level_coeffs(ds, var="ua"):
    """A_half [Pa], B_half [-] (length nlev+1, ordered top -> surface),
    vertical dimension name, and whether model levels must be reversed."""
    vdim = _vdim(ds[var])
    n = ds.sizes[vdim]
    has = lambda *names: all(nm in ds.variables for nm in names)
    p0 = float(_fixed(ds, "p0")) if has("p0") else 1e5
    A = B = None
    donor_used = False
    # files with invalid coefficients (all zero / denormal b): use the donor model right away
    if DONOR_FOR in VCOORD_DONOR and has("b") and not (np.nanmax(np.abs(_fixed(ds, "b"))) > 0.5):
        dA = _donor_coeffs(VCOORD_DONOR[DONOR_FOR], ds, vdim, n)
        if dA is not None:
            A, B = dA
            print(f"  note: vertical coefficients in the files are invalid; taken from "
                  f"{VCOORD_DONOR[DONOR_FOR]} (identical level coordinate and bounds)")
            flip = bool(A[0] + B[0] * 1e5 > A[-1] + B[-1] * 1e5)
            return (A[::-1].copy(), B[::-1].copy(), vdim, flip) if flip else (A, B, vdim, flip)

    bpair = (("ap_bnds", "b_bnds", 1.0) if has("ap_bnds", "b_bnds") else
             ("a_bnds", "b_bnds", p0) if has("a_bnds", "b_bnds") else None)
    # coefficients stored directly on the n+1 half levels (IPSL, dimension 'klevp1') take precedence:
    # there ap_bnds/b_bnds are bounds of the half levels (n+1 rows), not of the full levels
    for an, bn in (("ap", "b"), ("a", "b")):
        if has(an, bn) and ds[an].size == n + 1 and ds[bn].size == n + 1:
            A, B = _fixed(ds, an).ravel() * (p0 if an == "a" else 1.0), _fixed(ds, bn).ravel()
            bpair = None
            break
    if bpair is not None:
        an, bn, sc = bpair
        try:
            A, B = _pair_to_half(_fixed(ds, an), _fixed(ds, bn), sc, f"{an}/{bn}")
        except ValueError as err:
            A, B = _unique_interfaces(_fixed(ds, an), _fixed(ds, bn), sc)     # top -> surface
            if len(A) != n + 1:
                raise ValueError(f"{err}\n  rebuilding from unique interface values gave {len(A)} "
                                 f"interfaces for {n} levels") from None
            # level order of the data from the full-level coefficients (or the lev coordinate)
            pf = None
            if has("ap", "b") and ds["ap"].size == n:
                pf = _fixed(ds, "ap").ravel() + _fixed(ds, "b").ravel() * 1e5
            elif has("a", "b") and ds["a"].size == n:
                pf = _fixed(ds, "a").ravel() * p0 + _fixed(ds, "b").ravel() * 1e5
            elif vdim in ds.variables:
                pf = np.asarray(ds[vdim].values, float)
            if pf is None:
                raise ValueError(f"{err}\n  cannot determine level order to rebuild bounds") from None
            if pf[0] > pf[-1]:                                                # data surface -> top
                A, B = A[::-1], B[::-1]
            if has("ap", "b") or has("a", "b"):                               # sanity check
                ph = A + B * 1e5
                lo, hi = np.minimum(ph[:-1], ph[1:]), np.maximum(ph[:-1], ph[1:])
                bad = np.sum((pf < lo - 1.0) | (pf > hi + 1.0))
                if bad:
                    raise ValueError(f"{err}\n  rebuilt interfaces inconsistent with full levels "
                                     f"at {bad} levels") from None
            print(f"  note: {an}/{bn} in the files are scrambled; interfaces rebuilt from the "
                  f"{len(A)} unique values and checked against the full-level coefficients")
    elif A is None:  # half-level coefficients stored directly (e.g. IPSL on 'klevp1')
        for an, bn in (("ap", "b"), ("a", "b")):
            if has(an, bn) and ds[an].size == n + 1 and ds[bn].size == n + 1:
                scale = p0 if an == "a" else 1.0
                A, B = _fixed(ds, an).ravel() * scale, _fixed(ds, bn).ravel()
                break

    if A is None and vdim in ds.variables:     # pure sigma (e.g. FGOALS-g3)
        sn = ds[vdim].attrs.get("standard_name", "")
        bname = ds[vdim].attrs.get("bounds", f"{vdim}_bnds")
        if sn == "atmosphere_sigma_coordinate" and has("ptop", bname):
            S = _bnds_to_half(_fixed(ds, bname))
            ptop = float(_fixed(ds, "ptop"))
            A, B = ptop * (1.0 - S), S

    if A is None:
        raise NotImplementedError(
            f"cannot build half-level pressure for vertical dim {vdim!r} "
            f"(hybrid-height models such as ACCESS-CM2 are not supported); "
            f"variables: {sorted(ds.variables)}")
    if len(A) != n + 1:
        raise ValueError(f"{len(A)} half-level coefficients for {n} levels")

    # implausible coefficients (e.g. CNRM-CM6-1-HR: ap/b/ap_bnds/b_bnds stored as zeros or denormals):
    # take them from a model with the identical vertical grid (same level coordinate and bounds)
    ph_test = A + B * 1e5
    if not (np.all(np.isfinite(ph_test)) and np.ptp(ph_test) > 5e4 and np.all(np.abs(np.diff(ph_test)) > 0)):
        donor = VCOORD_DONOR.get(DONOR_FOR, None) if DONOR_FOR else None
        dA = _donor_coeffs(donor, ds, vdim, n) if donor else None
        if dA is None:
            raise ValueError(f"implausible vertical coefficients (half-level pressure range "
                             f"{np.nanmin(ph_test):.3g}..{np.nanmax(ph_test):.3g} Pa) and no donor model")
        A, B = dA
        donor_used = True
        print(f"  note: vertical coefficients in the files are invalid; taken from {donor} "
              f"(identical level coordinate and bounds)")

    # the coefficients must follow the order of the data. Most reliable: the full-level coefficients,
    # which are attached to the data dimension (CESM2 stores the bounds in the opposite order)
    pf_full = None
    if has("ap", "b") and ds["ap"].size == n and ds["b"].size == n:
        pf_full = _fixed(ds, "ap").ravel() + _fixed(ds, "b").ravel() * 1e5
    elif has("a", "b") and ds["a"].size == n and ds["b"].size == n:
        pf_full = _fixed(ds, "a").ravel() * p0 + _fixed(ds, "b").ravel() * 1e5
    # (not for donor coefficients: they are already in data order, and the file's own
    #  full-level coefficients are as corrupt as its bounds - CNRM-CM6-1-HR)
    pf_ok = (pf_full is not None and np.all(np.isfinite(pf_full)) and np.ptp(pf_full) > 5e4
             and (np.all(np.diff(pf_full) > 0) or np.all(np.diff(pf_full) < 0)))
    if donor_used:
        pass
    elif pf_ok:
        ph = A + B * 1e5
        if (pf_full[-1] > pf_full[0]) != (ph[-1] > ph[0]):
            print("  note: vertical coefficient bounds are stored in the opposite order to the data "
                  "levels; reversed to follow the full-level coefficients")
            A, B = A[::-1].copy(), B[::-1].copy()
    # otherwise: check against the level coordinate
    # (hybrid sigma value ~ p/p0 or a pressure: increasing = top -> surface)
    elif vdim in ds.variables and ds[vdim].size == n and n > 2:
        lv = np.asarray(ds[vdim].values, float)
        sn = ds[vdim].attrs.get("standard_name", "")
        if np.all(np.diff(lv) > 0) or np.all(np.diff(lv) < 0):
            data_down = lv[-1] > lv[0]                       # data ordered top -> surface
            if "height" in sn or "altitude" in sn:
                data_down = not data_down
            coef_down = A[-1] + B[-1] * 1e5 > A[0] + B[0] * 1e5
            if ("sigma" in sn or "pressure" in sn) and data_down != coef_down:
                print(f"  note: vertical coefficients are stored in the opposite order to '{vdim}' "
                      f"({sn}); reversed to follow the data")
                A, B = A[::-1].copy(), B[::-1].copy()

    flip = bool(A[0] + B[0] * 1e5 > A[-1] + B[-1] * 1e5)
    if flip:
        A, B = A[::-1], B[::-1]
    return A, B, vdim, flip


# ==================================================== spherical operators
LON_DERIV = "cd"     # "cd": centred differences, "fft": spectral (Fourier) derivative


def _ddlon(f):
    """d/dlambda per radian (regular longitudes assumed)."""
    if LON_DERIV == "fft":
        ax = f.get_axis_num("lon")
        n = f.sizes["lon"]
        F = np.fft.rfft(np.asarray(f.values), axis=ax)
        k = np.fft.rfftfreq(n, d=1.0 / n)          # integer zonal wavenumbers
        if n % 2 == 0:
            k[-1] = 0.0                            # no derivative at Nyquist
        shape = [1] * F.ndim
        shape[ax] = k.size
        d = np.fft.irfft(1j * k.reshape(shape) * F, n=n, axis=ax)
        return f.copy(data=d)
    dlon = np.deg2rad(np.mean(np.diff(f["lon"].values)))
    return (f.roll(lon=-1, roll_coords=False) - f.roll(lon=1, roll_coords=False)) / (2 * dlon)


def _ddlat(f):
    return f.differentiate("lat") * (180.0 / np.pi)


def _div(fx, fy):
    c = np.cos(np.deg2rad(fx["lat"]))
    return (_ddlon(fx) + _ddlat(fy * c)) / (A_EARTH * c)


# ================================================ spherical harmonics
class SH:
    """Spherical-harmonic divergence and gradient on a global Gaussian grid
    (triangular truncation N = min(nlat-1, nlon//2-1)); latitudes ascending."""
    def __init__(self, lat, lon):
        lat = np.asarray(lat, float); lon = np.asarray(lon, float)
        self.nlat, self.nlon = lat.size, lon.size
        mu, w = np.polynomial.legendre.leggauss(self.nlat)
        if not np.all(np.diff(lat) > 0):
            raise ValueError("latitudes must be ascending")
        if np.max(np.abs(np.sin(np.deg2rad(lat)) - mu)) > 1e-4:
            raise ValueError("not a Gaussian grid")
        dl = np.diff(lon)
        if not np.allclose(dl, 360.0 / self.nlon, atol=1e-3):
            raise ValueError("longitudes not regular and global")
        self.mu, self.w = mu, w
        self.cos = np.sqrt(1 - mu ** 2)
        self.N = N = min(self.nlat - 1, self.nlon // 2 - 1)
        self.P, self.H = [], []
        s = self.cos
        pmm = np.full_like(mu, np.sqrt(0.5))
        for m in range(N + 1):
            if m > 0:
                pmm = -np.sqrt((2 * m + 1) / (2 * m)) * s * pmm
            nn = N + 2 - m                     # n = m .. N+1
            P = np.zeros((nn, self.nlat))
            P[0] = pmm
            if nn > 1:
                P[1] = np.sqrt(2 * m + 3) * mu * pmm
            for k in range(2, nn):
                n = m + k
                a = np.sqrt((4 * n * n - 1) / (n * n - m * m))
                a1 = np.sqrt((4 * (n - 1) ** 2 - 1) / ((n - 1) ** 2 - m * m))
                P[k] = a * (mu * P[k - 1] - P[k - 2] / a1)
            n = np.arange(m, N + 1)
            eps = lambda nv: np.sqrt(np.maximum(nv * nv - m * m, 0) / (4.0 * nv * nv - 1))
            Pm1 = np.vstack([np.zeros((1, self.nlat)), P[:-2]])          # P_{n-1}
            H = (-n * eps(n + 1))[:, None] * P[1:] + ((n + 1) * eps(n))[:, None] * Pm1
            self.P.append(P[:-1])               # n = m .. N
            self.H.append(H)

    def _four(self, f):
        return np.fft.rfft(f, axis=-1) / self.nlon

    def _synth(self, G):
        return np.fft.irfft(G * self.nlon, n=self.nlon, axis=-1)

    def _nt(self, nt):
        return self.N if nt is None else int(min(nt, self.N))

    def divergence(self, u, v, nt=None):
        """horizontal divergence, triangular truncation at nt (default: N)."""
        nt = self._nt(nt)
        U = self._four(u * self.cos[:, None]); V = self._four(v * self.cos[:, None])
        wt = self.w / (1 - self.mu ** 2)
        G = np.zeros(U.shape, complex)
        for m in range(nt + 1):
            k = nt - m + 1
            P, H = self.P[m][:k], self.H[m][:k]
            coef = (1j * m * U[..., m] * wt) @ P.T - (V[..., m] * wt) @ H.T
            G[..., m] = coef @ P
        return self._synth(G) / A_EARTH

    def gradient(self, s, nt=None):
        nt = self._nt(nt)
        S = self._four(s)
        Gx = np.zeros(S.shape, complex); Gy = np.zeros(S.shape, complex)
        for m in range(nt + 1):
            k = nt - m + 1
            P, H = self.P[m][:k], self.H[m][:k]
            coef = (S[..., m] * self.w) @ P.T
            Gx[..., m] = 1j * m * (coef @ P)
            Gy[..., m] = coef @ H
        c = self.cos[:, None] * A_EARTH
        return self._synth(Gx) / c, self._synth(Gy) / c

    def filter(self, s, nt=None):
        """scalar field truncated to total wavenumber nt (analysis + synthesis)."""
        nt = self._nt(nt)
        S = self._four(s)
        G = np.zeros(S.shape, complex)
        for m in range(nt + 1):
            P = self.P[m][:nt - m + 1]
            G[..., m] = ((S[..., m] * self.w) @ P.T) @ P
        return self._synth(G)

    def spectrum(self, s):
        """power per total wavenumber n = 0..N (mean square contribution), leading dims kept."""
        S = self._four(s)
        E = np.zeros(S.shape[:-2] + (self.N + 1,))
        for m in range(self.N + 1):
            c = (S[..., m] * self.w) @ self.P[m].T          # n = m..N
            E[..., m:] += (1.0 if m == 0 else 2.0) * np.abs(c) ** 2
        return E / 2.0


def quadratic_truncation(nlat):
    """alias-free (quadratic) triangular truncation for a Gaussian grid with nlat latitudes:
    nlat >= (3N+1)/2. Equals the model truncation on quadratic grids (T63 96x192, T85 128x256,
    T106 160x320, T127 192x384) and is lower than it on linear grids (TL127 128x256 -> 85,
    TL159 160x320 -> 106, TL359 360x720 -> 239)."""
    return (2 * nlat - 1) // 3


SH_TRUNC = "none"      # "none" (N = min(nlat-1, nlon/2-1)), "auto" (quadratic) or an integer
ZONAL_FILTER = None    # None or factor c: remove zonal wavenumbers m > c * (nlon/2) * cos(lat) from ua, va, ps


def _trunc_for(sh, trunc=None):
    t = SH_TRUNC if trunc is None else trunc
    if t in ("auto", None):
        return min(quadratic_truncation(sh.nlat), sh.N)
    if t == "none":
        return sh.N
    return min(int(t), sh.N)


_SH_CACHE = {}


def _sh_for(lat, lon):
    key = (lat.size, lon.size, float(lat[0]), float(lon[0]))
    if key not in _SH_CACHE:
        _SH_CACHE[key] = SH(lat, lon)
    return _SH_CACHE[key]


# ===================================================== omega computation
def zonal_filter(x, c):
    """Remove zonal wavenumbers m > c * (nlon/2) * cos(lat) on every latitude circle.
    For data regridded from a reduced Gaussian grid (fewer points towards the poles) to a
    regular grid: wavenumbers the native latitude circle cannot carry are interpolation
    artefacts, which horizontal derivatives would amplify."""
    if c is None:
        return x
    nlon = x.sizes["lon"]
    mcut = c * (nlon / 2) * np.cos(np.deg2rad(x["lat"].values))
    m = np.arange(nlon // 2 + 1)
    keep = (m[None, :] <= mcut[:, None]).astype(float)                 # (lat, m)
    ax = x.get_axis_num("lon")
    F = np.fft.rfft(np.asarray(x.values, dtype=np.float64), axis=ax)
    shape = [1] * F.ndim
    shape[x.get_axis_num("lat")], shape[ax] = keep.shape
    F *= keep.reshape(shape)
    return x.copy(data=np.fft.irfft(F, n=nlon, axis=ax).astype(x.dtype))



# ======================================== undoing regridding from a reduced grid
UNREGRID = None        # None, or {lat_index_in_ascending_order: (scheme, n, offset)} from --regrid-check


def _interp_matrix(scheme, n, nlon, offset):
    """Matrix W (nlon x n) mapping values on a latitude circle with n points at
    lon_j = (j + offset) * 360/n to the regular circle lon_i = i * 360/nlon."""
    th = np.arange(nlon) * 360.0 / nlon
    d = 360.0 / n
    pos = th / d - offset                                   # position in source index units
    W = np.zeros((nlon, n))
    if scheme == "linear":
        j0 = np.floor(pos).astype(int)
        w = pos - j0
        W[np.arange(nlon), j0 % n] += 1 - w
        W[np.arange(nlon), (j0 + 1) % n] += w
    elif scheme == "nearest":
        W[np.arange(nlon), np.round(pos).astype(int) % n] = 1.0
    elif scheme == "conservative":
        h = 0.5 * (360.0 / nlon) / d                        # half target cell in source units
        for i in range(nlon):
            lo, hi = pos[i] - h, pos[i] + h
            for j in range(int(np.floor(lo - 0.5)) - 1, int(np.ceil(hi + 0.5)) + 1):
                ov = min(hi, j + 0.5) - max(lo, j - 0.5)
                if ov > 0:
                    W[i, j % n] += ov / (hi - lo)
    else:
        raise ValueError(scheme)
    return W


_PINV = {}


def _pinv(scheme, n, nlon, offset):
    key = (scheme, n, nlon, offset)
    if key not in _PINV:
        W = _interp_matrix(scheme, n, nlon, offset)
        _PINV[key] = (W, np.linalg.pinv(W))
    return _PINV[key]


def detect_reduced_grid(x, schemes=("conservative", "linear", "nearest"), offsets=(0.0, 0.5),
                        nmin=8, thr=1e-5, sharp=30.0, nmargin=4):
    """x: (samples, lat, lon) on a regular grid (lat ascending). For every latitude row, test
    interpolation from n equally spaced points (all n < nlon, schemes, offsets) and keep the best
    fit per n. A row counts as regridded from n* points only if the fit at n* is (a) accurate
    (relative rms residual < thr, i.e. at the float32 noise level) and (b) SHARP: at least `sharp`
    times better than at n*-1, and (c) the SMALLEST such n. (b) rejects false positives with n close
    to nlon, where interpolation is nearly the identity and fits any smooth row approximately
    (cnrm-1.1 accepted such rows at 23-34 deg); (c) rejects larger n that also fit staircase rows
    exactly near the poles (cnrm-1.2 picked e.g. 248 instead of 108 points at 73.5 deg)."""
    S, nlat, nlon = x.shape
    X = x.astype(np.float64)
    Xa = X - X.mean(2, keepdims=True)
    den = np.sqrt((Xa ** 2).mean((0, 2))) + 1e-30                      # per row
    ns = np.arange(nmin, nlon)
    best_r = np.full((nlat, ns.size), np.inf)
    best_c = [[None] * ns.size for _ in range(nlat)]
    for i, n in enumerate(ns):
        for sc in schemes:
            for of in offsets:
                W = _interp_matrix(sc, n, nlon, of)
                Wp = np.linalg.pinv(W)
                R = X - (X @ Wp.T) @ W.T
                r = np.sqrt((R ** 2).mean((0, 2))) / den
                better = r < best_r[:, i]
                best_r[better, i] = r[better]
                for j in np.where(better)[0]:
                    best_c[j][i] = (sc, int(n), of)
    out = []
    for j in range(nlat):
        hit = None
        for i in range(1, ns.size):                     # smallest n that fits and n-1 clearly does not
            if ns[i] > nlon - nmargin:                  # 253-255 of 256: indistinguishable from identity
                break
            r0 = best_r[j, i]
            if r0 < thr and best_r[j, i - 1] / max(r0, 1e-30) > sharp:
                hit = i
                break
        if hit is not None:
            sc, n, of = best_c[j][hit]
            out.append(dict(row=j, scheme=sc, n=n, offset=of, residual=float(best_r[j, hit]),
                            sharpness=float(best_r[j, hit - 1] / max(best_r[j, hit], 1e-30))))
        else:
            i = int(np.argmin(best_r[j]))
            out.append(dict(row=j, scheme="none", n=nlon, offset=0.0, residual=float(best_r[j, i]),
                            sharpness=float(best_r[j, i - 1] / max(best_r[j, i], 1e-30)) if i > 0 else 1.0))
    return out


def unregrid(x, table=None):
    """Replace every latitude row that was interpolated from a reduced circle by the exact
    band-limited (Fourier) interpolation of the recovered native values. x: DataArray with lat
    (ascending) and lon; table: {row: (scheme, n, offset)}."""
    table = UNREGRID if table is None else table
    if not table:
        return x
    xs = x.sortby("lat")
    nlon = xs.sizes["lon"]
    ax_lat, ax_lon = xs.get_axis_num("lat"), xs.get_axis_num("lon")
    v = np.moveaxis(np.asarray(xs.values, dtype=np.float64), (ax_lat, ax_lon), (-2, -1)).copy()
    for row, (sc, n, of) in table.items():
        W, Wp = _pinv(sc, n, nlon, of)
        a = v[..., row, :] @ Wp.T                           # native values on n points
        A = np.fft.rfft(a, axis=-1) / n
        if n % 2 == 0:
            A[..., -1] *= 0.5                               # split Nyquist term
        F = np.zeros(a.shape[:-1] + (nlon // 2 + 1,), complex)
        F[..., :A.shape[-1]] = A
        k = np.arange(A.shape[-1])                          # source longitudes start at offset*d
        F[..., :A.shape[-1]] *= np.exp(-2j * np.pi * k * of / n)
        v[..., row, :] = np.fft.irfft(F * nlon, n=nlon, axis=-1)
    v = np.moveaxis(v, (-2, -1), (ax_lat, ax_lon))
    return xs.copy(data=v.astype(x.dtype))


def reregrid(x, table=None):
    """Apply the model's own output regridding to a field on the regular grid: sample its band-limited
    Fourier series at the native reduced-grid longitudes of each row and remap with the detected
    scheme (inverse of unregrid). Used to compare a computed field with native CNRM output, which
    went through this (smoothing) remapping."""
    table = UNREGRID if table is None else table
    if not table:
        return x
    xs = x.sortby("lat")
    nlon = xs.sizes["lon"]
    ax_lat, ax_lon = xs.get_axis_num("lat"), xs.get_axis_num("lon")
    v = np.moveaxis(np.asarray(xs.values, dtype=np.float64), (ax_lat, ax_lon), (-2, -1)).copy()
    F = np.fft.rfft(v, axis=-1) / nlon
    for row, (sc, n, of) in table.items():
        W, _ = _pinv(sc, n, nlon, of)
        M = n // 2 + 1
        c = np.full(M, 2.0); c[0] = 1.0
        if n % 2 == 0:
            c[-1] = 1.0
        th = 2 * np.pi * (np.arange(n) + of) / n
        E = c[:, None] * np.exp(1j * np.arange(M)[:, None] * th[None, :])      # (M, n)
        a = np.real(F[..., row, :M] @ E)                                          # values on native points
        v[..., row, :] = a @ W.T
    v = np.moveaxis(v, (-2, -1), (ax_lat, ax_lon))
    return xs.copy(data=v.astype(x.dtype))


REREGRID = True       # with --unregrid: pass the computed omega through the same output regridding


def load_unregrid(path):
    import json
    with open(path) as fh:
        d = json.load(fh)
    return {int(r["row"]): (r["scheme"], int(r["n"]), float(r["offset"])) for r in d["rows"]
            if r["scheme"] != "none" and int(r["n"]) < d["nlon"]}


_BAD_ROWS_WARNED = False


def flux_terms(s, A, B, vdim, flip, deriv, trunc=None, zf="default"):
    """Mass-flux divergence div(v dp) per layer and v.grad(ps), for ONE time step.
    Returns (u, v, ps, dp, fd, vgrad) as DataArrays, levels ordered top -> surface.
    sh: divergence and grad(ln ps) truncated at T<trunc> (default SH_TRUNC: everything the grid allows).
    zf: zonal filter factor applied to ua, va, ps first (default ZONAL_FILTER; None = off)."""
    if deriv == "sh":
        s = s.sortby("lat")
    u, v, ps = s["ua"], s["va"], s["ps"]
    zf = ZONAL_FILTER if zf == "default" else zf
    if UNREGRID:
        u0, v0, ps0 = u, v, ps
        u, v, ps = unregrid(u), unregrid(v), unregrid(ps)
        # safety: rows whose reconstruction is not finite or gives unphysical ps keep the regridded values
        bad = (~np.isfinite(ps) | (ps < 20000.) | (ps > 120000.)).any("lon")
        for x in (u, v):
            bad = bad | (~np.isfinite(x)).any([d for d in x.dims if d != "lat"])
        if bool(bad.any()):
            global _BAD_ROWS_WARNED
            if not _BAD_ROWS_WARNED:
                print(f"  note: unregrid gave unphysical values on {int(bad.sum())} latitude rows; "
                      f"those rows keep the regridded values")
                _BAD_ROWS_WARNED = True
            u, v, ps = u.where(~bad, u0), v.where(~bad, v0), ps.where(~bad, ps0)
    if zf is not None:
        u, v, ps = zonal_filter(u, zf), zonal_filter(v, zf), zonal_filter(ps, zf)
    if flip:
        u = u.isel({vdim: slice(None, None, -1)})
        v = v.isel({vdim: slice(None, None, -1)})
    mk = lambda x: xr.DataArray(x, dims=vdim)
    dA, dB = mk(np.diff(A)), mk(np.diff(B))
    dp = dA + dB * ps
    if deriv == "sh":
        sh = _sh_for(ps["lat"].values, ps["lon"].values)
        u = u.transpose(vdim, "lat", "lon")
        v = v.transpose(vdim, "lat", "lon")
        nt = _trunc_for(sh, trunc)
        D = u.copy(data=sh.divergence(u.values.astype(np.float64), v.values.astype(np.float64), nt))
        gx, gy = sh.gradient(np.log(ps.values.astype(np.float64)), nt)
        gx, gy = ps.copy(data=gx) * ps, ps.copy(data=gy) * ps        # grad(ps) = ps * grad(ln ps)
        vgrad = u * gx + v * gy                                       # grid-point product, as in the model
        fd = dp * D + dB * vgrad
    else:
        fd = _div(u * dp, v * dp)
        c = np.cos(np.deg2rad(ps["lat"]))
        vgrad = u * (_ddlon(ps) / (A_EARTH * c)) + v * (_ddlat(ps) / A_EARTH)
    return u, v, ps, dp, fd, vgrad


def omega_from_terms(terms, A, B, vdim, fullevel, deriv):
    """omega [Pa/s] and full-level pressure on model levels from flux_terms()."""
    u, v, ps, dp, fd, vgrad = terms
    mk = lambda x: xr.DataArray(x, dims=vdim)
    above = fd.cumsum(vdim) - fd                       # layers strictly above
    p_minus = mk(A[:-1]) + mk(B[:-1]) * ps             # upper half level
    p_plus = mk(A[1:]) + mk(B[1:]) * ps                # lower half level
    p = 0.5 * (p_minus + p_plus)
    if fullevel == "sb":
        # Simmons & Burridge (1981), as in the ECMWF IFS documentation (Part III):
        # (w/p)_k = -[ln(p+/p-) * sum_{j<k} div(v dp)_j + alpha_k div(v dp)_k] / dp_k
        #           + [dB_k + C_k ln(p+/p-) / dp_k] / dp_k * v_k.grad(ps)
        top = p_minus <= 0
        dlnp = xr.where(top, 0.0, np.log(p_plus / p_minus.where(~top)))
        alpha = xr.where(top, np.log(2.0), 1.0 - p_minus / dp * dlnp)
        C = mk(A[1:] * B[:-1] - A[:-1] * B[1:])
        dB = mk(np.diff(B))
        om = p * (-(dlnp * above + alpha * fd) / dp + (dB + C * dlnp / dp) / dp * vgrad)
    else:
        Bf = mk(0.5 * (B[1:] + B[:-1]))
        om = -(above + 0.5 * fd) + Bf * vgrad
    if deriv != "sh":
        lat = om["lat"]
        om = om.where((lat != lat.max()) & (lat != lat.min()))   # drop polar rows
    return om, p


def to_plev(om, p, vdim, levels):
    """Linear interpolation in log(p) to `levels` [Pa]; NaN below ground / above top."""
    n = p.sizes[vdim]
    logp = np.log(p)
    out = []
    for pt in levels:
        k1 = (p < pt).sum(vdim).clip(1, n - 1)
        k0 = k1 - 1
        lp0, lp1 = logp.isel({vdim: k0}), logp.isel({vdim: k1})
        w = (np.log(pt) - lp0) / (lp1 - lp0)
        o0, o1 = om.isel({vdim: k0}), om.isel({vdim: k1})
        val = (o0 + w * (o1 - o0)).where((w >= 0) & (w <= 1))
        out.append(val.drop_vars([vdim], errors="ignore"))
    res = xr.concat(out, dim=xr.DataArray(np.asarray(levels, float), dims="plev", name="plev"))
    return res.transpose("plev", "lat", "lon")


def omega_on_plev(s, A, B, vdim, flip, levels, fullevel="mean", deriv=None):
    """omega [Pa/s] for ONE time step (Dataset with ua, va, ps) on pressure levels."""
    deriv = deriv or LON_DERIV
    terms = flux_terms(s, A, B, vdim, flip, deriv)
    om, p = omega_from_terms(terms, A, B, vdim, fullevel, deriv)
    return to_plev(om, p, vdim, levels)


# ============================================================== scoring
COLS = ["region", "grid", "plev_hPa", "agl_bin", "n", "r", "r_time", "std_ratio", "slope",
        "bias", "rmse", "sigma_ref", "rmse_pct", "asc_ratio", "desc_ratio",
        "mean_computed", "mean_reference"]
AGL_BINS = [(0, 50), (50, 100), (100, 200), (200, 300), (300, 500)]   # hPa above ground
NAN = float("nan")


def _q(v):
    """CSV field: quoted if it contains a comma (e.g. grid label 'daily, model levels')."""
    v = str(v)
    return '"' + v.replace('"', '""') + '"' if ("," in v or '"' in v) else v


def _row(**kw):
    r = {c: NAN for c in COLS}
    r.update(kw)
    return r


def _tails(ap, bp):
    """Ratio computed/native of the mean omega in the native lowest 5 % (strong ascent)
    and highest 5 % (strong descent), at the same points."""
    bv = bp.values
    vals = bv[np.isfinite(bv)]
    if vals.size < 200:
        return NAN, NAN
    q05, q95 = np.percentile(vals, [5, 95])
    av = ap.values
    asc = np.isfinite(av) & (bv < q05)
    desc = np.isfinite(av) & (bv > q95)
    ra = float(av[asc].mean() / bv[asc].mean()) if asc.any() else NAN
    rd = float(av[desc].mean() / bv[desc].mean()) if desc.any() else NAN
    return ra, rd


def _scores(a, b, label, mask, region="global", extremes=False):
    """a: computed, b: reference, dims (time, plev, lat, lon); mask (lat, lon) bool."""
    m = mask & a.notnull() & b.notnull()
    a, b = a.where(m), b.where(m)
    dims = ["time", "lat", "lon"]
    r = xr.corr(a, b, dim=["lat", "lon"]).mean("time")
    ratio = (a.std(["lat", "lon"]) / b.std(["lat", "lon"])).mean("time")
    rmse = np.sqrt(((a - b) ** 2).mean(dims))
    sigma = b.std(dims)
    bias = (a - b).mean(dims)
    slope = xr.cov(a, b, dim=dims, ddof=0) / b.var(dims)
    if a.sizes["time"] >= 10:
        wgt = (np.cos(np.deg2rad(a["lat"])) * mask.astype(float)).fillna(0.0)
        r_time = xr.corr(a, b, dim="time").weighted(wgt).mean(["lat", "lon"])
    else:
        r_time = None
    rows = []
    for p in a["plev"].values:
        ra, rd = _tails(a.sel(plev=p), b.sel(plev=p)) if extremes else (NAN, NAN)
        sg = float(sigma.sel(plev=p))
        rows.append(_row(region=region, grid=label, plev_hPa=int(round(p / 100)),
                         n=int(a.sizes["time"]), r=float(r.sel(plev=p)),
                         r_time=float(r_time.sel(plev=p)) if r_time is not None else NAN,
                         std_ratio=float(ratio.sel(plev=p)), slope=float(slope.sel(plev=p)),
                         bias=float(bias.sel(plev=p)), rmse=float(rmse.sel(plev=p)),
                         sigma_ref=sg, rmse_pct=100 * float(rmse.sel(plev=p)) / sg if sg > 0 else NAN,
                         asc_ratio=ra, desc_ratio=rd,
                         mean_computed=float(a.sel(plev=p).mean()),
                         mean_reference=float(b.sel(plev=p).mean())))
    return rows


def _scores_agl(a, b, ps, label, mask, region, ml=None):
    """Scores pooled by height above ground (ps - p, hPa) over all levels and times.
    Two passes per bin set: (1) tail thresholds from a subsample, (2) sums."""
    a = a.transpose("time", "plev", "lat", "lon")
    b = b.transpose("time", "plev", "lat", "lon")
    psv = ps.transpose("time", "lat", "lon").values.astype(np.float64)
    mk2 = mask.transpose("lat", "lon").values.astype(bool)[None]
    nb = len(AGL_BINS)
    rng = np.random.default_rng(0)
    samples = [[] for _ in range(nb)]

    def level_data(p):
        x = a.sel(plev=p).values.astype(np.float64)
        y = b.sel(plev=p).values.astype(np.float64)
        if ml is not None:                       # model levels: p = A_k + B_k * ps
            k = int(p)
            d = (psv - (ml[0][k] + ml[1][k] * psv)) / 100.0
        else:
            d = (psv - p) / 100.0
        ok = mk2 & np.isfinite(x) & np.isfinite(y)
        return x, y, d, ok

    for p in a["plev"].values:                                    # pass 1
        x, y, d, ok = level_data(p)
        for i, (lo, hi) in enumerate(AGL_BINS):
            yy = y[ok & (d >= lo) & (d < hi)]
            if yy.size > 200000:
                yy = rng.choice(yy, 200000, replace=False)
            samples[i].append(yy)
    thr = []
    for i in range(nb):
        v = np.concatenate(samples[i]) if samples[i] else np.array([])
        thr.append(np.percentile(v, [5, 95]) if v.size >= 200 else (NAN, NAN))

    S = np.zeros((nb, 11))   # n sx sy sxx syy sxy sd2 asc_x asc_y desc_x desc_y
    for p in a["plev"].values:                                    # pass 2
        x, y, d, ok = level_data(p)
        for i, (lo, hi) in enumerate(AGL_BINS):
            sel = ok & (d >= lo) & (d < hi)
            xs, ys = x[sel], y[sel]
            if xs.size == 0:
                continue
            q05, q95 = thr[i]
            S[i, :7] += [xs.size, xs.sum(), ys.sum(), (xs * xs).sum(), (ys * ys).sum(),
                         (xs * ys).sum(), ((xs - ys) ** 2).sum()]
            if np.isfinite(q05):
                S[i, 7] += xs[ys < q05].sum(); S[i, 8] += ys[ys < q05].sum()
                S[i, 9] += xs[ys > q95].sum(); S[i, 10] += ys[ys > q95].sum()
    rows = []
    for i, (lo, hi) in enumerate(AGL_BINS):
        n = S[i, 0]
        if n < 200:
            rows.append(_row(region=region, grid=label, agl_bin=f"{lo}-{hi}", n=int(n)))
            continue
        mx, my = S[i, 1] / n, S[i, 2] / n
        vx, vy = S[i, 3] / n - mx ** 2, S[i, 4] / n - my ** 2
        cov = S[i, 5] / n - mx * my
        rmse = np.sqrt(S[i, 6] / n)
        rows.append(_row(region=region, grid=label, agl_bin=f"{lo}-{hi}", n=int(n),
                         r=cov / np.sqrt(vx * vy), std_ratio=np.sqrt(vx / vy), slope=cov / vy,
                         bias=mx - my, rmse=rmse, sigma_ref=np.sqrt(vy),
                         rmse_pct=100 * rmse / np.sqrt(vy),
                         asc_ratio=S[i, 7] / S[i, 8] if S[i, 8] != 0 else NAN,
                         desc_ratio=S[i, 9] / S[i, 10] if S[i, 10] != 0 else NAN,
                         mean_computed=mx, mean_reference=my))
    return rows


def _verdict(slope):
    if not np.isfinite(slope):
        return "n/a"
    if slope < 0.95:
        return f"too weak by {100 * (1 - slope):.0f} %"
    if slope > 1.05:
        return f"too strong by {100 * (slope - 1):.0f} %"
    return "about right"


# ============================================================== printing
def _fmt(v, f):
    return f.format(v) if np.isfinite(v) else "--".rjust(len(f.format(0.0)))


def _print_scores(rows, title):
    print(f"\n{title}")
    print(f"  {'comparison':40s} {'hPa':>4s} {'n':>4s} {'r space':>7s} {'r time':>6s} {'slope':>6s}"
          f" {'std rat':>7s} {'bias':>8s} {'RMSE':>7s} {'sigma':>7s} {'RMSE/sig':>8s}  [Pa/s]")
    for x in rows:
        print(f"  {x['grid']:40s} {x['plev_hPa']:4d} {x['n']:4d} {_fmt(x['r'], '{:7.2f}')}"
              f" {_fmt(x['r_time'], '{:6.2f}')} {_fmt(x['slope'], '{:6.2f}')} {_fmt(x['std_ratio'], '{:7.2f}')}"
              f" {_fmt(x['bias'], '{:+8.4f}')} {_fmt(x['rmse'], '{:7.4f}')} {_fmt(x['sigma_ref'], '{:7.4f}')}"
              f" {_fmt(x['rmse_pct'], '{:7.0f}')}%")


def _print_regional(rows, comps, regions):
    print("\nREGIONS   cells: r time / RMSE in % of sigma / slope (computed vs native)")
    for comp in comps:
        rr = [x for x in rows if x["grid"] == comp]
        if not rr:
            continue
        levs = sorted({x["plev_hPa"] for x in rr}, reverse=True)
        print(f"\n  {comp}")
        print("  " + f"{'hPa':>5s}" + "".join(f"{g[:22]:>24s}" for g in regions))
        for lv in levs:
            cells = []
            for g in regions:
                x = next((y for y in rr if y["region"] == g and y["plev_hPa"] == lv), None)
                if x is None or not np.isfinite(x["rmse_pct"]):
                    cells.append(f"{'--':>24s}")
                else:
                    cells.append(f"{_fmt(x['r_time'], '{:.2f}'):>9s} /{x['rmse_pct']:4.0f}% / {x['slope']:.2f}")
            print("  " + f"{lv:5d}" + "".join(f"{c:>24s}" for c in cells))


def _print_agl(rows):
    print("\nNEAR THE SURFACE   scores pooled by height above ground (ps - p), all levels and times")
    print(f"  {'region / comparison':52s} {'hPa AGL':>8s} {'n [1e3]':>8s} {'r':>5s} {'RMSE/sig':>8s}"
          f" {'slope':>6s} {'bias':>8s} {'asc':>5s} {'desc':>5s}")
    for x in rows:
        print(f"  {(x['region'] + ': ' + x['grid'])[:52]:52s} {x['agl_bin']:>8s} {x['n'] / 1e3:8.0f}"
              f" {_fmt(x['r'], '{:5.2f}')} {_fmt(x['rmse_pct'], '{:7.0f}')}% {_fmt(x['slope'], '{:6.2f}')}"
              f" {_fmt(x['bias'], '{:+8.4f}')} {_fmt(100 * x['asc_ratio'], '{:4.0f}')}%"
              f" {_fmt(100 * x['desc_ratio'], '{:4.0f}')}%")
    print("  asc/desc: mean computed omega where the native omega is in its lowest/highest 5 %,"
          " as % of the native value (100 % = fully captured, < 100 % = underestimated)")


def _print_bias(rows, regions, comps=("6-hourly instants", "daily total")):
    print("\nUNDER- OR OVERESTIMATION   cells: slope / strongest ascent / strongest descent"
          "\n  slope < 1: computed omega too weak, > 1: too strong; ascent/descent: mean computed omega where the"
          "\n  native omega is in its lowest/highest 5 %, in % of the native value (< 100 %: underestimated)")
    for comp in comps:
        rr = [x for x in rows if x["grid"] == comp and isinstance(x["agl_bin"], float)
              and np.isfinite(x["asc_ratio"])]
        if not rr:
            continue
        levs = sorted({x["plev_hPa"] for x in rr}, reverse=True)
        print(f"\n  {comp}")
        print("  " + f"{'hPa':>5s}" + "".join(f"{g[:20]:>22s}" for g in regions))
        for lv in levs:
            cells = []
            for g in regions:
                x = next((y for y in rr if y["region"] == g and y["plev_hPa"] == lv), None)
                cells.append("--" if x is None else
                             f"{x['slope']:.2f} /{100 * x['asc_ratio']:4.0f}% /{100 * x['desc_ratio']:4.0f}%")
            print("  " + f"{lv:5d}" + "".join(f"{c:>22s}" for c in cells))


# ========================================================= daily helpers
def _dates(da):
    return np.asarray(da["time"].dt.strftime("%Y-%m-%d").values)


def _daily_from_windows(x, n_per_day):
    """Daily mean of window means labelled by time; days with fewer than
    n_per_day windows are dropped. Returns time coordinate as 'YYYY-MM-DD'."""
    d = _dates(x)
    x = x.assign_coords(day=("time", d))
    mean = x.groupby("day").mean("time")
    days, counts = np.unique(d, return_counts=True)
    full = [k for k, c in zip(days, counts) if c == n_per_day]
    return mean.sel(day=full).rename({"day": "time"})


def _daily_from_instants(x, n_per_day):
    """Trapezoidal daily mean over 00-24 UTC from equally spaced instants."""
    t = x["time"].values
    steps = {str(d) for d in np.diff(t)}
    if len(steps) != 1:
        raise ValueError(f"instants are not equally spaced: {sorted(steps)[:3]}")
    wm = (0.5 * (x + x.shift(time=-1))).isel(time=slice(0, -1))
    return _daily_from_windows(wm, n_per_day)


REF_SHIFT_DAYS = 0        # --ref-shift: add this many days to the reference's date labels


def _ref_dates(da, ds=None):
    """Date labels of a daily-mean reference. Uses the midpoint of time_bnds when the file has
    them (robust against time stamps at the end of the averaging interval), plus --ref-shift."""
    t = da["time"]
    d = None
    bn = t.attrs.get("bounds", "time_bnds")
    src = da.coords if bn in da.coords else (ds if ds is not None and bn in ds.variables else None)
    if src is not None and src[bn].shape[0] == t.size:
        try:
            b = src[bn].values
            mid = [b0 + (b1 - b0) / 2 for b0, b1 in b]
            d = np.asarray([f"{m.year:04d}-{m.month:02d}-{m.day:02d}" if hasattr(m, "year")
                            else str(np.datetime_as_string(m, unit="D")) for m in mid])
        except Exception:
            d = None
    if d is None:
        d = _dates(da)
    if REF_SHIFT_DAYS:
        d = np.asarray([str(np.datetime64(x) + np.timedelta64(REF_SHIFT_DAYS, "D")) for x in d])
    return d


def _time_info(ds, label):
    """Print the first time stamps (and bounds) of a dataset: shows time-label conventions."""
    t = ds["time"].values
    msg = f"  {label} times: " + ", ".join(str(x) for x in t[:3])
    bn = ds["time"].attrs.get("bounds", "time_bnds")
    if bn in ds.variables:
        b = ds[bn].values[0]
        msg += f"  (first bounds {b[0]} .. {b[1]})"
    print(msg)


def _lag_scan(comp, ref, dimlev, max_lag=3):
    """Correlation of daily computed vs reference omega when the reference is shifted by
    -max_lag..max_lag days (anomalies from the period mean, |lat| < 70, all compared levels).
    A maximum away from lag 0 means the two time axes are labelled differently."""
    days = [str(x) for x in comp["time"].values]
    try:
        dd = np.diff(np.asarray(days, dtype="datetime64[D]")).astype(int)
    except ValueError:
        return
    if len(days) < 10 or not np.all(dd == 1):
        print("  time-lag check skipped (days not consecutive)")
        return
    m = np.abs(comp["lat"].values) < 70
    c = np.asarray(comp.values, float)[..., m, :]
    r = np.asarray(ref.values, float)[..., m, :]
    c = c - np.nanmean(c, axis=0)
    r = r - np.nanmean(r, axis=0)
    out = {}
    for k in range(-max_lag, max_lag + 1):
        a = c[max(0, -k):len(days) - max(0, k)]
        b = r[max(0, k):max(0, k) + a.shape[0]]
        ok = np.isfinite(a) & np.isfinite(b)
        out[k] = float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 100 else np.nan
    cm = np.nanmean(np.asarray(comp.values, float), axis=0)[..., m, :].ravel()
    rm = np.nanmean(np.asarray(ref.values, float), axis=0)[..., m, :].ravel()
    ok = np.isfinite(cm) & np.isfinite(rm)
    rmean = float(np.corrcoef(cm[ok], rm[ok])[0, 1])
    best = max(out, key=lambda k: -np.inf if np.isnan(out[k]) else out[k])
    print("  time-lag check (computed day d vs reference day d+k; anomalies, |lat|<70): "
          + "  ".join(f"k={k:+d}: {v:.2f}" for k, v in out.items()))
    print(f"  correlation of the period-mean maps: {rmean:.2f}")
    if best != 0 and out[best] - out[0] > 0.1:
        print(f"  WARNING: the reference matches best when shifted by {best:+d} day(s) "
              f"(r {out[0]:.2f} -> {out[best]:.2f}): the time labels of the 6hrLev input and the "
              f"reference disagree. Check the time stamps printed above; to compare with the shift: "
              f"--ref-shift {-best}")


def _midpoint_rule(ds):
    """True if the 6-hourly instants lie at the centres of 6-h windows that tile the day
    (03/09/15/21 UTC, e.g. NorESM2-MM, SAM0-UNICON). Then the daily mean is the plain mean of the
    day's four instants (midpoint rule, covering 00-24 UTC). With instants at 00/06/12/18 UTC the
    trapezoidal rule over 00-24 UTC is used. Using the trapezoidal rule on 03..21 instants would
    average over 03 UTC to 03 UTC of the next day, i.e. 3 h off the reference day."""
    t = ds["time"].values
    if len(t) < 2:
        return False
    def secs(x):
        if hasattr(x, "hour"):
            return x.hour * 3600 + x.minute * 60 + x.second
        x = np.datetime64(x, "s")
        return int((x - x.astype("datetime64[D]")).astype(int))
    d = t[1] - t[0]
    step = d.total_seconds() if hasattr(d, "total_seconds") else float(d / np.timedelta64(1, "s"))
    off = secs(t[0]) % step
    if abs(off) < 60:
        return False
    if abs(off - step / 2) < 60:
        print(f"  note: 6hrLev instants at {off / 3600:g}, {(off + step) / 3600:g}, ... UTC (centres of the "
              f"6-h windows): daily means by the midpoint rule (mean of the four instants of each day)")
        return True
    print(f"  WARNING: 6hrLev instants offset by {off / 3600:g} h from 00 UTC; trapezoidal daily means "
          f"cover {off / 3600:g} UTC to {off / 3600:g} UTC of the next day")
    return False


def _on_dates(*arrays):
    common = sorted(set.intersection(*[set(a["time"].values) for a in arrays]))
    if not common:
        raise ValueError("no common days")
    return [a.sel(time=common) for a in arrays]


def _match_grid(comp, ref):
    same = comp.reindex(lat=ref["lat"].values, lon=ref["lon"].values,
                        method="nearest", tolerance=0.05)
    return same if bool(same.notnull().any()) else comp.interp(lat=ref["lat"], lon=ref["lon"])


# ============================================================== regions
def _load_land(model, exp, member):
    for e in dict.fromkeys((exp, "historical", "piControl")):
        d = version_dirs(model, "fx", "sftlf", e)
        if d:
            path = d.get(member) or d[sorted(d)[0]]
            files = sorted(glob.glob(f"{path}/*.nc"))
            if files:
                lf = xr.open_dataset(files[0])["sftlf"].load()
                return lf * 100 if float(lf.max()) <= 1.0 else lf
    return None


def _region_masks(ps_mean, land=None, box=None):
    lat, lon = ps_mean["lat"], ps_mean["lon"]
    one = xr.ones_like(ps_mean, dtype=bool)
    alat = np.abs(lat)
    m = {"global": (alat < 70) & one,
         "tropics": (alat < 23.5) & one,
         "subtropics": (alat >= 23.5) & (alat < 35) & one,
         "midlatitudes": (alat >= 35) & (alat < 65) & one}
    nh_mid = (lat >= 30) & (lat <= 60) & (ps_mean >= 90000.)
    if land is not None:
        m["NH midlat land"] = nh_mid & (land > 50)
    else:
        m["NH midlat lowland"] = nh_mid & one
    m["high terrain"] = (ps_mean < 85000.) & (alat < 70)
    def _inbox(lon0, lon1, lat0, lat1):
        L, b0, b1 = lon % 360, lon0 % 360, lon1 % 360
        inlon = ((L >= b0) & (L <= b1)) if b0 <= b1 else ((L >= b0) | (L <= b1))
        return inlon & (lat >= lat0) & (lat <= lat1) & one
    for name, bx in REGIONS.items():
        m[name] = _inbox(*bx)
    if box is not None:
        m["box"] = _inbox(*box)
    return {k: v.transpose("lat", "lon").fillna(False).astype(bool) for k, v in m.items()}


def _masks_on(masks, grid):
    return {k: (_match_grid(v.astype(float), grid).fillna(0) > 0.5) for k, v in masks.items()}


# ============================================================== summary
def _rng(vals, fmt):
    vals = [v for v in vals if np.isfinite(v)]
    if not vals:
        return "--"
    lo, hi = min(vals), max(vals)
    return fmt.format(lo) if fmt.format(lo) == fmt.format(hi) else fmt.format(lo) + "--" + fmt.format(hi)


def _summary(rows, model, member, table, kind, dtable, start, end, focus):
    G = lambda comp, reg="global": [x for x in rows if x["region"] == reg and x["grid"] == comp
                                    and isinstance(x["agl_bin"], float)]
    inst, meth, tot = G("6-hourly instants"), G("daily method"), G("daily total")
    agl = [x for x in rows if x["agl_bin"] == "0-50" and x["grid"] == "6-hourly instants"]
    six_h = 6 * 3600 / 100.0
    L = [f"omega_check {__version__}: {model} {member}, {start}..{end}, reference {table} ({kind})", ""]
    L.append("HOW BIG IS THE RMSE?  (6-hourly instants, native grid, |lat| < 70)")
    L.append("  sigma = standard deviation of the model's own omega at that level (all points and times).")
    L.append("  '6-h displacement' = vertical drift (hPa) if an omega error of that size persisted for 6 h.")
    for x in inst:
        L.append(f"  {x['plev_hPa']:4d} hPa: RMSE {x['rmse']:.3f} Pa/s = {x['rmse_pct']:.0f} % of sigma "
                 f"({x['sigma_ref']:.3f} Pa/s);  6-h displacement {x['rmse'] * six_h:4.1f} hPa "
                 f"vs typical {x['sigma_ref'] * six_h:4.1f} hPa;  r = {x['r']:.2f}")
    L.append("")
    L.append("UNDER- OR OVERESTIMATION  (6-hourly instants, global)")
    for x in inst:
        L.append(f"  {x['plev_hPa']:4d} hPa: computed omega {_verdict(x['slope'])} (slope {x['slope']:.2f});"
                 f" strongest ascent {100 * x['asc_ratio']:.0f} %, strongest descent"
                 f" {100 * x['desc_ratio']:.0f} % of native; mean bias {x['bias']:+.4f} Pa/s")
    zones = ["tropics", "subtropics", "midlatitudes", focus]
    for comp, title in (("6-hourly instants", "6-hourly instants"), ("daily total", "daily means incl. 6-hourly sampling")):
        if not G(comp):
            continue
        L.append("")
        L.append(f"BY ZONE: {title}   r time / RMSE in % of sigma / {'computed too weak (-) or too strong (+)'}")
        L.append("  " + f"{'hPa':>5s}" + "".join(f"{z[:20]:>26s}" for z in zones))
        for lv in [x["plev_hPa"] for x in G(comp)]:
            cells = []
            for z in zones:
                x = next((y for y in G(comp, z) if y["plev_hPa"] == lv), None)
                if x is None or not np.isfinite(x["rmse_pct"]):
                    cells.append("--")
                else:
                    cells.append(f"{_fmt(x['r_time'], '{:.2f}')} / {x['rmse_pct']:.0f}% / {100 * (x['slope'] - 1):+.0f}%")
            L.append("  " + f"{lv:5d}" + "".join(f"{c:>26s}" for c in cells))
    if agl:
        L.append("")
        L.append("WITHIN 50 hPa OF THE SURFACE  (6-hourly instants)")
        for x in agl:
            if np.isfinite(x["r"]):
                L.append(f"  {x['region']}: r = {x['r']:.2f}, RMSE = {x['rmse_pct']:.0f} % of sigma,"
                         f" {_verdict(x['slope'])} (slope {x['slope']:.2f}),"
                         f" strongest ascent {100 * x['asc_ratio']:.0f} %")
    if meth:
        L.append("")
        L.append("DAILY MEANS (method: computed vs native, both from 6-hourly instants)")
        for x in meth:
            L.append(f"  {x['plev_hPa']:4d} hPa: RMSE {x['rmse']:.3f} Pa/s = {x['rmse_pct']:.0f} % of sigma;"
                     f"  r time = {x['r_time']:.2f};  {_verdict(x['slope'])}")
    if tot:
        L.append("")
        L.append(f"DAILY MEANS incl. 6-hourly sampling (computed vs true daily mean, {dtable})")
        for x in tot:
            L.append(f"  {x['plev_hPa']:4d} hPa: RMSE {x['rmse']:.3f} Pa/s = {x['rmse_pct']:.0f} % of sigma;"
                     f"  r time = {x['r_time']:.2f};  {_verdict(x['slope'])}")

    hp = [x["plev_hPa"] for x in inst]
    what = "instantaneous" if kind == "instantaneous" else "6-hourly mean"
    sl = [x["slope"] for x in inst]
    tex = (f"This diagnostic approach was validated against the native {what} $\\omega$ of {model} "
           f"(table \\texttt{{{table}}}) for {start}--{end}, yielding pattern correlations of "
           f"{_rng([x['r'] for x in inst], '{:.2f}')} between {max(hp)} and {min(hp)}\\,hPa. "
           f"The corresponding root-mean-square errors are {_rng([x['rmse'] for x in inst], '{:.3f}')}\\,Pa\\,s$^{{-1}}$, "
           f"i.e.\\ {_rng([x['rmse_pct'] for x in inst], '{:.0f}')}\\,\\% of the standard deviation of $\\omega$; "
           f"regression slopes of {_rng(sl, '{:.2f}')} indicate that the diagnosed $\\omega$ "
           + ("slightly underestimates" if max(sl) < 0.98 else "slightly overestimates" if min(sl) > 1.02
              else "neither systematically under- nor overestimates")
           + " the magnitude of the native vertical velocity.")
    p_ref = min(hp, key=lambda h: abs(h - 500))
    zr = {z: next((x for x in G("6-hourly instants", z) if x["plev_hPa"] == p_ref), None)
          for z in ("tropics", "subtropics", "midlatitudes")}
    if all(v is not None and np.isfinite(v["rmse_pct"]) for v in zr.values()):
        tex += (f" At {p_ref}\\,hPa, the error amounts to {zr['tropics']['rmse_pct']:.0f}\\,\\% of the standard "
                f"deviation in the tropics, {zr['subtropics']['rmse_pct']:.0f}\\,\\% in the subtropics and "
                f"{zr['midlatitudes']['rmse_pct']:.0f}\\,\\% in the midlatitudes.")
    g_agl = next((x for x in agl if x["region"] == "global" and np.isfinite(x["r"])), None)
    if g_agl:
        tex += (f" Within 50\\,hPa of the surface, the correlation is {g_agl['r']:.2f} and the "
                f"root-mean-square error {g_agl['rmse_pct']:.0f}\\,\\% of the standard deviation.")
    if meth and all(np.isfinite(x["r_time"]) for x in meth):
        tex += (f" For daily means, the temporal correlation between computed and native $\\omega$ is "
                f"{_rng([x['r_time'] for x in meth], '{:.2f}')}, with root-mean-square errors of "
                f"{_rng([x['rmse'] for x in meth], '{:.3f}')}\\,Pa\\,s$^{{-1}}$ "
                f"({_rng([x['rmse_pct'] for x in meth], '{:.0f}')}\\,\\% of the standard deviation).")
    L += ["", "LATEX (numbers from this run):", tex]
    return "\n".join(L)


# ============================================================== figures
def _plt():
    import matplotlib
    if "ipykernel" not in sys.modules:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _figure(comp, ref, model, member, table, fname):
    plt = _plt()
    p_show = float(min(comp["plev"].values, key=lambda p: abs(p - 50000.)))
    fig, ax = plt.subplots(2, 2, figsize=(13, 6.5), constrained_layout=True)
    panels = [(comp.sel(plev=p_show).isel(time=0), "computed, first time", 1.0),
              (ref.sel(plev=p_show).isel(time=0), f"native {table}, first time", 1.0),
              (comp.sel(plev=p_show).mean("time"), "computed, period mean", 0.2),
              (ref.sel(plev=p_show).mean("time"), f"native {table}, period mean", 0.2)]
    for k, (fld, title, vmax) in enumerate(panels):
        a = ax.flat[k]
        im = fld.plot(ax=a, cmap="RdBu_r", vmin=-vmax, vmax=vmax, add_colorbar=False)
        a.set_title(title)
        if k % 2 == 1:
            fig.colorbar(im, ax=ax[k // 2, :], label="omega [Pa/s]")
    fig.suptitle(f"{model} {member}: omega at {p_show / 100:.0f} hPa")
    fig.savefig(fname, dpi=120)
    plt.close(fig)


def _figure_daily(comp_d, ref6_d, nat_d, model, member, dtable, fname):
    plt = _plt()
    p_show = float(min(comp_d["plev"].values, key=lambda p: abs(p - 50000.)))
    a, b, c = (x.sel(plev=p_show) for x in (comp_d, nat_d, ref6_d))
    fig, ax = plt.subplots(2, 2, figsize=(13, 6.5), constrained_layout=True)
    day0 = str(a["time"].values[0])
    im1 = a.isel(time=0).plot(ax=ax[0, 0], cmap="RdBu_r", vmin=-0.5, vmax=0.5, add_colorbar=False)
    ax[0, 0].set_title(f"computed daily mean, {day0}")
    b.isel(time=0).plot(ax=ax[0, 1], cmap="RdBu_r", vmin=-0.5, vmax=0.5, add_colorbar=False)
    ax[0, 1].set_title(f"native daily mean ({dtable}), {day0}")
    fig.colorbar(im1, ax=ax[0, :], label="omega [Pa/s]")
    if a.sizes["time"] >= 10:
        im2 = xr.corr(a, b, dim="time").plot(ax=ax[1, 0], cmap="viridis", vmin=0, vmax=1, add_colorbar=False)
        ax[1, 0].set_title("temporal r: computed vs native daily (total)")
        xr.corr(c, b, dim="time").plot(ax=ax[1, 1], cmap="viridis", vmin=0, vmax=1, add_colorbar=False)
        ax[1, 1].set_title("temporal r: native 6-hourly vs native daily (sampling)")
        fig.colorbar(im2, ax=ax[1, :], label="correlation of daily means")
    else:
        for k in (0, 1):
            ax[1, k].set_axis_off()
            ax[1, k].text(0.5, 0.5, "fewer than 10 days: no temporal correlation",
                          ha="center", va="center", transform=ax[1, k].transAxes)
    fig.suptitle(f"{model} {member}: daily-mean omega at {p_show / 100:.0f} hPa")
    fig.savefig(fname, dpi=120)
    plt.close(fig)


def _figure_agl(rows, model, member, regions, comp, fname):
    """Small multiples: RMSE (% of sigma) and regression slope vs height above ground,
    one line per region (surface at the bottom); global in dark ink, zones in categorical order."""
    plt = _plt()
    palette = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
    markers = ["o", "s", "^", "D", "v"]
    ink, muted, grid = "#2b2b29", "#6b6a64", "#e4e3df"
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True, sharey=True)
    k = 0
    for reg in regions:
        rr = [x for x in rows if x["region"] == reg and x["grid"] == comp
              and isinstance(x["agl_bin"], str) and np.isfinite(x["r"])]
        if not rr:
            continue
        if reg == "global":
            col, mk, lw = ink, "o", 2.5
        else:
            col, mk, lw = palette[k % len(palette)], markers[k % len(markers)], 2
            k += 1
        y = [np.mean([float(v) for v in x["agl_bin"].split("-")]) for x in rr]
        for ax, key in zip(axs, ("rmse_pct", "slope")):
            xv = [x[key] for x in rr]
            ax.plot(xv, y, color=col, lw=lw, marker=mk, ms=8, mec="white", mew=1.5, label=reg)
    axs[1].axvline(1.0, color=muted, lw=1, ls="--")
    axs[1].text(1.0, 1.0, " no bias", color=muted, fontsize=8, va="top",
                transform=axs[1].get_xaxis_transform())
    axs[0].set_xlabel("RMSE in % of the standard deviation of native omega", color=ink)
    axs[1].set_xlabel("regression slope computed vs native  (< 1: too weak)", color=ink)
    axs[0].set_ylabel("height above ground  (ps - p)  [hPa]", color=ink)
    axs[0].set_ylim(bottom=0)
    for ax in axs:
        ax.grid(color=grid, lw=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(muted)
        ax.tick_params(colors=muted)
    axs[1].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.suptitle(f"{model} {member}: computed vs native omega by height above ground - {comp}", color=ink)
    fig.savefig(fname, dpi=120)
    plt.close(fig)



def _figure_maps(comp, ref, land, model, member, table, kind, vlab, stem,
                 show=(92500., 85000., 50000.), snapshot=None, box=None):
    """Maps at 925/850/500 hPa (or nearest available):
    (1) snapshot: computed, native, computed - native;
    (2) skill over the period: temporal r, RMSE in % of local sigma, local regression slope."""
    plt = _plt()
    import matplotlib
    ink, muted, nan_col = "#2b2b29", "#6b6a64", "#d9d9d6"
    plevs = list(comp["plev"].values)
    levs = []
    for p in show:
        q = min(plevs, key=lambda x: abs(x - p))
        if abs(q - p) <= 2500 and q not in levs:
            levs.append(q)
    if not levs:
        return []
    m = np.abs(comp["lat"]) <= 85

    def cmap(name):
        c = matplotlib.colormaps[name].copy()
        c.set_bad(nan_col)
        return c

    def coast(ax):
        if land is not None:
            land.plot.contour(ax=ax, levels=[50], colors=[muted], linewidths=0.5, add_labels=False)

    def tidy(ax, title):
        ax.set_title(title, color=ink, fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(colors=muted, labelsize=8)

    # ---- snapshot time and optional zoom
    t0 = 0
    if snapshot:
        tt = [str(t)[:16].replace(" ", "T") for t in comp["time"].values]
        want = str(snapshot).replace(" ", "T")
        hits = [i for i, t in enumerate(tt) if t.startswith(want)]
        t0 = hits[0] if hits else 0
        if not hits:
            print(f"  snapshot {snapshot} not found; using first time {tt[0]}")
    tlabel = str(comp["time"].values[t0])[:16]
    if box is not None:
        lon0, lon1, lat0, lat1 = box
        if lon0 < 0 or lon1 < 0:                      # plot on -180..180
            relon = lambda x: x.assign_coords(lon=((x["lon"] + 180) % 360) - 180).sortby("lon")
            comp, ref = relon(comp), relon(ref)
            land = relon(land) if land is not None else None
        sel = dict(lon=slice(min(lon0, lon1), max(lon0, lon1)), lat=slice(min(lat0, lat1), max(lat0, lat1)))
        comp, ref = comp.sortby("lat").sel(**sel), ref.sortby("lat").sel(**sel)
        land = land.sortby("lat").sel(**sel) if land is not None else None
        m = np.abs(comp["lat"]) <= 85
    letters = "abcdefghijklmnopqrstuvwxyz"
    what = {"instantaneous": "instantaneous", "6-h mean": "6-h mean"}.get(kind, kind)
    out = []

    # ---- (1) snapshot
    fig, axs = plt.subplots(len(levs), 3, figsize=(16, 3.3 * len(levs) + 0.8), constrained_layout=True,
                            squeeze=False)
    for i, p in enumerate(levs):
        a = comp.sel(plev=p).isel(time=t0).where(m)
        b = ref.sel(plev=p).isel(time=t0).where(m & comp.sel(plev=p).isel(time=t0).notnull())
        vmax = float(np.nanpercentile(np.abs(b.values), 99)) or 1.0
        vmax = float(f"{vmax:.1g}")
        im = a.plot(ax=axs[i, 0], cmap=cmap("RdBu_r"), vmin=-vmax, vmax=vmax, add_colorbar=False)
        b.plot(ax=axs[i, 1], cmap=cmap("RdBu_r"), vmin=-vmax, vmax=vmax, add_colorbar=False)
        dv = vmax / 4
        imd = (a - b).plot(ax=axs[i, 2], cmap=cmap("RdBu_r"), vmin=-dv, vmax=dv, add_colorbar=False)
        both = a.notnull() & b.notnull()
        av, bv = a.where(both).values.ravel(), b.where(both).values.ravel()
        okv = np.isfinite(av) & np.isfinite(bv)
        r_s = np.corrcoef(av[okv], bv[okv])[0, 1] if okv.sum() > 10 else np.nan
        e_s = 100 * np.sqrt(np.mean((av[okv] - bv[okv]) ** 2)) / np.std(bv[okv]) if okv.sum() > 10 else np.nan
        for j, ttl in enumerate((f"computed ({vlab})", f"native {table}", "computed - native")):
            tidy(axs[i, j], f"({letters[3 * i + j]}) {p / 100:.0f} hPa: {ttl}")
            coast(axs[i, j])
        axs[i, 2].text(0.01, 0.02, f"r = {r_s:.3f}, RMSE = {e_s:.0f} % of sigma", transform=axs[i, 2].transAxes,
                       fontsize=8, color=ink, bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=2))
        fig.colorbar(im, ax=axs[i, :2], label="omega [Pa/s]", shrink=0.9, pad=0.01)
        fig.colorbar(imd, ax=axs[i, 2], label="difference [Pa/s]", shrink=0.9, pad=0.01)
    fig.suptitle(f"{model} {member}: omega ({what}) at {tlabel}   "
                 f"(negative = ascent; grey = below ground)", color=ink)
    fname = f"{stem}_maps.png"
    fig.savefig(fname, dpi=200)
    fig.savefig(fname.replace(".png", ".pdf"))
    plt.close(fig)
    out += [fname, fname.replace(".png", ".pdf")]

    # ---- (2) skill over the period
    fig, axs = plt.subplots(len(levs), 3, figsize=(16, 3.3 * len(levs) + 0.8), constrained_layout=True,
                            squeeze=False)
    for i, p in enumerate(levs):
        a, b = comp.sel(plev=p), ref.sel(plev=p)
        ok = a.notnull() & b.notnull()
        a, b = a.where(ok), b.where(ok)
        enough = ok.sum("time") >= 10
        r = xr.corr(a, b, dim="time").where(enough & m)
        pct = (100 * np.sqrt(((a - b) ** 2).mean("time")) / b.std("time")).where(enough & m)
        slope = (xr.cov(a, b, dim="time", ddof=0) / b.var("time")).where(enough & m)
        im1 = r.plot(ax=axs[i, 0], cmap=cmap("Blues"), vmin=0.8, vmax=1.0, add_colorbar=False)
        im2 = pct.plot(ax=axs[i, 1], cmap=cmap("Oranges"), vmin=0, vmax=50, add_colorbar=False)
        im3 = slope.plot(ax=axs[i, 2], cmap=cmap("PuOr"), vmin=0.8, vmax=1.2, add_colorbar=False)
        gm = lambda x: float(x.weighted(np.cos(np.deg2rad(x["lat"]))).mean())
        tidy(axs[i, 0], f"({letters[3 * i]}) {p / 100:.0f} hPa: temporal r (mean {gm(r):.3f})")
        tidy(axs[i, 1], f"({letters[3 * i + 1]}) {p / 100:.0f} hPa: RMSE in % of local sigma (mean {gm(pct):.0f} %)")
        tidy(axs[i, 2], f"({letters[3 * i + 2]}) {p / 100:.0f} hPa: slope, < 1 too weak (mean {gm(slope):.3f})")
        for j in range(3):
            coast(axs[i, j])
        for im, ax, lab in ((im1, axs[i, 0], "r"), (im2, axs[i, 1], "%"), (im3, axs[i, 2], "slope")):
            fig.colorbar(im, ax=ax, label=lab, shrink=0.9, pad=0.01)
    fig.suptitle(f"{model} {member}: computed ({vlab}) vs native {table} omega, "
                 f"{str(comp['time'].values[0])[:10]} to {str(comp['time'].values[-1])[:10]}, "
                 f"{comp.sizes['time']} times", color=ink)
    fname = f"{stem}_skillmaps.png"
    fig.savefig(fname, dpi=200)
    fig.savefig(fname.replace(".png", ".pdf"))
    plt.close(fig)
    out += [fname, fname.replace(".png", ".pdf")]
    return out


# ======================================================= CFday model levels
def _figure_zonal(zspec, zlats, klev, pnom, nat_spec, espec, variants, vl, kerr, model, member, fname):
    """(a, b) zonal power spectra of the INPUT ua by latitude; a cliff or a mirrored tail below
    m = nlon/2 at higher latitudes reveals regridding from a reduced grid. (c, d) zonal power
    spectra of the daily-mean omega error (computed - native) per variant, 40-60 and 0-20 deg."""
    plt = _plt()
    import matplotlib as mpl
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.2))
    cmap = mpl.colormaps["viridis"]
    for c, k in enumerate(klev):
        ax = axs[0, c]
        for j, L_ in enumerate(zlats):
            col = cmap(j / max(len(zlats) - 1, 1) * 0.9)
            P = zspec[k][j]
            m = np.arange(P.size)
            ax.semilogy(m[1:], P[1:], color=col, lw=1.3, label=f"{L_}\u00b0")
            mc = (P.size - 1) * np.cos(np.deg2rad(L_))
            ax.axvline(mc, color=col, lw=0.6, ls=":")
        ax.set_title(f"({'ab'[c]}) input ua, zonal power by latitude, level ~{pnom[k] / 100:.0f} hPa\n"
                     "(dotted: m = 128 cos(lat))", fontsize=10)
        ax.set_xlabel("zonal wavenumber m"); ax.set_ylabel("power [(m/s)$^2$]")
    axs[0, 1].legend(title="|lat|", fontsize=8, frameon=False, ncol=2)
    pal = ["#2b2b2b", "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#8c6bb1"]
    for c, nm in enumerate(("40-60", "0-20")):
        ax = axs[1, c]
        m = np.arange(nat_spec[nm].size)
        ax.semilogy(m[1:], nat_spec[nm][1:], color="#999999", lw=2.5, label="native omega (signal)")
        for j, v in enumerate(variants):
            ax.semilogy(m[1:], espec[(nm, v)][1:], color=pal[(j + 1) % len(pal)], lw=1.3, label=f"error, {vl(v)}")
        ax.set_title(f"({'cd'[c]}) daily-mean omega, level ~{pnom[kerr] / 100:.0f} hPa, |lat| {nm}\u00b0", fontsize=10)
        ax.set_xlabel("zonal wavenumber m"); ax.set_ylabel("power [(Pa/s)$^2$]")
    axs[1, 0].legend(fontsize=7.5, frameon=False)
    for ax in axs.flat:
        ax.grid(True, color="#e6e6e6", lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"{model} {member}: zonal spectra (input wind; omega error)", fontsize=11)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)


def run_cfday(model, start="2000-06-01", end="2000-08-31", exp="historical", member=None,
              outdir=".", make_figure=True, deriv="auto", box=None, fullevel="sb", snapshot=None,
              scan=False):
    """Daily means of computed omega on model levels vs CFday wap (daily mean, model levels).
    scan=True (spherical harmonics only): repeat for several truncations and compare."""
    t0 = time.time()
    lev = {v: version_dirs(model, "6hrLev", v, exp) for v in ("ua", "va", "ps")}
    cfd = version_dirs(model, "CFday", "wap", exp)
    common = sorted(set(lev["ua"]) & set(lev["va"]) & set(lev["ps"]) & set(cfd))
    if member is not None:
        common = [m for m in common if m == member]
    if not common:
        raise FileNotFoundError(f"{model}: no member with both 6hrLev (ua, va, ps) and CFday wap")
    member = common[0]
    print(f"omega_check {__version__} --cfday: {model}  member {member}  {start}..{end}")

    parts = _common_time([open_files(lev[v][member], start, end) for v in ("ua", "va", "ps")])
    ds = xr.merge([parts[0], parts[1][["va"]], parts[2][["ps"]]], compat="override", join="override")
    A, B, vdim, flip = half_level_coeffs(ds)
    n = len(A) - 1
    cf = open_files(cfd[member], start, end)
    print("  CFday cell_methods:", cf["wap"].attrs.get("cell_methods"))
    _time_info(ds, "6hrLev")
    _time_info(cf, "CFday")
    Ac, Bc, vdim_c, flip_c = half_level_coeffs(cf, var="wap")
    if len(Ac) != len(A):
        raise ValueError(f"CFday has {len(Ac) - 1} levels, 6hrLev has {n}")
    if not (np.allclose(Ac, A, rtol=1e-4, atol=1.0) and np.allclose(Bc, B, rtol=1e-4, atol=1e-6)):
        print("  WARNING: CFday and 6hrLev vertical coefficients differ; comparing level by level anyway")
    Af, Bf = 0.5 * (A[1:] + A[:-1]), 0.5 * (B[1:] + B[:-1])
    pnom = Af + Bf * 1e5
    print(f"  {n} model levels; {ds.sizes['lat']}x{ds.sizes['lon']} grid; {ds.sizes['time']} time steps")

    sh = None
    if deriv in ("auto", "sh"):
        try:
            sh = _sh_for(np.sort(ds["lat"].values), ds["lon"].values)
            deriv = "sh"
        except ValueError as e:
            if deriv == "sh":
                print(f"  spherical harmonics not possible ({e}); using centred differences")
            deriv = "cd"
    global LON_DERIV
    LON_DERIV = "fft" if deriv == "fft" else "cd"
    nt0 = _trunc_for(sh) if sh is not None else None
    # variants: (deriv, truncation, zonal-filter factor); the first one is the default
    variants = [(deriv, nt0, ZONAL_FILTER)]
    if scan:
        for zf in (None, 1.3, 1.15, 1.0, 0.9):
            variants.append((deriv, nt0, zf))
        variants.append(("cd", None, None))
        variants = list(dict.fromkeys(variants))
    vl = lambda v: (f"{fullevel}/{v[0]}" + (f"/T{v[1]}" if v[1] is not None and v[0] == "sh" else "")
                    + (f"/zf{v[2]:g}" if v[2] is not None else ""))
    v0 = variants[0]
    vlab = vl(v0)
    if sh is not None:
        print(f"  spherical harmonics: truncation T{nt0} (grid allows T{sh.N})")
    print(f"  grid_label {ds.attrs.get('grid_label', '?')}: {ds.attrs.get('grid', '')}")
    print("  variants:", ", ".join(vl(v) for v in variants), "  (default first)")

    # zonal spectra of the input wind (diagnostic for regridding from a reduced grid)
    k500 = int(np.argmin(np.abs(pnom - 50000.)))
    k850 = int(np.argmin(np.abs(pnom - 85000.)))
    latv = np.sort(ds["lat"].values)
    zlats = [0, 20, 30, 40, 45, 50, 55, 60, 70]
    zrows = {L_: [int(np.argmin(np.abs(latv - L_))), int(np.argmin(np.abs(latv + L_)))] for L_ in zlats}
    zspec = {k: np.zeros((len(zlats), ds.sizes["lon"] // 2 + 1)) for k in (k500, k850)}
    nz = 0

    # ---- streaming trapezoidal daily means on model levels (window mean -> day of window start)
    acc = {v: {} for v in variants}
    accp, cnt = {}, {}
    prev = {}
    mid = _midpoint_rule(ds)
    prevp = prevday = None
    lat_ref = lon_ref = None
    for i in range(ds.sizes["time"]):
        s = ds[["ua", "va", "ps"]].isel(time=i).load()
        day = str(s["time"].dt.strftime("%Y-%m-%d").values)
        cur = {}
        for v in variants:
            LON_DERIV = "fft" if v[0] == "fft" else "cd"
            terms = flux_terms(s, A, B, vdim, flip, v[0], trunc=v[1], zf=v[2])
            om, _ = omega_from_terms(terms, A, B, vdim, fullevel, v[0])
            om = om.sortby("lat").transpose(vdim, "lat", "lon")
            cur[v] = om.values.astype(np.float32)
        ps = s["ps"].sortby("lat")
        if lat_ref is None:
            lat_ref, lon_ref = om["lat"].values, om["lon"].values
        ua_ = s["ua"].sortby("lat").transpose(vdim, "lat", "lon")
        if flip:
            ua_ = ua_.isel({vdim: slice(None, None, -1)})
        for k in zspec:
            P = np.abs(np.fft.rfft(ua_.isel({vdim: k}).values.astype(np.float64), axis=-1) / ds.sizes["lon"]) ** 2
            for j, L_ in enumerate(zlats):
                zspec[k][j] += P[zrows[L_]].mean(0)
        nz += 1
        curp = ps.values.astype(np.float32)
        if mid:                       # instant = mean of its own 6-h window, which lies in its day
            prev, prevp, prevday = cur, curp, day
        if prevday is not None:
            if prevday not in cnt:
                accp[prevday], cnt[prevday] = np.zeros_like(curp), 0
                for v in variants:
                    acc[v][prevday] = np.zeros_like(cur[v])
            for v in variants:
                acc[v][prevday] += 0.5 * (prev[v] + cur[v])
            accp[prevday] += 0.5 * (prevp + curp)
            cnt[prevday] += 1
        prev, prevp, prevday = cur, curp, day
        if i % 40 == 0:
            print(f"    step {i + 1}/{ds.sizes['time']}  ({time.time() - t0:.0f} s)", flush=True)
    LON_DERIV = "fft" if deriv == "fft" else "cd"
    zspec = {k: z / nz for k, z in zspec.items()}
    days = sorted(d for d in cnt if cnt[d] == 4)
    if not days:
        raise ValueError("no complete days")
    crd = {"time": days, "lev": np.arange(n), "lat": lat_ref, "lon": lon_ref}
    comps = {v: xr.DataArray(np.stack([acc[v][d] / 4 for d in days]), dims=("time", "lev", "lat", "lon"),
                             coords=crd) for v in variants}
    del acc
    psd = xr.DataArray(np.stack([accp[d] / 4 for d in days]), dims=("time", "lat", "lon"),
                       coords={"time": days, "lat": lat_ref, "lon": lon_ref})
    print(f"  daily means of computed omega for {len(days)} days ({time.time() - t0:.0f} s)")

    # ---- native CFday on the same days, levels top -> surface
    ref = cf["wap"].transpose("time", vdim_c, "lat", "lon")
    if flip_c:
        ref = ref.isel({vdim_c: slice(None, None, -1)})
    ref = ref.astype(np.float32).load()
    ref = ref.rename({vdim_c: "lev"}).assign_coords(lev=np.arange(n), time=_ref_dates(ref, cf)).sortby("lat")
    ref = ref.drop_vars([c for c in ref.coords if c not in ("time", "lev", "lat", "lon")])
    comps = {v: _match_grid(c, ref) for v, c in comps.items()}
    psd = _match_grid(psd, ref)
    out = _on_dates(*comps.values(), psd, ref)
    comps = dict(zip(variants, out[:len(variants)]))
    psd, ref = out[-2], out[-1]
    if UNREGRID and REREGRID:
        # native CFday omega went through CNRM's conservative output regridding (smoothing on the
        # reduced rows); give the computed omega the same treatment so that like is compared with like
        comps = {v: reregrid(c) for v, c in comps.items()}
        print(f"  computed omega passed through the same output regridding as the native field "
              f"({len(UNREGRID)} rows; --no-reregrid to skip)")
    comp = comps[v0]
    print(f"  compared days: {comp.sizes['time']}")
    _lag_scan(comp, ref, None)

    # ---- regions
    land = _load_land(model, exp, member)
    if land is not None:
        land = _match_grid(land, ref)
    masks = _region_masks(psd.mean("time"), land, box)
    focus = next(k for k in masks if k.startswith("NH midlat"))
    regions = list(masks)
    key_regions = ["global", "tropics", "subtropics", "midlatitudes", focus] + (["box"] if box is not None else [])

    # ---- model levels nearest to nominal pressures
    targets = [92500., 85000., 70000., 60000., 50000., 25000.]
    ks = []
    for t in targets:
        k = int(np.argmin(np.abs(pnom - t)))
        if k not in ks:
            ks.append(k)
    sel = lambda x: x.isel(lev=ks).rename({"lev": "plev"}).assign_coords(plev=pnom[ks])
    allm = lambda x: x.rename({"lev": "plev"}).assign_coords(plev=np.arange(n).astype(float))
    rs, ra = sel(ref), allm(ref)
    lab = "daily, model levels"

    # ---- scan: compact scores for every variant
    scan_rows = {}
    if len(variants) > 1:
        for v in variants:
            r_ = _scores(sel(comps[v]), rs, lab, masks["global"], "global")
            for reg in ("global", "midlatitudes", focus, "tropics"):
                r_ += _scores_agl(allm(comps[v]), ra, psd, lab, masks[reg], reg, ml=(Af, Bf))[:1]
            scan_rows[v] = r_
            print(f"  scored {vl(v)} ({time.time() - t0:.0f} s)", flush=True)

    rows = []
    cs = sel(comp)
    for reg in regions:
        rows += _scores(cs, rs, lab, masks[reg], reg, extremes=reg in key_regions)
    ca = allm(comp)
    for reg in key_regions:
        rows += _scores_agl(ca, ra, psd, lab, masks[reg], reg, ml=(Af, Bf))

    print(f"\nCFDAY  {model} {member}: daily-mean computed omega ({vlab}) vs native CFday wap on model levels")
    print("  levels: model levels nearest to the nominal pressure (for ps = 1000 hPa):",
          ", ".join(f"{pnom[k] / 100:.0f}" for k in ks), "hPa")
    _print_scores([x for x in rows if x["region"] == "global" and isinstance(x["agl_bin"], float)],
                  "GLOBAL (|lat|<70)")
    _print_regional([x for x in rows if isinstance(x["agl_bin"], float)], [lab], regions)
    _print_agl([x for x in rows if isinstance(x["agl_bin"], str)])
    _print_bias(rows, key_regions, comps=(lab,))

    # ---- summary
    g = [x for x in rows if x["region"] == "global" and isinstance(x["agl_bin"], float)]
    agl = {x["region"]: x for x in rows if x["agl_bin"] == "0-50" and np.isfinite(x["r"])}
    L = [f"omega_check {__version__} --cfday: {model} {member}, {start}..{end}, variant {vlab}", ""]
    L += [f"INPUT GRID: grid_label {ds.attrs.get('grid_label', '?')} ({ds.attrs.get('grid', '')})", ""]
    L += ["DAILY MEANS ON MODEL LEVELS (computed vs CFday wap), global |lat|<70",
          "  nominal hPa   r time   r space   slope   std ratio   RMSE in % of sigma"]
    for x in g:
        L.append(f"  {x['plev_hPa']:8d}    {x['r_time']:6.3f}   {x['r']:6.3f}   {x['slope']:6.3f}"
                 f"   {x['std_ratio']:6.2f}      {x['rmse_pct']:5.0f} %")
    L.append("")
    L.append("WITHIN 50 hPa OF THE SURFACE (actual level pressure)")
    for reg, x in agl.items():
        L.append(f"  {reg:16s} r = {x['r']:.3f}, slope {x['slope']:.3f}, std ratio {x['std_ratio']:.2f}, "
                 f"RMSE {x['rmse_pct']:.0f} % of sigma")
    L += ["", "  std ratio > 1 with slope ~ 1 means extra variance in the computed field that is not in the",
          "  native one (noise), rather than a too weak or too strong signal."]
    if scan_rows:
        L += ["", "VARIANT SCAN  (r / RMSE in % of sigma; std ratio in brackets)",
              "  zf c: zonal wavenumbers m > c * 128 * cos(lat) removed from ua, va, ps before the derivatives",
              f"  {'variant':22s} {'~' + str(int(pnom[ks[0]] / 100)) + ' hPa':>18s} {'~' + str(int(pnom[ks[-2]] / 100)) + ' hPa':>18s}"
              f" {'0-50 global':>18s} {'0-50 midlat':>18s} {'0-50 ' + focus[:12]:>22s} {'0-50 tropics':>18s}"]
        hp_ = [x["plev_hPa"] for x in g]
        for v, r_ in scan_rows.items():
            lv = {x["plev_hPa"]: x for x in r_ if isinstance(x["agl_bin"], float)}
            ag = {x["region"]: x for x in r_ if x["agl_bin"] == "0-50"}
            cell = lambda x: f"{x['r']:.3f}/{x['rmse_pct']:3.0f}% ({x['std_ratio']:.2f})" if x else "--"
            L.append(f"  {vl(v):22s} {cell(lv.get(hp_[0])):>18s} {cell(lv.get(hp_[-2])):>18s}"
                     f" {cell(ag.get('global')):>18s} {cell(ag.get('midlatitudes')):>18s}"
                     f" {cell(ag.get(focus)):>22s} {cell(ag.get('tropics')):>18s}"
                     + ("   <- default" if v == v0 else ""))
    hp = [x["plev_hPa"] for x in g]
    tex = (f"Compared directly on the model levels with the native daily-mean $\\omega$ of {model} "
           f"(table \\texttt{{CFday}}), i.e.\\ without interpolation to pressure levels, daily means of the "
           f"diagnosed $\\omega$ yield temporal correlations of {_rng([x['r_time'] for x in g], '{:.2f}')} and "
           f"regression slopes of {_rng([x['slope'] for x in g], '{:.2f}')} between the model levels at "
           f"$\\sim${max(hp)} and $\\sim${min(hp)}\\,hPa")
    if "global" in agl:
        tex += (f", and a correlation of {agl['global']['r']:.2f} within 50\\,hPa of the surface")
    tex += "."
    L += ["", "LATEX (numbers from this run):", tex]
    text = "\n".join(L)
    os.makedirs(outdir, exist_ok=True)
    tag = f"{model}_cfday" + ("" if (SH_TRUNC == "none" and ZONAL_FILTER is None) else
                              "_" + vlab.replace("/", "_")) + ("_unregrid" + ("" if REREGRID else "_noreregrid") if UNREGRID else "")
    with open(os.path.join(outdir, f"omega_check_{tag}_summary.txt"), "w") as f:
        f.write(text + "\n")
    fmt = lambda x, c: ("" if isinstance(x[c], float) and not np.isfinite(x[c]) else
                        str(round(x[c], 5)) if isinstance(x[c], float) else _q(x[c]))
    with open(os.path.join(outdir, f"omega_check_{tag}.csv"), "w") as f:
        f.write("model,member,reference,start,end,variant," + ",".join(COLS) + "\n")
        for x in rows:
            f.write(f"{model},{member},CFday,{start},{end},{vlab}," + ",".join(fmt(x, c) for c in COLS) + "\n")
    if scan_rows:
        with open(os.path.join(outdir, f"omega_check_{tag}_scan.csv"), "w") as f:
            f.write("model,member,reference,start,end,variant," + ",".join(COLS) + "\n")
            for v, r_ in scan_rows.items():
                for x in r_:
                    f.write(f"{model},{member},CFday,{start},{end},{vl(v)}," + ",".join(fmt(x, c) for c in COLS) + "\n")
    print("\n" + text + "\n")
    if make_figure:
        for fname in _figure_maps(cs, rs, land, model, member, "CFday", "daily mean", vlab,
                                  os.path.join(outdir, f"omega_check_{tag}"), snapshot=snapshot):
            print(f"figure: {fname}")
        png = os.path.join(outdir, f"omega_check_{tag}_nearsurface.png")
        _figure_agl(rows, model, member, key_regions, lab, png)
        print(f"figure: {png}")
        # zonal spectra: input wind by latitude, and the error by latitude band
        lat_ = ref["lat"].values
        band = lambda lo, hi: (np.abs(lat_) >= lo) & (np.abs(lat_) < hi)
        espec = {}
        nat_spec = {}
        for nm, (lo, hi) in (("40-60", (40, 60)), ("0-20", (0, 20))):
            rv = ref.isel(lev=ks[0]).values[:, band(lo, hi)].astype(np.float64)
            nat_spec[nm] = (np.abs(np.fft.rfft(rv, axis=-1) / rv.shape[-1]) ** 2).mean((0, 1))
            for v in variants:
                d = comps[v].isel(lev=ks[0]).values[:, band(lo, hi)].astype(np.float64) - rv
                espec[(nm, v)] = (np.abs(np.fft.rfft(d, axis=-1) / d.shape[-1]) ** 2).mean((0, 1))
        png = os.path.join(outdir, f"omega_check_{tag}_zonalspectra.png")
        _figure_zonal(zspec, zlats, (k850, k500), pnom, nat_spec, espec, variants, vl, ks[0], model, member, png)
        print(f"figure: {png}")
    print(f"total time {time.time() - t0:.0f} s")
    return rows


def run_daily(model, start="2000-06-01", end="2000-08-31", exp="historical", member=None,
              outdir=".", make_figure=True, deriv="auto", box=None, fullevel="sb", snapshot=None):
    """Daily means of computed omega on pressure levels vs the model's daily-mean wap (day or Eday).
    For models without a sub-daily (E3hrPt, 6hrPlev) or model-level daily (CFday) reference.
    The reference is the true daily mean, so the scores include the 6-hourly sampling error."""
    t0 = time.time()
    lev, lev_ok = _lev_members(model, exp, start, end)
    ref_dirs = None
    for table in ("day", "Eday"):
        d = _covering(version_dirs(model, table, "wap", exp), start, end)
        common = sorted(lev_ok & set(d))
        if member is not None:
            common = [m for m in common if m == member]
        if common:
            ref_dirs, dtable, member = d, table, common[0]
            break
    if ref_dirs is None:
        raise FileNotFoundError(f"{model}: no member with both 6hrLev (ua, va, ps) and day/Eday wap for {start}..{end}")
    print(f"omega_check {__version__} --daily: {model}  member {member}  {start}..{end}  reference {dtable} wap")

    parts = _common_time([open_files(lev[v][member], start, end) for v in ("ua", "va", "ps")])
    ds = xr.merge([parts[0], parts[1][["va"]], parts[2][["ps"]]], compat="override", join="override")
    A, B, vdim, flip = half_level_coeffs(ds)
    dnat = open_files(ref_dirs[member], start, end)
    print("  reference cell_methods:", dnat["wap"].attrs.get("cell_methods"))
    _time_info(ds, "6hrLev")
    _time_info(dnat, "reference")
    levels = sorted((float(p) for p in dnat["plev"].values if 10000. - 1 <= p <= 100000. + 1), reverse=True)
    print("  compared levels [hPa]:", [int(round(p / 100)) for p in levels])
    print(f"  {ds.sizes[vdim]} model levels; {ds.sizes['lat']}x{ds.sizes['lon']} grid; {ds.sizes['time']} time steps")

    if deriv in ("auto", "sh"):
        try:
            _sh_for(np.sort(ds["lat"].values), ds["lon"].values)
            deriv = "sh"
        except ValueError as e:
            if deriv == "sh":
                print(f"  spherical harmonics not possible ({e}); using centred differences")
            deriv = "cd"
    global LON_DERIV
    LON_DERIV = "fft" if deriv == "fft" else "cd"
    vlab = f"{fullevel}/{deriv}"
    print(f"  variant {vlab}")

    # ---- streaming trapezoidal daily means on pressure levels (window mean -> day of window start)
    acc, accp, cnt = {}, {}, {}
    mid = _midpoint_rule(ds)
    prev = prevp = prevday = None
    lat_ref = lon_ref = None
    for i in range(ds.sizes["time"]):
        s = ds[["ua", "va", "ps"]].isel(time=i).load()
        terms = flux_terms(s, A, B, vdim, flip, deriv)
        om, p = omega_from_terms(terms, A, B, vdim, fullevel, deriv)
        op = to_plev(om, p, vdim, levels).sortby("lat")
        ps = terms[2].sortby("lat")
        if lat_ref is None:
            lat_ref, lon_ref = op["lat"].values, op["lon"].values
        cur, curp = op.values.astype(np.float32), ps.values.astype(np.float32)
        day = str(s["time"].dt.strftime("%Y-%m-%d").values)
        if mid:                       # instant = mean of its own 6-h window, which lies in its day
            prev, prevp, prevday = cur, curp, day
        if prev is not None:
            if prevday not in acc:
                acc[prevday], accp[prevday], cnt[prevday] = np.zeros_like(cur), np.zeros_like(curp), 0
            # NaN (below ground) in either instant -> NaN in the window mean
            acc[prevday] += 0.5 * (prev + cur)
            accp[prevday] += 0.5 * (prevp + curp)
            cnt[prevday] += 1
        prev, prevp, prevday = cur, curp, day
        if i % 80 == 0:
            print(f"    step {i + 1}/{ds.sizes['time']}  ({time.time() - t0:.0f} s)", flush=True)
    days = sorted(d for d in cnt if cnt[d] == 4)
    if not days:
        raise ValueError("no complete days")
    comp = xr.DataArray(np.stack([acc[d] / 4 for d in days]), dims=("time", "plev", "lat", "lon"),
                        coords={"time": days, "plev": np.asarray(levels, float), "lat": lat_ref, "lon": lon_ref})
    psd = xr.DataArray(np.stack([accp[d] / 4 for d in days]), dims=("time", "lat", "lon"),
                       coords={"time": days, "lat": lat_ref, "lon": lon_ref})
    print(f"  daily means of computed omega for {len(days)} days ({time.time() - t0:.0f} s)")

    ref = dnat["wap"].sel(plev=levels, method="nearest").astype(np.float32).load()
    ref = ref.assign_coords(plev=np.asarray(levels, float), time=_ref_dates(ref, dnat)).transpose("time", "plev", "lat", "lon")
    ref = ref.sortby("lat")
    ref = ref.drop_vars([c for c in ref.coords if c not in ("time", "plev", "lat", "lon")])
    comp, psd = _match_grid(comp, ref), _match_grid(psd, ref)
    comp, psd, ref = _on_dates(comp, psd, ref)
    if UNREGRID and REREGRID:
        # native daily wap went through CNRM's output regridding: same treatment for the computed field
        comp = reregrid(comp.fillna(0.0)).where(comp.notnull())
        print(f"  computed omega passed through the same output regridding as the native field ({len(UNREGRID)} rows)")
    print(f"  compared days: {comp.sizes['time']}")
    _lag_scan(comp, ref, None)

    land = _load_land(model, exp, member)
    if land is not None:
        land = _match_grid(land, ref)
    masks = _region_masks(psd.mean("time"), land, box)
    focus = next(k for k in masks if k.startswith("NH midlat"))
    regions = list(masks)
    key_regions = ["global", "tropics", "subtropics", "midlatitudes", focus] + (["box"] if box is not None else [])

    lab = "daily total"
    rows = []
    for reg in regions:
        rows += _scores(comp, ref, lab, masks[reg], reg, extremes=reg in key_regions)
    for reg in key_regions:
        rows += _scores_agl(comp, ref, psd, lab, masks[reg], reg)

    print(f"\nDAILY  {model} {member}: daily-mean computed omega ({vlab}) vs native {dtable} wap (true daily mean)")
    _print_scores([x for x in rows if x["region"] == "global" and isinstance(x["agl_bin"], float)],
                  "GLOBAL (|lat|<70)")
    _print_regional([x for x in rows if isinstance(x["agl_bin"], float)], [lab], regions)
    _print_agl([x for x in rows if isinstance(x["agl_bin"], str)])
    _print_bias(rows, key_regions, comps=(lab,))

    g = [x for x in rows if x["region"] == "global" and isinstance(x["agl_bin"], float) and np.isfinite(x["r"])]
    agl = {x["region"]: x for x in rows if x["agl_bin"] == "0-50" and np.isfinite(x["r"])}
    L = [f"omega_check {__version__} --daily: {model} {member}, {start}..{end}, reference {dtable} wap, variant {vlab}", "",
         "DAILY MEANS ON PRESSURE LEVELS (computed from 6-hourly instants vs the model's true daily mean), global |lat|<70",
         "  includes the error of estimating a daily mean from four 6-hourly snapshots (sampling)",
         "  hPa    r time   r space   slope   RMSE in % of sigma"]
    for x in g:
        L.append(f"  {x['plev_hPa']:4d}    {x['r_time']:6.3f}   {x['r']:6.3f}   {x['slope']:6.3f}   {x['rmse_pct']:5.0f} %")
    L += ["", "WITHIN 50 hPa OF THE SURFACE"]
    for reg, x in agl.items():
        L.append(f"  {reg:16s} r = {x['r']:.3f}, slope {x['slope']:.3f}, RMSE {x['rmse_pct']:.0f} % of sigma")
    hp = [x["plev_hPa"] for x in g]
    tex = (f"Daily means of the diagnosed $\\omega$ of {model}, computed from the 6-hourly instants, agree with the "
           f"model's native daily-mean $\\omega$ (table \\texttt{{{dtable}}}) with temporal correlations of "
           f"{_rng([x['r_time'] for x in g], '{:.2f}')} and regression slopes of {_rng([x['slope'] for x in g], '{:.2f}')} "
           f"between {max(hp)} and {min(hp)}\\,hPa." if g else "")
    L += ["", "LATEX (numbers from this run):", tex]
    text = "\n".join(L)
    os.makedirs(outdir, exist_ok=True)
    tag = f"{model}_daily" + ("_unregrid" + ("" if REREGRID else "_noreregrid") if UNREGRID else "")
    with open(os.path.join(outdir, f"omega_check_{tag}_summary.txt"), "w") as f:
        f.write(text + "\n")
    with open(os.path.join(outdir, f"omega_check_{tag}.csv"), "w") as f:
        f.write("model,member,reference,start,end,variant," + ",".join(COLS) + "\n")
        for x in rows:
            vals = [("" if isinstance(x[c], float) and not np.isfinite(x[c]) else
                     str(round(x[c], 5)) if isinstance(x[c], float) else _q(x[c])) for c in COLS]
            f.write(f"{model},{member},{dtable},{start},{end},{vlab}," + ",".join(vals) + "\n")
    print("\n" + text + "\n")
    if make_figure:
        for fname in _figure_maps(comp, ref, land, model, member, dtable, "daily mean", vlab,
                                  os.path.join(outdir, f"omega_check_{tag}"), snapshot=snapshot):
            print(f"figure: {fname}")
        png = os.path.join(outdir, f"omega_check_{tag}_nearsurface.png")
        _figure_agl(rows, model, member, key_regions, lab, png)
        print(f"figure: {png}")
    print(f"total time {time.time() - t0:.0f} s")
    return rows


# ================================================================ main
VARIANTS_COMPARE = [("mean", "cd"), ("sb", "cd"), ("mean", "sh"), ("sb", "sh")]


def _vtag(v):
    return f"{v[0]}/{v[1]}"


def _tkey(t):
    """(year, month, day, minutes of day rounded to TOL) for numpy datetime64 or cftime."""
    if isinstance(t, np.datetime64):
        t = datetime.datetime.utcfromtimestamp((t - np.datetime64("1970-01-01T00:00:00")) / np.timedelta64(1, "s"))
    step = TOL.total_seconds() / 60
    m = int(round((t.hour * 60 + t.minute + t.second / 60) / step) * step)
    return (t.year, t.month, t.day, m)


def _match_time(om, ps_all, nat, nat_full, kind):
    """Align computed omega (and ps) with the reference in time and space."""
    if kind == "instantaneous":
        ref6 = nat_full.sel(time=om["time"].values, method="nearest", tolerance=TOL)
        comp = om.assign_coords(time=ref6["time"].values)
        ps6 = ps_all.assign_coords(time=ref6["time"].values)
    else:
        bname = nat["time"].attrs.get("bounds")
        if bname not in nat.variables:
            bname = next(n for n in nat.variables if "time" in n and "bnd" in n.replace("bounds", "bnd"))
        # match by calendar date and time of day (rounded to TOL), so that cftime and numpy
        # datetimes, or 'noleap' and 'proleptic_gregorian' of the same date, compare equal
        idx = {_tkey(t): k for k, t in enumerate(om["time"].values)}
        keep, means, pmeans = [], [], []
        for i, (b0, b1) in enumerate(nat[bname].values):
            k0, k1 = idx.get(_tkey(b0)), idx.get(_tkey(b1))
            if k0 is None or k1 is None:
                continue
            means.append(om.isel(time=[k0, k1]).mean("time"))
            pmeans.append(ps_all.isel(time=[k0, k1]).mean("time"))
            keep.append(i)
        if not keep:
            # instants at the window CENTRES (e.g. NorESM2-MM: 6hrLev at 03, 09, 15, 21 UTC, 6hrPlev
            # windows 00-06, ...): the instant itself is the (midpoint-rule) estimate of the window mean
            for i, (b0, b1) in enumerate(nat[bname].values):
                k = idx.get(_tkey(b0 + (b1 - b0) / 2))
                if k is None:
                    continue
                means.append(om.isel(time=k))
                pmeans.append(ps_all.isel(time=k))
                keep.append(i)
            if keep:
                print(f"  note: 6hrLev instants lie at the centres of the 6-h mean windows; each window mean "
                      f"is compared with the single instant at its centre ({len(keep)} windows)")
        if not keep:
            b = nat[bname].values
            raise ValueError("no reference window is bounded by two 6hrLev instants; first window "
                             f"{b[0][0]} .. {b[0][1]}, 6hrLev times {om['time'].values[0]} .. "
                             f"{om['time'].values[-1]}")
        ref6 = nat_full.isel(time=keep)
        comp = xr.concat(means, dim="time").assign_coords(time=ref6["time"].values)
        ps6 = xr.concat(pmeans, dim="time").assign_coords(time=ref6["time"].values)
    comp = _match_grid(comp, ref6).transpose("time", "plev", "lat", "lon")
    ps6 = _match_grid(ps6, ref6).transpose("time", "lat", "lon")
    return comp, ps6, ref6.transpose("time", "plev", "lat", "lon")


def _figure_compare(rows, variants, model, member, fname):
    """Small multiples: slope and RMSE (% of sigma) vs pressure, one line per variant."""
    plt = _plt()
    palette = ["#2b2b29", "#2a78d6", "#eb6834", "#1baf7a"]
    markers = ["o", "s", "^", "D"]
    ink, muted, grid = "#2b2b29", "#6b6a64", "#e4e3df"
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True, sharey=True)
    for k, v in enumerate(variants):
        rr = [x for x in rows if x["grid"] == f"{_vtag(v)} | 6-hourly instants" and x["region"] == "global"
              and isinstance(x["agl_bin"], float)]
        if not rr:
            continue
        y = [x["plev_hPa"] for x in rr]
        for ax, key in zip(axs, ("slope", "rmse_pct")):
            xv = [x[key] for x in rr]
            ax.plot(xv, y, color=palette[k], lw=2, marker=markers[k], ms=8, mec="white", mew=1.5,
                    label=_vtag(v))
    axs[0].axvline(1.0, color=muted, lw=1, ls="--")
    axs[0].text(1.0, 1.0, " no bias", color=muted, fontsize=8, va="top",
                transform=axs[0].get_xaxis_transform())
    axs[0].set_xlabel("regression slope computed vs native  (< 1: too weak)", color=ink)
    axs[1].set_xlabel("RMSE in % of the standard deviation of native omega", color=ink)
    axs[0].set_ylabel("pressure [hPa]", color=ink)
    axs[0].invert_yaxis()
    for ax in axs:
        ax.grid(color=grid, lw=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(muted)
        ax.tick_params(colors=muted)
    axs[1].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  title="full level / derivative", title_fontsize=8)
    fig.suptitle(f"{model} {member}: effect of full-level formula and derivative method (global, 6-hourly)",
                 color=ink)
    fig.savefig(fname, dpi=120)
    plt.close(fig)


def _compare_report(oms, variants, ps_all, nat, nat_full, kind, masks, focus, box,
                    model, member, table, start, end, outdir, make_figure):
    regions = ["global", "tropics", "subtropics", "midlatitudes", focus] + (["box"] if box is not None else [])
    rows, ref6, ps6 = [], None, None
    for v in variants:
        comp, ps6, ref6 = _match_time(oms[v], ps_all, nat, nat_full, kind)
        lab = _vtag(v)
        for reg in regions:
            rows += _scores(comp, ref6, f"{lab} | 6-hourly instants", masks[reg], reg,
                            extremes=reg in ("global", focus))
        for reg in ("global", focus):
            rows += _scores_agl(comp, ref6, ps6, f"{lab} | 6-hourly instants", masks[reg], reg)
        if kind == "instantaneous":
            cd_, rd_ = _on_dates(_daily_from_instants(comp, 4), _daily_from_instants(ref6, 4))
        else:
            cd_, rd_ = _on_dates(_daily_from_windows(comp, 4), _daily_from_windows(ref6, 4))
        rows += _scores(cd_, rd_, f"{lab} | daily method", masks["global"], "global")
        print(f"  scored {lab}")

    get = lambda v, comp, reg="global", lev=None, agl=None: next(
        (x for x in rows if x["grid"] == f"{_vtag(v)} | {comp}" and x["region"] == reg
         and ((agl is None and isinstance(x["agl_bin"], float) and x["plev_hPa"] == lev)
              or (agl is not None and x["agl_bin"] == agl))), None)
    levs = sorted({x["plev_hPa"] for x in rows if isinstance(x["agl_bin"], float)}, reverse=True)
    W = 24
    head = "  " + f"{'':>12s}" + "".join(f"{_vtag(v):>{W}s}" for v in variants)

    def table_(title, keyrows, cell):
        print(f"\n{title}")
        print(head)
        for label, kw in keyrows:
            cells = []
            for v in variants:
                x = get(v, **kw)
                cells.append(cell(x) if x is not None and np.isfinite(x["r"]) else "--")
            print("  " + f"{label:>12s}" + "".join(f"{c:>{W}s}" for c in cells))

    print(f"\nCOMPARISON  {model} {member} {start}..{end}, reference {table} ({kind})")
    print("  variants: full level = mean (half-level average) or sb (Simmons & Burridge);"
          " derivative = cd (centred differences) or sh (spherical harmonics)")
    rs = lambda x: f"{x['r']:.3f} / {x['slope']:.3f} / {x['rmse_pct']:.0f}%"
    table_("GLOBAL, 6-hourly instants   cells: r / slope / RMSE in % of sigma",
           [(f"{lv} hPa", dict(comp="6-hourly instants", lev=lv)) for lv in levs], rs)
    table_("GLOBAL, strongest 5 % ascent / descent reproduced (100 % = fully)",
           [(f"{lv} hPa", dict(comp="6-hourly instants", lev=lv)) for lv in levs],
           lambda x: f"{100 * x['asc_ratio']:.0f}% / {100 * x['desc_ratio']:.0f}%")
    table_("NEAR THE SURFACE (height above ground)   cells: r / slope / RMSE in % of sigma",
           [(f"{reg[:7]} {b}", dict(comp="6-hourly instants", reg=reg, agl=b))
            for reg in ("global", focus) for b in ("0-50", "50-100", "100-200")], rs)
    show = [lv for lv in levs if lv in (925, 850, 700, 500, 250)] or levs
    table_("BY REGION, slope (< 1: too weak)",
           [(f"{reg[:7]} {lv}", dict(comp="6-hourly instants", reg=reg, lev=lv))
            for reg in regions[1:] for lv in show], lambda x: f"{x['slope']:.3f}")
    table_("DAILY MEANS (method), global   cells: r time / slope / RMSE in % of sigma",
           [(f"{lv} hPa", dict(comp="daily method", lev=lv)) for lv in levs],
           lambda x: f"{x['r_time']:.3f} / {x['slope']:.3f} / {x['rmse_pct']:.0f}%")

    # ---- attribution of the underestimation
    base = ("mean", "cd")
    lines = ["", "WHAT EACH FIX DOES  (change in slope; positive = less underestimation)"]
    for v, name in ((("sb", "cd"), "issue 2: Simmons-Burridge full levels"),
                    (("mean", "sh"), "issue 3: spherical-harmonic derivatives"),
                    (("sb", "sh"), "both together")):
        if v not in variants or base not in variants:
            continue
        parts = []
        for lab, kw in [(f"{lv} hPa", dict(comp="6-hourly instants", lev=lv)) for lv in show] + \
                       [("0-50 hPa AGL", dict(comp="6-hourly instants", agl="0-50"))]:
            a, b = get(v, **kw), get(base, **kw)
            if a is not None and b is not None and np.isfinite(a["slope"]) and np.isfinite(b["slope"]):
                parts.append(f"{lab} {a['slope'] - b['slope']:+.3f}")
        lines.append(f"  {name:42s} " + ", ".join(parts))
    print("\n".join(lines))

    os.makedirs(outdir, exist_ok=True)
    csv = os.path.join(outdir, f"omega_check_{model}_compare.csv")
    with open(csv, "w") as f:
        f.write("model,member,reference,start,end," + ",".join(COLS) + "\n")
        for x in rows:
            vals = [("" if isinstance(x[c], float) and not np.isfinite(x[c]) else
                     str(round(x[c], 5)) if isinstance(x[c], float) else _q(x[c])) for c in COLS]
            f.write(f"{model},{member},{table},{start},{end}," + ",".join(vals) + "\n")
    print(f"\nscores: {csv}")
    if make_figure:
        png = os.path.join(outdir, f"omega_check_{model}_compare.png")
        _figure_compare(rows, variants, model, member, png)
        print(f"figure: {png}")
    return rows


def run(model, start="2000-06-01", end="2000-08-31", levels=None, exp="historical",
        member=None, outdir=".", make_figure=True, deriv="auto", box=None, fullevel="sb",
        compare=False, snapshot=None):
    global LON_DERIV
    is_default = (deriv == "auto" and fullevel == "sb")
    LON_DERIV = "fft" if (deriv == "fft" and not compare) else "cd"
    t0 = time.time()

    # ---- member, sub-daily reference and daily reference
    lev = {v: version_dirs(model, "6hrLev", v, exp) for v in ("ua", "va", "ps")}
    lev_members = set(lev["ua"]) & set(lev["va"]) & set(lev["ps"])
    if not lev_members:
        raise FileNotFoundError(f"{model}: no member with ua, va, ps in 6hrLev/{exp} under {POOL}")
    daily_dirs = {t: version_dirs(model, t, "wap", exp) for t in ("day", "Eday")}
    has_daily = set(daily_dirs["day"]) | set(daily_dirs["Eday"])
    ref = None
    for table, kind in (("E3hrPt", "instantaneous"), ("6hrPlev", "6-h mean")):
        w = version_dirs(model, table, "wap", exp)
        common = sorted(lev_members & set(w))
        if member is not None:
            common = [m for m in common if m == member]
        if common:
            pref = [m for m in common if m in has_daily] or common
            ref = (table, kind, pref[0], w[pref[0]])
            break
    if ref is None:
        raise FileNotFoundError(f"{model}: no member with both 6hrLev and E3hrPt/6hrPlev wap")
    table, kind, member, wapdir = ref
    dtable = next((t for t in ("day", "Eday") if member in daily_dirs[t]), None)
    variants = list(VARIANTS_COMPARE) if compare else [(fullevel, deriv)]
    print(f"omega_check {__version__}: {model}  member {member}  {start}..{end}")
    print(f"reference: {table} wap ({kind});  daily reference: "
          f"{dtable + ' wap' if dtable else 'none for this member'}")

    # ---- sub-daily reference and levels
    nat = open_files(wapdir, start, end)
    print("  reference cell_methods:", nat["wap"].attrs.get("cell_methods"))
    plev_ref = nat["plev"].values
    if levels is None:
        levels = sorted((float(p) for p in plev_ref if 10000. - 1 <= p <= 100000. + 1), reverse=True)
    else:
        levels = [p for p in levels if np.min(np.abs(plev_ref - p)) < 100]
    if not levels:
        raise ValueError(f"no usable levels in reference levels {plev_ref}")
    print("  compared levels [hPa]:", [int(round(p / 100)) for p in levels])
    nat = nat.sel(plev=levels, method="nearest")

    # ---- 6hrLev input
    parts = _common_time([open_files(lev[v][member], start, end) for v in ("ua", "va", "ps")])
    ds = xr.merge([parts[0], parts[1][["va"]], parts[2][["ps"]]],
                  compat="override", join="override")
    print("  6hrLev ua cell_methods:", ds["ua"].attrs.get("cell_methods"))
    A, B, vdim, flip = half_level_coeffs(ds)
    print(f"  {ds.sizes[vdim]} model levels, top {A[0] / 100 + B[0] * 1000:.2f} hPa,"
          f" {'flipped' if flip else 'top->surface'};  {ds.sizes['lat']}x{ds.sizes['lon']} grid,"
          f" {ds.sizes['time']} time steps")

    # ---- derivative method: auto -> spherical harmonics on Gaussian grids, else centred differences
    if any(v[1] == "auto" for v in variants):
        try:
            _sh_for(np.sort(ds["lat"].values), ds["lon"].values)
            auto = "sh"
            why = "Gaussian grid (spectral model)"
        except ValueError as e:
            auto, why = "cd", f"{e} (grid-point model)"
        print(f"  derivatives: auto -> {auto}  [{why}]")
        variants = [(f, auto if d == "auto" else d) for f, d in variants]
    print("  variants:", ", ".join(_vtag(v) for v in variants))
    # ---- spherical harmonics only on Gaussian grids
    if any(v[1] == "sh" for v in variants):
        try:
            _sh_for(np.sort(ds["lat"].values), ds["lon"].values)
            _sh = _sh_for(np.sort(ds['lat'].values), ds['lon'].values)
            print(f"  spherical harmonics: Gaussian grid {_sh.nlat}x{_sh.nlon}, derivatives truncated at "
                  f"T{_trunc_for(_sh)} (grid allows T{_sh.N}; --trunc {SH_TRUNC})")
        except ValueError as e:
            print(f"  spherical harmonics not possible ({e}); using centred differences instead")
            variants = list(dict.fromkeys((f, "cd") for f, d in variants))

    # ---- computed omega, one time step at a time
    steps = {v: [] for v in variants}
    ps_steps = []
    derivs = sorted({v[1] for v in variants})
    for i in range(ds.sizes["time"]):
        s = ds[["ua", "va", "ps"]].isel(time=i).load()
        terms = {d: flux_terms(s, A, B, vdim, flip, d) for d in derivs}
        for v in variants:
            om_i, p_i = omega_from_terms(terms[v[1]], A, B, vdim, v[0], v[1])
            steps[v].append(to_plev(om_i, p_i, vdim, levels).astype(np.float32).sortby("lat"))
        ps_steps.append(s["ps"].astype(np.float32).sortby("lat"))
    tvals = ds["time"].values
    oms = {v: xr.concat(steps[v], dim="time").assign_coords(time=tvals, plev=np.asarray(levels, float))
           for v in variants}
    ps_all = xr.concat(ps_steps, dim="time").assign_coords(time=tvals)
    print(f"  omega computed ({time.time() - t0:.0f} s)")

    nat_full = nat["wap"].astype(np.float32).load().assign_coords(plev=np.asarray(levels, float))
    comp, ps6, ref6 = _match_time(oms[variants[0]], ps_all, nat, nat_full, kind)

    # ---- regions
    land = _load_land(model, exp, member)
    if land is not None:
        land = _match_grid(land, ref6)
    print("  land fraction:", "sftlf found" if land is not None else "not found (NH midlat lowland includes ocean)")
    masks = _region_masks(ps6.mean("time"), land, box)
    regions = list(masks)
    focus = next(k for k in regions if k.startswith("NH midlat"))
    key_regions = ["global", "tropics", "subtropics", "midlatitudes", focus] + (["box"] if box is not None else [])
    for k in regions:
        print(f"    region {k:24s}: {int(masks[k].sum()):6d} grid points")

    if compare:
        rows = _compare_report(oms, variants, ps_all, nat, nat_full, kind, masks, focus, box,
                               model, member, table, start, end, outdir, make_figure)
        print(f"total time {time.time() - t0:.0f} s")
        return rows

    rows = []
    # ---- sub-daily scores
    for reg in regions:
        rows += _scores(comp, ref6, "6-hourly instants", masks[reg], reg,
                        extremes=reg in key_regions)
    cg = dict(lat=4, lon=4, boundary="trim")
    both = comp.notnull() & ref6.notnull()                     # same points before averaging
    cc, rc = comp.where(both).coarsen(**cg).mean(), ref6.where(both).coarsen(**cg).mean()
    gm = (np.abs(cc["lat"]) < 70) & xr.ones_like(cc.isel(time=0, plev=0), dtype=bool)
    rows += _scores(cc, rc, "6-hourly, 4x4 coarse", gm, "global")
    agl_regions = key_regions
    for reg in agl_regions:
        rows += _scores_agl(comp, ref6, ps6, "6-hourly instants", masks[reg], reg)

    # ---- daily scores
    daily = None
    if kind == "instantaneous":
        comp_d, ref6_d, ps_d6 = (_daily_from_instants(x, 4) for x in (comp, ref6, ps6))
        lab6 = "native@6h"
    else:
        comp_d, ref6_d, ps_d6 = (_daily_from_windows(x, 4) for x in (comp, ref6, ps6))
        lab6 = "native 6h-means"
    comp_d, ref6_d, ps_d6 = _on_dates(comp_d, ref6_d, ps_d6)
    for reg in regions:
        rows += _scores(comp_d, ref6_d, "daily method", masks[reg], reg,
                        extremes=reg in key_regions)
    for reg in agl_regions:
        rows += _scores_agl(comp_d, ref6_d, ps_d6, "daily method", masks[reg], reg)
    if dtable:
        try:
            dnat = open_files(daily_dirs[dtable][member], start, end)
            dlev = [p for p in levels if np.min(np.abs(dnat["plev"].values - p)) < 100]
            if len(dlev) < len(levels):
                print(f"  note: {dtable} lacks some levels; daily-total comparison on "
                      f"{[int(round(p / 100)) for p in dlev]} hPa")
            nat_d = dnat["wap"].sel(plev=dlev, method="nearest").astype(np.float32).load()
            nat_d = nat_d.assign_coords(plev=np.asarray(dlev, float), time=_ref_dates(nat_d, dnat))
            nat_d = nat_d.transpose("time", "plev", "lat", "lon")
            c_d = _match_grid(comp_d.sel(plev=dlev), nat_d)
            r6_d = _match_grid(ref6_d.sel(plev=dlev), nat_d)
            p_d = _match_grid(ps_d6, nat_d)
            c_d, r6_d, p_d, nat_d = _on_dates(c_d, r6_d, p_d, nat_d)
            dm = _masks_on(masks, nat_d)
            for reg in regions:
                ex = reg in key_regions
                rows += _scores(r6_d, nat_d, f"sampling {lab6} vs daily", dm[reg], reg)
                rows += _scores(c_d, nat_d, "daily total", dm[reg], reg, extremes=ex)
            if kind == "instantaneous":
                r3_d = _match_grid(_daily_from_instants(nat_full.sel(plev=dlev), 8), nat_d)
                r3_d, nd2 = _on_dates(r3_d, nat_d)
                rows += _scores(r3_d, nd2, "sampling native@3h vs daily", dm["global"], "global")
            for reg in agl_regions:
                rows += _scores_agl(c_d, nat_d, p_d, "daily total", dm[reg], reg)
            daily = (c_d, r6_d, nat_d)
        except (FileNotFoundError, ValueError, KeyError) as e:
            print(f"  daily reference skipped: {e}")

    # ---- print
    vlab = _vtag(variants[0])
    glob_rows = [x for x in rows if x["region"] == "global" and isinstance(x["agl_bin"], float)]
    sub = [x for x in glob_rows if x["grid"].startswith("6-hourly")]
    _print_scores(sub, f"SUB-DAILY (global, |lat|<70)  {model} {member}: computed ({vlab}) vs {table} ({kind})")
    day_rows = [x for x in glob_rows if not x["grid"].startswith("6-hourly")]
    if day_rows:
        _print_scores(day_rows, "DAILY MEANS (global)  method = computed vs native, both from 6-hourly; "
                                "sampling = native 6-hourly vs true daily; total = computed vs true daily")
    _print_regional([x for x in rows if isinstance(x["agl_bin"], float)],
                    ["6-hourly instants", "daily method", "daily total"], regions)
    _print_agl([x for x in rows if isinstance(x["agl_bin"], str)])
    _print_bias(rows, key_regions)

    # ---- output files
    os.makedirs(outdir, exist_ok=True)
    fl, dv = variants[0]
    tag = model + ("" if is_default else f"_{fl}_{dv}")
    csv = os.path.join(outdir, f"omega_check_{tag}.csv")
    with open(csv, "w") as f:
        f.write("model,member,reference,daily_reference,start,end,variant," + ",".join(COLS) + "\n")
        for x in rows:
            vals = []
            for c in COLS:
                v = x[c]
                if isinstance(v, float):
                    vals.append("" if not np.isfinite(v) else str(round(v, 5)))
                else:
                    vals.append(_q(v))
            f.write(f"{model},{member},{table},{dtable or ''},{start},{end},{vlab}," + ",".join(vals) + "\n")
    print(f"\nscores: {csv}")
    text = _summary(rows, model, member, table, kind, dtable, start, end, focus)
    text = text.replace(f"reference {table} ({kind})", f"reference {table} ({kind}); variant {vlab}", 1)
    txt = os.path.join(outdir, f"omega_check_{tag}_summary.txt")
    with open(txt, "w") as f:
        f.write(text + "\n")
    print("\n" + text + "\n")
    print(f"summary: {txt}")
    if make_figure:
        png = os.path.join(outdir, f"omega_check_{tag}.png")
        _figure(comp, ref6, model, member, table, png)
        print(f"figure: {png}")
        if daily is not None:
            png_d = os.path.join(outdir, f"omega_check_{tag}_daily.png")
            _figure_daily(*daily, model, member, dtable, png_d)
            print(f"figure: {png_d}")
        for comp_name, suffix in (("6-hourly instants", "6h"), ("daily total", "daily")):
            if any(x["grid"] == comp_name and isinstance(x["agl_bin"], str) for x in rows):
                png_a = os.path.join(outdir, f"omega_check_{tag}_nearsurface_{suffix}.png")
                _figure_agl(rows, model, member, agl_regions, comp_name, png_a)
                print(f"figure: {png_a}")
        for fname in _figure_maps(comp, ref6, land, model, member, table, kind, vlab,
                                  os.path.join(outdir, f"omega_check_{tag}"), snapshot=snapshot):
            print(f"figure: {fname}")
        if box is not None:
            for fname in _figure_maps(comp, ref6, land, model, member, table, kind, vlab,
                                      os.path.join(outdir, f"omega_check_{tag}_box"),
                                      snapshot=snapshot, box=box):
                print(f"figure: {fname}")
    print(f"total time {time.time() - t0:.0f} s")
    return rows


def regrid_check(model, exp="historical", member=None, outdir="."):
    """How was the 6hrLev output regridded? Per latitude row: distinct values, equal neighbours
    (staircase = conservative / nearest-neighbour), vanishing second differences (linear), plus
    zonal power spectra by latitude."""
    lev = version_dirs(model, "6hrLev", "ua", exp)
    if not lev:
        raise FileNotFoundError(f"{model}: no 6hrLev ua")
    member = member if member in lev else sorted(lev)[0]
    f = sorted(glob.glob(os.path.join(lev[member], "*.nc")))[0]
    ds = xr.open_dataset(f)
    L = [f"omega_check {__version__} --regrid-check: {f}"]
    for k in ("grid", "grid_label", "nominal_resolution"):
        L.append(f"  {k}: {ds.attrs.get(k)}")
    u = ds["ua"].isel(time=0).load()
    vd = [d for d in u.dims if d not in ("lat", "lon")][0]
    nlon = u.sizes["lon"]
    L += ["", f"  per latitude row, all {u.sizes[vd]} levels pooled ({nlon} longitudes):",
          "    lat   nlon*cos  distinct values  equal neighbours  2nd diff ~ 0",
          "    (many equal neighbours: staircase regridding; 2nd diff ~ 0: linear interpolation)"]
    for lat in (0, 20, 30, 40, 45, 50, 55, 60, 70, 80):
        x = u.sel(lat=lat, method="nearest").transpose(vd, "lon").values.astype(np.float64)
        sc = np.abs(x).max(axis=1, keepdims=True) + 1e-12
        d1 = np.diff(np.concatenate([x, x[:, :1]], 1), axis=1) / sc
        d2 = np.diff(np.concatenate([x[:, -1:], x, x[:, :1]], 1), n=2, axis=1) / sc
        uniq = np.median([np.unique(np.round(r / s_, 7)).size for r, s_ in zip(x, sc[:, 0])])
        L.append(f"   {lat:4d}   {nlon * np.cos(np.deg2rad(lat)):6.0f}   {uniq:12.0f}      "
                 f"{np.mean(np.abs(d1) < 1e-6):10.2f}      {np.mean(np.abs(d2) < 1e-5):10.2f}")
    os.makedirs(outdir, exist_ok=True)
    # ---- which latitude rows were interpolated from fewer points, and how?
    print("  detecting the native points per latitude circle (1-3 minutes) ...", flush=True)
    t1 = time.time()
    us = ds["ua"].isel(time=slice(0, 2)).sortby("lat").load()
    vdn = [d for d in us.dims if d not in ("time", "lat", "lon")][0]
    nl = us.sizes[vdn]
    ks_ = sorted({0, nl // 4, nl // 2, 3 * nl // 4, nl - 1})
    x = us.isel({vdn: ks_}).transpose("time", vdn, "lat", "lon").values.reshape(-1, us.sizes["lat"], nlon)
    rows = detect_reduced_grid(x)
    latv = us["lat"].values
    L += ["", f"  NATIVE POINTS PER LATITUDE CIRCLE (smallest n with fit < 1e-5 relative rms and >30x better than n-1; {time.time() - t1:.0f} s)",
          "    |lat| range        n points   scheme         offset   max residual   min sharpness"]
    i = 0
    half = len(rows) // 2
    rs_ = rows[half:]                                        # northern rows (ascending lat)
    while i < len(rs_):
        j = i
        while j + 1 < len(rs_) and (rs_[j + 1]["n"], rs_[j + 1]["scheme"], rs_[j + 1]["offset"]) == (rs_[i]["n"], rs_[i]["scheme"], rs_[i]["offset"]):
            j += 1
        grp = rs_[i:j + 1]
        L.append(f"    {latv[grp[0]['row']]:6.1f} - {latv[grp[-1]['row']]:5.1f}     {grp[0]['n']:6d}     "
                 f"{(grp[0]['scheme'] if grp[0]['scheme'] != 'none' else 'full circle'):12s}   {grp[0]['offset']:5.2f}   "
                 f"{max(g['residual'] for g in grp):.1e}      {min(g['sharpness'] for g in grp):.0f}")
        i = j + 1
    sym = all((rows[k]["n"], rows[k]["scheme"], rows[k]["offset"]) == (rows[-1 - k]["n"], rows[-1 - k]["scheme"], rows[-1 - k]["offset"]) for k in range(half))
    L.append(f"    (northern hemisphere shown; southern hemisphere {'identical' if sym else 'DIFFERENT - see json'})")
    nred = sum(1 for r in rows if r["scheme"] != "none" and r["n"] < nlon)
    if nred:
        L += ["", f"  -> {nred} of {len(rows)} rows were interpolated from fewer points. The native values can be",
              "     recovered exactly:  python omega_check_cnrm.py " + model + " --cfday --unregrid " +
              os.path.join(outdir, f"regrid_{model}.json")]
    else:
        L += ["", "  -> no row is reproduced by interpolation from fewer points (other regridding, or native grid)."]
    import json
    with open(os.path.join(outdir, f"regrid_{model}.json"), "w") as fh:
        json.dump({"model": model, "file": f, "nlon": nlon, "lat": [float(v_) for v_ in latv], "rows": rows}, fh, indent=1)
    text = "\n".join(L)
    print(text)
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, f"regrid_check_{model}.txt"), "w") as fh:
        fh.write(text + "\n")
    plt = _plt()
    uu = ds["ua"].isel(time=slice(0, 8)).load()
    lats = (0, 30, 45, 50, 55, 60, 70)
    fig, ax = plt.subplots(figsize=(8, 5))
    for j, lat in enumerate(lats):
        x = uu.sel(lat=[lat, -lat], method="nearest").values.astype(np.float64)
        P = (np.abs(np.fft.rfft(x, axis=-1) / nlon) ** 2).reshape(-1, nlon // 2 + 1).mean(0)
        col = plt.cm.viridis(j / (len(lats) - 1) * 0.9)
        ax.semilogy(np.arange(1, P.size), P[1:], color=col, lw=1.3, label=f"{lat}\u00b0")
        ax.axvline(nlon / 2 * np.cos(np.deg2rad(lat)), color=col, ls=":", lw=0.8)
    ax.set_xlabel("zonal wavenumber m"); ax.set_ylabel("power [(m/s)$^2$]")
    ax.set_title(f"{model} 6hrLev ua: zonal spectra by |lat| (dotted: m = nlon/2 cos lat)", fontsize=10)
    ax.legend(title="|lat|", frameon=False, fontsize=8)
    ax.grid(True, color="#e6e6e6", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    png = os.path.join(outdir, f"regrid_check_{model}.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    print(f"\nfigure: {png}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="experimental omega check for regridded (gr) models")
    ap.add_argument("model")
    ap.add_argument("--start", default="2000-06-01")
    ap.add_argument("--end", default="2000-08-31")
    ap.add_argument("--levels", type=float, nargs="+", default=None,
                    help="pressure levels in hPa (default: all reference levels 1000-100 hPa)")
    ap.add_argument("--box", type=float, nargs=4, default=None, metavar=("LON0", "LON1", "LAT0", "LAT1"),
                    help="extra region, e.g. --box -10 30 40 55")
    ap.add_argument("--member", default=None)
    ap.add_argument("--exp", default="historical")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--no-figure", action="store_true")
    ap.add_argument("--fullevel", choices=("mean", "sb"), default="sb",
                    help="full-level omega: Simmons & Burridge (default) or half-level mean")
    ap.add_argument("--deriv", choices=("auto", "cd", "fft", "sh"), default="auto",
                    help="horizontal derivatives: auto (default: spherical harmonics on Gaussian grids, "
                         "else centred differences), cd, fft, sh")
    ap.add_argument("--fft", action="store_true", help="same as --deriv fft")
    ap.add_argument("--ref-shift", type=int, default=0, metavar="DAYS",
                    help="add DAYS to the date labels of daily-mean references (diagnosis of time-label errors)")
    ap.add_argument("--snapshot", default=None,
                    help="time for the snapshot maps, e.g. 2000-06-15T12 (default: first time)")
    ap.add_argument("--cfday", action="store_true",
                    help="compare daily means on model levels with CFday wap (no pressure-level interpolation)")
    ap.add_argument("--trunc", default="none",
                    help="spherical harmonics: triangular truncation of the derivatives. none (default): "
                         "everything the grid allows; auto: quadratic truncation (2*nlat-1)//3; or an integer")
    ap.add_argument("--zonal-filter", type=float, default=None, metavar="C",
                    help="remove zonal wavenumbers m > C*(nlon/2)*cos(lat) from ua, va, ps before the "
                         "derivatives (for data regridded from a reduced Gaussian grid, e.g. CNRM 'gr')")
    ap.add_argument("--scan", action="store_true",
                    help="with --cfday: compare zonal-filter factors 1.3/1.15/1.0/0.9, no filter and centred "
                         "differences in one pass (+ _scan.csv, _zonalspectra.png)")
    ap.add_argument("--unregrid", default=None, metavar="JSON",
                    help="recover the native reduced-grid values (table from --regrid-check) and re-interpolate "
                         "them spectrally before computing omega")
    ap.add_argument("--daily", action="store_true",
                    help="compare daily means on pressure levels with day/Eday wap (models without 6hrPlev/CFday wap)")
    ap.add_argument("--no-reregrid", action="store_true",
                    help="with --unregrid: compare the computed omega WITHOUT applying CNRM's output regridding to it")
    ap.add_argument("--regrid-check", action="store_true",
                    help="inspect how the 6hrLev output was regridded (1 minute), then exit")
    ap.add_argument("--compare", action="store_true",
                    help="compute mean/cd, sb/cd, mean/sh, sb/sh and compare them")
    a = ap.parse_args(argv)
    global DONOR_FOR, REF_SHIFT_DAYS
    DONOR_FOR = a.model
    REF_SHIFT_DAYS = a.ref_shift
    if REF_SHIFT_DAYS:
        print(f"  note: reference date labels shifted by {REF_SHIFT_DAYS:+d} day(s) (--ref-shift)")
    levels = [x * 100 for x in a.levels] if a.levels else None
    global SH_TRUNC, ZONAL_FILTER
    SH_TRUNC = a.trunc if a.trunc in ("auto", "none") else int(a.trunc)
    ZONAL_FILTER = a.zonal_filter
    global UNREGRID, REREGRID
    REREGRID = not a.no_reregrid
    if a.unregrid:
        UNREGRID = load_unregrid(a.unregrid)
        print(f"unregrid: {len(UNREGRID)} latitude rows reconstructed from their native points ({a.unregrid})")
    if a.regrid_check:
        regrid_check(a.model, a.exp, a.member, a.outdir)
        return
    if a.daily:
        run_daily(a.model, a.start, a.end, a.exp, a.member, a.outdir, not a.no_figure,
                  "fft" if a.fft else a.deriv, a.box, a.fullevel, a.snapshot)
        return
    if a.cfday:
        run_cfday(a.model, a.start, a.end, a.exp, a.member, a.outdir, not a.no_figure,
                  "fft" if a.fft else a.deriv, a.box, a.fullevel, a.snapshot, a.scan)
        return
    run(a.model, a.start, a.end, levels, a.exp, a.member, a.outdir, not a.no_figure,
        "fft" if a.fft else a.deriv, a.box, a.fullevel, a.compare, a.snapshot)


if __name__ == "__main__":
    main()
