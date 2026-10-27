#!/usr/bin/env python3
"""LAGRANTO P- and S-files from CMIP6 6-hourly model-level output (6hrLev) on spiritx.

Replaces stages 01-03 of the ERA5 pipeline (case_studies/June_2026) for one CMIP6
model / experiment / member. Stages 04-06 are then the unchanged ERA5 ones.

  P-files  T, U, V, OMEGA [Pa s-1], PS [hPa], hyam/hybm/hyai/hybi
           OMEGA is diagnosed from 6hrLev ua, va, ps with the continuity equation on the
           model's native grid, with the validated functions of tools/omega/omega_check.py
           (spherical-harmonic derivatives on Gaussian grids, Simmons & Burridge full levels).
  S-files  THETA, TCLIM, DTDT, PS, hyam/hybm
           TCLIM = the model's own climatology by the paper rule: 21-day running-mean T at the
           same 6-hourly slot, averaged over a centred CLIM_NYEARS window around the event year,
           shifted to stay inside CLIM_FIRST_YEAR..CLIM_LAST_YEAR (historical: 1950-2014).
           DTDT = centred difference of TCLIM over +/-6 h (as stage 03).

Everything is computed on the native grid and then interpolated linearly to a regular GRID_DEG
grid, LAT_MIN..90N, lon 0..360-GRID_DEG (as the ERA5 P/S files). Levels top -> surface,
lev = 1..N, with the coefficients of exactly these N levels (LAGRANTO takes the level pressures
from the bottom-most N coefficients). Written with lagranto_io.write_lagranto_netcdf.

Not yet supported (refused with a message):
  * models whose 6hrLev output is regridded from a reduced Gaussian grid (CNRM family,
    EC-Earth3): they need the unregrid step of tools/omega/omega_check_cnrm.py first;
  * non-standard calendars (noleap/360_day): time stamps must decode to datetime64.

Usage (config.sh sourced, see submit_cmip6.sh):
  python prepare_cmip6.py            # all P/S files of the run
  python prepare_cmip6.py --test     # first time step only, with plausibility checks
  python prepare_cmip6.py --fix-attrs  # only set the variable attributes in existing P/S files
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "tools", "omega"))
sys.path.insert(0, HERE)
import omega_check as oc  # noqa: E402
from lagranto_io import write_lagranto_netcdf  # noqa: E402

ENV = os.environ
MODEL = ENV["CMIP6_MODEL"]
EXP = ENV.get("CMIP6_EXP", "historical")
MEMBER = ENV.get("CMIP6_MEMBER", "")
DATA_DIR = ENV["DATA_DIR"]
TRAJ_START = pd.Timestamp(ENV["TRAJ_START"])
TRAJ_END = pd.Timestamp(ENV["TRAJ_END"]) + pd.Timedelta(hours=18)
BACK_DAYS = int(ENV.get("BACK_DAYS", "15"))
GRID_DEG = float(ENV.get("GRID_DEG", "1.0"))
LAT_MIN = float(ENV.get("LAT_MIN", "0"))
CLIM_NYEARS = int(ENV.get("CLIM_NYEARS", "9"))
CLIM_FIRST_YEAR = int(ENV.get("CLIM_FIRST_YEAR", "1950"))
CLIM_LAST_YEAR = int(ENV.get("CLIM_LAST_YEAR", "2014"))
HALF_WIN = 10                            # days: 21-day running mean
DT = pd.Timedelta(hours=6)
KAPPA = 287.0 / 1004.0                   # as stage 03
REGRIDDED = ("CNRM-CM6-1", "CNRM-ESM2-1", "CNRM-CM6-1-HR", "EC-Earth3", "EC-Earth3-Veg",
             "EC-Earth3-Veg-LR", "EC-Earth3-AerChem", "EC-Earth3-CC")
LEV_ATTRS = {"standard_name": "hybrid_sigma_pressure", "long_name": "hybrid level at layer midpoints",
             "formula": "hyam hybm (mlev=hyam+hybm*aps)", "formula_terms": "ap: hyam b: hybm ps: aps",
             "units": "level", "positive": "down"}
LAT_ATTRS = {"standard_name": "latitude", "long_name": "latitude", "units": "degrees_north", "axis": "Y"}
LON_ATTRS = {"standard_name": "longitude", "long_name": "longitude", "units": "degrees_east", "axis": "X"}
# variable attributes as in the ERA5 P/S files. LAGRANTO reads the units of the hybrid coefficients:
# without them ("NetCDF: Attribute not found") it puts every level at the wrong pressure.
VAR_ATTRS = {
    "hyai": {"long_name": "hybrid A coefficient at layer interfaces", "units": "Pa"},
    "hybi": {"long_name": "hybrid B coefficient at layer interfaces", "units": "1"},
    "hyam": {"long_name": "hybrid A coefficient at layer midpoints", "units": "Pa"},
    "hybm": {"long_name": "hybrid B coefficient at layer midpoints", "units": "1"},
    "T": {"standard_name": "air_temperature", "long_name": "Temperature", "units": "K"},
    "U": {"standard_name": "eastward_wind", "long_name": "U component of wind", "units": "m s**-1"},
    "V": {"standard_name": "northward_wind", "long_name": "V component of wind", "units": "m s**-1"},
    "OMEGA": {"standard_name": "lagrangian_tendency_of_air_pressure", "long_name": "Vertical velocity",
              "units": "Pa s**-1"},
    "PS": {"units": "hPa", "long_name": "Surface Pressure"},
    "THETA": {"long_name": "Potential Temperature", "units": "K"},
    "TCLIM": {"standard_name": "air_temperature", "long_name": "Temperature", "units": "K"},
    "DTDT": {"units": "K s-1", "long_name": "Time derivative of T_clim"},
}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ------------------------------------------------------------------ input
def covering(var, start, end):
    """{member: version dir} of the 6hrLev files of `var` overlapping [start, end]."""
    return oc._covering(oc.version_dirs(MODEL, "6hrLev", var, EXP),
                        start.strftime("%Y-%m-%d"), end.strftime("%Y-%m-%d"))


def open_var(path, start, end):
    ds = oc.open_files(path, start.strftime("%Y-%m-%d"), (end + pd.Timedelta(days=1)).strftime("%Y-%m-%d"))
    if not np.issubdtype(ds["time"].dtype, np.datetime64):
        raise SystemExit(f"{MODEL}: calendar '{ds['time'].encoding.get('calendar')}' not supported yet")
    return ds.sel(time=slice(start, end))


def on_levels(da, vdim, flip):
    """top -> surface, vertical dim renamed to 'lev'"""
    if flip:
        da = da.isel({vdim: slice(None, None, -1)})
    return da.rename({vdim: "lev"}) if vdim != "lev" else da


# ------------------------------------------------------------------ grid
LAT_T = np.round(np.arange(LAT_MIN, 90.0 + 1e-9, GRID_DEG), 6)
LON_T = np.round(np.arange(0.0, 360.0 - 1e-9, GRID_DEG), 6)


def to_regular(da, name, units, lev=None):
    """Native (lev, lat, lon) -> regular grid, linear, lon periodic, extrapolated at the pole rows."""
    da = da.sortby("lat")
    lon = da["lon"].values
    da = xr.concat([da.isel(lon=[-1]).assign_coords(lon=[lon[-1] - 360.0]), da,
                    da.isel(lon=[0]).assign_coords(lon=[lon[0] + 360.0])], dim="lon")
    out = da.interp(lat=LAT_T, lon=LON_T, method="linear", kwargs={"fill_value": None})
    out = out.drop_vars([c for c in out.coords if c not in ("lev", "lat", "lon")])
    if "lev" in out.dims:
        out = out.transpose("lev", "lat", "lon").assign_coords(lev=lev)
    out = out.astype("float32")
    out.attrs = {"units": units}
    out.name = name
    return out


# ------------------------------------------------------------------ climatology
def clim_years(year):
    """Centred CLIM_NYEARS window, shifted (not shortened) at the ends of the run."""
    y0 = min(max(year - CLIM_NYEARS // 2, CLIM_FIRST_YEAR), CLIM_LAST_YEAR - CLIM_NYEARS + 1)
    return list(range(y0, y0 + CLIM_NYEARS))


def climatology(targets, member, vdim, flip):
    """TCLIM on the native grid at the target times: {time: DataArray (lev, lat, lon)}."""
    years = clim_years(targets[len(targets) // 2].year)
    log(f"climatology {years[0]}-{years[-1]} ({len(years)} years), 21-day running mean per 6-h slot")
    acc = {}
    cnt = {t: 0 for t in targets}
    for y in years:
        sh = lambda t: t.replace(year=y)
        w0 = sh(targets[0]).normalize() - pd.Timedelta(days=HALF_WIN)
        w1 = sh(targets[-1]).normalize() + pd.Timedelta(days=HALF_WIN, hours=18)
        dirs = covering("ta", w0, w1)
        if member not in dirs:
            raise SystemExit(f"{MODEL} {EXP} {member}: no 6hrLev ta covering {w0:%Y-%m-%d}..{w1:%Y-%m-%d}")
        ta = on_levels(open_var(dirs[member], w0, w1)["ta"], vdim, flip)
        for hh in (0, 6, 12, 18):
            mine = [t for t in targets if t.hour == hh]
            if not mine:
                continue
            x = ta.isel(time=(ta["time"].dt.hour == hh).values).load().astype("float32")
            days = pd.DatetimeIndex(x["time"].values).normalize()
            if len(days) != len(pd.date_range(w0, w1.normalize(), freq="D")):
                raise SystemExit(f"{MODEL} ta {y} {hh:02d} UTC: {len(days)} days in "
                                 f"{w0:%Y-%m-%d}..{w1:%Y-%m-%d} (missing or repeated time steps)")
            r = x.rolling(time=2 * HALF_WIN + 1, center=True, min_periods=2 * HALF_WIN + 1).mean()
            for t in mine:
                k = np.flatnonzero(days == sh(t).normalize())
                v = r.isel(time=int(k[0])).drop_vars("time")
                if bool(np.isnan(v.values).any()):
                    raise SystemExit(f"{MODEL} ta {y}: running mean undefined at {sh(t)}")
                acc[t] = v if t not in acc else acc[t] + v
                cnt[t] += 1
        del ta
        log(f"  {y} done")
    bad = [t for t in targets if cnt[t] != len(years)]
    if bad:
        raise SystemExit(f"climatology incomplete for {len(bad)} times, e.g. {bad[:3]}")
    return {t: acc[t] / len(years) for t in targets}


# ------------------------------------------------------------------ attributes of existing files
def fix_attrs():
    """Set VAR_ATTRS in place in the P/S files of DATA_DIR (files written before the attributes)."""
    import glob
    import netCDF4
    files = sorted(glob.glob(os.path.join(DATA_DIR, "[PS][0-9]*_[0-9][0-9]")))
    if not files:
        raise SystemExit(f"no P/S files in {DATA_DIR}")
    for f in files:
        with netCDF4.Dataset(f, "a") as nc:
            for v, att in VAR_ATTRS.items():
                if v in nc.variables:
                    nc.variables[v].setncatts(att)
    log(f"attributes set in {len(files)} files in {DATA_DIR}")


# ------------------------------------------------------------------ main
def main():
    if "--fix-attrs" in sys.argv[1:]:
        return fix_attrs()
    test = "--test" in sys.argv[1:]
    if MODEL in REGRIDDED:
        raise SystemExit(f"{MODEL}: 6hrLev output is regridded from a reduced Gaussian grid; the unregrid "
                         "step (tools/omega/omega_check_cnrm.py --unregrid) is not wired in here yet")
    oc.DONOR_FOR = MODEL
    oc.LON_DERIV = "cd"
    os.makedirs(DATA_DIR, exist_ok=True)
    times = pd.date_range(TRAJ_START - pd.Timedelta(days=BACK_DAYS), TRAJ_END, freq=DT)
    if test:
        times = times[:1]
    targets = list(pd.date_range(times[0] - DT, times[-1] + DT, freq=DT))   # TCLIM incl. +/-6 h
    log(f"{MODEL} {EXP}: P/S files {times[0]} .. {times[-1]} ({len(times)} steps), "
        f"{GRID_DEG} deg grid {LAT_MIN}-90N -> {DATA_DIR}")

    # one member that has ua, va, ps and ta for the whole period
    dirs = {v: covering(v, times[0], times[-1]) for v in ("ua", "va", "ps", "ta")}
    common = sorted(set.intersection(*(set(d) for d in dirs.values())))
    member = MEMBER or (common[0] if common else "")
    if member not in common:
        raise SystemExit(f"{MODEL} {EXP}: member '{member}' lacks 6hrLev ua/va/ps/ta for "
                         f"{times[0]:%Y-%m-%d}..{times[-1]:%Y-%m-%d} under {oc.POOL} (covering: {common})")
    log(f"member {member}: " + ", ".join(f"{v} {dirs[v][member]}" for v in dirs))
    ds = {v: open_var(dirs[v][member], times[0], times[-1]) for v in dirs}
    uvps = xr.merge([ds["ua"][["ua"]], ds["va"][["va"]], ds["ps"][["ps"]]], compat="override", join="inner")

    A, B, vdim, flip = oc.half_level_coeffs(ds["ua"], "ua")
    At, Bt, vdim_t, flip_t = oc.half_level_coeffs(ds["ta"], "ta")
    if not (len(At) == len(A) and np.allclose(A, At) and np.allclose(B, Bt)):
        raise SystemExit("ua and ta have different vertical coefficients")
    n = len(A) - 1
    deriv = ENV.get("OMEGA_DERIV", "auto")
    if deriv == "auto":
        try:
            oc._sh_for(np.sort(uvps["lat"].values), uvps["lon"].values)
            deriv = "sh"
        except ValueError:
            deriv = "cd"
    pm = 0.5 * (A[1:] + A[:-1]) + 0.5 * (B[1:] + B[:-1]) * 1e5
    log(f"{n} levels (top {pm[0] / 100:.3g} hPa, lowest {pm[-1] / 100:.1f} hPa at ps 1000 hPa); "
        f"native grid {uvps.sizes['lat']} x {uvps.sizes['lon']}; omega derivatives '{deriv}'")
    if n > 120:
        log("WARNING: more than 120 levels - the LAGRANTO trace on spiritx crashed above ~120 levels")

    clim = climatology(targets, member, vdim_t, flip_t)
    lev = np.arange(1, n + 1, dtype="float64")
    hyam, hybm = 0.5 * (A[1:] + A[:-1]), 0.5 * (B[1:] + B[:-1])
    coords = {"lev": ("lev", lev, LEV_ATTRS)}

    def finish(dsx, t):
        dsx = dsx.assign_coords(coords)
        dsx["lat"].attrs, dsx["lon"].attrs = LAT_ATTRS, LON_ATTRS
        for v in dsx.data_vars:
            dsx[v].attrs = dict(VAR_ATTRS[v])
        return dsx.expand_dims(time=[np.datetime64(t, "ns")])

    for t in times:
        tag = t.strftime("%Y%m%d_%H")
        terms = oc.flux_terms(uvps.sel(time=t).load(), A, B, vdim, flip, deriv)
        om, _ = oc.omega_from_terms(terms, A, B, vdim, "sb", deriv)
        u, v, ps = terms[0], terms[1], terms[2]
        ta = on_levels(ds["ta"]["ta"].sel(time=t), vdim_t, flip_t).load()
        rn = lambda x: x.rename({vdim: "lev"}) if vdim != "lev" else x
        T = to_regular(ta, "T", "K", lev)
        U = to_regular(rn(u), "U", "m s-1", lev)
        V = to_regular(rn(v), "V", "m s-1", lev)
        W = to_regular(rn(om), "OMEGA", "Pa s-1", lev)
        PS = (to_regular(ps, "PS", "hPa") / 100.0).astype("float32")
        PS.attrs = {"units": "hPa", "long_name": "Surface Pressure"}
        P = xr.Dataset({"T": T, "U": U, "V": V, "OMEGA": W, "PS": PS,
                        "hyai": ("nhyi", A), "hybi": ("nhyi", B),
                        "hyam": ("nhym", hyam), "hybm": ("nhym", hybm)})
        write_lagranto_netcdf(finish(P, t), os.path.join(DATA_DIR, f"P{tag}"))

        TC = to_regular(clim[t], "TCLIM", "K", lev)
        dTdt = ((to_regular(clim[t + DT], "x", "K", lev) - to_regular(clim[t - DT], "x", "K", lev))
                / (2 * 6 * 3600.0)).astype("float32")
        dTdt.attrs = {"units": "K s-1", "long_name": "Time derivative of T_clim"}
        pfull = xr.DataArray(hyam, dims="lev", coords={"lev": lev}) \
            + xr.DataArray(hybm, dims="lev", coords={"lev": lev}) * PS * 100.0
        TH = (T * (1e5 / pfull) ** KAPPA).astype("float32").transpose("lev", "lat", "lon")
        TH.attrs = {"units": "K", "long_name": "Potential Temperature"}
        S = xr.Dataset({"THETA": TH, "TCLIM": TC, "DTDT": dTdt, "PS": PS,
                        "hyam": ("nhym", hyam), "hybm": ("nhym", hybm)})
        write_lagranto_netcdf(finish(S, t), os.path.join(DATA_DIR, f"S{tag}"))
        log(f"P/S {tag}")

        if test:
            k5 = int(np.argmin(np.abs(pm - 5e4)))
            kb = n - 1
            print(f"  checks at {tag} (regular grid):\n"
                  f"    PS {float(PS.min()):.0f}..{float(PS.max()):.0f} hPa\n"
                  f"    T  lowest level {float(T[kb].min()):.1f}..{float(T[kb].max()):.1f} K, "
                  f"~500 hPa (lev {k5 + 1}) {float(T[k5].mean()):.1f} K mean\n"
                  f"    OMEGA ~500 hPa: mean {float(W[k5].mean()):+.4f}, std {float(W[k5].std()):.3f}, "
                  f"|max| {float(abs(W[k5]).max()):.2f} Pa/s (ERA5 1 deg typically std 0.1-0.2)\n"
                  f"    T - TCLIM lowest level: mean {float((T[kb] - TC[kb]).mean()):+.2f}, "
                  f"std {float((T[kb] - TC[kb]).std()):.2f} K\n"
                  f"    DTDT ~500 hPa |max| {float(abs(dTdt[k5]).max()) * 86400:.3f} K/day\n"
                  f"    THETA lowest level {float(TH[kb].min()):.1f}..{float(TH[kb].max()):.1f} K",
                  flush=True)
    log("done")


if __name__ == "__main__":
    main()
