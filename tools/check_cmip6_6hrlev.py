#!/usr/bin/env python3
"""
Read-only integrity check of CMIP6 6hrLev model-level data on the ESPRI
mirror (/bdd/CMIP6), for every model/member/grid that has hus, ta, ua and va.

Per model/member it reports resolution, vertical grid, calendar and time
coverage, and flags anything that would break a trajectory run:
  files      unreadable, truncated (last timestep unreadable), dates in the
             filename != time axis, time dim != variable length
  time       not strictly 6-hourly: gaps, duplicates, overlapping files;
             the 4 variables do not share one identical time axis
  vertical   lev / hybrid coefficients differ between files or variables;
             a formula_terms variable (ps, ap or a/p0, b, ...) is missing
  horizontal lat/lon differ between files or variables
  values     fill/NaN or out-of-range values in sampled timesteps (first,
             middle and last of every file), for the variable and for ps

Default: every activity and experiment, every member fully checked.
  --quick  only lists directories and parses filenames (no file opened)
  --best   per model/experiment fully check the longest member only, and fall
           back to the next longest if it fails (up to --max-tries)

Also reports the atmosphere model (global 'source' attribute), spectral vs
grid-point dynamical core (DYCORE table below), native vs regridded
horizontal grid (grid_label + 'grid' attribute) and the number of model
levels between 200 and 1000 hPa for a 1013.25 hPa surface pressure.

Writes nothing to the data. Output in --out: ensembles_*.tsv (one row per
model/experiment: longest clean member + model/grid info), matrix_*.tsv
(model x experiment), summary_*.tsv (one row per checked member) and
issues_*.tsv.

  python check_cmip6_6hrlev.py [--quick|--best] [--activity CMIP] [--exp historical]
                               [--njobs 32] [--out DIR] [--only M1,M2]
"""
import argparse
import glob
import hashlib
import os
import re
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

import cftime
import netCDF4
import numpy as np

ROOT = "/bdd/CMIP6"
VARS = ("hus", "ta", "ua", "va")
STEP = 0.25                                    # 6 h in days
REF = "days since 1850-01-01 00:00:00"
CAL_ALIAS = {"gregorian": "standard", "365_day": "noleap", "366_day": "all_leap"}
RANGES = {"ta": (150., 350.), "hus": (-1e-4, 0.06), "ua": (-250., 250.),
          "va": (-250., 250.), "ps": (3.0e4, 1.15e5)}
UNITS = {"ta": {"K"}, "hus": {"1", "kg kg-1", "kg/kg"}, "ua": {"m s-1", "m/s"},
         "va": {"m s-1", "m/s"}, "ps": {"Pa"}}

# Dynamical core of the atmosphere model. Not in the CMIP6 metadata, so taken
# from the model documentation; "?" = not in this list, look it up.
#   spectral   = spectral transform model, fields saved on its Gaussian grid
#   grid point = finite-volume / finite-difference / semi-Lagrangian grid model
_SP, _GP = "spectral", "grid point"
DYCORE = {
    "ACCESS-CM2": (_GP, "UM, N96 lat-lon"), "ACCESS-ESM1-5": (_GP, "UM HadGAM2, N96 lat-lon"),
    "AWI-CM-1-1-MR": (_SP, "ECHAM6.3 T127"), "AWI-ESM-1-1-LR": (_SP, "ECHAM6.3 T63"),
    "BCC-CSM2-MR": (_SP, "BCC-AGCM3 T106"), "BCC-ESM1": (_SP, "BCC-AGCM3 T42"),
    "CanESM5": (_SP, "CanAM5 T63"), "CanESM5-CanOE": (_SP, "CanAM5 T63"),
    "CESM2": (_GP, "CAM6 finite volume 0.9x1.25"), "CESM2-WACCM": (_GP, "WACCM6 finite volume 0.9x1.25"),
    "CESM2-FV2": (_GP, "CAM6 finite volume 1.9x2.5"), "CESM2-WACCM-FV2": (_GP, "WACCM6 finite volume 1.9x2.5"),
    "CMCC-CM2-SR5": (_GP, "CAM5 finite volume"), "CMCC-ESM2": (_GP, "CAM5 finite volume"),
    "CMCC-CM2-HR4": (_GP, "CAM4 finite volume"),
    "CNRM-CM6-1": (_SP, "ARPEGE-Climat 6.3 T127, reduced Gaussian"),
    "CNRM-CM6-1-HR": (_SP, "ARPEGE-Climat 6.3 T359, reduced Gaussian"),
    "CNRM-ESM2-1": (_SP, "ARPEGE-Climat 6.3 T127, reduced Gaussian"),
    "E3SM-1-0": (_GP, "EAM spectral element (ne30)"), "E3SM-1-1": (_GP, "EAM spectral element (ne30)"),
    "EC-Earth3": (_SP, "IFS cy36r4 T255, reduced Gaussian"), "EC-Earth3-Veg": (_SP, "IFS cy36r4 T255, reduced Gaussian"),
    "EC-Earth3-AerChem": (_SP, "IFS cy36r4 T255, reduced Gaussian"), "EC-Earth3-CC": (_SP, "IFS cy36r4 T255, reduced Gaussian"),
    "EC-Earth3-Veg-LR": (_SP, "IFS cy36r4 T159, reduced Gaussian"),
    "FGOALS-g3": (_GP, "GAMIL3 finite difference"), "FGOALS-f3-L": (_GP, "FAMIL2 cubed sphere"),
    "GFDL-CM4": (_GP, "AM4 cubed sphere C96"), "GFDL-ESM4": (_GP, "AM4.1 cubed sphere C96"),
    "GISS-E2-1-G": (_GP, "GISS ModelE2.1, 2x2.5 lat-lon"), "GISS-E2-1-H": (_GP, "GISS ModelE2.1, 2x2.5 lat-lon"),
    "GISS-E2-2-G": (_GP, "GISS ModelE2.2, 2x2.5 lat-lon"), "GISS-E3-G": (_GP, "GISS ModelE3 cubed sphere C90"),
    "HadGEM3-GC31-LL": (_GP, "UM ENDGame N96"), "HadGEM3-GC31-MM": (_GP, "UM ENDGame N216"),
    "UKESM1-0-LL": (_GP, "UM ENDGame N96"), "KACE-1-0-G": (_GP, "UM N96"),
    "INM-CM4-8": (_GP, "INM finite difference 2x1.5"), "INM-CM5-0": (_GP, "INM finite difference 2x1.5"),
    "IPSL-CM6A-LR": (_GP, "LMDZ6 finite difference lat-lon 144x143"),
    "IPSL-CM6A-LR-INCA": (_GP, "LMDZ6 finite difference lat-lon 144x143"),
    "IPSL-CM5A2-INCA": (_GP, "LMDZ5 finite difference lat-lon 96x95"),
    "MIROC6": (_SP, "CCSR-AGCM T85"), "MIROC-ES2L": (_SP, "CCSR-AGCM T42"), "MIROC-ES2H": (_SP, "CCSR-AGCM T85"),
    "MPI-ESM1-2-LR": (_SP, "ECHAM6.3 T63"), "MPI-ESM1-2-HR": (_SP, "ECHAM6.3 T127"),
    "MPI-ESM-1-2-HAM": (_SP, "ECHAM6.3-HAM T63"), "NESM3": (_SP, "ECHAM6.3 T63"),
    "MRI-ESM2-0": (_SP, "MRI-AGCM3.5 TL159"),
    "NorESM2-LM": (_GP, "CAM-OSLO finite volume 1.9x2.5"), "NorESM2-MM": (_GP, "CAM-OSLO finite volume 0.9x1.25"),
    "SAM0-UNICON": (_GP, "CAM5 finite volume"), "TaiESM1": (_GP, "CAM5 finite volume"),
}


def hsh(*arrs):
    m = hashlib.md5()
    for a in arrs:
        m.update(np.round(np.ma.getdata(np.asarray(a, "f8")), 6).tobytes())
    return m.hexdigest()[:10]


def natural(s):
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", s)]


def discover(activity, exp, only):
    """key = (model, experiment, member, grid, institution, activity)"""
    pat = os.path.join(ROOT, activity, "*", "*", exp, "*", "6hr[lL]ev")
    for d in sorted(glob.glob(pat)):
        act, inst, model, ex, member = d.split(os.sep)[-6:-1]
        if only and model not in only:
            continue
        grids = None
        for v in VARS:
            g = {os.path.basename(p) for p in glob.glob(os.path.join(d, v, "*")) if os.path.isdir(p)}
            grids = g if grids is None else grids & g
        for g in sorted(grids or ()):
            yield (model, ex, member, g, inst, act), d


def quick_var(key, vdir, var):
    """Filenames only: number of files and the period they claim to cover."""
    d, ver = version_dir(os.path.join(vdir, var, key[3]))
    files = sorted(glob.glob(os.path.join(d, "*.nc"))) if d else []
    spans = [m.groups() for f in files if (m := re.search(r"_(\d{8,12})-(\d{8,12})\.nc$", f))]
    return key, dict(var=var, files=len(files), version=ver,
                     start=min((s for s, _ in spans), default=None),
                     end=max((e for _, e in spans), default=None))


def version_dir(vg):
    latest = os.path.join(vg, "latest")
    if os.path.isdir(latest):
        return latest, os.path.basename(os.path.realpath(latest))
    vs = sorted(glob.glob(os.path.join(vg, "v*")))
    return (vs[-1], os.path.basename(vs[-1])) if vs else (None, None)


def pressure_profile(nc, terms, ps0=101325.):
    """Pressure (hPa) of every level for a reference surface pressure."""
    g = lambda k: np.asarray(nc.variables[terms[k]][...], "f8")
    try:
        if "ap" in terms and "b" in terms:
            p = g("ap") + g("b") * ps0
        elif {"a", "b", "p0"} <= terms.keys():
            p = g("a") * float(g("p0")) + g("b") * ps0
        elif {"sigma", "ptop"} <= terms.keys():
            pt = float(g("ptop"))
            p = pt + g("sigma") * (ps0 - pt)
        elif {"a", "orog"} <= terms.keys():
            # hybrid height (UM): height over sea level (orog = 0) -> US standard atmosphere
            z = g("a")
            p = np.where(z <= 11000., ps0 * (1. - 2.25577e-5 * z) ** 5.25588,
                         22632. * np.exp(-(z - 11000.) / 6341.6))
        else:
            return None
    except Exception:
        return None
    return p / 100.


def atmos_component(source):
    """'atmos: LMDZ (NPv6, N96; 144 x 143 ...)' from the global 'source' attribute."""
    m = re.search(r"(?:^|\n)\s*atmos:\s*([^\n]+)", source)
    return " ".join(m[1].split()) if m else "?"


def stamp(d, n):
    return d.strftime("%Y%m%d%H%M")[:n]


def scan_var(key, vdir, var):
    out = dict(var=var, issues=[], files=0, version=None, info={}, cal=set(),
               levh=set(), coefh=set(), gridh=set(), t=None)
    issue = lambda f, msg: out["issues"].append((var, os.path.basename(f) if f else "-", msg))
    d, out["version"] = version_dir(os.path.join(vdir, var, key[3]))
    files = sorted(glob.glob(os.path.join(d, "*.nc"))) if d else []
    out["files"] = len(files)
    if not files:
        issue(None, "no .nc files")
        return key, out
    stats = {n: [np.inf, -np.inf, 0, 0] for n in (var, "ps")}   # min max nmasked nnan
    tchunks = []
    for f in files:
        try:
            with netCDF4.Dataset(f) as nc:
                v = nc.variables[var]
                tv = nc.variables["time"]
                cal = CAL_ALIAS.get(getattr(tv, "calendar", "standard").lower(),
                                    getattr(tv, "calendar", "standard").lower())
                dates = cftime.num2date(np.asarray(tv[:], "f8"), tv.units, cal)
                days = np.asarray(cftime.date2num(dates, REF, cal), "f8")
                nt = len(days)
                out["cal"].add(cal)
                if v.shape[0] != nt:
                    issue(f, f"{var} has {v.shape[0]} steps, time axis {nt}")
                m = re.search(r"_(\d{8,12})-(\d{8,12})\.nc$", f)
                if m and nt and (stamp(dates[0], len(m[1])) != m[1] or stamp(dates[-1], len(m[2])) != m[2]):
                    issue(f, f"filename {m[1]}-{m[2]} but time axis {dates[0]} .. {dates[-1]}")
                if getattr(v, "units", "?") not in UNITS[var]:
                    issue(f, f"units '{getattr(v, 'units', '?')}'")

                # horizontal grid
                dims = v.dimensions
                hdims = dims[2:]
                hvars = [nc.variables[n][:] for n in ("lat", "lon") if n in nc.variables]
                out["gridh"].add(hsh(*hvars))
                # vertical grid
                levvar = nc.variables.get(dims[1])
                terms = dict(re.findall(r"(\w+):\s*(\w+)", getattr(levvar, "formula_terms", ""))) if levvar is not None else {}
                if levvar is None:
                    issue(f, f"no coordinate variable for vertical dim '{dims[1]}'")
                else:
                    out["levh"].add(hsh(levvar[:]))
                if not terms:
                    issue(f, "lev has no formula_terms (pressure not reconstructable)")
                coefs = []
                for k, n in sorted(terms.items()):
                    if n not in nc.variables:
                        issue(f, f"formula term {k}: '{n}' missing from file")
                    elif k != "ps" and n != dims[1]:
                        coefs.append(nc.variables[n][...])
                out["coefh"].add(hsh(*coefs))
                psname = terms.get("ps")
                ps = nc.variables.get(psname) if psname else None
                if ps is not None and ps.shape[0] != nt:
                    issue(f, f"ps has {ps.shape[0]} steps, time axis {nt}")

                if not out["info"]:
                    lat = np.asarray(nc.variables["lat"][:]) if "lat" in nc.variables else None
                    lon = np.asarray(nc.variables["lon"][:]) if "lon" in nc.variables else None
                    p = pressure_profile(nc, terms)
                    out["info"] = dict(
                        shape="x".join(str(nc.dimensions[h].size) for h in hdims),
                        dlat=float(np.median(np.abs(np.diff(lat)))) if lat is not None and lat.ndim == 1 else np.nan,
                        dlon=float(np.median(np.abs(np.diff(lon)))) if lon is not None and lon.ndim == 1 else np.nan,
                        nlev=nc.dimensions[dims[1]].size,
                        vcoord=getattr(levvar, "standard_name", "?") if levvar is not None else "?",
                        ptop=float(p.min()) if p is not None else np.nan,
                        pbot=float(p.max()) if p is not None else np.nan,
                        order=("top->bottom" if p[0] < p[-1] else "bottom->top") if p is not None else "?",
                        n200_1000=int(((p >= 200.) & (p <= 1000.)).sum()) if p is not None else "?",
                        nominal=getattr(nc, "nominal_resolution", "?"),
                        atmos=atmos_component(getattr(nc, "source", "")),
                        grid_label=getattr(nc, "grid_label", key[3]),
                        grid_desc=" ".join(str(getattr(nc, "grid", "?")).split()),
                        hour0=f"{dates[0].hour:02d}" if nt else "?")

                # sampled values: first/middle/last step (also catches truncated files)
                for i in sorted({0, nt // 2, nt - 1} - {-1}):
                    for n, vv in ((var, v), ("ps", ps)):
                        if vv is None:
                            continue
                        try:
                            a = vv[i]
                        except Exception as e:
                            issue(f, f"{n} unreadable at {dates[i]}: {type(e).__name__}: {e}")
                            continue
                        a = np.ma.masked_invalid(a) if not np.ma.isMaskedArray(a) else a
                        s = stats[n]
                        s[2] += int(np.ma.count_masked(a))
                        s[3] += int(np.isnan(np.ma.getdata(a)).sum())
                        if a.count():
                            s[0], s[1] = min(s[0], float(a.min())), max(s[1], float(a.max()))
                tchunks.append(days)
        except Exception as e:
            issue(f, f"cannot open/read: {type(e).__name__}: {e}")

    for name, label in (("cal", "calendars"), ("levh", "lev values"),
                        ("coefh", "hybrid coefficients"), ("gridh", "lat/lon")):
        if len(out[name]) > 1:
            issue(None, f"{label} differ between files ({len(out[name])} variants)")
    for n, (lo, hi, nmask, nnan) in stats.items():
        if n == "ps" and lo == np.inf:
            continue
        if nmask or nnan:
            issue(None, f"{n}: {nmask} fill/masked and {nnan} NaN values in sampled steps")
        if lo < RANGES[n][0] or hi > RANGES[n][1]:
            issue(None, f"{n}: sampled range {lo:.4g} .. {hi:.4g} outside {RANGES[n]}")
    out["range"] = {n: (s[0], s[1]) for n, s in stats.items()}

    if tchunks:
        t = np.concatenate(tchunks)
        out["t"] = t
        cal = next(iter(out["cal"]))
        dt = np.diff(t)
        gaps = np.where(dt > STEP + 1e-4)[0]
        back = np.where(dt < STEP - 1e-4)[0]
        fmt = lambda i: str(cftime.num2date(t[i], REF, cal))[:16]
        if len(gaps):
            miss = int(np.round((dt[gaps] / STEP)).sum() - len(gaps))
            issue(None, f"{len(gaps)} time gaps, {miss} steps missing, first after "
                        + ", ".join(fmt(i) for i in gaps[:5]))
        if len(back):
            issue(None, f"{len(back)} duplicate/overlapping/irregular steps, first at "
                        + ", ".join(fmt(i + 1) for i in back[:5]))
        out["start"], out["end"] = fmt(0), fmt(len(t) - 1)
        out["nsteps"] = len(t)
        out["missing"] = int(np.round((t[-1] - t[0]) / STEP)) + 1 - len(np.unique(t))
        out["years"] = (t[-1] - t[0] + STEP) / {"360_day": 360.}.get(cal, 365.2425)
    return key, out


def regrid_status(model, label):
    if label == "gn":
        return "native"
    if model.startswith("IPSL"):          # LMDZ writes its own lon-lat grid but labels it gr
        return "native (LMDZ lon-lat, labelled gr)"
    return "regridded"


def stamp_years(start, end):
    """Length in years of a filename period YYYYMMDD[hhmm]-YYYYMMDD[hhmm]."""
    f = lambda s: int(s[:4]) + (int(s[4:6]) - 1) / 12. + (int(s[6:8]) - 1) / 365.
    return f(end) - f(start) + 0.25 / 365. if start and end else -1.


def print_table(rows, show):
    w = {c: max([len(c)] + [len(str(r[c])) for r in rows]) for c in show}
    print("\n" + "  ".join(c.ljust(w[c]) for c in show))
    for r in rows:
        print("  ".join(str(r[c]).ljust(w[c]) for c in show))


def write_tsv(path, rows):
    cols = list(rows[0]) if rows else []
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.writelines("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)


def exp_order(e):
    return (e not in ("piControl", "historical"), natural(e))


def matrix_table(rows, out, tag, cell):
    """Model x experiment matrix, one cell(rows) per entry."""
    models = sorted({r["model"] for r in rows}, key=str.lower)
    exps = sorted({r["experiment"] for r in rows}, key=exp_order)
    by = {}
    for r in rows:
        by.setdefault((r["model"], r["experiment"]), []).append(r)
    matrix = [dict(model=m, **{e: cell(by[m, e]) if (m, e) in by else "" for e in exps}) for m in models]
    write_tsv(os.path.join(out, f"matrix_{tag}.tsv"), matrix)
    print_table(matrix, ["model"] + exps)


def run(func, combos, njobs, thread=False, verbose=True):
    res = {k: {} for k in combos}
    with (ThreadPoolExecutor if thread else ProcessPoolExecutor)(njobs) as ex:
        futs = [ex.submit(func, k, d, v) for k, d in combos.items() for v in VARS]
        for n, fu in enumerate(as_completed(futs), 1):
            k, o = fu.result()
            res[k][o["var"]] = o
            if verbose:
                print(f"  [{n}/{len(futs)}] {k[0]} {k[1]} {k[2]} {o['var']}: {o['files']} files, "
                      f"{len(o['issues'])} issues", flush=True)
    return res


def assess(k, o):
    """One summary row + all issues for one model/experiment/member/grid."""
    model, exp, member, grid, inst, act = k
    iss = [i for v in VARS for i in o[v]["issues"]]
    # the 4 variables must share grid, levels, coefficients, calendar and time axis
    for name, label in (("gridh", "lat/lon"), ("levh", "lev values"),
                        ("coefh", "hybrid coefficients"), ("cal", "calendar")):
        if len({frozenset(o[v][name]) for v in VARS}) > 1:
            iss.append(("all", "-", f"{label} differ between variables"))
    ts = [o[v]["t"] for v in VARS]
    if all(t is not None for t in ts) and len({hsh(t) for t in ts}) > 1:
        iss.append(("all", "-", "time axes differ between variables: " + "; ".join(
            f"{v} {o[v]['start']}..{o[v]['end']} ({o[v]['nsteps']})" for v in VARS)))

    r = o["ta"] if o["ta"]["t"] is not None else next((o[v] for v in VARS if o[v]["t"] is not None), o["ta"])
    info = r["info"] or {}
    dyc = DYCORE.get(model, ("?", "?"))
    label = info.get("grid_label", grid)
    row = dict(
        model=model, experiment=exp, member=member, grid=grid, activity=act, institution=inst,
        atmos_model=info.get("atmos", "?"), dycore=dyc[0], dycore_detail=dyc[1],
        nominal=info.get("nominal", "?"), shape=info.get("shape", "?"),
        dlat=f"{info.get('dlat', np.nan):.3g}", dlon=f"{info.get('dlon', np.nan):.3g}",
        regridded=regrid_status(model, label), grid_desc=info.get("grid_desc", "?"),
        nlev=info.get("nlev", "?"), n200_1000=info.get("n200_1000", "?"), vcoord=info.get("vcoord", "?"),
        ptop_hPa=f"{info.get('ptop', np.nan):.3g}", pbot_hPa=f"{info.get('pbot', np.nan):.4g}",
        lev_order=info.get("order", "?"), calendar=",".join(sorted(r["cal"])) or "?",
        first_hour=info.get("hour0", "?"),
        start=r.get("start", "-"), end=r.get("end", "-"),
        years=f"{r.get('years', 0):.1f}", nsteps=r.get("nsteps", 0),
        missing=max((o[v].get("missing", 0) for v in VARS), default="?"),
        files="/".join(str(o[v]["files"]) for v in VARS),
        versions=",".join(sorted({str(o[v]["version"]) for v in VARS})),
        ta_range="{:.4g}..{:.4g}".format(*o["ta"].get("range", {}).get("ta", (np.nan, np.nan))),
        hus_min="{:.3g}".format(o["hus"].get("range", {}).get("hus", (np.nan,))[0]),
        status="OK" if not iss else f"{len(iss)} issues")
    return row, [(model, exp, member, grid, *i) for i in iss]


def ensemble_table(rows, members):
    """Per model/experiment: the longest member that passed every check."""
    by = {}
    for r in rows:
        by.setdefault((r["model"], r["experiment"]), []).append(r)
    out = []
    for (m, e) in sorted(members, key=lambda me: (me[0].lower(), exp_order(me[1]))):
        rs = by.get((m, e), [])
        ok = [r for r in rs if r["status"] == "OK"]
        b = max(ok or rs, key=lambda r: float(r["years"]))
        out.append(dict(
            model=m, experiment=e, n_members=len(members[m, e]), checked=len(rs), clean=len(ok),
            member=b["member"] if ok else "-", start=b["start"][:10] if ok else "-",
            end=b["end"][:10] if ok else "-", years=b["years"] if ok else "0",
            atmos_model=b["atmos_model"], dycore=b["dycore"], dycore_detail=b["dycore_detail"],
            nominal=b["nominal"], shape=b["shape"], dlon=b["dlon"], dlat=b["dlat"],
            grid=b["grid"], regridded=b["regridded"], grid_desc=b["grid_desc"],
            nlev=b["nlev"], n200_1000=b["n200_1000"], ptop_hPa=b["ptop_hPa"], calendar=b["calendar"],
            status="OK" if ok else f"no clean member ({len(rs)} checked)",
            all_members=",".join(sorted(members[m, e], key=natural))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--activity", default="*", help="CMIP, ScenarioMIP, ... (glob, default all)")
    ap.add_argument("--exp", default="*", help="historical, ssp585, ... (glob, default all)")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--quick", action="store_true",
                      help="inventory from directory and file names only (minutes, no file opened)")
    mode.add_argument("--best", action="store_true",
                      help="per model/experiment fully check only the longest member, and the "
                           "next longest if it fails (up to --max-tries)")
    ap.add_argument("--max-tries", type=int, default=3)
    ap.add_argument("--njobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--out", default=f"/homedata/{os.environ.get('USER', 'mmayer')}/cmip6_check")
    ap.add_argument("--only", default="", help="comma-separated model names")
    args = ap.parse_args()
    only = set(filter(None, args.only.split(",")))
    os.makedirs(args.out, exist_ok=True)
    tag = "_".join(s.replace("*", "all") for s in (args.activity, args.exp)) + \
        ("_quick" if args.quick else "_best" if args.best else "")

    combos = dict(discover(args.activity, args.exp, only))
    print(f"{len(combos)} model/experiment/member/grid combos with {'+'.join(VARS)} in "
          f"{ROOT}/{args.activity}/*/*/{args.exp}/*/6hrLev", flush=True)
    members = {}
    for k in combos:
        members.setdefault((k[0], k[1]), set()).add(k[2])
    order = lambda ks: sorted(ks, key=lambda k: (k[0].lower(), exp_order(k[1]), natural(k[2]), k[3]))

    if args.quick or args.best:
        q = run(quick_var, combos, max(args.njobs, 16), thread=True, verbose=False)
        qspan = {}
        for k, o in q.items():
            starts, ends = [o[v]["start"] for v in VARS], [o[v]["end"] for v in VARS]
            s, e = (max(starts), min(ends)) if all(starts) and all(ends) else (None, None)
            qspan[k] = (s, e, stamp_years(s, e), starts, ends)

    if args.quick:
        rows = []
        for k in order(combos):
            model, exp, member, grid, inst, act = k
            o, (s, e, _, starts, ends) = q[k], qspan[k]
            rows.append(dict(
                model=model, experiment=exp, member=member, grid=grid, activity=act, institution=inst,
                start=s or "-", end=e or "-", files="/".join(str(o[v]["files"]) for v in VARS),
                versions=",".join(sorted({str(o[v]["version"]) for v in VARS})),
                note="" if len(set(starts)) == 1 and len(set(ends)) == 1 and all(o[v]["files"] for v in VARS)
                     else "vars differ: " + "; ".join(f"{v} {o[v]['start']}-{o[v]['end']} ({o[v]['files']})" for v in VARS)))
        write_tsv(os.path.join(args.out, f"summary_{tag}.tsv"), rows)
        yrs = lambda rs: "{}-{}".format(min(r["start"][:4] for r in rs), max(r["end"][:4] for r in rs))
        print("\nmembers with all 4 variables (common period from filenames):")
        matrix_table(rows, args.out, tag, lambda rs: f"{len(rs)} ({yrs(rs)})")
        print_table(rows, ["model", "experiment", "member", "grid", "start", "end", "files", "versions", "note"])
        print(f"\n-> {args.out}/{{matrix,summary}}_{tag}.tsv")
        return

    rows, issues = [], []
    if args.best:
        cand = {}
        for k in combos:
            cand.setdefault((k[0], k[1]), []).append(k)
        pending = {e: sorted(ks, key=lambda k: (-qspan[k][2], natural(k[2]), k[3]))[:args.max_tries]
                   for e, ks in cand.items()}
        rnd = 0
        while pending:
            rnd += 1
            batch = {e: ks.pop(0) for e, ks in pending.items()}
            print(f"\nround {rnd}: checking {len(batch)} members", flush=True)
            res = run(scan_var, {k: combos[k] for k in batch.values()}, args.njobs)
            for e, k in batch.items():
                row, iss = assess(k, res[k])
                rows.append(row)
                issues += iss
                if row["status"] == "OK" or not pending[e]:
                    del pending[e]
        rows.sort(key=lambda r: (r["model"].lower(), exp_order(r["experiment"]), natural(r["member"])))
    else:
        res = run(scan_var, combos, args.njobs)
        for k in order(res):
            row, iss = assess(k, res[k])
            rows.append(row)
            issues += iss

    write_tsv(os.path.join(args.out, f"summary_{tag}.tsv"), rows)
    with open(os.path.join(args.out, f"issues_{tag}.tsv"), "w") as fh:
        fh.write("model\texperiment\tmember\tgrid\tvar\tfile\tissue\n")
        fh.writelines("\t".join(map(str, i)) + "\n" for i in issues)
    ens = ensemble_table(rows, members)
    write_tsv(os.path.join(args.out, f"ensembles_{tag}.tsv"), ens)

    print("\nlongest fully checked period per model/experiment:")
    matrix_table(ens, args.out, tag, lambda rs: f"{rs[0]['start'][:4]}-{rs[0]['end'][:4]}" if rs[0]["status"] == "OK" else "FAIL")
    print("\nmodels (grid/levels of the checked member):")
    seen, models = set(), []
    for r in ens:
        if r["model"] not in seen and r["status"] == "OK":
            seen.add(r["model"])
            models.append(r)
    print_table(models, ["model", "atmos_model", "dycore", "dycore_detail", "shape", "dlon", "dlat", "nominal",
                         "grid", "regridded", "nlev", "n200_1000", "ptop_hPa", "calendar"])
    print_table(ens, ["model", "experiment", "n_members", "clean", "member", "start", "end", "years", "status"])
    print(f"\n{sum(r['status'] == 'OK' for r in rows)}/{len(rows)} checked members clean; "
          f"{len(issues)} issues -> {args.out}/issues_{tag}.tsv")
    for i in issues[:60]:
        print("  " + " | ".join(map(str, i)))
    if len(issues) > 60:
        print(f"  ... {len(issues) - 60} more in the issues file")
    print(f"\n-> {args.out}/{{ensembles,matrix,summary,issues}}_{tag}.tsv")


if __name__ == "__main__":
    main()
