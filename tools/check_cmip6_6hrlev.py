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

Writes nothing to the data. Output: summary.tsv + issues.tsv in --out.

  python check_cmip6_6hrlev.py [--exp historical] [--activity CMIP]
                               [--njobs 32] [--out DIR] [--only M1,M2]
"""
import argparse
import glob
import hashlib
import os
import re
from concurrent.futures import ProcessPoolExecutor, as_completed

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


def hsh(*arrs):
    m = hashlib.md5()
    for a in arrs:
        m.update(np.round(np.ma.getdata(np.asarray(a, "f8")), 6).tobytes())
    return m.hexdigest()[:10]


def discover(activity, exp, only):
    pat = os.path.join(ROOT, activity, "*", "*", exp, "*", "6hr[lL]ev")
    for d in sorted(glob.glob(pat)):
        inst, model, member = d.split(os.sep)[-5], d.split(os.sep)[-4], d.split(os.sep)[-2]
        if only and model not in only:
            continue
        grids = None
        for v in VARS:
            g = {os.path.basename(p) for p in glob.glob(os.path.join(d, v, "*")) if os.path.isdir(p)}
            grids = g if grids is None else grids & g
        for g in sorted(grids or ()):
            yield (model, member, g, inst), d


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
        else:
            return None
    except Exception:
        return None
    return p / 100.


def stamp(d, n):
    return d.strftime("%Y%m%d%H%M")[:n]


def scan_var(key, vdir, var):
    out = dict(var=var, issues=[], files=0, version=None, info={}, cal=set(),
               levh=set(), coefh=set(), gridh=set(), t=None)
    issue = lambda f, msg: out["issues"].append((var, os.path.basename(f) if f else "-", msg))
    d, out["version"] = version_dir(os.path.join(vdir, var, key[2]))
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
                        nominal=getattr(nc, "nominal_resolution", "?"),
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--activity", default="CMIP")
    ap.add_argument("--exp", default="historical")
    ap.add_argument("--njobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--out", default=f"/homedata/{os.environ.get('USER', 'mmayer')}/cmip6_check")
    ap.add_argument("--only", default="", help="comma-separated model names")
    args = ap.parse_args()
    only = set(filter(None, args.only.split(",")))
    os.makedirs(args.out, exist_ok=True)

    combos = dict(discover(args.activity, args.exp, only))
    print(f"{len(combos)} model/member/grid combos with {'+'.join(VARS)} in "
          f"{ROOT}/{args.activity}/*/*/{args.exp}/*/6hrLev", flush=True)
    res = {k: {} for k in combos}
    with ProcessPoolExecutor(args.njobs) as ex:
        futs = [ex.submit(scan_var, k, d, v) for k, d in combos.items() for v in VARS]
        for n, fu in enumerate(as_completed(futs), 1):
            k, o = fu.result()
            res[k][o["var"]] = o
            print(f"  [{n}/{len(futs)}] {k[0]} {k[1]} {o['var']}: {o['files']} files, "
                  f"{len(o['issues'])} issues", flush=True)

    rows, issues = [], []
    for k in sorted(res, key=lambda k: (k[0].lower(), k[1], k[2])):
        model, member, grid, inst = k
        o = res[k]
        iss = [i for v in VARS for i in o[v]["issues"]]
        # the 4 variables must share grid, levels, coefficients, calendar and time axis
        for name, label in (("gridh", "lat/lon"), ("levh", "lev values"),
                            ("coefh", "hybrid coefficients"), ("cal", "calendar")):
            if len({frozenset(o[v][name]) for v in VARS}) > 1:
                iss.append(("all", "-", f"{label} differ between variables"))
        ts = [o[v]["t"] for v in VARS]
        if any(t is None for t in ts):
            pass
        elif len({hsh(t) for t in ts}) > 1:
            iss.append(("all", "-", "time axes differ between variables: " + "; ".join(
                f"{v} {o[v]['start']}..{o[v]['end']} ({o[v]['nsteps']})" for v in VARS)))
        issues += [(model, member, grid, *i) for i in iss]

        r = o["ta"] if o["ta"]["t"] is not None else next((o[v] for v in VARS if o[v]["t"] is not None), o["ta"])
        info = r["info"] or {}
        rows.append(dict(
            model=model, member=member, grid=grid, institution=inst,
            nominal=info.get("nominal", "?"), shape=info.get("shape", "?"),
            dlat=f"{info.get('dlat', np.nan):.3g}", dlon=f"{info.get('dlon', np.nan):.3g}",
            nlev=info.get("nlev", "?"), vcoord=info.get("vcoord", "?"),
            ptop_hPa=f"{info.get('ptop', np.nan):.3g}", pbot_hPa=f"{info.get('pbot', np.nan):.4g}",
            lev_order=info.get("order", "?"), calendar=",".join(sorted(r["cal"])) or "?",
            first_hour=info.get("hour0", "?"),
            start=r.get("start", "-"), end=r.get("end", "-"),
            years=f"{r.get('years', 0):.1f}", nsteps=r.get("nsteps", 0), missing=max((o[v].get("missing", 0) for v in VARS), default="?"),
            files="/".join(str(o[v]["files"]) for v in VARS),
            versions=",".join(sorted({str(o[v]["version"]) for v in VARS})),
            ta_range="{:.4g}..{:.4g}".format(*o["ta"].get("range", {}).get("ta", (np.nan, np.nan))),
            hus_min="{:.3g}".format(o["hus"].get("range", {}).get("hus", (np.nan,))[0]),
            status="OK" if not iss else f"{len(iss)} issues"))

    cols = list(rows[0]) if rows else []
    with open(os.path.join(args.out, f"summary_{args.activity}_{args.exp}.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.writelines("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)
    with open(os.path.join(args.out, f"issues_{args.activity}_{args.exp}.tsv"), "w") as fh:
        fh.write("model\tmember\tgrid\tvar\tfile\tissue\n")
        fh.writelines("\t".join(map(str, i)) + "\n" for i in issues)

    show = ["model", "member", "grid", "shape", "dlon", "dlat", "nlev", "ptop_hPa", "calendar",
            "start", "end", "years", "missing", "files", "status"]
    w = {c: max([len(c)] + [len(str(r[c])) for r in rows]) for c in show}
    print("\n" + "  ".join(c.ljust(w[c]) for c in show))
    for r in rows:
        print("  ".join(str(r[c]).ljust(w[c]) for c in show))
    print(f"\n{sum(r['status'] == 'OK' for r in rows)}/{len(rows)} combos clean; "
          f"{len(issues)} issues -> {args.out}/issues_{args.activity}_{args.exp}.tsv")
    for i in issues[:60]:
        print("  " + " | ".join(map(str, i)))
    if len(issues) > 60:
        print(f"  ... {len(issues) - 60} more in the issues file")


if __name__ == "__main__":
    main()
