#!/usr/bin/env python
"""Flag models whose diagnosed omega does not reach the scores of the site's reference model.

Scores come from the omega_check CSVs (results/<model>/**/omega_check_<model>*.csv). Instantaneous comparisons
always come first; only a model without one is judged on daily means, and its status says so (PASS_DAILY):
    instantaneous (E3hrPt, CFsubhr), benchmark MPI-ESM1-2-HR on every site (no other model has one on E3hrPt:
    on spiritx keep its Levante CSV in results/MPI-ESM1-2-HR/levante/)
    > daily on model levels (CFday) > daily on pressure levels (day/Eday), benchmark: the site's reference
    model, spiritx MRI-ESM2-0 (which has no instantaneous reference at all), Levante MPI-ESM1-2-HR.
in the N Atlantic-Europe region at 850-250 hPa (global if a file has no such region), plus the layer
within 50 hPa of the surface in the midlatitudes. It fails a level if
    r_time   < reference r_time - TOL_R
    |slope-1| > |reference slope-1| + TOL_SLOPE
    RMSE %   > reference RMSE % + TOL_RMSE
FLAG means "not shown to be as good as the reference", not "wrong": daily scores contain the sampling
error of four instants, which grows with resolution (MPI-ESM1-2-HR: r 0.99-1.00 instantaneous, 0.90-0.95
daily). Where a day table is a plain mean of the same four instants (EC-Earth3) that error is absent, so
such a model and one with a true daily mean do not face the same test. omega_check --daily finds the
definition from ua and scores against it; the table records it ("no definition check" for older files).

The separate column "criterion" is a fixed limit, not a comparison with the reference: it decides whether a model
is left out when the Lagrangian analysis is reviewed. Flagged models still run. On the model's best comparison, in
the same region, at any of 850-250 hPa it fails if
    r_time < MIN_R  or  RMSE % > MAX_RMSE  or  slope outside SLOPE_RANGE
which gives EXCLUDE. A failure on daily means only is REVIEW while within DAILY_NEAR of those limits, as the
sampling error alone takes MPI-ESM1-2-HR from r 0.997-1.000 and RMSE 2-8 % to r 0.90-0.96 and RMSE 23-31 %.
A model in REFERENCE_IN_DOUBT is at most REVIEW. The layer near the surface is not part of it: EC-Earth3 has
RMSE 42 % at 1000 hPa with r 0.992-0.999 at 850-500 hPa.

    python omega_benchmark.py --site spiritx --results /homedata/$USER/omega_check/results
    python omega_benchmark.py --chain origin_contributions_zerocross.nc     # quality of a finished run
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

__version__ = "1.1"
REFERENCE = {"spiritx": "MRI-ESM2-0", "levante": "MPI-ESM1-2-HR"}
TOL_R, TOL_SLOPE, TOL_RMSE = 0.02, 0.03, 10.0
LEVELS = (850, 700, 500, 250)
INSTANT_REFERENCE = "MPI-ESM1-2-HR"
TIERS = [("instantaneous", ("E3hrPt", "CFsubhr"), "6-hourly instants"),
         ("daily, model levels", ("CFday",), "daily, model levels"),
         ("daily, pressure levels", ("day", "Eday", "E3hrPt", "6hrPlev"), "daily total")]   # with a sub-daily run in
         # the same file the reference column holds that table, "daily total" is still the comparison with day wap
REGION = "N Atlantic-Europe"
# no benchmark possible, with the reason
NOT_COMPUTABLE = {"ACCESS-ESM1-5": "hybrid height coordinate, omega not implemented",
                  "IPSL-CM6A-MR1": "piControl only, no historical JJA 2000 run"}
# fixed limits of the criterion column; daily scores near them (margins as in the docstring) stay under review
MIN_R, MAX_RMSE, SLOPE_RANGE = 0.9, 50.0, (0.8, 1.2)
DAILY_NEAR = dict(r=0.8, rmse=75.0, slope=(0.7, 1.3))
# scored against a wap that may not be the true omega
REFERENCE_IN_DOUBT = {"IPSL-CM6A-LR": "wap looks like the mass flux through hybrid surfaces, not omega (925 hPa r 0.26 vs omega, "
                                       "0.83 vs that flux), but the flux from 6hrLev winds is itself too weak (slope 0.6), open",
                      "IPSL-CM6A-LR-INCA": "same model code as IPSL-CM6A-LR, wap in doubt"}
# chain benchmark of a finished run, set by MRI-ESM2-0 historical 1953: valid 100 %, gap +0.004 K, p99 0.51 K
CHAIN_VALID, CHAIN_GAP, CHAIN_P99 = 99.5, 0.05, 1.0


def load_scores(results, model):
    """All score rows of a model; with and without the reduced-grid reconstruction, the former wins."""
    files = [f for f in glob.glob(os.path.join(results, model, "**", f"omega_check_{model}*.csv"), recursive=True)
             if "window_scan" not in f and "_scan" not in os.path.basename(f)]
    rows = {}
    for f in sorted(files, key=lambda x: ("unregrid" in x, os.path.getmtime(x))):
        d = pd.read_csv(f)
        for (ref, grid), g in d.groupby(["reference", "grid"]):
            rows[(ref, grid)] = g.assign(file=f)      # later files replace earlier ones
    return rows


def pick(rows, tier):
    _, refs, grid = tier
    for r in refs:
        if (r, grid) in rows:
            return r, rows[(r, grid)]
    return None, None


def at_levels(g, region):
    g = g[(g["region"] == region) & g["agl_bin"].isna()]
    out = {}
    for lev in LEVELS:
        if len(g):
            k = (g["plev_hPa"] - lev).abs().idxmin()
            if abs(g.loc[k, "plev_hPa"] - lev) <= 0.15 * lev:
                out[lev] = g.loc[k]
    return out


def near_surface(g):
    g = g[(g["region"] == "midlatitudes") & (g["agl_bin"].astype(str) == "0-50")]
    return {"0-50": g.iloc[0]} if len(g) else {}


def failures(cand, ref):
    bad = []
    for lev in cand:
        if lev not in ref:
            continue
        c, r = cand[lev], ref[lev]
        rc, rr = c["r_time"], r["r_time"]
        if not np.isfinite(rc) or not np.isfinite(rr):      # near-surface rows carry r, not r_time
            rc, rr = c["r"], r["r"]
        if rc < rr - TOL_R:
            bad.append(f"{lev}: r {rc:.2f} (ref {rr:.2f})")
        if abs(c["slope"] - 1) > abs(r["slope"] - 1) + TOL_SLOPE:
            bad.append(f"{lev}: slope {c['slope']:.2f} (ref {r['slope']:.2f})")
        if c["rmse_pct"] > r["rmse_pct"] + TOL_RMSE:
            bad.append(f"{lev}: RMSE {c['rmse_pct']:.0f}% (ref {r['rmse_pct']:.0f}%)")
    return bad


def criterion(model, rows):
    """Fixed limits on the model's best comparison: (PASS, REVIEW, EXCLUDE or -, details)."""
    for tier in TIERS:
        ct, cg = pick(rows, tier)
        if cg is None:
            continue
        levels = at_levels(cg, REGION if (cg["region"] == REGION).any() else "global")
        if not levels:
            continue
        instant = tier[0] == "instantaneous"
        bad, far = [], False
        for lev, c in levels.items():
            r = c["r_time"] if np.isfinite(c["r_time"]) else c["r"]
            if r < MIN_R:
                bad.append(f"{lev}: r {r:.2f}")
            if c["rmse_pct"] > MAX_RMSE:
                bad.append(f"{lev}: RMSE {c['rmse_pct']:.0f}%")
            if not SLOPE_RANGE[0] <= c["slope"] <= SLOPE_RANGE[1]:
                bad.append(f"{lev}: slope {c['slope']:.2f}")
            far = far or r < DAILY_NEAR["r"] or c["rmse_pct"] > DAILY_NEAR["rmse"] \
                or not DAILY_NEAR["slope"][0] <= c["slope"] <= DAILY_NEAR["slope"][1]
        worst = f"{tier[0]} vs {ct}: r >= {min(c['r_time'] if np.isfinite(c['r_time']) else c['r'] for c in levels.values()):.3f}, " \
                f"RMSE <= {max(c['rmse_pct'] for c in levels.values()):.0f}%"
        if not bad:
            return "PASS", worst
        status = "EXCLUDE" if instant or far else "REVIEW"
        note = [REFERENCE_IN_DOUBT[model]] if model in REFERENCE_IN_DOUBT else []
        return ("REVIEW" if note else status), f"{tier[0]} vs {ct}: " + "; ".join(bad + note)
    return "-", "no score at 850-250 hPa"


def benchmark(model, ref_rows, rows):
    if model in NOT_COMPUTABLE:
        return "NOT_BENCHMARKED", "-", NOT_COMPUTABLE[model]
    if not rows:
        return "NOT_BENCHMARKED", "-", "no score file"
    for tier in TIERS:
        rt, rg = pick(ref_rows, tier)
        ct, cg = pick(rows, tier)
        if rg is None or cg is None:
            continue
        region = REGION if (rg["region"] == REGION).any() and (cg["region"] == REGION).any() else "global"
        bad = failures(at_levels(cg, region), at_levels(rg, region)) + failures(near_surface(cg), near_surface(rg))
        where = f"{tier[0]} vs {ct}, {region}" + (", reconstructed" if "unregrid" in cg["file"].iloc[0] else "")
        if tier[0].startswith("daily, pressure"):
            where += ", " + (str(cg["day_definition"].iloc[0]) if "day_definition" in cg else "no definition check")
        if tier[0] != "instantaneous":
            where += ", NO INSTANTANEOUS CHECK"
        return ("FLAG", where, "; ".join(bad)) if bad else ("PASS" if tier[0] == "instantaneous" else "PASS_DAILY", where, "")
    return "NOT_BENCHMARKED", "-", "no comparison shared with the reference: " + ", ".join(sorted(f"{a}/{b}" for a, b in rows))


def chain(path):
    import xarray as xr
    ds = xr.open_dataset(path)
    ta, ap = ds["TA"].values, ds["arrival_p"].values
    s = sum(ds[v].values for v in ("advection", "adiabatic", "diabatic", "seasonal"))
    m = np.isfinite(ap)
    ok = m & np.isfinite(ta) & np.isfinite(s)
    gap = (ta - s)[ok]
    valid, mean, p99 = 100 * ok.sum() / m.sum(), gap.mean(), np.percentile(np.abs(gap), 99)
    bad = [x for x, b in ((f"valid {valid:.1f} %", valid < CHAIN_VALID), (f"mean gap {mean:+.3f} K", abs(mean) > CHAIN_GAP),
                          (f"p99 |gap| {p99:.2f} K", p99 > CHAIN_P99)) if b]
    print(f"{path}: valid {valid:.2f} %, mean TA - sum of terms {mean:+.3f} K, p99 |gap| {p99:.2f} K -> "
          + ("FLAG: " + "; ".join(bad) if bad else "PASS"))
    return bool(bad)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", choices=sorted(REFERENCE))
    ap.add_argument("--results", default=None)
    ap.add_argument("--models", nargs="*", default=None, help="default: every folder in --results")
    ap.add_argument("--out", default=None, help="default: <results>/omega_benchmark.tsv")
    ap.add_argument("--chain", nargs="+", default=None, help="origin_contributions*.nc of finished runs")
    a = ap.parse_args()
    if a.chain:
        return 1 if any([chain(f) for f in a.chain]) else 0
    if not (a.site and a.results):
        ap.error("--site and --results are needed")
    refm = REFERENCE[a.site]
    ref_rows = {**load_scores(a.results, INSTANT_REFERENCE), **load_scores(a.results, refm)}
    if not ref_rows:
        sys.exit(f"no scores of the reference model {refm} in {a.results}")
    models = a.models or sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(a.results, "*", ""))
                                 if os.path.basename(os.path.dirname(p)) != "logs")
    models = sorted(set(models) | set(NOT_COMPUTABLE))
    lines = [f"model\tstatus\tcriterion\tcompared in\tdetails (reference model {refm}, omega_benchmark {__version__})\tcriterion details"]
    out = a.out or os.path.join(a.results, "omega_benchmark.tsv")
    for m in models:
        rows = load_scores(a.results, m)
        if m == refm:
            has_instant = pick(rows, TIERS[0])[1] is not None
            st, where, why = "REFERENCE", "-", "" if has_instant else "NO INSTANTANEOUS CHECK"
        else:
            st, where, why = benchmark(m, ref_rows, rows)
        crit, cwhy = criterion(m, rows) if rows else ("-", "")
        lines.append(f"{m}\t{st}\t{crit}\t{where}\t{why}\t{cwhy}")
        print(f"{m:20s} {st:16s} {crit:8s} {where:45s} {why}" + (f"  | {cwhy}" if cwhy else ""))
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("->", out)


if __name__ == "__main__":
    sys.exit(main())
