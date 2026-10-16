#!/usr/bin/env python3
"""Reduced-resolution copy of a run's P/S files (resolution experiments).

Two steps per file, both with cdo:
  1. vertical : keep NLEV model levels below ~100 hPa (cdo sellevel; lists below)
  2. horizontal (optional): conservative regrid to a coarser regular grid,
     0-90N (cdo remapcon to RES x RES degrees)
Both are linear and act per level / per column, so the order does not change
the result; levels are reduced first.

The output is written back in the layout LAGRANTO reads (lagranto_io), with
the hybrid coefficients of the source file, and p/dp (P-files only)
recomputed from the reduced/regridded surface pressure with the stage 01
formula (level k: p = hyam[k] + hybm[k]*ps, dp = interfaces k+1 minus k).

Usage (source config.sh first, or run via reduce_resolution.slurm):
    python reduce_resolution.py SRC_DATA_DIR DST_DATA_DIR --nlev 60|30 [--res 1.5|2] [--jobs 16]
    python reduce_resolution.py SRC DST --nlev 30 --res 2 --test    # one P + one S file, with checks
"""
import argparse
import glob
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lagranto_io import write_lagranto_netcdf  # noqa: E402

# ERA5 L137 model levels below ~100 hPa (levels 61..137; level 61 = 104 hPa at a
# 1013.25 hPa surface), thinned evenly in level index, so the L137 spacing
# (dense near the ground) is kept. Lowest level 137 and level 61 always kept.
LEVEL_SETS = {
    60: [61, 62, 64, 65, 66, 67, 69, 70, 71, 73, 74, 75, 76, 78, 79, 80, 82, 83, 84, 85,
         87, 88, 89, 91, 92, 93, 94, 96, 97, 98, 100, 101, 102, 104, 105, 106, 107, 109, 110, 111,
         113, 114, 115, 116, 118, 119, 120, 122, 123, 124, 125, 127, 128, 129, 131, 132, 133, 134, 136, 137],
    30: [61, 64, 66, 69, 71, 74, 77, 79, 82, 85, 87, 90, 92, 95, 98, 100, 103, 106, 108, 111,
         113, 116, 119, 121, 124, 127, 129, 132, 134, 137],
}
CDO = os.environ.get("CDO", "cdo")
COEFFS = ("hyai", "hybi", "hyam", "hybm")


def gridfile(res, path):
    nlon, nlat = int(round(360 / res)), int(round(90 / res)) + 1      # 0..360-res, 0..90N
    with open(path, "w") as f:
        f.write(f"gridtype = lonlat\nxsize = {nlon}\nysize = {nlat}\n"
                f"xfirst = 0\nxinc = {res}\nyfirst = 0\nyinc = {res}\n")
    return path


def p_dp(src, ps_hpa):
    """Stage 01 formula on the kept levels; ps in hPa (P-file PS)."""
    lev = src["lev"].values.astype(int); k = lev - 1
    ps = ps_hpa * 100.0
    c = {n: np.asarray(src[n].values).reshape(-1, src[n].shape[-1])[0] for n in COEFFS}   # 1-D, also if (time, n)
    am, bm = c["hyam"][k], c["hybm"][k]
    at, bt = c["hyai"][k], c["hybi"][k]
    ab, bb = c["hyai"][k + 1], c["hybi"][k + 1]
    ps4 = ps.values[..., None]                                       # (time, lat, lon, 1)
    p = am + bm * ps4
    dp = (ab + bb * ps4) - (at + bt * ps4)
    dims = ("time", "lat", "lon", "lev")
    return (xr.DataArray(p, dims=dims, attrs={"units": "Pa"}),
            xr.DataArray(dp, dims=dims, attrs={"units": "Pa"}))


def reduce_file(args):
    src_path, dst_path, levels, grid = args
    with xr.open_dataset(src_path) as src0:
        src = src0.load()
    has_p = "p" in src
    # cdo gets the gridded fields only (p/dp have the LAGRANTO dim order, which
    # cdo does not read; they are recomputed below). sellevel drops every field
    # without the selected levels, so 3-D fields (sellevel [+ remapcon]) and
    # 2-D fields like PS (remapcon only) go through cdo separately.
    cin = src.drop_vars([v for v in ("p", "dp", *COEFFS) if v in src])
    v3 = [v for v in cin.data_vars if "lev" in cin[v].dims]
    v2 = [v for v in cin.data_vars if "lev" not in cin[v].dims]
    parts = []
    with tempfile.TemporaryDirectory(dir=os.path.dirname(dst_path)) as tmp:
        for names, sel in ((v3, True), (v2, False)):
            if not names:
                continue
            part = cin[names]
            if sel:                                   # plain level axis: the model level number (selected by value)
                part["lev"].attrs = {"axis": "Z", "long_name": "model level number"}
            else:
                part = part.drop_vars("lev", errors="ignore")
            if not sel and not grid:                  # 2-D and no regrid: unchanged
                parts.append(part)
                continue
            f0, f1 = os.path.join(tmp, f"in{int(sel)}.nc"), os.path.join(tmp, f"out{int(sel)}.nc")
            part.to_netcdf(f0)
            ops = ([f"-remapcon,{grid}"] if grid else []) + (["-sellevel," + ",".join(map(str, levels))] if sel else [])
            subprocess.run([CDO, "-s", "-O", "-f", "nc4", *ops, f0, f1], check=True)
            with xr.open_dataset(f1) as r0:
                parts.append(r0.load())
    parts = [x.drop_vars([v for v in x.variables if v.endswith("_bnds")], errors="ignore") for x in parts]
    red = xr.merge(parts, join="override", compat="override")
    red = red.drop_vars([v for v in red.variables if v.endswith("_bnds") or v in COEFFS], errors="ignore")
    keep_lev = red["lev"].values.astype(int)
    if list(keep_lev) != list(levels):
        raise ValueError(f"{src_path}: cdo kept levels {list(keep_lev)}")
    out = red.assign_coords(time=src["time"].values)
    missing = [v for v in out.data_vars if v not in src]
    if missing:
        raise ValueError(f"{src_path}: cdo output has unexpected variables {missing}")
    for v in list(out.data_vars):                                    # original attributes + dtype
        out[v] = out[v].astype(src[v].dtype)
        out[v].attrs = src[v].attrs
    for c in ("lat", "lon", "lev"):
        out[c].attrs = src[c].attrs
    for v in COEFFS:                                                 # full coefficient arrays, as in the source
        if v in src:
            out[v] = src[v]
    if has_p:
        out = out.assign_coords(lev=red["lev"].values.astype(src["lev"].dtype))
        p, dp = p_dp(src.sel(lev=out["lev"].values), out["PS"])
        out["p"] = p.assign_coords(time=out.time, lat=out.lat, lon=out.lon, lev=out.lev).astype(src["p"].dtype)
        out["dp"] = dp.assign_coords(time=out.time, lat=out.lat, lon=out.lon, lev=out.lev).astype(src["dp"].dtype)
        out["p"].attrs, out["dp"].attrs = src["p"].attrs, src["dp"].attrs
    out = out[[v for v in src.data_vars if v in out]]                # source variable order
    out.attrs = dict(src.attrs, reduced_resolution=f"model levels {','.join(map(str, levels))}"
                     + (f"; cdo remapcon to {os.path.basename(grid)}" if grid else ""))
    tmp_out = dst_path + ".tmp"
    write_lagranto_netcdf(out, tmp_out)
    os.replace(tmp_out, dst_path)
    return os.path.basename(dst_path)


def check(src_path, dst_path, levels, res):
    """Sanity checks of one reduced file against its source."""
    a, b = xr.open_dataset(src_path), xr.open_dataset(dst_path)
    print(f"\n{os.path.basename(dst_path)}: {dict(b.sizes)}")
    print(f"  variables: {list(b.data_vars)}  (source: {list(a.data_vars)})")
    assert list(b["lev"].values.astype(int)) == list(levels), "levels"
    if res is None:
        for v in b.data_vars:
            if "lev" in b[v].dims and v in a:
                d = float(np.abs(b[v] - a[v].sel(lev=b["lev"])).max())
                print(f"  {v:8s} max |reduced - source| on kept levels: {d:.3g}")
    else:
        assert abs(float(b["lat"].diff("lat").mean()) - res) < 1e-6, "grid"
        for v in [x for x in b.data_vars if "lat" in b[x].dims]:
            sa = a[v].sel(lev=b["lev"]) if "lev" in a[v].dims else a[v]
            wa, wb = np.cos(np.deg2rad(a["lat"])), np.cos(np.deg2rad(b["lat"]))
            ma = float(sa.weighted(wa).mean()); mb = float(b[v].weighted(wb).mean())
            print(f"  {v:8s} area mean source {ma:12.5g}   reduced {mb:12.5g}")
    if "p" in b:
        lev = b["lev"].values.astype(int) - 1
        hyam = np.asarray(b["hyam"].values).reshape(-1, b["hyam"].shape[-1])[0]
        hybm = np.asarray(b["hybm"].values).reshape(-1, b["hybm"].shape[-1])[0]
        pp = hyam[lev] + hybm[lev] * b["PS"].values[..., None] * 100
        print(f"  p consistent with hyam/hybm and PS: max diff {float(np.abs(b['p'].values - pp).max()):.3g} Pa")
    print(f"  time {b['time'].values}  units {b['time'].encoding.get('units')}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src"); ap.add_argument("dst")
    ap.add_argument("--nlev", type=int, required=True, choices=sorted(LEVEL_SETS))
    ap.add_argument("--res", type=float, help="target grid spacing in degrees (default: keep the grid)")
    ap.add_argument("--jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--test", action="store_true", help="only the first P and S file, with checks")
    a = ap.parse_args()
    levels = LEVEL_SETS[a.nlev]
    files = sorted(glob.glob(os.path.join(a.src, "P[0-9]*")) + glob.glob(os.path.join(a.src, "S[0-9]*")))
    files = [f for f in files if not f.endswith(".tmp")]
    if not files:
        sys.exit(f"no P/S files in {a.src}")
    if a.test:
        files = [f for f in files if os.path.basename(f).startswith("P")][:1] + \
                [f for f in files if os.path.basename(f).startswith("S")][:1]
    os.makedirs(a.dst, exist_ok=True)
    grid = gridfile(a.res, os.path.join(a.dst, f"grid_{a.res}deg.txt")) if a.res else None
    with open(os.path.join(a.dst, "RESOLUTION.txt"), "w") as fh:
        fh.write(f"source {a.src}\nmodel levels ({len(levels)}): {','.join(map(str, levels))}\n"
                 f"grid: {'cdo remapcon to %g deg, 0-90N' % a.res if a.res else 'unchanged'}\n")
    todo = [(f, os.path.join(a.dst, os.path.basename(f)), levels, grid) for f in files
            if a.test or not os.path.exists(os.path.join(a.dst, os.path.basename(f)))]
    print(f"{len(levels)} levels, grid {a.res or 'unchanged'}: {len(files)} files, {len(todo)} to write -> {a.dst}", flush=True)
    with ProcessPoolExecutor(a.jobs) as ex:
        for i, name in enumerate(ex.map(reduce_file, todo), 1):
            if i % 20 == 0 or i == len(todo):
                print(f"  {i}/{len(todo)} {name}", flush=True)
    if a.test:
        for f, d, lv, _ in todo:
            check(f, d, lv, a.res)
    print("done")


if __name__ == "__main__":
    main()
