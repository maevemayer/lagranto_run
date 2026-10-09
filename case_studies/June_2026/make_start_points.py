#!/usr/bin/env python3
"""
Write the LAGRANTO start-point list (one "lon lat" per line) for stage 04.

The arrival box comes from config.sh (ARRIVAL_LON_MIN/MAX, ARRIVAL_LAT_MIN/MAX,
ARRIVAL_DLON/DLAT). Longitudes are given in -180..180 so a box across 0E is a
plain range; they are written in 0..360 like the P/S-file grid. Points run
lat-major, lon ascending - identical to the old hard-coded list for the
default box (-5.5..5.0E, 44.5..55.0N, 0.5 deg).

Usage: python make_start_points.py OUTFILE
"""
import os
import sys


def axis(lo, hi, step):
    n = int(round((hi - lo) / step)) + 1
    if n < 1:
        sys.exit(f"empty arrival axis: {lo}..{hi} step {step}")
    return [lo + step * i for i in range(n)]


def main():
    env = os.environ
    lon_min, lon_max = float(env["ARRIVAL_LON_MIN"]), float(env["ARRIVAL_LON_MAX"])
    lat_min, lat_max = float(env["ARRIVAL_LAT_MIN"]), float(env["ARRIVAL_LAT_MAX"])
    dlon, dlat = float(env["ARRIVAL_DLON"]), float(env["ARRIVAL_DLAT"])
    if lon_max < lon_min:                     # e.g. 350 .. 10 given in 0..360
        lon_max += 360.0
    lons = [lon % 360.0 for lon in axis(lon_min, lon_max, dlon)]
    lats = axis(lat_min, lat_max, dlat)
    # 1 decimal as before for 0.5/0.1 grids, more only if the step needs it
    nd = 1 if all(abs(round(x, 1) - x) < 1e-9 for x in lons + lats) else 3
    with open(sys.argv[1], "w") as f:
        for lat in lats:
            for lon in lons:
                f.write(f"{lon:.{nd}f} {lat:.{nd}f}\n")
    print(f"{len(lons)} x {len(lats)} = {len(lons) * len(lats)} start points "
          f"({lon_min}..{lon_max}E, {lat_min}..{lat_max}N)", flush=True)


if __name__ == "__main__":
    main()
