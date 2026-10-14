#!/usr/bin/env python3
"""
Check every climatology (tbar) file against the paper's definition:
a centred 9-year window (year +/- 4) of 21-day running means, clamped to
FIRST..LAST (1940..2025) so it never extends past 2025 - years 2021-2026 all
get 2017-2025, as in Hotz et al. (2024).

The years each file actually averages are read from its cdo history (the
'cdo cat ... t_YYYY-MM-DD_6hHH_mean.nc ...' step that built it).

Usage: python check_climatology_years.py TBAR_DIR [TBAR_DIR ...]
(reads only the history attribute, 8 processes in parallel; prints progress)
"""
import collections
import glob
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

try:                                    # fast: read only the global attribute
    import netCDF4

    def read_history(path):
        with netCDF4.Dataset(path) as nc:
            return getattr(nc, "history", "")
except ImportError:                     # fallback (e.g. a laptop without netCDF4)
    import h5netcdf

    def read_history(path):
        with h5netcdf.File(path, "r") as nc:
            h = nc.attrs.get("history", "")
            return h.decode() if isinstance(h, bytes) else h

FIRST, LAST, N = 1940, 2025, 9


def paper_years(year, n=N, first=FIRST, last=LAST):
    start, end = year - n // 2, year + n // 2
    if end > last:
        start, end = last - n + 1, last
    elif start < first:
        start, end = first, first + n - 1
    return list(range(start, end + 1))


def main():
    bad = 0
    for folder in sys.argv[1:]:
        files = sorted(glob.glob(os.path.join(folder, "t_*_tbar.nc")))
        print(f"\n{folder}: {len(files)} climatology files")
        summary = collections.Counter()
        # processes, not threads: the HDF5 library under netCDF4 is not thread-safe
        with ProcessPoolExecutor(max_workers=8) as pool:
            for i, (f, hist) in enumerate(zip(files, pool.map(read_history, files, chunksize=32)), 1):
                year = int(re.search(r"t_(\d{4})-\d\d-\d\d_6h\d\d_tbar\.nc", os.path.basename(f)).group(1))
                found = sorted({int(y) for y in re.findall(r"t_(\d{4})-\d\d-\d\d_6h\d\d_mean\.nc", hist)})
                want = paper_years(year)
                ok = found == want
                summary[(year, found[0] if found else None, found[-1] if found else None, len(found), ok)] += 1
                if not ok:
                    bad += 1
                    if bad <= 10:
                        print(f"  MISMATCH {os.path.basename(f)}: averages {found}, paper rule {want[0]}-{want[-1]}")
                if i % 250 == 0 or i == len(files):
                    print(f"  ... {i}/{len(files)} checked", flush=True)
        for (year, a, b, n, ok), cnt in sorted(summary.items()):
            print(f"  {year}: {cnt:4d} files average {a}-{b} ({n} years)  {'OK' if ok else 'NOT the paper rule'}")
    print(f"\n{'ALL files follow the paper definition' if bad == 0 else f'{bad} files do NOT follow it'}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
