"""
ERA5 input reader for spiritx (shared by stages 01 and 02).

On DKRZ every variable sat in its own daily GRIB file (spectral T/U/V on model
levels). On spiritx the FLEXPART-input pool has ONE daily GRIB file holding all
fields:

    ${ERA5_DIR}/YYYY/ERA5ENYYYYMMDD.grb

so each variable is picked out with `cdo -selname` before the (unchanged)
setgridtype -> sp2gp -> remapbil -> sellonlatbox chain. Fields that are already
on a regular grid simply fall through sp2gp (it raises and is skipped, exactly
as on DKRZ for non-spectral fields).
"""
import os
from functools import lru_cache

import xarray as xr
from cdo import Cdo

if os.environ.get("CDO_BIN_DIR"):
    os.environ["PATH"] = os.environ["CDO_BIN_DIR"] + ":" + os.environ["PATH"]

ERA5_DIR = os.environ["ERA5_DIR"]
ERA5_FILE_PATTERN = os.environ.get("ERA5_FILE_PATTERN",
                                   "{year}/ERA5EN{year}{month:02d}{day:02d}.grb")
HOURS_CSV = os.environ.get("HOURS_CSV", "0,6,12,18")
LEVELS = os.environ.get("LEVELS", "").strip()
REMAP = os.environ.get("REMAP_GRID", "r720x360")
# file-name tag of the intermediate files: the grid they are on
GRID_TAG = "native" if REMAP.lower() in ("", "none", "native") else REMAP
FILE_TAG = f"{GRID_TAG}_6h"                 # e.g. t_2022-07-18_native_6h.nc
LATBOX = os.environ.get("DOMAIN_LATBOX", "0,360,0,90").strip()   # NH cut from the start

# pipeline variable -> GRIB shortName in the daily file
GRIB_NAMES = {
    "t":  os.environ.get("ERA5_NAME_T", "t"),
    "u":  os.environ.get("ERA5_NAME_U", "u"),
    "v":  os.environ.get("ERA5_NAME_V", "v"),
    "w":  os.environ.get("ERA5_NAME_W", "w"),
    "sp": os.environ.get("ERA5_NAME_SP", "sp"),
}
CDO_GRIB_OPTS = "--eccodes"


def era5_file(year, month, day):
    return os.path.join(ERA5_DIR, ERA5_FILE_PATTERN.format(year=year, month=month, day=day))


@lru_cache(maxsize=None)
def _names_in(path):
    return frozenset(Cdo().showname(input=path, options=CDO_GRIB_OPTS)[0].split())


def _select(var, path):
    """cdo operator chain that extracts `var` (named as the pipeline expects)."""
    grib_name = GRIB_NAMES[var]
    names = _names_in(path)
    if grib_name in names:
        sel = f"-selname,{grib_name}"
        if grib_name != var:
            sel = f"-chname,{grib_name},{var} {sel}"
        return sel
    if var == "sp" and "lnsp" in names:
        # surface pressure from log(sp) on model level 1
        return "-chname,lnsp,sp -exp -selname,lnsp"
    raise KeyError(f"{grib_name!r} not found in {path} (has: {' '.join(sorted(names))})")


def _check_hours(path, var_file, var):
    """Every day must provide all HOURS_CSV snapshots (00/06/12/18): a day with
    fewer would silently give a climatology/P-file without some slots."""
    want = sorted(int(h) for h in HOURS_CSV.split(","))
    with xr.open_dataset(path) as ds:
        got = sorted(int(h) for h in ds["time"].dt.hour.values) if "time" in ds.dims else []
    if got != want:
        os.remove(path)
        raise RuntimeError(f"{var_file}: {var} has hours {got} at {HOURS_CSV} UTC, "
                           f"expected {want} - check the time steps of this file "
                           f"(cdo -s showtime {var_file})")


def regrid_daily(year, month, day, var, out_dir, tag):
    """Regrid one day of `var` (REMAP_GRID, or native), keeping ONLY hours 00/06/12/18.

    Same chain as the DKRZ notebooks (setgridtype -> sp2gp -> remapbil) except:
      * -selname picks the field out of the all-variable daily file,
      * -selhour is applied first (fewer timesteps through the transforms),
      * -sellevidx is inserted when LEVELS is set (3-D vars only),
      * NO timselmean afterwards.
    Returns the output path, or None if the day's input file is missing.
    """
    cdo = Cdo()
    date_str = f"{year}-{month:02d}-{day:02d}"
    var_file = era5_file(year, month, day)
    if not os.path.exists(var_file):
        print(f"Skipped {var_file}: file does not exist", flush=True)
        return None

    reg_file = os.path.join(out_dir, f"reg_{var}_{date_str}.nc")
    gp_file = os.path.join(out_dir, f"gp_{var}_{date_str}.nc")
    regrid_file = os.path.join(out_dir, f"{var}_{date_str}_{tag}.nc")

    if os.path.exists(regrid_file):
        print(f"{regrid_file} already exists - skipping.", flush=True)
        return regrid_file

    # Step 1: select the field + 6-hourly snapshots (+ optional levels) FIRST,
    #         then convert reduced->regular Gaussian and GRIB->netCDF.
    sel = f"-selhour,{HOURS_CSV} {_select(var, var_file)}"
    if LEVELS and var != "sp":
        sel = f"-sellevidx,{LEVELS} {sel}"     # model levels by index, e.g. 18/137
    # --reduce_dim drops the size-1 level axis of sp derived from lnsp; only
    # for that field: it would also drop a time axis that has a single step
    opts = f"-f nc {CDO_GRIB_OPTS}"
    if var == "sp" and "-selname,lnsp" in sel:
        opts += " --reduce_dim"
    cdo.setgridtype("regular", input=f"{sel} {var_file}", output=reg_file, options=opts)
    _check_hours(reg_file, var_file, var)

    # Step 2: spectral -> gridpoint (pass-through if not spectral)
    try:
        cdo.sp2gp(input=reg_file, output=gp_file, options="-P 16")
        input_for_remap = gp_file
    except Exception as e:
        print(f"sp2gp failed or not needed for {var}: {e}", flush=True)
        input_for_remap = reg_file

    # Step 3: bilinear remap to REMAP_GRID - or, with REMAP_GRID=none, keep the
    #         native grid and only put latitudes south->north (the order the
    #         remap produced on Levante and LAGRANTO is known to read) - then
    #         cut to the domain box (default NH), all in one chained call so
    #         nothing outside the box is ever written.
    if REMAP.lower() in ("", "none", "native"):
        with xr.open_dataset(input_for_remap) as ds:
            lat = ds["lat"].values
        grid_op = "-invertlat" if lat[0] > lat[-1] else ""
    else:
        grid_op = f"-remapbil,{REMAP}"
    if LATBOX:
        cdo.sellonlatbox(LATBOX, input=f"{grid_op} {input_for_remap}", output=regrid_file)
    elif grid_op:
        cdo.copy(input=f"{grid_op} {input_for_remap}", output=regrid_file)
    else:
        cdo.copy(input=input_for_remap, output=regrid_file)

    for f in (reg_file, gp_file):
        if os.path.exists(f):
            os.remove(f)
    return regrid_file
