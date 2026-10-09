"""
Write P- and S-files in the netCDF layout LAGRANTO reads (shared by stages 01, 03).

Plain xarray.to_netcdf gives a layout that differs from the P-files LAGRANTO is
known to read on spiritx (the run/ test case):
  * hyai/hybi/hyam/hybm pick up a time dimension when a dataset is split per
    time step and re-expanded -> written back as 1-D coefficient arrays,
  * time is encoded as int64 "days since ..." -> double "hours since ...",
    unlimited time dimension,
  * every variable gets _FillValue = NaN -> no fill values.
Data values are unchanged.
"""
import pandas as pd

HYBRID_COEFFS = ("hyai", "hybi", "hyam", "hybm")


def write_lagranto_netcdf(ds, path, **kwargs):
    ds = ds.copy()
    for name in HYBRID_COEFFS:
        if name in ds and "time" in ds[name].dims:
            ds[name] = ds[name].isel(time=0, drop=True)
    if "time" in ds.coords:
        # encode time by hand: xarray would write "hours since 2022-07-18T18:00:00"
        # (ISO 'T'); the working P-files have "hours since 1999-12-31 18:00:00"
        times = pd.DatetimeIndex(ds["time"].values)
        t0 = times[0]
        hours = ((times - t0) / pd.Timedelta(hours=1)).to_numpy(dtype="float64")
        ds = ds.assign_coords(time=("time", hours, {
            "standard_name": "time", "axis": "T", "calendar": "proleptic_gregorian",
            "units": f"hours since {t0:%Y-%m-%d %H:%M:%S}"}))
    encoding = {name: {"_FillValue": None} for name in ds.variables}
    ds.to_netcdf(path, encoding=encoding, unlimited_dims=["time"], **kwargs)
