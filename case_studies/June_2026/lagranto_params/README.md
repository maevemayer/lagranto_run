# LAGRANTO parameter files for stage 04

`tracevars.shellscript` is copied into each run directory. Its order and scale
factors are what stage 05 reads and undoes:

| var | factor | file | stage 05 |
|---|---|---|---|
| T | 1 | P | |
| OMEGA | 100 | P | `/100` |
| THETA | 1 | S | |
| TCLIM | 1 | S | |
| TCLIM:+10HPA | 1 | S | T_clim 10 hPa below the parcel (p + 10 hPa) |
| TCLIM:-10HPA | 1 | S | T_clim 10 hPa above the parcel (p - 10 hPa) |
| DTDT | 1000000 | S | `/1e6` |

The two shifted TCLIM columns replace the old precomputed DTDP field: stage 05
computes d(T_clim)/dp along the trajectory from them (see the case README).
Change the order or a factor and you must change `05_traj_postprocess_6h.py`
with it.

`create_startf.param` is not needed: `create_startf.sh` writes it on every call.
