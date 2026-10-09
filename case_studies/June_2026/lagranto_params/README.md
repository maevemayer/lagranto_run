# LAGRANTO parameter files for stage 04

`tracevars.shellscript` (from Levante) is copied into each run directory. Its
order and scale factors are what stage 05 reads and undoes:

| var | factor | file | stage 05 |
|---|---|---|---|
| T | 1 | P | |
| OMEGA | 100 | P | `/100` |
| THETA | 1 | S | |
| TCLIM | 1 | S | |
| DTDP | 1000 | S | `/1e3` |
| DTDT | 1000000 | S | `/1e6` |

Change one and you must change `05_traj_postprocess_6h.py` with it.

`create_startf.param` is not needed: `create_startf.sh` writes it on every call.
