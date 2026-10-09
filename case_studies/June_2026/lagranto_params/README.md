# LAGRANTO parameter files for stage 04

Stage 04 copies these into each run directory. They are not here yet: copy them
from Levante (the paths the old `04_run_parallel.csh` used):

```bash
scp levante:/work/uc1275/u301827/01_LAGRANTO/2026_May/create_startf.param .
scp levante:/work/uc1275/u301827/01_LAGRANTO/tracevars.shellscript .
```

`submit_all.sh` refuses to submit while either file is missing.
