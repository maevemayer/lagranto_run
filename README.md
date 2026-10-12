# lagranto_run

LAGRANTO backward-trajectory pipeline (ERA5, 6-hourly) on IPSL **spiritx**:
ERA5 -> P-files, 9-year T climatology -> S-files, 15-day back-trajectories,
process decomposition (advection / adiabatic / diabatic / seasonal) and origin
contributions.

| path | what |
|---|---|
| `case_studies/June_2026/` | the pipeline (stages 01-06, `config.sh`, `submit_all.sh`); see its README |
| `environment.yml` | the conda env (pipeline + colindex2 + analysis) |
| `tools/check_spiritx_env.sh` | read-only check of Slurm, env, cdo, storage, LAGRANTO and the ERA5 input |

## Setting up (new user on spiritx)

Everything below runs on spiritx. Paths use `$USER`, so nothing needs editing
for a new user except where noted.

### 1. Access
- **This repository** is private on GitHub: ask for access, then add an SSH key
  on spiritx (`ssh-keygen -t ed25519`, paste `~/.ssh/id_ed25519.pub` into
  GitHub -> Settings -> SSH keys).
- **ERA5 input**: read access to `/proju/flexpart/flexpart_in/ERA5/EN-true/`
  (1-deg, 3-hourly, all fields per day; check with `ls` on that folder).
- Storage: `/scratchx/$USER` (netCDF output) and `/homedata/$USER` (env, logs,
  trajectories). `$HOME` only holds code.

### 2. Code
```bash
mkdir -p ~/code && cd ~/code && git clone git@github.com:maevemayer/lagranto_run.git
```

### 3. Conda env (once, ~30 min)
On a compute node, with micromamba (conda's own solver hangs on this env):
```bash
srun --partition=zen4 --cpus-per-task=4 --mem=16G --time=02:00:00 --pty bash
mkdir -p /homedata/$USER/bin && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /homedata/$USER bin/micromamba
export MAMBA_ROOT_PREFIX=/homedata/$USER/conda CONDA_PKGS_DIRS=/homedata/$USER/conda/pkgs
/homedata/$USER/bin/micromamba create -y -p /homedata/$USER/conda/thesis -f ~/code/lagranto_run/environment.yml
conda activate /homedata/$USER/conda/thesis && python -m ipykernel install --user --name thesis --display-name "thesis (spiritx)"
exit
```
`config.sh` activates `/homedata/$USER/conda/thesis` in every job.

### 4. LAGRANTO
LAGRANTO is not part of this repo. You need a compiled `lagranto.shellscript`
tree (with `caltra/caltra`, `trace/trace`, `create_startf/create_startf`);
`config.sh` expects it at `~/code/lagranto.era5.test/lagranto.shellscript`, or
set `export LAGRANTO_ROOT=/path/to/lagranto.shellscript` before submitting.
Two settings matter:
- `goodies/datelist.sh` must have `set interval = "6"` (6-hourly input files).
- `trace` segfaults on P/S files with more than ~120 model levels (tested:
  120 ok, 128 crash). The pipeline therefore keeps levels 18-137 (`LEVELS` in
  `config.sh`); do not change that unless your `trace` handles 137 levels.

### 5. Check
```bash
bash ~/code/lagranto_run/tools/check_spiritx_env.sh 20220718
```
Every line should read `ok`; the ERA5 section should show t, u, v, w on 137
hybrid levels plus lnsp, 3-hourly, on a global 1-deg grid.

### 6. Run an event
Edit the **EVENT block** at the top of `case_studies/June_2026/config.sh`
(run name, arrival days, arrival box, pressure levels, climatology years), then:
```bash
cd ~/code/lagranto_run/case_studies/June_2026 && ./submit_all.sh
squeue -u $USER
```
`submit_all.sh` derives the ERA5 period (arrival days - 15 days), checks the
LAGRANTO settings, refuses to start while jobs of the same run are queued, and
chains stages 01-06 with Slurm dependencies. A full event (12 days of
arrivals) takes about an hour on `zen4`.

### Where the output goes
| what | where |
|---|---|
| `origin_contributions.nc` (final result) | `/homedata/$USER/lagranto_run/<RUN_NAME>/results_6h/` (+ `/scratchx` original) |
| LAGRANTO trajectories + per-start logs | `/homedata/$USER/lagranto_run/<RUN_NAME>/results_6h/` |
| P/S files, climatology, decomposition netCDFs | `/scratchx/$USER/lagranto_run/<RUN_NAME>/` (may be purged) |
| Slurm logs | `/homedata/$USER/lagranto_run/slurm_logs/` |

Details on each stage, the climatology rule, levels and cleanup:
`case_studies/June_2026/README.md`.
