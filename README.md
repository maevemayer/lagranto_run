# lagranto_run

LAGRANTO backward-trajectory pipeline (ERA5, 6-hourly) on IPSL **spiritx**:
ERA5 → P-files, 9-year T climatology → S-files, 15-day back-trajectories,
process decomposition (advection / adiabatic / diabatic / seasonal) and origin
contributions.

| path | what |
|---|---|
| `case_studies/June_2026/` | the May/June 2026 heatwave event + analogues; see its README |
| `tools/check_spiritx_env.sh` | read-only check of Slurm, conda/cdo, storage, LAGRANTO and the ERA5 input |

LAGRANTO itself is not in this repo; it is compiled on spiritx at
`/home/mmayer/code/lagranto.era5.test/lagranto.shellscript`.

## First time on spiritx

```bash
git clone <this repo> && cd lagranto_run
bash tools/check_spiritx_env.sh          # then fix the [CHECK] items in config.sh
cd case_studies/June_2026
./submit_all.sh                          # event set in the EVENT block of config.sh
```
