# Gadi (NCI, project os85) tools

SU2 v8.5.0 install and PBS drivers for running the srp_M2.5 cases on Gadi. Everything lives under
`/scratch/os85/lp3540/su2/` (no `/g/data`). Rules: `.claude/rules/HPC_Guidelines_Gadi_M3.md` in the Simulation folder.

| File | What it does |
|---|---|
| `su2_env.sh` | Environment: `gcc/13.2.0`, `intel-mpi/2021.13.1`, `python3/3.11.7`, the venv and SU2 paths. Copied to `SU2-v8.5.0/su2_env.sh`. |
| `su2_setup_login.sh` | Login-node part of the install (needs internet): venv (mpi4py built against this Intel MPI), `git clone` of v8.5.0, meson configure with the M3 options. |
| `build_su2.pbs` | Compile and install (`ninja install`), 12 cores, 1 h, at most 24 SU. |
| `su2_stage.pbs` | One SU2_CFD stage: `qsub -l walltime=.. -v STAGE=<folder>,CFG=<cfg> su2_stage.pbs`. 48 ranks, `-t 1`. Writes `<cfg>.log` ending in `SU2_CFD exit=<rc>`, like the M3 drivers. |
| `setup_restart_test.sh` | Builds `SST_p0j7.4945_..._meshFineJet2/test_gadi/`: the first N iterations of t3 from a copy of the t2 restart. |
| `compare_history.py` | Row-by-row comparison of two `history.csv` files (M3 vs Gadi) and the time per iteration. |

Cost in `normal`: 2 SU per core-hour, charged on max(ncpus, mem/4 GB) x walltime used; estimate from the requested
walltime. A full 48-core node is 96 SU per hour.

Not ported yet: the watchers (`watch_fj.sh`, box snapshots, live convergence plot) and the chain drivers, which have
M3 paths and ParaView 5.13 hard-coded. Gadi has `paraview/5.12.0-mesa` and `6.1.1-mesa`.
