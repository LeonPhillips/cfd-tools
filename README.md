# cfd-tools

Tool and utility scripts for my SU2 simulations on MASSIVE M3: meshing, run drivers, monitoring and analysis.
Code only: no meshes, solutions, logs or figures (`.gitignore` blocks the common data types).

## Where it lives

| Copy | Path | Use |
|---|---|---|
| Edit here | `Simulation/scripts/` on the PC (OneDrive) | Edit, commit, push |
| GitHub | `cfd-tools` (private) | The master copy |
| M3 clone | `/scratch/qk75/LeonPhillips/PhD/cfd/su2/cfd-tools` | `git pull` before use; never edit on M3 |

Workflow: edit on the PC, commit, push, then `git pull` in `/scratch/qk75/LeonPhillips/PhD/cfd/su2/cfd-tools` on M3. Run folders on scratch call the scripts
from the M3 clone (or copy a pinned version in, noting the commit hash), and keep their outputs on scratch.

## Layout

```
scripts/
├── m3/                 small helpers piped to M3 with  ssh ... 'bash -s' < scripts/m3/<name>.sh
└── srp_M2.5/           supersonic retropropulsion, M 2.5 (case folder on M3:
    │                   /scratch/qk75/LeonPhillips/PhD/cfd/su2/cases/srp_M2.5/)
    ├── mesh/           Gmsh mesh generation and mesh comparison
    └── runs/<run>/     run drivers, watchers and analysis, at the same relative paths as on M3
```

`srp_M2.5/` was imported as-is from M3 on 2026-11-02 (115 files), so each script sits where it was used. Some tools
exist in several run folders: `box_change.py`, `plot_conv_roi.py`, `resample_box.py` and `interp.sh` differ between
copies, and the rest are identical. Most reusable versions (as of the import):
- line sampling and plots: `runs/SA_p0j8.015_p0inf3.247_meshCoarse/{line_profiles.py,plot_lines.py,check_run.py}`
- live monitoring: `runs/SSTsarkar_p0j8.015_p0inf3.247_tools/{watch_fj.sh,box_snap.py,box_change.py,plot_conv_roi.py}`
- box resampling and mirror contours: `runs/grid_compare/`
- synthetic schlieren vs experiment: `runs/SST_p0j8.015_p0inf3.247_meshFineJet2/analysis/schlieren/`

New shared tools should go into a common folder rather than being copied into each run.
