# Live convergence monitor (srp_M2.5, user 2026-11-04)

Default convergence check for new runs (see the Simulation folder's `.claude/rules/convergence.md`). It replaces
`watch_fj.sh` / `box_snap.py` / `box_change.py` for new chains. Nothing here stops a run.

## What it makes, per stage, every 1000 iterations

| File | Content |
|---|---|
| `convergence.png` | Residuals, probe Mach, box change, mass flows, CFL, Mach difference map (x 0-16, r 0-6 mm) |
| `conv_lines.png` | Mach along each of the six lines for every snapshot of the stage, coloured by iteration; latest in black (latest mean / p99 / max \|dMach\| in the titles) |
| `conv_features.png` | 2x2 against iteration: Mach disk x, Mach disk r, bow shock x, stagnation x (all logged in features.csv) |
| `box_change.csv`, `line_change.csv`, `features.csv` | One row per snapshot (a few kB in all) |
| `line_history.npz` | Mach on the six lines of every snapshot (~34 kB per snapshot, ~1.7 MB per 50k iterations) |

Each PNG is also copied to the run folder as `*_latest.png`.

- **Lines:** r = 0.01, 1.3, 2.3, 3 and 5 mm, x 1-20 mm, plus the jet shear layer from (x, r) = (2.0, 2.34) to
  (7.84, 4.70) mm. Each has 2000 points and is sampled for Mach only.
- **Features:** as defined in `../grid_compare/feature_positions.py`. The Mach disk is the normal jet shock at about
  8.6 mm, and its r is the subsonic-behind edge. The shock crossing is where the barrel shocks meet on the axis (about
  4.5 mm).
- **Storage:** each `flow_<iter>.vtu` is deleted after sampling. Only the latest two `snaps/flow_<iter>.npz` are kept
  (about 2.3 MB each: a 0.025 mm box over x 0-20, r 0-9.5 mm, plus the line Mach). SU2's own `flow.vtu`/`restart.dat` stay.

## Using it in a chain

1. **Stage config:**
   - `OUTPUT_WRT_FREQ= <restart>, 1000, <surface>`
   - `WRT_VOLUME_OVERWRITE= NO`
   - add `CUSTOM` to `HISTORY_OUTPUT`
   - `CUSTOM_OUTPUTS= '<line printed by probes_md.py>'`
2. **Probes:** run `pvbatch --mesa probes_md.py <starting solution flow.vtu>` (on a compute node). It places probes
   1 mm upstream and 1 mm behind the Mach disk at r = 0, 1 and 2 mm (user, 2026-11-04).
3. **Watcher:** start it next to SU2: `(sleep 30; PREV=<previous stage or -> bash <clone>/srp_M2.5/runs/monitor/watch_conv.sh <stage>) &`

## Files

- `watch_conv.sh`: the watcher (`ONCE=1` for one pass).
- `snap.py` (pvbatch): box + lines for each snapshot.
- `conv_update.py`: the three CSVs, plus pruning.
- `plot_conv_roi.py`, `plot_conv_lines.py`, `plot_conv_features.py`: the three figures.
- `probes_md.py` (pvbatch): the probe line.
- `plot_chain.py <run folder> "<title>"`: the three figures for a whole chain (t1 + t2 + t3), on one cumulative iteration
  axis (each stage offset by the ITER of the stages before it). Writes `<run>/{convergence,conv_lines,conv_features}_chain.png`.
- `plot_grid_lines.py` / `plot_grid_features.py <out.png> "<title>" label=<t3 stage folder> ...`: the line and feature
  plots with one curve or point per mesh (each mesh's final snapshot), coarsest first.
- **Default report (user, 2026-11-05):** a request for convergence plots means `plot_chain.py` for each mesh plus the two
  grid plots, not the per-stage figures.
- `test_offline.sh`: tests the post-processing on existing solutions (no SU2); `MANY=1` uses seven. Both passed on
  2026-11-04 on m3f102.
