# SA_p0j8.015_p0inf3.247 on meshFineJet2

Confirmed by the user on 2026-10-29.

- **Physics and boundary conditions:** as the meshFineJet / meshVeryFine runs of this case. The stage configs are made
  from the very-fine `v3_cap3.cfg`, with only the solution, CFL, iterations, mesh and outputs changed.
- **Mesh:** `mesh/srp_M2.5_meshFineJet2.su2` (1,017,931 cells). It is meshFineJet with the jet box (x 1.5-8.5 mm,
  r 0-3.6 mm) at 0.00884 mm instead of 0.0125 mm.
- **Start:** meshFineJet's final solution (`meshFineJet/2_jet_staged/j3_cap3/restart.dat`, end of j3), interpolated with
  SU2's restart interpolation (`0_interp`).
- **Stages:** each runs its full length. There is no early exit and no stop on convergence.

  | Stage | CFL cap | Starting CFL | Iterations |
  |---|---|---|---|
  | j1 | 0.5 | 0.1 | 5000 |
  | j2 | 1 | 0.5 | 25000 |
  | j3 | 3 | 1 | 20000 |

  All stages use a CFL floor of 0.1 and ramp x1.01 per iteration.
- **Monitoring:**
  - Probes p1-p7, the same as meshFineJet.
  - A jet-box snapshot every 1000 iterations, plotted as the latest minus the previous one.
  - A live convergence plot copied into figures/.
  - The watcher logs only.
- **Analysis at the end:**
  - Line profiles at r = 0.01, 1.3, 2.3, 3, 6 and 9 mm, x 1-20 mm.
  - End of j2 vs end of j3.
  - fineJet j3 / fineJet2 j3 / very fine v3.
