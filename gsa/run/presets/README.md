# GSA run presets

Extra `run_gsa_srp` arguments, read by `run_gsa_gadi.pbs` (`PRESET=<name>`). Each file is one line of MATLAB
that is appended to `run_gsa_srp(pwd, h, m, shifts ...)`.

- `gridW3_u`, `gridW3_med`, `gridW3_fine` (grid study, 2026-11-12): W3 shock smoothing (threshold 0.05,
  filter 0.3 mm, blend 0.35 mm, lip region x <= 3 mm untouched), alpha_rho 3e-5, no adjoint, mu_t filter 0.25 mm,
  and the same base flow for every grid (built on a uniform 0.0177 mm helper grid, then interpolated).
  `_u` is the uniform grid at h; `_med` and `_fine` are the meshMedium / meshFine cell sizes in the SU2 box
  (0.025 / 0.0177 mm in x 0.5-20, r 0-9.5 mm) with h outside (0.05 mm).
