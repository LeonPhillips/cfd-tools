function p = srp_case_defaults()
% SRP_CASE_DEFAULTS  Parameters for build_srp_case (srp_M2.5, SST p0j 6.6 bar meshFine t3), 2026-11-09.
% Confirmed by the user: sigma_visc = 1, domain x -9..20 mm, r 0..27 mm (choice B), uniform masked grid,
% local shock smoothing (method), D and freestream scaling. Smoothing widths and sponge sizes below are
% first choices, to be confirmed from the check plots.

root = getenv('GSA_ROOT'); if isempty(root), root = '/scratch/qk75/LeonPhillips/PhD/GSA'; end   % M3 default; Gadi sets GSA_ROOT
p.baseflow_prefix = fullfile(root, 'srp_p0j6.6_meshFine_t3', 'baseflow', 'su2L');   % x -30..23 mm export

% reference (freestream, from t3_cap3.cfg; see gsa/README.md)
p.D = 0.018;
p.ref = struct('U', 567.162, 'T', 128.067, 'rho', 0.51693, 'mu', 8.8606e-6, ...
               'Re', 5.9559e5, 'Ma', 2.5, 'Pr', 0.72, 'kappa', 1.4);

% geometry (mm), from the wall_noslip surface of the mesh (x axial, jet towards +x)
p.x_exit_mm    = 1.876;   % nozzle lip tip / exit plane
p.x_face_mm    = 0.496;   % body front face
p.r_body_mm    = 9.0;     % body radius
p.r_lip_out_mm = 2.88;    % lip (nozzle tube) outer radius, x 0.496..1.876
p.r_exit_mm    = 2.27;    % nozzle exit radius (inner contour before the lip fillet)
p.x_inlet_mm   = -8.124;  % plenum inlet face (SU2 inlet_jet, total conditions)
p.r_plenum_mm  = 2.879;   % plenum radius

% nozzle interior (user, 2026-11-10): model the plenum + converging-diverging nozzle up to the inlet face.
% The inner wall is staircased from the mesh's wall contour. Inlet perturbation BC: Dirichlet on all
% variables (method's inlet) - to be confirmed.
p.nozzle_interior = true;

% body-side slip wall (user, 2026-11-10): the cylinder r = 9 mm behind the face is a slip wall, and the base-flow
% boundary layer there (within bodyside_bl_mm of the wall) is replaced by the values at that distance
p.bodyside_slip   = true;
p.bodyside_bl_mm  = 0.5;

% domain and grid
p.x_range_mm = [-27 20];  % outflow moved from -9 to -27 mm (user, 2026-11-10)
p.r_max_mm   = 27;
p.h_mm       = 0.1;       % 0.1 smoke test, 0.07 sweep, 0.05 production, 0.04 grid check
% stretched grid (user, 2026-11-12): spacing h_box_mm inside the SU2 refinement box (x 0.5-20, r 0-9.5 mm), h_mm
% outside, with smooth transitions (cell-to-cell growth <= max_growth). h_box_mm = [] gives the uniform grid.
% meshMedium: h_box 0.025, h 0.05; meshFine: h_box 0.0177, h 0.05.
p.grid = struct('h_box_mm', [], 'box_x_mm', [0.5 20], 'box_r_mm', 9.5, 'max_growth', 0.05);
% base flow built on a separate uniform grid of spacing bf_h_mm (smoothing, mu_t filter, sensor there), then
% interpolated (linear) to the GSA grid. [] = on the GSA grid itself (uniform grids only).
p.bf_h_mm    = [];

% effective viscosity
p.sigma_visc      = 1;    % MU = MU_Sutherland + sigma_visc * mu_t/mu_inf
p.mut_sigma_cells = 5;    % Gaussian smoothing of mu_t (method: 5 cells); 0 = none
p.mut_sigma_mm    = [];   % if set, the mu_t filter width in mm (overrides mut_sigma_cells; 5 cells at 0.05 = 0.25)

% local shock smoothing (method's smooth_shock; sigmas in mm here).
% Variant B, chosen by the user 2026-11-09: narrow filter, sensor off near the lip (x <= 3 mm) so the
% lip expansion and the start of the shear layer stay as in the CFD.
p.smooth = struct('sensor_sigma_mm', 0.15, 'threshold', 0.15, 'filter_sigma_mm', 0.1, ...
                  'blend_sigma_mm', 0.15, 'x_min_mm', 3, 'r_min_mm', -Inf);

% sponges (method's GetSponge shape; widths/transitions in mm)
p.sponge = struct('amp_in', 5,  'width_in_mm', 3,  'trans_in_mm', 1, ...   % freestream inflow, x max
                  'amp_out', 1.5, 'width_out_mm', 12, 'trans_out_mm', 4, ... % outflow, x min: long and gentle (method-like)
                  'amp_top', 5, 'width_top_mm', 3, 'trans_top_mm', 1);      % top, r max
end
