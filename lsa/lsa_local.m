%% lsa_local.m - locally parallel, spatial, inviscid linear stability analysis (LSA) of the srp_M2.5 jet
%
% WHAT THIS DOES
%   At each axial station x, it takes the radial profile of the RANS mean flow, treats it as if it did not change
%   in x (the "locally parallel" assumption) and finds which small waves of the form
%           q'(x, r, theta, t) = q_hat(r) * exp( i*(alpha*x + m*theta - omega*t) )
%   the profile supports. The frequency omega and azimuthal order m are chosen (real); the axial wavenumber
%   alpha = alpha_r + i*alpha_i is the unknown (complex). That is a SPATIAL analysis: alpha_i < 0 means the wave
%   grows as it travels in +x (because exp(i*alpha*x) = exp(i*alpha_r*x) * exp(-alpha_i*x)).
%
% METHOD (follows Tran, Nogueira, Weightman & Edgington-Mitchell 2026, JFM 1046 A14, section 2.3 and App. D)
%   - Inviscid compressible linearised Euler equations in cylindrical coordinates.
%   - State vector q = [nu, u_x, u_r, u_theta, p], nu = 1/rho (specific volume). Using nu instead of rho makes
%     every equation at most bilinear in the unknowns, which keeps the operator simple.
%   - The mean radial and azimuthal velocities are DROPPED (only U(r), nu(r), p(r) are kept). See the caveat below.
%   - Chebyshev collocation in r, with a mapping that puts many points in the shear layers.
%   - The result is a generalised eigenvalue problem  L*q_hat = alpha*R*q_hat , solved with MATLAB eig.
%
% CONFIRMED ASSUMPTIONS (user, 2026-11-10; also in README.md next to this file)
%   1. Base flow: SST_p0j6.6_p0inf3.247_meshFine, stage t3_cap3 (the same as the GSA), exported by
%      cfd-tools gsa/baseflow/export_baseflow.py to su2L_flow.bin/.json.
%   2. Inviscid, no eddy viscosity; mean u_r dropped (as Tran 2026).
%   3. Scales: length D = 18 mm (body diameter), velocity a_inf = 226.865 m/s, density rho_inf = 0.51693 kg/m3,
%      pressure rho_inf*a_inf^2 (so the free-stream pressure is 1/gamma). Frequencies reported as
%      St_D = f*D/U_inf (U_inf = 567.162 m/s) and in Hz.
%   4. Stations from just downstream of the nozzle exit (x_e = 1.876 mm, r_e = 2.389 mm, from the exported wall)
%      to x = 8.0 mm (the axis is subsonic by 8.5 mm: Mach disk), every 0.1 mm.
%   5. r_max just inside the bow shock at each station, Dirichlet conditions there.
%   6. Base-flow profiles smoothed (penalised least squares, no toolbox); Tran's tanh fit is NOT used here.
%   7. m = 0 and 1; blind scan St_D 0.05-2.
%   8. N = 250, clustering at the jet and outer shear layers; results only trusted if they do not move with N.
%   9. KH and upstream (guided-jet-like) wave spectra only: no closure / STIM frequency prediction.
%
% CAVEAT (found 2026-11-10 when looking at the profiles): the mean radial velocity is NOT small here
%   (|V|/U_max ~ 0.6-0.7 near the jet boundary, against "small" in Tran's flow), because the underexpanded jet is
%   still spreading and the counterflow turns outwards. Dropping V is the confirmed choice, but it is a strong
%   approximation for this flow; treat growth rates as indicative and compare with the GSA.
%
% HOW TO RUN (MATLAB r2023b or later, no toolboxes needed; Parallel Computing Toolbox used if present)
%   Set MODE in the workspace before calling the script, e.g. on M3 (on a compute node, not the login node):
%       module load matlab/r2023b
%       matlab -batch "MODE='selfcheck'; lsa_local"
%   MODE = 'selfcheck'   : thin tanh jet against the exact compressible vortex sheet (checks the operator)
%          'single'      : one station, one St: full spectrum, mode labels, eigenfunctions (for looking)
%          'sensitivity' : a few stations/St with N, r_max and smoothing varied (is the answer converged?)
%          'scan'        : all stations x all St x m = 0, 1 (the main result)
%          'profiles'    : all stations, base-flow profiles and fits only (no eigenvalues): check the tanh fits
%   Optional workspace variables: BASEFLOW_DIR, OUT_DIR (default to the M3 paths below), NWORKERS,
%   PROFILE ('tanh', the default: Tran's fitted profiles; or 'smooth': the smoothed RANS profile),
%   SINGLE_X [mm] and SINGLE_ST (station and St_D of MODE='single'; defaults 3.0 and 1.271).
%
% OUTPUT: <OUT_DIR>/<MODE>/ : a .mat file with everything, PNG figures, and a text log.
%
% Notation in comments: "Tran eq. D5" = equation D5 of Tran et al. 2026 (App. D). Overbars (mean quantities) are
% written as Ub, nub, pb; their r-derivatives as dUb, dnub, dpb.

%% 1. Settings ---------------------------------------------------------------------------------------------------
% Everything a user might change is here, in one struct S, so the rest of the file has no hidden constants.

if ~exist('MODE', 'var'), MODE = 'single'; end          % default mode when run interactively
S = struct();                                           % S collects all settings
S.mode = MODE;                                          % which part of the analysis to run

% --- paths (M3 defaults; override by setting BASEFLOW_DIR / OUT_DIR before calling the script)
if exist('BASEFLOW_DIR', 'var'), S.baseflow_dir = BASEFLOW_DIR;
else, S.baseflow_dir = '/scratch/qk75/LeonPhillips/PhD/GSA/srp_p0j6.6_meshFine_t3/baseflow'; end
if exist('OUT_DIR', 'var'), S.out_root = OUT_DIR;
else, S.out_root = '/scratch/qk75/LeonPhillips/PhD/LSA'; end
S.flow_prefix = 'su2L';                                 % su2L_flow.bin/.json and su2L_wall.bin/.json

% --- reference scales (assumption 3), all SI
S.D_ref   = 0.018;                                      % length scale: body diameter [m]
S.a_inf   = 226.865;                                    % free-stream sound speed [m/s]  (T = 128.067 K)
S.rho_inf = 0.51693;                                    % free-stream density [kg/m3]
S.U_inf   = 567.162;                                    % free-stream speed [m/s] (M = 2.5), used only for St_D
S.gamma   = 1.4;                                        % ratio of specific heats (air, ideal gas)
S.M_inf   = S.U_inf / S.a_inf;                          % = 2.5; converts St_D to omega (see section 6)

% --- geometry found from the exported wall (su2L_wall.bin), 2026-11-10
S.x_exit_mm = 1.876;                                    % nozzle exit plane [mm]
S.r_exit_mm = 2.389;                                    % nozzle exit (inner lip) radius [mm]

% --- stations and frequencies (assumptions 4 and 7)
S.x_mm  = 1.9 : 0.1 : 8.0;                              % axial stations [mm]
S.St    = 0.05 : 0.05 : 2.0;                            % Strouhal numbers St_D = f D / U_inf
S.m     = [0 1];                                        % azimuthal orders

% --- numerics (assumption 8)
S.N        = 250;                                       % Chebyshev polynomial order: N+1 collocation points
S.N_check  = 300;                                       % second resolution: modes that move are spurious
S.conv_tol = 1e-3;                                      % relative |alpha_N - alpha_Ncheck| / |alpha| to count as converged
S.dr_fine_mm  = 0.005;                                  % spacing of the fine uniform line used to sample the RANS data [mm]
S.smooth_mm   = 0.05;                                   % smoothing length l [mm] (about one RANS cell in the plume box)

% --- which base-flow profile the operator uses (user, 2026-11-10: Tran's tanh fit, kept as a switch)
%   'tanh'   : Tran et al. 2026 eqs. 3.1-3.3: two tanh shear layers per variable, fitted to the RANS profile
%              (removes the barrel-shock jumps inside the jet, which smoothing cannot remove without also
%              blurring the thin shear layer; the sensitivity run of 2026-11-10 showed that problem)
%   'smooth' : the Whittaker-smoothed RANS profile itself
if exist('PROFILE', 'var'), S.profile = PROFILE; else, S.profile = 'tanh'; end
S.fit_delta0_mm = [0.1 0.3 1.0];                        % starting half-thicknesses tried for each shear layer [mm]
S.fit_delta_min_mm = 0.05;                              % thinnest allowed layer [mm]: one RANS cell; thinner is a
                                                        % fit snapping onto a jump, not a resolved shear layer
S.rmax_gap_mm = 0.5;                                    % r_max = (bow-shock radius) - gap - 3*l   [mm]
S.clu_amp  = [30 8];                                    % mapping: extra point density at [jet, outer] shear layer
S.clu_wid_mm = [0.4 1.0];                               % mapping: width of each cluster [mm]

% --- which eigenvalues to keep and how to label them
S.alpha_max_fac = 5;                                    % keep |alpha_r| < fac*omega + 20 (removes the far, unresolved spectrum)
S.alpha_i_max   = 50;                                   % keep |alpha_i| < this (1/D units)
S.neutral_tol   = 0.02;                                 % |alpha_i| < tol*|alpha_r| counts as (nearly) neutral
S.sl_window_mm  = 0.6;                                  % a KH mode's |p| peak must lie this close to its shear layer [mm]

% --- 'single' mode
if exist('SINGLE_X', 'var'), S.single_x_mm = SINGLE_X; else, S.single_x_mm = 3.0; end    % station [mm]
if exist('SINGLE_ST', 'var'), S.single_St = SINGLE_ST; else, S.single_St = 1.271; end   % St_D (user, 2026-11-12)

% --- 'sensitivity' mode
S.sens_x_mm = [2.5 4.0 6.0 7.5];
S.sens_St   = [0.25 1.0 2.0];

% --- parallel workers (0 = serial)
if exist('NWORKERS', 'var'), S.nworkers = NWORKERS; else, S.nworkers = 0; end

% --- output folder for this mode, and a log file that records the settings
S.out_dir = fullfile(S.out_root, S.mode);               % e.g. .../LSA/scan
if ~exist(S.out_dir, 'dir'), mkdir(S.out_dir); end      % make it if needed
diary(fullfile(S.out_dir, 'lsa_local.log'));            % copy everything printed to a log file
fprintf('lsa_local.m  mode=%s  %s\n', S.mode, datestr(now)); %#ok<TNOW1,DATST>
disp(S);                                                % record the settings used

%% 2. Self-check against the exact vortex sheet (MODE='selfcheck') ----------------------------------------------
% Before trusting the operator on the real flow, check it on a case with a known answer. A round jet with a very
% thin tanh shear layer must approach the compressible vortex-sheet dispersion relation as the layer gets thinner.
if strcmp(S.mode, 'selfcheck')
    run_selfcheck(S);                                   % all the work is in the local function at the bottom
    diary off; return;                                  % nothing else to do in this mode
end

%% 3. Load the RANS base flow -------------------------------------------------------------------------------------
% The exporter wrote float64 columns, column-major, so MATLAB can read the whole block with one fread.
meta = jsondecode(fileread(fullfile(S.baseflow_dir, [S.flow_prefix '_flow.json'])));  % N (rows) and column names
fid  = fopen(fullfile(S.baseflow_dir, [S.flow_prefix '_flow.bin']), 'r');            % open the binary file
A    = fread(fid, [meta.N, numel(meta.names)], 'double');                             % N x K matrix of node data
fclose(fid);
col  = @(name) A(:, strcmp(meta.names, name));          % helper: pick a column by its name
fprintf('base flow: %d nodes from %s\n', meta.N, meta.source);

%% 4. Non-dimensionalise ------------------------------------------------------------------------------------------
% Lengths by D, velocity by a_inf, density by rho_inf, pressure by rho_inf*a_inf^2. Then the far free stream has
% nu = 1, p = 1/gamma, U = -M_inf (it flows in -x). The inviscid equations keep exactly the same form after this
% scaling, so no extra factors appear in the operator.
xn  = col('x') / S.D_ref;                               % x / D
rn  = col('r') / S.D_ref;                               % r / D
Un  = col('Velocity_x') / S.a_inf;                      % U / a_inf (local "Mach-like" velocity)
nun = S.rho_inf ./ col('Density');                      % nu * rho_inf = rho_inf / rho
pn  = col('Pressure') / (S.rho_inf * S.a_inf^2);        % p / (rho_inf a_inf^2)
Mn  = col('Mach');                                      % local Mach number (only used to find the bow shock)
clear A                                                 % the raw block is no longer needed (saves memory)

%% 5. Radial profiles at each station -----------------------------------------------------------------------------
% For each station x we (a) sample the unstructured RANS nodes on a fine uniform radial line, (b) find the two
% shear layers and the bow shock, (c) smooth the profiles and take their r-derivatives. This is done serially
% once, so the eigenvalue loop (section 6) only needs small arrays.
switch S.mode                                           % which stations does this mode need?
    case 'single',      xs_mm = S.single_x_mm;
    case 'sensitivity', xs_mm = S.sens_x_mm;
    case 'scan',        xs_mm = S.x_mm;
    case 'profiles',    xs_mm = S.x_mm;                 % profiles only: check the fits, no eigenvalues
    otherwise, error('unknown MODE %s', S.mode);
end
keep = xn > (min(xs_mm) - 0.5)*1e-3/S.D_ref & xn < (max(xs_mm) + 0.5)*1e-3/S.D_ref;  % only nodes near the stations
FU  = scatteredInterpolant(xn(keep), rn(keep), Un(keep),  'linear', 'none');  % linear interpolation on the nodes'
Fnu = scatteredInterpolant(xn(keep), rn(keep), nun(keep), 'linear', 'none');  % Delaunay triangulation; one
Fp  = scatteredInterpolant(xn(keep), rn(keep), pn(keep),  'linear', 'none');  % interpolant per variable
FM  = scatteredInterpolant(xn(keep), rn(keep), Mn(keep),  'linear', 'none');
clear xn rn Un nun pn Mn

prof = struct([]);                                      % one entry per station
for k = 1:numel(xs_mm)
    prof_k = radial_profile(xs_mm(k), FU, Fnu, Fp, FM, S);   % sample, find layers, smooth (local function)
    prof = [prof, prof_k]; %#ok<AGROW>
    fprintf('x = %.2f mm: jet layer r1 = %.3f mm, outer layer r2 = %.3f mm, bow shock %.2f mm, r_max = %.2f mm\n', ...
        xs_mm(k), prof_k.r1*S.D_ref*1e3, prof_k.r2*S.D_ref*1e3, prof_k.r_bs*S.D_ref*1e3, prof_k.rmax*S.D_ref*1e3);
    if strcmp(S.profile, 'tanh')                        % fit quality: Tran reports R^2 > 0.97 for every variable
        fprintf('            tanh fit R^2: U %.4f, rho %.4f, p %.4f\n', prof_k.fit.U.R2, prof_k.fit.rho.R2, prof_k.fit.p.R2);
    end
end
plot_profiles(prof, S);                                 % raw vs used profiles, with r1, r2 marked
if strcmp(S.mode, 'profiles')                           % this mode stops here
    save(fullfile(S.out_dir, 'lsa_profiles.mat'), 'prof', 'S');
    diary off; return;
end

%% 6. Solve the eigenvalue problems -------------------------------------------------------------------------------
% omega in our units: omega = 2*pi*f*D/a_inf = 2*pi*St_D*(U_inf/a_inf) = 2*pi*St_D*M_inf.
St2omega = @(St) 2*pi*St*S.M_inf;

switch S.mode
    case 'single'
        % One station and one St, both m: keep every eigenvalue and eigenvector so we can look at them.
        out = struct();
        for im = 1:numel(S.m)
            sol = solve_case(prof(1), S.m(im), St2omega(S.single_St), S.N, S, true);   % true = keep eigenvectors
            chk = solve_case(prof(1), S.m(im), St2omega(S.single_St), S.N_check, S, -1);   % -1: eigenvalues only
            sol = mark_converged(sol, chk, S);          % flag eigenvalues that also appear at N_check
            sol = label_modes(sol, prof(1), S);         % KH_jet / KH_outer / upstream_jet / other
            out.sol{im} = sol;
            print_modes(sol, S);
        end
        plot_single(out, prof(1), S);
        save(fullfile(S.out_dir, sprintf('x_%1.1f_St_%2.2f_lsa_single.mat', S.single_x_mm, S.single_St)), ...
             'out', 'prof', 'S', '-v7.3');                % file name carries x and St (user, 2026-11-12)

    case 'sensitivity'
        % Repeat a few cases with one numerical choice changed at a time, and tabulate the KH and upstream alphas.
        variants = { 'base',        struct();
                     'N150',        struct('N', 150, 'N_check', 200);
                     'N350',        struct('N', 350, 'N_check', 400);
                     'rmax-2.5mm',  struct('rmax_gap_mm', S.rmax_gap_mm + 2.5);
                     'smooth0.025', struct('smooth_mm', 0.025);
                     'smooth0.10',  struct('smooth_mm', 0.10) };
        if strcmp(S.profile, 'tanh'), variants = variants(1:4, :); end   % smoothing does not enter a tanh fit
        T = {};                                         % rows of the results table
        for iv = 1:size(variants, 1)
            Sv = S; f = fieldnames(variants{iv, 2});    % copy the settings and apply this variant's changes
            for j = 1:numel(f), Sv.(f{j}) = variants{iv, 2}.(f{j}); end
            for k = 1:numel(xs_mm)
                pk = prof(k);
                if ~strcmp(variants{iv,1}, 'base') && (Sv.smooth_mm ~= S.smooth_mm || Sv.rmax_gap_mm ~= S.rmax_gap_mm)
                    pk = radial_profile(xs_mm(k), FU, Fnu, Fp, FM, Sv);   % profile depends on these two settings
                end
                for im = 1:numel(S.m)
                    for is = 1:numel(S.sens_St)
                        sol = solve_case(pk, S.m(im), St2omega(S.sens_St(is)), Sv.N, Sv, false);
                        chk = solve_case(pk, S.m(im), St2omega(S.sens_St(is)), Sv.N_check, Sv, -1);   % -1: eigenvalues only
                        sol = label_modes(mark_converged(sol, chk, Sv), pk, Sv);
                        [kh, up] = pick_branches(sol);  % most unstable KH_jet, least damped upstream_jet
                        T(end+1, :) = {variants{iv,1}, xs_mm(k), S.m(im), S.sens_St(is), ...
                                       real(kh), imag(kh), real(up), imag(up)}; %#ok<SAGROW>
                    end
                end
            end
            fprintf('variant %s done\n', variants{iv,1});
        end
        T = cell2table(T, 'VariableNames', {'variant','x_mm','m','St','KH_ar','KH_ai','UP_ar','UP_ai'});
        writetable(T, fullfile(S.out_dir, 'lsa_sensitivity.csv'));
        print_sensitivity(T);
        save(fullfile(S.out_dir, 'lsa_sensitivity.mat'), 'T', 'prof', 'S');

    case 'scan'
        % All stations x St x m. Each (station, m, St) is independent, so loop stations in parallel.
        nx = numel(prof); nm = numel(S.m); nS = numel(S.St);
        R = cell(nx, 1);                                % results per station (parfor needs a sliced output)
        if S.nworkers > 0 && isempty(gcp('nocreate')), parpool(S.nworkers); end
        parfor (k = 1:nx, S.nworkers)                   % runs as a normal for-loop when nworkers = 0
            Rk = struct('alpha', {cell(nm, nS)}, 'label', {cell(nm, nS)}, 'cg', {cell(nm, nS)}, ...
                        'rpk', {cell(nm, nS)}, 'KH', nan(nm, nS), 'UP', nan(nm, nS));
            for im = 1:nm
                for is = 1:nS
                    om  = 2*pi*S.St(is)*S.M_inf;
                    sol = solve_case(prof(k), S.m(im), om, S.N, S, false);
                    chk = solve_case(prof(k), S.m(im), om, S.N_check, S, -1);   % -1: eigenvalues only
                    sol = label_modes(mark_converged(sol, chk, S), prof(k), S);
                    c = sol.conv;                       % store only converged eigenvalues
                    Rk.alpha{im,is} = sol.alpha(c);  Rk.label{im,is} = sol.label(c);
                    Rk.cg{im,is}    = sol.cg(c);     Rk.rpk{im,is}   = sol.rpk(c);
                    [Rk.KH(im,is), Rk.UP(im,is)] = pick_branches(sol);
                end
            end
            R{k} = Rk;
            fprintf('station x = %.2f mm done\n', prof(k).x_mm);
        end
        save(fullfile(S.out_dir, 'lsa_scan.mat'), 'R', 'prof', 'S', '-v7.3');
        plot_scan(R, prof, S);
end
fprintf('finished %s\n', datestr(now)); %#ok<TNOW1,DATST>
diary off;


%% ===============================================================================================================
%% Local functions
%% ===============================================================================================================

function P = radial_profile(x_mm, FU, Fnu, Fp, FM, S)
% Sample the mean flow along r at station x, find the shear layers and the bow shock, smooth and differentiate.
    L  = S.D_ref*1e3;                                   % D in mm (converts mm <-> 1/D units)
    xq = x_mm / L;                                      % station in units of D
    dr = S.dr_fine_mm / L;                              % fine spacing in units of D
    r  = (0:dr:28/L)';                                  % fine line out to 28 mm (the export box ends at 30 mm)
    X  = xq*ones(size(r));                              % same x for every point on the line
    U = FU(X, r); nu = Fnu(X, r); p = Fp(X, r); M = FM(X, r);   % linear interpolation of the RANS nodes

    % --- jet shear layer r1: where U first falls below half of its jet maximum, moving outwards from the peak
    inner = r < 10/L;                                   % the jet is well inside 10 mm at every station
    [Umax, iU] = max(U .* inner);                       % jet maximum and where it is
    i1 = iU - 1 + find(U(iU:end) < 0.5*Umax, 1);        % first point beyond the peak below half of it
    r1 = r(i1);

    % --- bow shock: first radius beyond the jet where the flow is back to free stream (U = -M_inf, Mach = 2.5)
    ibs = find(r > r1 + 1/L & U < -0.98*S.M_inf & M > 0.98*S.M_inf, 1);
    r_bs = r(ibs);

    % --- r_max: stay inside the bow shock, by a gap plus three smoothing lengths (so the smoothing does not
    %     pull the shock jump into the domain)
    rmax = r_bs - (S.rmax_gap_mm + 3*S.smooth_mm)/L;

    % --- cut the line at r_max and smooth each variable (penalised least squares, see whittaker below)
    in  = r <= rmax;  r = r(in);
    lam = (S.smooth_mm / S.dr_fine_mm)^4;               % smoothing weight for a cut-off length of smooth_mm
    Us = whittaker(U(in), lam); nus = whittaker(nu(in), lam); ps = whittaker(p(in), lam);

    % --- derivatives of the smoothed profiles (second-order central differences on the fine line)
    dUs = gradient(Us, r); dnus = gradient(nus, r); dps = gradient(ps, r);

    % --- outer shear layer r2: steepest |dU/dr| between the jet layer and r_max
    win = r > r1 + 0.5/L & r < rmax - 0.5/L;            % search window (avoid the jet layer and the edge)
    [~, j2] = max(abs(dUs) .* win);
    r2 = r(j2);

    P = struct('x_mm', x_mm, 'r', r, 'U', Us, 'nu', nus, 'p', ps, 'dU', dUs, 'dnu', dnus, 'dp', dps, ...
               'U_raw', U(in), 'nu_raw', nu(in), 'p_raw', p(in), 'r1', r1, 'r2', r2, 'r_bs', r_bs, 'rmax', rmax);

    % --- Tran's analytic profiles (S.profile = 'tanh'): replace the smoothed profiles by two-tanh fits
    if strcmp(S.profile, 'tanh')
        d0 = S.fit_delta0_mm / L;                       % starting half-thicknesses in units of D
        dmin = S.fit_delta_min_mm / L;                  % minimum half-thickness in units of D
        fU = fit_tanh2(r, U(in),       r1, r2, d0, dmin);     % U(r)          (Tran eq. 3.1)
        fR = fit_tanh2(r, 1./nu(in),   r1, r2, d0, dmin);     % rho(r) = 1/nu (Tran eq. 3.2 fits 1/nu, i.e. density)
        fP = fit_tanh2(r, p(in),       r1, r2, d0, dmin);     % p(r)          (Tran eq. 3.3)
        [P.U, P.dU] = eval_tanh2(fU, r);                % analytic values and exact r-derivatives
        [rho, drho] = eval_tanh2(fR, r);
        P.nu = 1 ./ rho; P.dnu = -drho ./ rho.^2;       % nu = 1/rho, so dnu/dr = -(drho/dr)/rho^2
        [P.p, P.dp] = eval_tanh2(fP, r);
        P.fit = struct('U', fU, 'rho', fR, 'p', fP);
        P.r1 = fU.R(1); P.r2 = fU.R(2);                 % shear-layer radii of the U fit drive the grid mapping
    end
end

function F = fit_tanh2(r, y, Rg1, Rg2, d0, dmin)
% Fit Tran et al. (2026) eqs. 3.1-3.3, two tanh shear layers, to a profile y(r):
%   y(r) = A1*T(r; R1, b1) + A2*T(r; R2, b2) + c,   T(r; R, b) = 0.5*(1 + tanh(b*(R/r - r/R)))
% T is 1 on the axis, 1/2 at r = R and 0 far out, so A1 + A2 + c is the core value, A2 + c the value between the
% layers and c the value outside. Tran writes A1 = (1-h)*w1 and A2 = h*w2 (h: inner/outer ratio, w: weight) and
% fixes c at the ambient value (0 for U, 1 for 1/nu and p); here c is FITTED because outside the jet the srp
% flow is the moving, compressed shock layer, not still ambient air. (Only A1, A2 are identifiable, not h and w.)
% The model is linear in (A1, A2, c), so those come from linear least squares for each trial (R1, b1, R2, b2),
% and fminsearch (base MATLAB) searches only the four layer parameters ("variable projection"). They are searched
% as q = log([R1, delta1 - dmin, R2, delta2 - dmin]), so radii stay positive and no layer gets thinner than dmin.
% A half-thickness delta corresponds to b = R/(2*delta): near r = R, b*(R/r - r/R) ~ -2*b*(r - R)/R.
    ok = isfinite(y); r = r(ok); y = y(ok);
    best = struct('cost', Inf);
    for da = d0                                         % a few starting thicknesses (multi-start)
        for db = d0
            q0 = log([Rg1, da - dmin, Rg2, db - dmin]);
            opt = optimset('MaxFunEvals', 4000, 'MaxIter', 4000, 'TolX', 1e-7, 'TolFun', 1e-12, 'Display', 'off');
            q = fminsearch(@(q) tanh2_cost(q, r, y, dmin), q0, opt);
            c = tanh2_cost(q, r, y, dmin);
            if c < best.cost, best = struct('cost', c, 'q', q); end
        end
    end
    [~, coef, R, b] = tanh2_cost(best.q, r, y, dmin);
    F = struct('R', R, 'b', b, 'delta', R./(2*b), 'A', coef(1:2)', 'c', coef(3));   % delta: half-thicknesses
    yf = eval_tanh2(F, r);
    F.R2 = 1 - sum((y - yf).^2) / sum((y - mean(y)).^2);   % coefficient of determination
end

function [cost, coef, R, b] = tanh2_cost(q, r, y, dmin)
% Sum of squared residuals of the best linear (A1, A2, c) for q = log([R1, delta1 - dmin, R2, delta2 - dmin]).
    v = exp(q); dl = dmin + v([2 4]);                   % half-thicknesses, never below dmin
    bb = v([1 3]) ./ (2*dl);                            % b = R/(2*delta)
    [R, k] = sort(v([1 3])); b = bb(k);                 % keep R1 < R2 (inner layer first)
    B = [tanhT(r, R(1), b(1)), tanhT(r, R(2), b(2)), ones(size(r))];
    coef = B \ y;                                       % linear least squares for A1, A2, c
    cost = sum((B*coef - y).^2);
    if R(2) > r(end) || R(1) < 1e-4 || R(2) - R(1) < 1e-3, cost = cost * 10 + 1; end   % keep layers in the domain
end

function [y, dy] = eval_tanh2(F, r)
% Value and exact derivative of a fitted two-tanh profile.
    [T1, dT1] = tanhT(r, F.R(1), F.b(1));
    [T2, dT2] = tanhT(r, F.R(2), F.b(2));
    y  = F.A(1)*T1 + F.A(2)*T2 + F.c;
    dy = F.A(1)*dT1 + F.A(2)*dT2;
end

function [T, dT] = tanhT(r, R, b)
% T(r) = 0.5*(1 + tanh(b*(R/r - r/R))) and dT/dr = 0.5*sech^2(.)*b*(-R/r^2 - 1/R). At r = 0 the argument is
% +infinity, T = 1 and dT = 0 (sech^2 underflows to exactly 0 before R/r^2 overflows, so no 0*Inf occurs).
    rr = max(r, eps);
    t  = tanh(b*(R./rr - rr/R));
    T  = 0.5*(1 + t);
    dT = 0.5*(1 - t.^2) .* b .* (-R./rr.^2 - 1/R);
end

function z = whittaker(y, lam)
% Whittaker smoother: z minimises  sum (z - y)^2 + lam * sum (second difference of z)^2.
% Its transfer function is 1/(1 + lam*(k*dr)^4), so it removes wiggles shorter than about 2*pi*lam^(1/4)*dr
% while leaving smooth trends alone. It needs no toolbox and has no end effects of a moving average.
    n  = numel(y);
    D2 = diff(speye(n), 2);                             % (n-2) x n second-difference matrix
    z  = (speye(n) + lam*(D2'*D2)) \ y;                 % sparse solve of the normal equations
end

function [x, D] = cheb(N)
% Chebyshev points and differentiation matrix (Trefethen 2000, Spectral Methods in MATLAB, p. 54).
% Returns x(j) = cos(pi*j/N), j = 0..N (from +1 down to -1), and D with (D*f)(j) = f'(x_j).
    x = cos(pi*(0:N)'/N);
    c = [2; ones(N-1, 1); 2] .* (-1).^(0:N)';
    X = repmat(x, 1, N+1);
    dX = X - X';
    D = (c*(1./c)') ./ (dX + eye(N+1));                 % off-diagonal entries
    D = D - diag(sum(D, 2));                            % diagonal: rows of a differentiation matrix sum to zero
end

function G = make_grid(P, N, S)
% Collocation points r_j in [0, r_max] and the radial derivative matrix Dr.
% Mapping: choose a point density g(r) = 1 + A1*exp(-((r-r1)/w1)^2) + A2*exp(-((r-r2)/w2)^2) (high in the two
% shear layers). Let s(r) = int_0^r g / int_0^rmax g, which runs from 0 to 1. Put Chebyshev points uniformly in
% the Chebyshev sense in s and map back to r. Then dr/ds = Gtot/g(r) exactly, and by the chain rule
% d/dr = (ds/dr) d/ds = (g/Gtot) d/ds.
    L = S.D_ref*1e3;
    rf = linspace(0, P.rmax, 20001)';                   % very fine line for the integral
    g  = @(r) 1 + S.clu_amp(1)*exp(-((r - P.r1)/(S.clu_wid_mm(1)/L)).^2) ...
                + S.clu_amp(2)*exp(-((r - P.r2)/(S.clu_wid_mm(2)/L)).^2);
    Gc = cumtrapz(rf, g(rf));                           % running integral of g
    Gtot = Gc(end);
    [xi, Dxi] = cheb(N);                                % xi from +1 to -1
    xi = flipud(xi); Dxi = rot90(Dxi, 2);               % reorder so index 1 is the axis (xi = -1, r = 0)
    s  = (xi + 1)/2;                                    % s in [0, 1]
    r  = interp1(Gc/Gtot, rf, s);                       % invert s(r) to get the collocation radii
    r(1) = 0; r(end) = P.rmax;                          % make the ends exact
    Ds = 2*Dxi;                                         % d/ds = (dxi/ds) d/dxi = 2 d/dxi
    Dr = diag(g(r)/Gtot) * Ds;                          % d/dr = (ds/dr) d/ds
    G = struct('r', r, 'Dr', Dr);
end

function sol = solve_case(P, m, omega, N, S, keepvec)
% keepvec: true = also return eigenvectors V, false = labels only, -1 = eigenvalues only (convergence check)
% Build L and R at this station, apply boundary conditions and solve L q = alpha R q.
    G  = make_grid(P, N, S);
    r  = G.r; n = numel(r);
    % mean flow and its derivatives on the collocation points (pchip: shape-preserving, no overshoot)
    Ub  = interp1(P.r, P.U,  r, 'pchip');  dUb  = interp1(P.r, P.dU,  r, 'pchip');
    nub = interp1(P.r, P.nu, r, 'pchip');  dnub = interp1(P.r, P.dnu, r, 'pchip');
    pb  = interp1(P.r, P.p,  r, 'pchip');  dpb  = interp1(P.r, P.dp,  r, 'pchip');
    [L, R] = build_LR(r, G.Dr, Ub, dUb, nub, dnub, pb, dpb, m, omega, S.gamma);
    [L, R, bc] = apply_bcs(L, R, G.Dr, n, m);
    sol = eig_solve(L, R, bc, omega, keepvec, S);
    sol.r = r; sol.m = m; sol.omega = omega; sol.N = N; sol.Ub = Ub;   % Ub kept for the continuum test
end

function [L, R] = build_LR(r, Dr, Ub, dUb, nub, dnub, pb, dpb, m, omega, gam)
% The linearised Euler operator for q = [nu; u; v; w; p] (each a block of n values on the r_j).
% Derivation (inserting q = qbar(r) + q_hat(r) exp(i(alpha x + m theta - omega t)), keeping first order, mean
% V = W = 0, d/dx -> i*alpha, d/dtheta -> i*m, d/dt -> -i*omega, and putting every alpha term on the right):
%   specific volume  D(nu)/Dt = nu div(u):
%       -i w nu + dnub v - nub v' - nub v/r - i m nub w/r      = alpha ( -i Ub nu + i nub u )    (Tran D1-D3)
%   axial momentum   Du/Dt = -nu dp/dx:
%       -i w u + dUb v                                         = alpha ( -i Ub u  - i nub p )    (Tran D4-D5)
%   radial momentum  Dv/Dt = -nu dp/dr:
%       -i w v + dpb nu + nub p'                               = alpha ( -i Ub v )               (Tran D6-D7)
%   azimuthal mom.   Dw/Dt = -nu (1/r) dp/dtheta:
%       -i w w + i m nub p / r                                 = alpha ( -i Ub w )               (Tran D8)
%   energy (as p)    Dp/Dt = -gamma p div(u):
%       -i w p + (dpb + gam pb/r) v + gam pb v' + i m gam pb w/r = alpha ( -i Ub p - i gam pb u ) (Tran D9-D11)
% (here w on the left means omega and the w in "w/r" is the azimuthal velocity; v' = dv/dr.)
    n  = numel(r);
    I  = eye(n);  Z = zeros(n);
    ir = [0; 1./r(2:end)];                              % 1/r, with 0 on the axis: that row is replaced by a BC
    d  = @(f) diag(f);                                  % shorthand for a diagonal matrix
    iw = -1i*omega*I;                                   % the -i*omega term, the same in every equation
    % --- left-hand operator L (rows: equations; columns: nu, u, v, w, p)
    L = [ iw,     Z,      d(dnub) - d(nub)*Dr - d(nub.*ir),          -1i*m*d(nub.*ir),       Z;
          Z,      iw,     d(dUb),                                    Z,                      Z;
          d(dpb), Z,      iw,                                        Z,                      d(nub)*Dr;
          Z,      Z,      Z,                                         iw,                     1i*m*d(nub.*ir);
          Z,      Z,      d(dpb) + gam*d(pb.*ir) + gam*d(pb)*Dr,     1i*m*gam*d(pb.*ir),     iw ];
    % --- right-hand operator R (multiplies alpha)
    R = [ -1i*d(Ub), 1i*d(nub),       Z,          Z,          Z;
          Z,         -1i*d(Ub),       Z,          Z,          -1i*d(nub);
          Z,         Z,               -1i*d(Ub),  Z,          Z;
          Z,         Z,               Z,          -1i*d(Ub),  Z;
          Z,         -1i*gam*d(pb),   Z,          Z,          -1i*d(Ub) ];
end

function [L, R, bc] = apply_bcs(L, R, Dr, n, m)
% Replace the equation rows at the axis (j = 1) and at r_max (j = n) by boundary conditions. Setting the R rows to
% zero there turns those rows into constraints; they give infinite eigenvalues, which eig_solve removes.
% bc returns which rows are boundary rows (needed for the group velocity in eig_solve).
    blk = @(k) (k-1)*n + (1:n);                         % columns of variable k (1 nu, 2 u, 3 v, 4 w, 5 p)
    row = @(k, j) (k-1)*n + j;                          % the row of equation k at point j
    e1  = zeros(1, n); e1(1) = 1;                       % row vector that picks the value at the axis
    en  = zeros(1, n); en(n) = 1;                       % row vector that picks the value at r_max
    bc  = false(5*n, 1);
    % --- axis (Lesshafft & Huerre 2007; Tran eqs. 2.5-2.9): the perturbation must be a smooth Cartesian field
    for k = 1:5, L(row(k,1), :) = 0; R(row(k,1), :) = 0; bc(row(k,1)) = true; end   % clear the 5 axis rows
    if m == 0                                           % axisymmetric: no cross-flow, flat profiles at r = 0
        L(row(1,1), blk(1)) = Dr(1,:);                  % d(nu)/dr = 0
        L(row(2,1), blk(2)) = Dr(1,:);                  % du/dr = 0
        L(row(3,1), blk(3)) = e1;                       % v = 0
        L(row(4,1), blk(4)) = e1;                       % w = 0
        L(row(5,1), blk(5)) = Dr(1,:);                  % dp/dr = 0
    elseif abs(m) == 1                                  % |m| = 1: scalars vanish, a uniform cross-flow is allowed
        L(row(1,1), blk(1)) = e1;                       % nu = 0
        L(row(2,1), blk(2)) = e1;                       % u = 0
        L(row(3,1), blk(3)) = Dr(1,:);                  % dv/dr = 0
        L(row(4,1), blk(3)) = e1;                       % v + i*sign(m)*w = 0: a cross-flow U_y gives
        L(row(4,1), blk(4)) = 1i*sign(m)*e1;            %   v = U_y cos(theta), w = -U_y sin(theta)
        L(row(5,1), blk(5)) = e1;                       % p = 0
    else                                                % |m| >= 2: everything zero at the axis
        for k = 1:5, L(row(k,1), blk(k)) = e1; end
    end
    % --- r_max: Dirichlet, all perturbations zero (Tran eq. 2.10)
    for k = 1:5
        L(row(k,n), :) = 0; R(row(k,n), :) = 0; bc(row(k,n)) = true;
        L(row(k,n), blk(k)) = en;
    end
end

function sol = eig_solve(L, R, bc, omega, keepvec, S)
% Solve L q = alpha R q (dense QZ). Also return left eigenvectors, used for the group velocity.
% keepvec = -1 means "eigenvalues only" (the N_check solve: it only has to say where the eigenvalues are, and
% skipping the eigenvectors makes it several times cheaper).
    inwin = @(a) isfinite(a) & abs(real(a)) < S.alpha_max_fac*omega + 20 & abs(imag(a)) < S.alpha_i_max;
    if isequal(keepvec, -1)
        a = eig(L, R); sol = struct('alpha', a(inwin(a))); return;
    end
    [V, Dg, W] = eig(L, R);                             % right vectors V, eigenvalues Dg, left vectors W
    a = diag(Dg);
    ok = inwin(a);                                      % finite, and inside the window we can resolve
    a = a(ok); V = V(:, ok); W = W(:, ok);
    % --- group velocity from first-order perturbation theory: dL/domega = -i*E, E = identity on equation rows
    %     d(alpha)/d(omega) = y' (dL/domega) q / (y' R q), with y the left and q the right eigenvector
    n5 = size(L, 1);
    E = double(~bc);                                    % 1 on equation rows, 0 on boundary rows (no omega there)
    dadw = zeros(size(a));
    for j = 1:numel(a)
        q = V(:, j); y = W(:, j);
        dadw(j) = (y' * (-1i*E.*q)) / (y' * (R*q));
    end
    cg = 1 ./ real(dadw);                               % group velocity d(omega)/d(alpha_r), units a_inf
    n = n5/5;
    pabs = abs(V(4*n+1:5*n, :));                        % |p_hat| on the r_j (labelling, plots)
    uabs = abs(V(  n+1:2*n, :));                        % |u_hat| (axial velocity): a KH wave peaks at its shear layer
    sol = struct('alpha', a, 'cg', cg, 'pabs', pabs, 'uabs', uabs);
    if keepvec, sol.V = V; end
end

function sol = mark_converged(sol, chk, S)
% An eigenvalue counts as converged if (1) the N_check solution has one within conv_tol (relative), and
% (2) that match is much closer (10x) than the nearest OTHER eigenvalue of the same solve. Test (2) matters:
% a discretised continuous spectrum is a dense cloud of eigenvalues, so test (1) alone is always passed by some
% neighbour in the cloud even though no individual eigenvalue is converged (found 2026-11-10 on x = 3 mm).
    sol.conv = false(size(sol.alpha));
    for j = 1:numel(sol.alpha)
        dmin  = min(abs(chk.alpha - sol.alpha(j)));                      % distance to the N_check spectrum
        dself = min(abs(sol.alpha([1:j-1, j+1:end]) - sol.alpha(j)));    % distance to its own nearest neighbour
        if isempty(dself), dself = Inf; end
        sol.conv(j) = ~isempty(dmin) && dmin < S.conv_tol*max(1, abs(sol.alpha(j))) && dmin < 0.1*dself;
    end
end

function sol = label_modes(sol, P, S)
% Label each eigenvalue by its physics (heuristics, not a Briggs-Bers analysis; see README):
%   KH_jet       : growing (alpha_i < 0), downstream group velocity, |u| peak in the jet shear layer
%   KH_outer     : the same, peak in the outer shear layer
%   upstream_jet : nearly neutral, upstream group velocity (cg < 0), most |p|^2 energy inside the jet
%   continuum    : alpha = omega/U(r_j) at some collocation point: a disturbance simply carried by the mean flow
%                  (entropy/vorticity); the inviscid problem has a continuous spectrum there, discretised as a
%                  dense set of nearly neutral eigenvalues. Not a wave; never counted as KH or upstream.
%   other        : everything else (acoustic/duct-like waves, ...)
    L = S.D_ref*1e3; w = S.sl_window_mm/L;
    r = sol.r; na = numel(sol.alpha);
    sol.label = repmat({'other'}, na, 1); sol.rpk = nan(na, 1); sol.rupk = nan(na, 1);
    for j = 1:na
        pa = sol.pabs(:, j); [~, ip] = max(pa); sol.rpk(j) = r(ip);   % where |p| peaks
        [~, iu] = max(sol.uabs(:, j)); sol.rupk(j) = r(iu);                % where |u| peaks
        a = sol.alpha(j); cg = sol.cg(j);
        cph = sol.omega/real(a);                                 % phase speed (units a_inf)
        if min(abs(a*sol.Ub - sol.omega)) < 1e-3*sol.omega        % (i) omega - alpha*U = 0 at some r_j, or
            sol.label{j} = 'continuum'; continue;
        end
        if abs(imag(a)) < S.neutral_tol*abs(real(a)) && abs(cg - cph) < 0.02*abs(cph) ...
                && cph >= min(sol.Ub) && cph <= max(sol.Ub)       % (ii) neutral, non-dispersive (cg = c) and
            sol.label{j} = 'continuum'; continue;                 %      moving at a speed the mean flow has
        end
        grows = imag(a) < -1e-3*abs(a);                           % clearly growing, not a neutral mode at -1e-13
        Ein = trapz(r(r <= P.r1), pa(r <= P.r1).^2 .* r(r <= P.r1));   % |p|^2 r dr inside the jet
        Eall = trapz(r, pa.^2 .* r);
        if grows && cg > 0 && abs(r(iu) - P.r1) < w
            sol.label{j} = 'KH_jet';
        elseif grows && cg > 0 && abs(r(iu) - P.r2) < 2*w
            sol.label{j} = 'KH_outer';
        elseif abs(imag(a)) < S.neutral_tol*abs(real(a)) && cg < 0 && Ein/Eall > 0.5
            sol.label{j} = 'upstream_jet';
        end
    end
end

function [kh, up] = pick_branches(sol)
% Most unstable converged KH_jet mode, and least damped converged upstream_jet mode (NaN if none).
    c = sol.conv;
    i = find(c & strcmp(sol.label, 'KH_jet'));
    if isempty(i), kh = NaN; else, [~, j] = min(imag(sol.alpha(i))); kh = sol.alpha(i(j)); end
    i = find(c & strcmp(sol.label, 'upstream_jet'));
    if isempty(i), up = NaN; else, [~, j] = min(abs(imag(sol.alpha(i)))); up = sol.alpha(i(j)); end
end

function print_modes(sol, S)
    fprintf('\nm = %d, omega = %.3f (St_D = %.3f): %d eigenvalues kept, %d converged\n', sol.m, sol.omega, ...
        sol.omega/(2*pi*S.M_inf), numel(sol.alpha), sum(sol.conv));
    fprintf('  of the converged: %d continuum, %d other\n', sum(sol.conv & strcmp(sol.label, 'continuum')), sum(sol.conv & strcmp(sol.label, 'other')));
    i = find(sol.conv & ~ismember(sol.label, {'other', 'continuum'}));   % print only the labelled waves
    for j = i'
        fprintf('  %-13s alpha = %8.3f %+8.3fi   cg = %7.3f a_inf   |p| peak r = %.3f mm\n', sol.label{j}, ...
            real(sol.alpha(j)), imag(sol.alpha(j)), sol.cg(j), sol.rpk(j)*S.D_ref*1e3);
    end
end

function print_sensitivity(T)
    fprintf('\n%-12s %5s %2s %5s %18s %18s\n', 'variant', 'x', 'm', 'St', 'KH alpha', 'upstream alpha');
    for j = 1:height(T)
        fprintf('%-12s %5.2f %2d %5.2f %8.3f%+8.3fi %8.3f%+8.3fi\n', T.variant{j}, T.x_mm(j), T.m(j), T.St(j), ...
            T.KH_ar(j), T.KH_ai(j), T.UP_ar(j), T.UP_ai(j));
    end
end

% ---------------------------------------------------------------------------------------------- self-check ------
function run_selfcheck(S)
% Round jet, D = 1 (R = 0.5), velocity in a_inf, ambient at rest: U = Uj*f(r), nu = 1 + (nuj - 1)*f(r),
% p = 1/gamma, with Michalke's profile f = 0.5*(1 + tanh(0.25*(R/theta)*(R/r - r/R))).
% As R/theta grows, the KH alpha must approach the vortex-sheet root. The finite-thickness correction is first
% order in alpha*theta, so the error should roughly halve each time theta halves, and the Richardson extrapolation
% 2*alpha(R/theta = 80) - alpha(R/theta = 40) (an estimate of the zero-thickness limit) should sit on the root.
    gam = S.gamma; R0 = 0.5; rmax = 10;                 % domain 0..20 jet radii
    Sx = S; Sx.D_ref = 1e-3;                            % so the "mm" mapping settings below are in units of D
    Sx.clu_amp = [60 0]; Sx.clu_wid_mm = [0.05 1];      % one tight cluster at r = R (the layer is only R/80 thick)
    cases = struct('Mj', {1.5, 1.5}, 'nuj', {1.0, 1.5});% isothermal and heated (nu ratio = T ratio)
    Stj = [0.2 0.4 0.6]; thetas = [80 40 20]; ms = [0 1];   % thinnest layer first: it gives the Newton start
    T = {}; X = {};
    for c = 1:numel(cases)
        Uj = cases(c).Mj*sqrt(cases(c).nuj);            % U_j/a_inf = M_j * a_j/a_inf, a_j^2/a_inf^2 = nu_j
        for m = ms
            for St = Stj
                om = 2*pi*St*Uj;                        % St_j = f D / U_j, D = 1
                aVS = NaN; aKH = nan(size(thetas));     % vortex-sheet root and the LSA KH value per R/theta
                for it = 1:numel(thetas)
                    P = selfcheck_profile(thetas(it), R0, Uj, cases(c).nuj, gam, rmax);
                    sol = solve_case(P, m, om, S.N, Sx, false);
                    chk = solve_case(P, m, om, S.N_check, Sx, -1);   % -1: eigenvalues only
                    sol = mark_converged(sol, chk, S);
                    % KH candidates: converged, clearly growing (not a neutral mode with alpha_i ~ -1e-13),
                    % downstream group velocity, |u| (axial velocity) peak near r = R (|p| need not peak there)
                    i = find(sol.conv & imag(sol.alpha) < -1e-3*abs(sol.alpha) & sol.cg > 0);
                    pk = arrayfun(@(j) sol.r(find(sol.uabs(:,j) == max(sol.uabs(:,j)), 1)), i);
                    i = i(abs(pk - R0) < 0.1);
                    if isempty(i), continue; end
                    if isnan(aVS)                       % thinnest layer: Newton for the root from its most
                        [~, j] = min(imag(sol.alpha(i)));   % unstable KH candidate
                        aVS = vortex_sheet_root(sol.alpha(i(j)), om, m, Uj, cases(c).nuj, R0);
                    end
                    [~, j] = min(abs(sol.alpha(i) - aVS));  % the candidate on the branch that ends at the root
                    aKH(it) = sol.alpha(i(j));
                end
                for it = 1:numel(thetas)
                    T(end+1, :) = {cases(c).Mj, cases(c).nuj, m, St, thetas(it), real(aKH(it)), imag(aKH(it)), ...
                                   real(aVS), imag(aVS), abs(aKH(it) - aVS)/abs(aVS)}; %#ok<AGROW>
                end
                aRE = 2*aKH(1) - aKH(2);                % Richardson: zero-thickness estimate from R/theta 80, 40
                X(end+1, :) = {cases(c).nuj, m, St, real(aRE), imag(aRE), abs(aRE - aVS)/abs(aVS), ...
                               abs(aKH(2) - aVS)/abs(aKH(1) - aVS)}; %#ok<AGROW>
            end
        end
    end
    T = cell2table(T, 'VariableNames', {'Mj','nuj','m','St_j','R_theta','LSA_ar','LSA_ai','VS_ar','VS_ai','rel_err'});
    X = cell2table(X, 'VariableNames', {'nuj','m','St_j','extrap_ar','extrap_ai','extrap_rel_err','err40_over_err80'});
    disp(T); disp(X);
    writetable(T, fullfile(S.out_dir, 'selfcheck.csv'));
    writetable(X, fullfile(S.out_dir, 'selfcheck_extrapolated.csv'));
    % pass (revised 2026-11-10 after the first runs): (1) a KH mode found in every case; (2) the error falls at
    % first order, err(40)/err(80) in 1.8-2.5; (3) the two-point extrapolation within 10 %. What is left after a
    % two-point extrapolation is O((alpha theta)^2), so it grows with St: about 2 % at St_j 0.4 and 5-7 % at 0.6.
    ok = all(isfinite(X.extrap_rel_err)) && all(X.err40_over_err80 > 1.8 & X.err40_over_err80 < 2.5) ...
         && all(X.extrap_rel_err < 0.10);
    fprintf('\nSELF-CHECK: max rel. error of the extrapolated alpha: %.3g; err(40)/err(80) %.2f-%.2f (2 = first order)\n', ...
        max(X.extrap_rel_err), min(X.err40_over_err80), max(X.err40_over_err80));
    if ok, fprintf('SELF-CHECK PASSED\n'); else, fprintf('SELF-CHECK FAILED\n'); end
end

function P = selfcheck_profile(Rt, R0, Uj, nuj, gam, rmax)
% Analytic profile on a fine line, in the same struct format as radial_profile (D_ref is set to 1 mm by the
% caller so the mm-based mapping widths below are in units of D).
    r = linspace(0, rmax, 40001)';
    b = 0.25*Rt;                                        % Michalke's profile 2 gradient parameter
    arg = b*(R0./r - r/R0); arg(1) = 50;                % r = 0: tanh -> 1
    f = 0.5*(1 + tanh(arg));
    df = 0.5*(1 - tanh(arg).^2) .* b .* (-R0./r.^2 - 1/R0); df(1) = 0;
    P = struct('x_mm', 0, 'r', r, 'U', Uj*f, 'nu', 1 + (nuj - 1)*f, 'p', ones(size(r))/gam, ...
               'dU', Uj*df, 'dnu', (nuj - 1)*df, 'dp', zeros(size(r)), 'r1', R0, 'r2', 3*R0, ...
               'r_bs', rmax, 'rmax', rmax);
end

function a = vortex_sheet_root(a, om, m, Uj, nuj, R0)
% Compressible vortex sheet at r = R0 between a jet (U = Uj, nu = nuj) and still ambient (nu = 1).
% Pressure: p = A I_m(lj r) inside, B K_m(lo r) outside, with lj^2 = a^2 - (om - a Uj)^2/nuj (a_j^2 = nuj in our
% units), lo^2 = a^2 - om^2. Matching pressure and particle displacement eta = nu p'/(om - a U)^2 gives
%   nuj lj I_m'(lj R)/((om - a Uj)^2 I_m(lj R)) - lo K_m'(lo R)/(om^2 K_m(lo R)) = 0.
    Ip = @(z) 0.5*(besseli(m-1, z) + besseli(m+1, z));
    Kp = @(z) -0.5*(besselk(m-1, z) + besselk(m+1, z));
    F = @(a) nuj*sqrt(a^2 - (om - a*Uj)^2/nuj)*Ip(sqrt(a^2 - (om - a*Uj)^2/nuj)*R0) ...
             / ((om - a*Uj)^2*besseli(m, sqrt(a^2 - (om - a*Uj)^2/nuj)*R0)) ...
           - sqrt(a^2 - om^2)*Kp(sqrt(a^2 - om^2)*R0)/(om^2*besselk(m, sqrt(a^2 - om^2)*R0));
    for it = 1:50                                       % secant/Newton with a numerical derivative
        h = 1e-6*max(1, abs(a));
        da = -F(a)/((F(a + h) - F(a - h))/(2*h));
        a = a + da;
        if abs(da) < 1e-12*max(1, abs(a)), break; end
    end
end

% ---------------------------------------------------------------------------------------------- plots -----------
function plot_profiles(prof, S)
    L = S.D_ref*1e3; f = figure('Visible', 'off', 'Position', [0 0 1500 900]);
    names = {'U', 'nu', 'p'}; ylab = {'U / a_\infty', '\nu \rho_\infty', 'p / \rho_\infty a_\infty^2'};
    for v = 1:3
        subplot(2, 3, v); hold on;
        for k = 1:numel(prof)
            plot(prof(k).r*L, prof(k).([names{v} '_raw']), ':', 'Color', [.6 .6 .6]);
            plot(prof(k).r*L, prof(k).(names{v}));
        end
        xlabel('r [mm]'); ylabel(ylab{v}); title(sprintf('%s: RANS (grey dotted) and used (%s)', names{v}, S.profile));
        subplot(2, 3, v + 3); hold on;
        for k = 1:numel(prof), plot(prof(k).r*L, prof(k).(['d' names{v}])); end
        xlabel('r [mm]'); ylabel(['d' names{v} '/dr (1/D units)']);
    end
    sgtitle(sprintf('srp\\_M2.5 LSA base-flow profiles, %d station(s), profile = %s', numel(prof), S.profile));
    exportgraphics(f, fullfile(S.out_dir, 'lsa_profiles.png'), 'Resolution', 130); close(f);
end

function plot_single(out, P, S)
    L = S.D_ref*1e3; f = figure('Visible', 'off', 'Position', [0 0 1500 900]);
    cols = struct('KH_jet', [0.85 0.1 0.1], 'KH_outer', [0.9 0.6 0], 'upstream_jet', [0 0.4 0.9], ...
                  'continuum', [0.3 0.7 0.3], 'other', [.4 .4 .4]);
    for im = 1:numel(out.sol)
        s = out.sol{im};
        subplot(2, 2, im); hold on;
        % every series carries its own DisplayName, so the legend cannot get out of step with the plotted data
        plot(real(s.alpha(~s.conv)), imag(s.alpha(~s.conv)), '.', 'Color', [.8 .8 .8], 'DisplayName', 'not converged');
        for lab = fieldnames(cols)'
            i = s.conv & strcmp(s.label, lab{1});
            if ~any(i), continue; end
            mk = 'o'; if strcmp(lab{1}, 'continuum'), mk = '.'; end
            plot(real(s.alpha(i)), imag(s.alpha(i)), mk, 'Color', cols.(lab{1}), 'MarkerFaceColor', cols.(lab{1}), ...
                'DisplayName', sprintf('%s (%d)', lab{1}, sum(i)));
        end
        xlabel('\alpha_r D'); ylabel('\alpha_i D'); grid on; yline(0, 'k', 'HandleVisibility', 'off');
        title(sprintf('m = %d, St_D = %.2f, x = %.2f mm', s.m, S.single_St, P.x_mm));
        legend('Interpreter', 'none', 'Location', 'best');
        subplot(2, 2, im + 2); hold on;
        n = numel(s.r);
        for lab = {'KH_jet', 'KH_outer', 'upstream_jet'}
            i = find(s.conv & strcmp(s.label, lab{1}));
            if isempty(i), continue; end
            if strcmp(lab{1}, 'upstream_jet'), [~, j] = min(abs(imag(s.alpha(i)))); else, [~, j] = min(imag(s.alpha(i))); end
            pv = abs(s.V(4*n+1:5*n, i(j))); plot(s.r*L, pv/max(pv), 'Color', cols.(lab{1}), 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%s  \\alpha = %.2f%+.2fi', lab{1}, real(s.alpha(i(j))), imag(s.alpha(i(j)))));
        end
        xline(P.r1*L, 'k--', 'r1', 'HandleVisibility', 'off'); xline(P.r2*L, 'k:', 'r2', 'HandleVisibility', 'off');
        xlabel('r [mm]'); ylabel('|p| (normalised)'); legend('Interpreter', 'none'); xlim([0 P.rmax*L]);
    end
    exportgraphics(f, fullfile(S.out_dir, sprintf('x_%1.1f_St_%2.2f_lsa_single_spectrum.png', S.single_x_mm, ...
        S.single_St)), 'Resolution', 130); close(f);
end

function plot_scan(R, prof, S)
    L = S.D_ref*1e3; nx = numel(prof); x = [prof.x_mm];
    for im = 1:numel(S.m)
        KH = cell2mat(cellfun(@(r) r.KH(im, :), R, 'UniformOutput', false));   % nx x nS
        UP = cell2mat(cellfun(@(r) r.UP(im, :), R, 'UniformOutput', false));
        f = figure('Visible', 'off', 'Position', [0 0 1500 900]);
        subplot(2, 2, 1); imagesc(S.St, x, -imag(KH)); axis xy; colorbar;
        xlabel('St_D'); ylabel('x [mm]'); title(sprintf('m = %d: KH_{jet} growth rate -\\alpha_i D', S.m(im)));
        subplot(2, 2, 2); hold on; cm = parula(nx);
        for k = 1:nx, plot(S.St, -imag(KH(k, :)), 'Color', cm(k, :)); end
        xlabel('St_D'); ylabel('-\alpha_i D'); title('KH_{jet} growth rate (colour = x, blue upstream)'); grid on;
        subplot(2, 2, 3); hold on;
        for k = 1:nx, plot(S.St, real(KH(k, :)), '-', 'Color', cm(k, :)); plot(S.St, real(UP(k, :)), '--', 'Color', cm(k, :)); end
        xlabel('St_D'); ylabel('\alpha_r D'); title('\alpha_r: KH_{jet} (solid), upstream_{jet} (dashed)'); grid on;
        subplot(2, 2, 4); imagesc(S.St, x, double(~isnan(UP))); axis xy; colorbar;
        xlabel('St_D'); ylabel('x [mm]'); title('upstream_{jet} wave found (1) or not (0)');
        sgtitle(sprintf('srp\\_M2.5 local LSA, m = %d (inviscid, mean V dropped; smoothing %.3f mm, N = %d)', ...
            S.m(im), S.smooth_mm, S.N));
        exportgraphics(f, fullfile(S.out_dir, sprintf('lsa_scan_m%d.png', S.m(im))), 'Resolution', 130); close(f);
    end
    % jet and outer shear-layer radii and r_max along x, for reference
    f = figure('Visible', 'off'); plot(x, [prof.r1]*L, x, [prof.r2]*L, x, [prof.rmax]*L, x, [prof.r_bs]*L);
    legend('r1 jet layer', 'r2 outer layer', 'r_{max}', 'bow shock'); xlabel('x [mm]'); ylabel('r [mm]'); grid on;
    exportgraphics(f, fullfile(S.out_dir, 'lsa_scan_radii.png'), 'Resolution', 130); close(f);
end
