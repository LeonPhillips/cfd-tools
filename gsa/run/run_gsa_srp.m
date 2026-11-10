function R = run_gsa_srp(outdir, h_mm, m_list, shifts, varargin)
% RUN_GSA_SRP  Global stability analysis of the srp_M2.5 base flow (supervisor's method; see gsa/README.md).
%
%   R = run_gsa_srp(outdir, h_mm, m_list, shifts, 'name', value, ...)
%     h_mm    grid spacing (0.1 smoke, 0.07 sweep, 0.05 production)
%     m_list  azimuthal wavenumbers, e.g. [0 1]
%     shifts  complex eigenvalue targets omega (omega = 2 pi f D / U_inf; Im > 0 = growth)
%   options: 'nEig' (50), 'p' (300), 'tol' (1e-8), 'adjoint' (true), 'nsave' (8 modes saved/plotted),
%            'solver' ('mumps'), 'case' (struct of overrides for srp_case_defaults, e.g. struct('sigma_visc',0.5)),
%            'outflow_bc' ('neumann'), 'alpha_rho' (0; density damping on rho', x and r)
%
% For each m: build the case (build_srp_case), the operator (GetLinProblemVarMu, full linearisation),
% the boundary conditions, then for each shift solve with eigs_shift_invert (one MUMPS factorisation shared
% by direct and adjoint). Writes <outdir>/gsa_h<h>_m<m>_s<k>.mat and spectrum / mode plots.
%
% Boundary conditions (method): Dirichlet row-only (1e8) on all variables at the freestream inflow (r) and
% the jet exit (j), on u,v,w at walls (w); Neumann (1e4) on T at walls and on all variables at x min (l) and
% r max (t). Axis: LinStab2D symmetric/antisymmetric derivatives for the given m.

o = struct('nEig', 50, 'p', 300, 'tol', 1e-8, 'adjoint', true, 'nsave', 8, 'solver', 'mumps', 'case', struct(), ...
           'outflow_bc', 'neumann', ...   % x-min outflow: 'neumann' (method) or 'dirichlet' (LinStab2D examples)
           'alpha_rho', 0);               % density damping (method: 3e-5 or 1e-4), see GetLinProblemVarMu
for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k+1}; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
fls_mumps_available();

pc = srp_case_defaults();
pc.h_mm = h_mm;
for f = fieldnames(o.case)', pc.(f{1}) = o.case.(f{1}); end
[mesh, BF, info] = build_srp_case(pc);
U = pc.ref.U; D = pc.D;
f_Hz = @(w) real(w) * U / (2*pi*D);
htag = strrep(sprintf('%g', h_mm), '.', 'p');

R = struct('m', {}, 'shift', {}, 'omega', {}, 'omega_adj', {}, 'info', {});
for m = m_list
    t0 = tic;
    [L0, idx] = GetLinProblemVarMu(mesh, BF, m, struct('alpha_rho', o.alpha_rho));
    [L0, dD1] = BC_DirichletRowOnly(L0, idx, {'r', mesh.inflow_set}, {'rho','u','v','w','T'}, 1e8);
    [L0, dD2] = BC_DirichletRowOnly(L0, idx, {'w'}, {'u','v','w'}, 1e8);
    [L0, dN1] = BC_NeumannRows(L0, mesh, idx, {'w'}, {'T'}, 1e4);
    if isfield(idx, 'si')                    % body-side slip wall: no normal (radial) velocity, rest zero gradient
        [L0, dS1] = BC_DirichletRowOnly(L0, idx, {'s'}, {'v'}, 1e8);
        [L0, dS2] = BC_NeumannRows(L0, mesh, idx, {'s'}, {'rho','u','w','T'}, 1e4);
        dN1 = [dN1; dS1; dS2];
    end
    if strcmpi(o.outflow_bc, 'dirichlet')
        [L0, dN2] = BC_NeumannRows(L0, mesh, idx, {'t'}, {'rho','u','v','w','T'}, 1e4);
        [L0, dD3] = BC_DirichletRowOnly(L0, idx, {'l'}, {'rho','u','v','w','T'}, 1e8);
        dN2 = [dN2; dD3];
    else
        [L0, dN2] = BC_NeumannRows(L0, mesh, idx, {'l','t'}, {'rho','u','v','w','T'}, 1e4);
    end
    L0 = L0 / -1i;                           % eigenvalues as omega: q ~ exp(-1i omega t)
    [W, invW] = GetCompEnergyNorm(mesh, BF, 'axy');
    bcdofs = unique([dD1; dD2; dN1; dN2]);
    t_build = toc(t0);
    fprintf('m = %d: operator n = %d, nnz = %d, built in %.0f s\n', m, size(L0,1), nnz(L0), t_build);

    for ks = 1:numel(shifts)
        sigma = shifts(ks);
        opts = struct('tol', o.tol, 'p', o.p, 'disp', 0);
        if o.adjoint
            [V, w, Va, wa, si] = eigs_shift_invert(L0, o.nEig, sigma, opts, struct('type', o.solver), W, invW);
        else
            [V, w, ~, ~, si] = eigs_shift_invert(L0, o.nEig, sigma, opts, struct('type', o.solver));
            Va = []; wa = [];
        end
        % diagnostics per mode
        ng = mesh.ngp;
        sp = repmat(mesh.sponge(mesh.usedInd) > 0.5, 5, 1);
        Wd = diag(W);
        E  = sum(abs(V).^2 .* Wd, 1);
        f_sponge = sum(abs(V(sp,:)).^2 .* Wd(sp), 1) ./ E;   % energy fraction in the sponges
        [f_grid, f_shock] = mode_metrics(V, Wd, mesh, info.shock);
        pair = nan(numel(w), 1);
        if ~isempty(wa), for j = 1:numel(w), pair(j) = min(abs(w(j) - conj(wa))); end, end
        [~, ord] = sort(imag(w), 'descend');
        T = table((1:numel(w))', w(ord), real(w(ord))/(2*pi), f_Hz(w(ord))/1e3, imag(w(ord)), ...
                  si.resid_direct(ord).', f_sponge(ord).', f_grid(ord).', f_shock(ord).', pair(ord), ...
                  'VariableNames', {'k','omega','St_D','f_kHz','growth','resid','sponge_frac', ...
                                    'grid_frac','shock_frac','adj_pair_dist'});
        fprintf('m = %d, shift %s: %d modes, max residual %.1e, solve %.0f s\n', m, num2str(sigma), ...
                numel(w), max(si.resid_direct), si.t_total_s);
        disp(T(1:min(20, height(T)), :));

        % save
        keep = ord(1:min(o.nsave, numel(ord)));
        modes = V(:, keep);
        P = pressure_modes(modes, BF, mesh, pc.ref);
        S = struct('m', m, 'shift', sigma, 'omega', w, 'omega_adj', wa, 'table', T, 'solver', si, ...
                   'modes', modes, 'p_modes', P, 'keep', keep, 'h_mm', h_mm, 'case_params', pc, ...
                   'n', size(L0,1), 'nnz', nnz(L0), 't_build_s', t_build, 'bc_dofs', numel(bcdofs), ...
                   'alpha_rho', o.alpha_rho, 'outflow_bc', o.outflow_bc);
        if ~isempty(Va), S.modes_adj = Va(:, keep); end
        tag = sprintf('gsa_h%s_m%d_s%d', htag, m, ks);
        save(fullfile(outdir, [tag '.mat']), '-struct', 'S', '-v7.3');
        plot_spectrum_modes(S, mesh, pc, fullfile(outdir, tag));
        R(end+1) = struct('m', m, 'shift', sigma, 'omega', w, 'omega_adj', wa, 'info', si); %#ok<AGROW>
    end
    clear L0 W invW
end
end

% ------------------------------------------------------------------------------------------------
function [f_grid, f_shock] = mode_metrics(V, Wd, mesh, shock)
% Diagnostics to separate numerical from physical modes (energy-norm weighted, all 5 variables):
%  f_grid  : fraction of the energy at grid scale, |q - S q|^2 / |q|^2, with S the masked 3x3 [1 2 1]x[1 2 1]/16
%            filter (a pure 2h checkerboard gives ~1; a smooth mode gives O(h^2))
%  f_shock : fraction of the energy within 3 cells of the shock-sensor points
ng = mesh.ngp; u = mesh.usedInd; sz = size(mesh.X);
fluid = false(sz); fluid(u) = true;
k1 = [1 2 1]; K = (k1' * k1) / 16;
den = conv2(double(fluid), K, 'same');
band = imdilate(shock, strel('square', 7)) & fluid;
inband = repmat(band(u), 5, 1);
nm = size(V, 2);
f_grid = zeros(1, nm); f_shock = zeros(1, nm);
for k = 1:nm
    q = V(:, k);
    Etot = sum(abs(q).^2 .* Wd);
    Eg = 0;
    for v = 1:5
        g = zeros(sz); g(u) = q((v-1)*ng + (1:ng));
        s = conv2(g, K, 'same') ./ max(den, eps);
        d = g(u) - s(u);
        Eg = Eg + sum(abs(d).^2 .* Wd((v-1)*ng + (1:ng)));
    end
    f_grid(k)  = Eg / Etot;
    f_shock(k) = sum(abs(q(inband)).^2 .* Wd(inband)) / Etot;
end
end

function P = pressure_modes(V, BF, mesh, ref)
% p' = (rho' T + rho T') / (kappa Ma^2), q = [rho, u, v, w, T]
ng = mesh.ngp; u = mesh.usedInd;
rho = BF.RHO(u); T = BF.T(u);
P = (V(1:ng,:) .* T + rho .* V(4*ng+1:5*ng,:)) / (ref.kappa * ref.Ma^2);
end

function plot_spectrum_modes(S, mesh, pc, base)
set(0, 'DefaultFigureVisible', 'off');
mm = 1e-3 / pc.D;
f = figure('Position', [0 0 900 600]); hold on; grid on;
plot(real(S.omega)/(2*pi), imag(S.omega), 'bo', 'DisplayName', 'direct');
if ~isempty(S.omega_adj), plot(real(S.omega_adj)/(2*pi), -imag(S.omega_adj), 'rx', 'DisplayName', 'conj(adjoint)'); end
plot(real(S.shift)/(2*pi), imag(S.shift), 'k+', 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', 'shift');
yline(0, 'k-', 'HandleVisibility', 'off');
xlabel('St_D = \omega_r / 2\pi'); ylabel('\omega_i (growth rate)');
ax = gca; xt = ax.XTick;
title(sprintf('m = %d, h = %g mm, shift %s; f (kHz) = St_D x %.2f', S.m, S.h_mm, num2str(S.shift), pc.ref.U/pc.D/1e3));
legend('Location', 'best');
exportgraphics(f, [base '_spectrum.png'], 'Resolution', 110); close(f);

n = min(4, size(S.p_modes, 2));
f = figure('Position', [0 0 1600 900]);
for k = 1:n
    g = nan(size(mesh.X)); g(mesh.usedInd) = real(S.p_modes(:,k));
    g = g / max(abs(g(:)), [], 'omitnan');
    subplot(2, 2, k);
    pcolor(mesh.X/mm, mesh.Y/mm, g); shading flat; axis equal; xlim([-3 18]); ylim([0 12]); caxis([-1 1]);
    colorbar; colormap(redblue());
    w = S.omega(S.keep(k));
    title(sprintf('Re(p''), \\omega = %.3f %+.3fi, St_D %.3f, f %.1f kHz', real(w), imag(w), real(w)/(2*pi), ...
          real(w)*pc.ref.U/(2*pi*pc.D)/1e3));
end
exportgraphics(f, [base '_modes.png'], 'Resolution', 100); close(f);
end

function c = redblue()
n = 128; t = linspace(0, 1, n)';
c = [[t; ones(n,1)], [t; flipud(t)], [ones(n,1); flipud(t)]];
end
