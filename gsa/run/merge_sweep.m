function merge_sweep(case_dir, out_base, varargin)
% MERGE_SWEEP  Collect the eigenvalues of a shift sweep (many run_gsa_srp outputs), remove repeats found from
% several shifts, match modes between grids, write CSVs and a spectrum plot.
%
%   merge_sweep(case_dir, out_base, 'u05', 'sweep_u05_*', 'fine', 'sweep_fine_*', ...)
%
% Each label/pattern pair names one grid; pattern is a folder glob under case_dir whose gsa_*_m<m>_s<k>.mat
% files are read (only the 'table', 'm' and 'shift' variables, via matfile). St = omega_r / 2 pi (D = 18 mm,
% U_inf), f in kHz from the run table.
% Writes:
%   <out_base>_modes.csv    every unique mode per grid and m (St, f, growth, grid/shock/sponge fractions,
%                           residual, times found, file and table row of the best copy, match on the other grid)
%   <out_base>_matches.csv  modes found on the first grid and on another grid within the match tolerance
%   <out_base>_spectrum.png growth vs St, one panel per m, colour = grid_frac, matched modes ringed
% Options (name/value after the pairs, names start with '-'): '-tol_same' (1e-4, relative |d omega| for the
% same eigenvalue from two shifts), '-tol_match' (0.02, absolute |d omega| between grids), '-growth_min' (-1).

o = struct('tol_same', 1e-4, 'tol_match', 0.02, 'growth_min', -1);
k = find(cellfun(@(a) ischar(a) && startsWith(a, '-'), varargin), 1);
if ~isempty(k)
    opt = varargin(k:end); varargin = varargin(1:k-1);
    for j = 1:2:numel(opt), o.(opt{j}(2:end)) = opt{j+1}; end
end
labels = varargin(1:2:end); pats = varargin(2:2:end);

A = table();
for g = 1:numel(labels)
    dirs = dir(fullfile(case_dir, pats{g}));
    dirs = dirs([dirs.isdir]);
    for d = dirs'
        fs = dir(fullfile(d.folder, d.name, 'gsa_*_m*_s*.mat'));
        for f = fs'
            fn = fullfile(f.folder, f.name);
            M = matfile(fn);
            T = M.table; m = M.m; sh = M.shift;
            T.grid  = repmat(string(labels{g}), height(T), 1);
            T.m     = repmat(m, height(T), 1);
            T.shift = repmat(sh, height(T), 1);
            T.file  = repmat(string(fullfile(d.name, f.name)), height(T), 1);
            T.row   = T.k;
            A = [A; T(:, {'grid','m','shift','omega','St_D','f_kHz','growth','resid','sponge_frac', ...
                         'grid_frac','shock_frac','file','row'})]; %#ok<AGROW>
        end
        fprintf('%s: %s, %d files\n', labels{g}, d.name, numel(fs));
    end
end
if isempty(A), error('merge_sweep: no gsa_*.mat files found'); end

% unique eigenvalues per grid and m: same omega (relative tol) from several shifts -> keep the lowest residual
U = table();
for g = unique(A.grid)'
    for m = unique(A.m(A.grid == g))'
        B = A(A.grid == g & A.m == m, :);
        B = sortrows(B, 'resid');
        used = false(height(B), 1); keep = []; nfound = [];
        for i = 1:height(B)
            if used(i), continue; end
            same = ~used & abs(B.omega - B.omega(i)) <= o.tol_same * max(abs(B.omega(i)), 1);
            used(same) = true; keep(end+1) = i; nfound(end+1) = nnz(same); %#ok<AGROW>
        end
        C = B(keep, :); C.n_found = nfound(:);
        U = [U; C]; %#ok<AGROW>
    end
end
U = sortrows(U, {'grid', 'm', 'St_D'});

% match between grids (to the nearest mode of each other grid, same m)
grids = string(labels(:));
U.match_grid = strings(height(U), 1); U.match_dist = nan(height(U), 1); U.match_omega = nan(height(U), 1);
for i = 1:height(U)
    best = inf;
    for g = grids'
        if g == U.grid(i), continue; end
        o2 = U.omega(U.grid == g & U.m == U.m(i));
        if isempty(o2), continue; end
        [dmin, j] = min(abs(o2 - U.omega(i)));
        if dmin < best, best = dmin; U.match_grid(i) = g; U.match_omega(i) = o2(j); end
    end
    U.match_dist(i) = best;
end
U.matched = U.match_dist <= o.tol_match;

W = U; W.omega_r = real(W.omega); W.omega_i = imag(W.omega); W.match_omega_r = real(W.match_omega);
W.match_omega_i = imag(W.match_omega); W.omega = []; W.match_omega = []; W.shift = [];
writetable(W, [out_base '_modes.csv']);
Mt = W(W.matched & W.grid == grids(1), :);
writetable(Mt, [out_base '_matches.csv']);
fprintf('unique modes: %s; matched within %.3g: %d (first grid)\n', ...
        strjoin(arrayfun(@(g) sprintf('%s %d', g, nnz(U.grid == g)), grids', 'UniformOutput', false), ', '), ...
        o.tol_match, height(Mt));

% plot
set(0, 'DefaultFigureVisible', 'off');
ms = unique(U.m)';
f = figure('Position', [0 0 1400 520*numel(ms)]);
mk = {'o', 's', '^', 'd'};
for im = 1:numel(ms)
    subplot(numel(ms), 1, im); hold on; grid on; box on;
    for ig = 1:numel(grids)
        S = U(U.grid == grids(ig) & U.m == ms(im) & U.growth >= o.growth_min, :);
        scatter(S.St_D, S.growth, 28, S.grid_frac, mk{1 + mod(ig-1, 4)}, 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.3, 'DisplayName', char(grids(ig)));
        Sm = S(S.matched, :);
        plot(Sm.St_D, Sm.growth, 'o', 'MarkerSize', 13, 'Color', [0.85 0.1 0.1], 'LineWidth', 1.2, ...
             'HandleVisibility', 'off');
    end
    yline(0, 'k-', 'HandleVisibility', 'off');
    colormap(parula); c = colorbar; c.Label.String = 'grid\_frac'; clim([0 1]);
    xlabel('St = f D / U_\infty  (D = 18 mm)'); ylabel('\omega_i (growth rate, U_\infty / D)');
    ax = gca; ax2 = axes('Position', ax.Position, 'XAxisLocation', 'top', 'Color', 'none', 'YTick', [], ...
        'XLim', ax.XLim * U_over_D_kHz(A), 'Box', 'off'); xlabel(ax2, 'f (kHz)'); axes(ax); %#ok<LAXES>
    title(ax2, sprintf('m = %d: %s; red ring = matched between grids (|d\\omega| <= %.3g)', ms(im), ...
          strjoin(cellstr(grids), ' / '), o.tol_match));
    legend('Location', 'southwest');
end
exportgraphics(f, [out_base '_spectrum.png'], 'Resolution', 120); close(f);
end

function c = U_over_D_kHz(A)
% f_kHz / St_D from the run tables
k = find(abs(A.St_D) > 0, 1);
c = A.f_kHz(k) / A.St_D(k);
end
