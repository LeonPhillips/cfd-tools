% freejet_summary.m: free-jet LSA summary at one station (user, 2026-11-24): KH_jet vs St_j, the N / r_max checks
% (if run at that station), free jet vs counterflow at ~40 kHz, and KH vs St_j for every station scanned so far.
% Run on a compute node in /scratch/qk75/LeonPhillips/PhD/LSA, e.g.  matlab -batch "FJ_X=4; run('.../freejet_summary.m')"
% FJ_X [mm] picks the station (default 3.0). Writes PNG + CSV into single_freeJet/.
if ~exist('FJ_X', 'var'), FJ_X = 3.0; end
W = '/scratch/qk75/LeonPhillips/PhD/LSA'; D = fullfile(W, 'single_freeJet');
xs = sprintf('x%.1f', FJ_X); xf = sprintf('x_%.1f', FJ_X);
pick = @(s, r1) pickKH(s, r1);

% --- 1. the scan at this station
f = dir(fullfile(D, [xf '_St_*_n91_freeJet_lsa_single.mat'])); A = load(fullfile(D, f(1).name));
S = A.S; P = A.prof(1); St = S.single_St; nm = numel(S.m);
fprintf('scan file %s, r1 = %.3f mm, r2 = %.3f mm, fit R^2 U %.4f rho %.4f p %.4f\n', f(1).name, P.r1*4.6, P.r2*4.6, ...
    P.fit.U.R2, P.fit.rho.R2, P.fit.p.R2);
T = [];
for im = 1:nm
    for is = 1:numel(St)
        s = A.out.sol{is, im}; a = pick(s, P.r1);
        T = [T; St(is), St(is)*S.U_j/S.D_j/1e3, S.m(im), real(a), imag(a), (s.omega/real(a))*S.a_inf/S.U_j]; %#ok<AGROW>
    end
end
Tt = array2table(T, 'VariableNames', {'St_j', 'f_kHz', 'm', 'alpha_r_De', 'alpha_i_De', 'c_over_Uj'});
writetable(Tt, fullfile(D, ['fj_' xs '_KH_vs_Stj.csv']));
fg = figure('Visible', 'off', 'Position', [0 0 1400 450]); cols = lines(2);
labs = {'\alpha_r D_e', '-\alpha_i D_e', 'c_r / U_j'};
for k = 1:3
    subplot(1, 3, k); hold on;
    for im = 1:nm
        q = T(:, 3) == S.m(im); y = {T(q,4), -T(q,5), T(q,6)};
        plot(T(q, 1), y{k}, 'o-', 'Color', cols(im,:), 'MarkerSize', 3, 'DisplayName', sprintf('m = %d', S.m(im)));
    end
    xlabel('St_j = f D_j / U_j'); grid on; legend('Location', 'best'); ylabel(labs{k});
end
sgtitle(sprintf('Free jet, x = %.1f mm: most unstable converged KH\\_jet mode (r_1 = %.2f mm)', FJ_X, P.r1*4.6));
exportgraphics(fg, fullfile(D, ['fj_' xs '_KH_vs_Stj.png']), 'Resolution', 130); close(fg);
for sj = [0.1 0.2 0.3 0.405 0.455]
    for im = 1:nm
        r = T(abs(T(:,1) - sj) < 1e-9 & T(:,3) == S.m(im), :);
        fprintf('St_j %.3f (%.1f kHz) m %d: KH %.4f%+.4fi  c/U_j %.3f\n', sj, r(2), S.m(im), r(4), r(5), r(6));
    end
end

% --- 2. checks at five St_j: base (N 250, r_max 20), N 350, r_max 15 and 25 (only where they were run)
tags = {'', '_N350', '_rmax15', '_rmax25'}; names = {'base', 'N350', 'rmax15', 'rmax25'};
if ~isempty(dir(fullfile(D, [xf '_St_*_n5_freeJet_lsa_single.mat'])))
    fid = fopen(fullfile(D, ['fj_' xs '_checks.csv']), 'w'); fprintf(fid, 'case,St_j,m,alpha_r,alpha_i,rel_diff_vs_base\n');
    B = [];
    for it = 1:numel(tags)
        if it == 2, pat = [xf '_St_*_n5_N350_freeJet_lsa_single.mat'];
        else, pat = sprintf('%s_St_*_n5_freeJet%s_lsa_single.mat', xf, tags{it}); end
        g = dir(fullfile(D, pat)); if isempty(g), fprintf('missing %s\n', pat); continue; end
        C = load(fullfile(D, g(1).name)); Pc = C.prof(1);
        for im = 1:nm
            for is = 1:numel(C.S.single_St)
                a = pick(C.out.sol{is, im}, Pc.r1);
                if it == 1, B(is, im) = a; end %#ok<SAGROW>
                fprintf(fid, '%s,%.3f,%d,%.5f,%.5f,%.2e\n', names{it}, C.S.single_St(is), C.S.m(im), real(a), imag(a), ...
                    abs(a - B(is, im))/abs(B(is, im)));
            end
        end
    end
    fclose(fid);
end

% --- 3. free jet (scan entry St_j 0.405 = 40.07 kHz) vs counterflow (St_D 1.271 = 40.05 kHz) at this station
cf = fullfile(W, 'single', sprintf('%s_St_1.27_lsa_single.mat', xf));
if exist(cf, 'file')
    G = load(cf); [~, i405] = min(abs(St - 0.405));
    fg = figure('Visible', 'off', 'Position', [0 0 1300 800]);
    for im = 1:nm
        for c = 1:2
            if c == 1, s = A.out.sol{i405, im}; Pp = P; lab = sprintf('free jet, St_j %.3f', St(i405));
            else, s = G.out.sol{1, im}; Pp = G.prof(1); lab = 'counterflow, St_D 1.271'; end
            [a, j] = pick(s, Pp.r1);
            if isnan(j), fprintf('%s m %d: no converged KH_jet\n', lab, S.m(im)); continue; end
            if isfield(s, 'Vidx'), k = find(s.Vidx == j); else, k = j; end   % older files keep every vector
            n = numel(s.r); v = s.V(:, k);
            u = abs(v(n+1:2*n)); p = abs(v(4*n+1:5*n));
            subplot(2, 2, im); hold on; plot(s.r*4.6, u/max(u), 'LineWidth', 1.2, 'DisplayName', sprintf('%s: \\alpha = %.2f%+.2fi', lab, real(a), imag(a)));
            subplot(2, 2, 2 + im); hold on; plot(s.r*4.6, p/max(p), 'LineWidth', 1.2, 'DisplayName', lab);
            fprintf('%s m %d: KH %.4f%+.4fi, r1 %.3f mm\n', lab, S.m(im), real(a), imag(a), Pp.r1*4.6);
        end
        subplot(2, 2, im); xlim([0 8]); grid on; xlabel('r [mm]'); ylabel('|u| / max'); title(sprintf('m = %d: KH |u|', S.m(im))); legend('Location', 'best');
        subplot(2, 2, 2 + im); xlim([0 8]); grid on; xlabel('r [mm]'); ylabel('|p| / max'); title(sprintf('m = %d: KH |p|', S.m(im))); legend('Location', 'best');
    end
    sgtitle(sprintf('x = %.1f mm, ~40 kHz: KH\\_jet eigenfunctions, free jet vs counterflow (both N = 250)', FJ_X));
    exportgraphics(fg, fullfile(D, ['fj_vs_cf_' xs '_40kHz_KH.png']), 'Resolution', 130); close(fg);
else
    fprintf('no counterflow file %s\n', cf);
end

% --- 4. KH vs St_j at every station scanned so far
g = dir(fullfile(D, 'fj_x*_KH_vs_Stj.csv'));
if numel(g) > 1
    fg = figure('Visible', 'off', 'Position', [0 0 1400 800]); cx = lines(numel(g));
    for ig = 1:numel(g)
        R = readmatrix(fullfile(D, g(ig).name)); xl = regexp(g(ig).name, 'x[0-9.]+(?=_KH)', 'match', 'once');
        for im = 1:nm
            q = R(:, 3) == S.m(im) & R(:, 1) >= 0.04;           % below St_j ~ 0.04 the pick is not the KH branch
            y = {R(q,4), -R(q,5), R(q,6)};
            for k = 1:3
                subplot(2, 3, 3*(im-1) + k); hold on;
                plot(R(q, 1), y{k}, '-', 'Color', cx(ig,:), 'LineWidth', 1.2, 'DisplayName', sprintf('x = %s mm', xl(2:end)));
                grid on; xlabel('St_j'); ylabel(labs{k}); title(sprintf('m = %d', S.m(im))); legend('Location', 'best');
            end
        end
    end
    sgtitle('Free jet: most unstable converged KH\_jet mode against St_j at each station (St_j \geq 0.04)');
    exportgraphics(fg, fullfile(D, 'fj_KH_vs_Stj_stations.png'), 'Resolution', 130); close(fg);
end
disp('ANALYSE_DONE');

function [a, j] = pickKH(s, r1) %#ok<INUSD>
% most unstable converged KH_jet mode with alpha_r > 0 (the label alone also catches an alpha_r < 0 mode)
    i = find(s.conv & strcmp(s.label, 'KH_jet') & real(s.alpha) > 0);
    if isempty(i), a = NaN; j = NaN; return; end
    [~, k] = min(imag(s.alpha(i))); j = i(k); a = s.alpha(j);
end
