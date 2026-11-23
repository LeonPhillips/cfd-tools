% freejet_summary.m: free-jet LSA at x = 3 mm (user, 2026-11-24): KH_jet vs St_j, the N / r_max checks, and free jet vs counterflow at 40 kHz.
% Run on a compute node in /scratch/qk75/LeonPhillips/PhD/LSA. Writes PNG + CSV into single_freeJet/.
W = '/scratch/qk75/LeonPhillips/PhD/LSA'; D = fullfile(W, 'single_freeJet');
pick = @(s, r1) pickKH(s, r1);

% --- 1. the scan
f = dir(fullfile(D, '*_n91_freeJet_lsa_single.mat')); A = load(fullfile(D, f(1).name));
S = A.S; P = A.prof(1); St = S.single_St; nm = numel(S.m);
fprintf('scan file %s, r1 = %.3f mm, r2 = %.3f mm\n', f(1).name, P.r1*4.6, P.r2*4.6);
T = [];
for im = 1:nm
    for is = 1:numel(St)
        s = A.out.sol{is, im}; a = pick(s, P.r1);
        T = [T; St(is), St(is)*S.U_j/S.D_j/1e3, S.m(im), real(a), imag(a), (s.omega/real(a))*S.a_inf/S.U_j]; %#ok<AGROW>
    end
end
Tt = array2table(T, 'VariableNames', {'St_j', 'f_kHz', 'm', 'alpha_r_De', 'alpha_i_De', 'c_over_Uj'});
writetable(Tt, fullfile(D, 'fj_x3.0_KH_vs_Stj.csv'));
fg = figure('Visible', 'off', 'Position', [0 0 1400 450]); cols = lines(2);
for k = 1:3
    subplot(1, 3, k); hold on;
    for im = 1:nm
        q = T(:, 3) == S.m(im); y = {T(q,4), -T(q,5), T(q,6)};
        plot(T(q, 1), y{k}, 'o-', 'Color', cols(im,:), 'MarkerSize', 3, 'DisplayName', sprintf('m = %d', S.m(im)));
    end
    xlabel('St_j = f D_j / U_j'); grid on; legend('Location', 'best');
    labs = {'\alpha_r D_e', '-\alpha_i D_e', 'c_r / U_j'}; ylabel(labs{k});
end
sgtitle(sprintf('Free jet, x = 3.0 mm: most unstable converged KH\\_jet mode (r_1 = %.2f mm)', P.r1*4.6));
exportgraphics(fg, fullfile(D, 'fj_x3.0_KH_vs_Stj.png'), 'Resolution', 130); close(fg);

% --- 2. checks at five St_j: base (N 250, r_max 20), N 350, r_max 15 and 25
tags = {'', '_N350', '_rmax15', '_rmax25'}; names = {'base', 'N350', 'rmax15', 'rmax25'};
fid = fopen(fullfile(D, 'fj_x3.0_checks.csv'), 'w'); fprintf(fid, 'case,St_j,m,alpha_r,alpha_i,rel_diff_vs_base\n');
B = [];
for it = 1:numel(tags)
    if it == 2, pat = '*_n5_N350_freeJet_lsa_single.mat'; else, pat = sprintf('*_n5_freeJet%s_lsa_single.mat', tags{it}); end
    g = dir(fullfile(D, pat)); if isempty(g), fprintf('missing %s\n', pat); continue; end
    C = load(fullfile(D, g(1).name)); Pc = C.prof(1);
    for im = 1:nm
        for is = 1:numel(C.S.single_St)
            a = pick(C.out.sol{is, im}, Pc.r1);
            if it == 1, B(is, im) = a; end %#ok<SAGROW>
            fprintf(fid, '%s,%.3f,%d,%.5f,%.5f,%.2e\n', names{it}, C.S.single_St(is), C.S.m(im), real(a), imag(a), ...
                abs(a - B(is, im))/abs(B(is, im)));
            fprintf('%-7s St_j %.3f m %d: KH %.4f%+.4fi  rel diff %.2e\n', names{it}, C.S.single_St(is), C.S.m(im), ...
                real(a), imag(a), abs(a - B(is, im))/abs(B(is, im)));
        end
    end
end
fclose(fid);

% --- 3. free jet (St_j 0.405 = 40.07 kHz) vs counterflow (St_D 1.271 = 40.05 kHz) at x = 3 mm
F = load(fullfile(D, 'x_3.0_St_0.41_freeJet_lsa_single.mat')); G = load(fullfile(W, 'single', 'x_3.0_St_1.27_lsa_single.mat'));
fg = figure('Visible', 'off', 'Position', [0 0 1300 800]);
for im = 1:nm
    for c = 1:2
        if c == 1, X = F; lab = 'free jet, St_j 0.405'; else, X = G; lab = 'counterflow, St_D 1.271'; end
        s = X.out.sol{1, im}; Pp = X.prof(1); [a, j] = pick(s, Pp.r1);
        k = find(s.Vidx == j); n = numel(s.r); v = s.V(:, k);
        u = abs(v(n+1:2*n)); p = abs(v(4*n+1:5*n));
        subplot(2, 2, im); hold on; plot(s.r*4.6, u/max(u), 'LineWidth', 1.2, 'DisplayName', sprintf('%s: \\alpha = %.2f%+.2fi', lab, real(a), imag(a)));
        subplot(2, 2, 2 + im); hold on; plot(s.r*4.6, p/max(p), 'LineWidth', 1.2, 'DisplayName', lab);
        fprintf('%s m %d: KH %.4f%+.4fi, r1 %.3f mm\n', lab, S.m(im), real(a), imag(a), Pp.r1*4.6);
    end
    subplot(2, 2, im); xlim([0 8]); grid on; xlabel('r [mm]'); ylabel('|u| / max'); title(sprintf('m = %d: KH |u|', S.m(im))); legend('Location', 'best');
    subplot(2, 2, 2 + im); xlim([0 8]); grid on; xlabel('r [mm]'); ylabel('|p| / max'); title(sprintf('m = %d: KH |p|', S.m(im))); legend('Location', 'best');
end
sgtitle('x = 3.0 mm, ~40 kHz: KH\_jet eigenfunctions, free jet vs counterflow (both N = 250)');
exportgraphics(fg, fullfile(D, 'fj_vs_cf_x3.0_40kHz_KH.png'), 'Resolution', 130); close(fg);
disp('ANALYSE_DONE');

function [a, j] = pickKH(s, r1)
% most unstable converged KH_jet mode with alpha_r > 0 (the label alone also catches an alpha_r < 0 mode)
    i = find(s.conv & strcmp(s.label, 'KH_jet') & real(s.alpha) > 0);
    if isempty(i), a = NaN; j = NaN; return; end
    [~, k] = min(imag(s.alpha(i))); j = i(k); a = s.alpha(j);
end
