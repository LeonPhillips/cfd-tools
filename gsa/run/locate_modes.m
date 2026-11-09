function locate_modes(matfile)
% LOCATE_MODES  Where the saved modes of a run_gsa_srp result live: energy in each sponge (x max, x min, r max),
% near each boundary set, and the location of max |p'|. Also writes <matfile>_where.png (full-domain |p'|).
S = load(matfile);
f = fieldnames(S); pc = S.(f{find(startsWith(f, 'case'), 1)});   % saved as 'case' (a keyword: load renames it)
[mesh, BF] = build_srp_case(pc); %#ok<ASGLU>
mm = 1e-3 / pc.D; u = mesh.usedInd; ng = mesh.ngp;
X = mesh.X(u)/mm; Y = mesh.Y(u)/mm;
x0 = min(X); x1 = max(X); R = max(Y); w = pc.sponge;
zone = struct('x_max_sponge', X > x1 - w.width_in_mm, 'x_min_sponge', X < x0 + w.width_out_mm, ...
              'r_max_sponge', Y > R - w.width_top_mm);
[W] = GetCompEnergyNorm(mesh, BF, 'axy'); Wd = full(diag(W));
fprintf('%s\n', matfile);
for k = 1:numel(S.keep)
    q = S.modes(:, k);
    e = reshape(abs(q).^2 .* Wd, ng, 5); e = sum(e, 2); E = sum(e);
    p = abs(S.p_modes(:, k)); [~, i] = max(p);
    s = sprintf('mode %d omega %.3f%+.3fi: max|p''| at (%.2f, %.2f) mm;', k, real(S.omega(S.keep(k))), ...
                imag(S.omega(S.keep(k))), X(i), Y(i));
    for f = fieldnames(zone)', s = [s sprintf(' %s %.2f', f{1}, sum(e(zone.(f{1})))/E)]; end %#ok<AGROW>
    for b = {'ri','ji','wi','li','ti'}
        near = false(ng, 1); near(mesh.idx.(b{1})) = true;
        s = [s sprintf(' | at %s %.2f', b{1}, sum(e(near))/E)]; %#ok<AGROW>
    end
    disp(s);
end
set(0, 'DefaultFigureVisible', 'off');
f = figure('Position', [0 0 1500 900]);
for k = 1:min(4, numel(S.keep))
    g = nan(size(mesh.X)); g(u) = abs(S.p_modes(:, k)); g = g / max(g(:));
    subplot(2,2,k); pcolor(mesh.X/mm, mesh.Y/mm, log10(max(g, 1e-6))); shading flat; axis equal tight;
    caxis([-6 0]); colorbar; title(sprintf('log10 |p''|, mode %d, \\omega = %.3f%+.3fi', k, ...
        real(S.omega(S.keep(k))), imag(S.omega(S.keep(k)))));
end
exportgraphics(f, [erase(matfile, '.mat') '_where.png'], 'Resolution', 100); close(f);
end
