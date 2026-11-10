function smoothing_check(outdir, h_mm, variants)
% SMOOTHING_CHECK  Compare shock-smoothing settings on the srp case without solving anything.
%   smoothing_check(outdir, 0.05, {struct('tag','B'), struct('tag','W2','threshold',0.05,'filter_sigma_mm',0.2, ...
%                   'blend_sigma_mm',0.25)})
% Each variant overrides fields of srp_case_defaults().smooth. For each it prints the number of sensor points,
% whether each spurious-mode hotspot is inside the smoothing zone (blend weight), and how sharp the base flow is
% there: the largest density change per cell, max(h |grad rho| / rho), within 0.3 mm of the hotspot.
% Writes <outdir>/smoothing_check_h<h>.png (blend weight and Mach per variant, and rho along the axis).
hot = [11.7 8.8; 7.3 4.5; 4.8 0.03; 8.5 0.7; 8.2 1.8; 13.8 0.3];
hname = {'bow shock off axis', 'triple point', 'axis crossing', 'Mach disk axis', 'Mach disk r 1.8', ...
         'bow shock axis'};
if ~exist(outdir, 'dir'), mkdir(outdir); end
set(0, 'DefaultFigureVisible', 'off');
nv = numel(variants);
f = figure('Position', [0 0 1700 420*(nv+1)]);
for k = 1:nv
    v = variants{k};
    pc = srp_case_defaults(); pc.h_mm = h_mm;
    for fn = fieldnames(v)', if ~strcmp(fn{1}, 'tag'), pc.smooth.(fn{1}) = v.(fn{1}); end, end
    [mesh, BF, info] = build_srp_case(pc);
    mm = 1e-3 / pc.D; u = mesh.usedInd;
    X = mesh.X / mm; Y = mesh.Y / mm;
    rho = BF.RHO(u);
    g = zeros(size(mesh.X)); g(u) = hypot(mesh.Dx*rho, mesh.Dy*rho) * (h_mm*mm) ./ rho;
    fprintf('\n=== %s: threshold %.3f, filter %.2f mm, blend %.2f mm: %d sensor points\n', v.tag, ...
            pc.smooth.threshold, pc.smooth.filter_sigma_mm, pc.smooth.blend_sigma_mm, nnz(info.shock));
    for j = 1:size(hot, 1)
        near = hypot(X - hot(j,1), Y - hot(j,2)) < 0.3 & ~isnan(g) & g > 0;
        [~, i] = min(hypot(X(:) - hot(j,1), Y(:) - hot(j,2)));
        fprintf('  %-20s (%5.1f,%4.1f): blend weight %.2f, max drho/rho per cell %.3f\n', hname{j}, ...
                hot(j,1), hot(j,2), info.wblend(i), max(g(near)));
    end
    gf = g(u); gf(X(u) <= 3) = 0; [~, order] = sort(gf, 'descend');
    fprintf('  whole domain x > 3 mm: max drho/rho per cell %.3f; top 5 at', gf(order(1)));
    Xu = X(u); Yu = Y(u); wb = info.wblend(u);
    for i = order(1:5)', fprintf(' (%.2f,%.2f) %.3f w%.2f', Xu(i), Yu(i), gf(i), wb(i)); end
    s = info.shock; fprintf('\n  sensor points with r > 3 mm and x < 7 mm (shear layer?): %d\n', ...
                            nnz(s & X < 7 & Y > 3 & Y < 6));
    M = nan(size(mesh.X)); M(u) = hypot(BF.U(u), BF.V(u)) ./ sqrt(BF.T(u)) * pc.ref.Ma;
    subplot(nv+1, 2, 2*k-1); pcolor(X, Y, info.wblend); shading flat; axis equal; xlim([-2 16]); ylim([0 11]);
    colorbar; hold on; plot(hot(:,1), hot(:,2), 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    title(sprintf('%s: blend weight (x = spurious-mode hotspots)', v.tag));
    subplot(nv+1, 2, 2*k); pcolor(X, Y, M); shading flat; axis equal; xlim([-2 16]); ylim([0 11]); colorbar;
    caxis([0 5]); title(sprintf('%s: Mach', v.tag));
    subplot(nv+1, 2, [2*nv+1 2*nv+2]); hold on;
    r0 = mesh.Y(:,1) / mm; [~, i0] = min(abs(r0 - 0)); x0 = X(i0, :); rr = BF.RHO(i0, :); rr(~ismember(sub2ind(size(X), i0*ones(size(x0)), 1:numel(x0)), u)) = nan;
    plot(x0, rr, 'DisplayName', v.tag);
end
subplot(nv+1, 2, [2*nv+1 2*nv+2]); xlim([2 16]); grid on; legend('Location', 'best');
xlabel('x (mm)'); ylabel('\rho / \rho_\infty'); title('density on the first row (r = h/2)');
exportgraphics(f, fullfile(outdir, sprintf('smoothing_check_h%s.png', strrep(sprintf('%g', h_mm), '.', 'p'))), ...
               'Resolution', 90);
close(f);
end
