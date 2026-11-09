function check_srp_case(mesh, BF, info, outdir, tag)
% CHECK_SRP_CASE  Check plots for build_srp_case: grid/mask/boundaries/sponge, Mach raw vs smoothed with the
% shock sensor, mu_t, and Mach on the monitor lines (CFD nodes vs grid raw vs grid smoothed).
if ~exist(outdir, 'dir'), mkdir(outdir); end
p = info.p; mm = 1e-3 / p.D;
X = mesh.X / mm; Y = mesh.Y / mm;           % mm for plotting
solid = mesh.mask;
M  = @(U, V, T) p.ref.Ma * hypot(U, V) ./ sqrt(T);
Mr = M(info.raw.U, info.raw.V, info.raw.T);  Mr(solid) = NaN;
Ms = M(BF.U, BF.V, BF.T);                    Ms(solid) = NaN;
set(0, 'DefaultFigureVisible', 'off');

%% 1. grid, mask, boundary sets, sponge
f = figure('Position', [0 0 1500 650]);
subplot(1,2,1); hold on;
Xu = mesh.X(mesh.usedInd)/mm; Yu = mesh.Y(mesh.usedInd)/mm;
plot(X(solid), Y(solid), 's', 'Color', [.75 .75 .75], 'MarkerSize', 2, 'MarkerFaceColor', [.75 .75 .75]);
c = struct('ji', 'r', 'wi', 'k', 'li', 'b', 'ri', 'm', 'ti', 'g');
for b = fieldnames(c)'
    k = mesh.idx.(b{1}); plot(Xu(k), Yu(k), '.', 'Color', c.(b{1}), 'MarkerSize', 8, 'DisplayName', b{1});
end
k = mesh.idx.wi; quiver(Xu(k), Yu(k), mesh.wall_normal(:,1), mesh.wall_normal(:,2), 0.3, 'k');
axis equal; xlim([-1.5 3.5]); ylim([0 4]); grid on; xlabel('x (mm)'); ylabel('r (mm)');
title(sprintf('mask (grey), ji red, wi black + normals; h = %.3f mm', p.h_mm));
subplot(1,2,2);
s = mesh.sponge; s(solid) = NaN;
pcolor(X, Y, s); shading flat; axis equal tight; colorbar; xlabel('x (mm)'); ylabel('r (mm)');
title(sprintf('sponge (in %g, out %g, top %g)', p.sponge.amp_in, p.sponge.amp_out, p.sponge.amp_top));
exportgraphics(f, fullfile(outdir, [tag '_grid_sponge.png']), 'Resolution', 110); close(f);

%% 2. Mach raw / smoothed / sensor / mu_t
f = figure('Position', [0 0 1600 1000]);
lim = {[-3 18], [0 12]};
ax = subplot(2,2,1); pcolor(X, Y, Mr); shading flat; axis equal; xlim(lim{1}); ylim(lim{2}); colorbar; caxis([0 4]);
title('Mach, interpolated (raw)'); colormap(ax, turbo);
ax = subplot(2,2,2); pcolor(X, Y, Ms); shading flat; axis equal; xlim(lim{1}); ylim(lim{2}); colorbar; caxis([0 4]);
hold on; contour(X, Y, double(info.shock), [0.5 0.5], 'w', 'LineWidth', 0.5);
title(sprintf('Mach, smoothed (white: shock sensor, thr %.2f; filter %.2f mm, blend %.2f mm)', ...
      p.smooth.threshold, p.smooth.filter_sigma_mm, p.smooth.blend_sigma_mm)); colormap(ax, turbo);
ax = subplot(2,2,3); d = Ms - Mr; pcolor(X, Y, d); shading flat; axis equal; xlim(lim{1}); ylim(lim{2}); colorbar;
caxis([-0.5 0.5]); title('Mach smoothed - raw'); colormap(ax, redblue_map());
ax = subplot(2,2,4); mt = BF.MU_t; mt(solid) = NaN; pcolor(X, Y, log10(max(mt, 1e-2))); shading flat; axis equal;
xlim(lim{1}); ylim(lim{2}); colorbar; title(sprintf('log_{10}(\\mu_t/\\mu_\\infty), max %.0f (sigma_{visc} = %g)', ...
      max(mt(:)), p.sigma_visc)); colormap(ax, parula);
exportgraphics(f, fullfile(outdir, [tag '_mach_mut.png']), 'Resolution', 110); close(f);

%% 3. Mach on the monitor lines: CFD nodes vs grid
S = read_flow(p.baseflow_prefix);
F = scatteredInterpolant(S.x*1e3, S.r*1e3, S.Mach, 'natural', 'none');
lines = {[1 0.01 20 0.01], [1 1.3 20 1.3], [1 2.3 20 2.3], [1 3 20 3], [1 5 20 5], [2.0 2.34 7.84 4.70]};
f = figure('Position', [0 0 1500 1000]);
Mr(isnan(Mr)) = 0; Ms(isnan(Ms)) = 0;
for k = 1:numel(lines)
    L = lines{k}; t = linspace(0, 1, 2000)';
    xl = L(1) + t*(L(3)-L(1)); rl = L(2) + t*(L(4)-L(2));
    rg = max(rl, Y(1,1));                   % the grid's first row is at dr/2 (axis line r = 0.01 mm)
    subplot(3,2,k); hold on;
    plot(xl, F(xl, rl), 'k-', 'DisplayName', 'CFD (SU2 nodes)');
    plot(xl, interp2(X, Y, Mr, xl, rg), 'b--', 'DisplayName', sprintf('grid, raw (r >= %.3f mm)', Y(1,1)));
    plot(xl, interp2(X, Y, Ms, xl, rg), 'r-', 'DisplayName', 'grid, smoothed');
    grid on; xlabel('x (mm)'); ylabel('Mach');
    title(sprintf('(%.2f, %.2f) -> (%.2f, %.2f) mm', L)); if k == 1, legend('Location', 'best'); end
end
exportgraphics(f, fullfile(outdir, [tag '_mach_lines.png']), 'Resolution', 110); close(f);
fprintf('check plots written to %s\n', outdir);
end

function S = read_flow(prefix)
m = jsondecode(fileread([prefix '_flow.json']));
fid = fopen([prefix '_flow.bin'], 'r'); a = fread(fid, [m.N, numel(m.names)], 'double'); fclose(fid);
S = struct(); for k = 1:numel(m.names), S.(m.names{k}) = a(:, k); end
end

function c = redblue_map()
n = 128; t = linspace(0, 1, n)';
c = [[t; ones(n,1)], [t; flipud(t)], [ones(n,1); flipud(t)]];
end
