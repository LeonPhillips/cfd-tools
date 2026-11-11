function test_stretched_grid(h_out, h_box)
% TEST_STRETCHED_GRID  Checks of the stretched srp grid (build_srp_case with p.grid.h_box_mm), small sizes.
%   test_stretched_grid(0.2, 0.1)
% 1. Grid: spacing range, max cell-to-cell growth, box spacing, exit plane on a node, first row at h_box/2.
% 2. Derivatives on the masked, deformed mesh against analytic ones (fluid points, away from walls):
%    f = sin(3x) cos(2r) (even in r) and g = sin(3x) sin(2r) (odd in r), with x, r in units of D;
%    Dx, D2x, Dy_symm / D2y_symm on f, Dy_asymm / D2y_asymm on g, Dxy_symm on f. Compared with the uniform
%    mesh at h_box: errors should be of the same order (or smaller) in the box.
% 3. Weights: sum(W r) over the fluid vs the uniform mesh (same region, should agree to O(h)).
% 4. Base flow on the stretched grid vs on the uniform grid (same helper grid): Mach difference at shared x.
if nargin < 1, h_out = 0.2; end
if nargin < 2, h_box = 0.1; end
pc = srp_case_defaults(); pc.h_mm = h_out; pc.grid.h_box_mm = h_box; pc.mut_sigma_mm = 0.25;
[ms, Bs] = build_srp_case(pc);
pu = srp_case_defaults(); pu.h_mm = h_box; pu.mut_sigma_mm = 0.25;
[mu, Bu] = build_srp_case(pu);
mm = 1e-3 / pc.D;

% 1. grid
xs = ms.X(1, :) / mm; rs = ms.Y(:, 1) / mm;
dx = diff(xs); dr = diff(rs);
gx = max(abs(dx(2:end) ./ dx(1:end-1) - 1)); gr = max(abs(dr(2:end) ./ dr(1:end-1) - 1));
fprintf('\n1. grid: dx %.4f..%.4f mm, dr %.4f..%.4f mm, max growth x %.3f r %.3f\n', min(dx), max(dx), ...
        min(dr), max(dr), gx, gr);
inb = xs(1:end-1) >= pc.grid.box_x_mm(1) & xs(2:end) <= pc.grid.box_x_mm(2);
fprintf('   box: dx %.6f..%.6f (target %.4f); exit plane node offset %.1e mm; first row r = %.6f mm\n', ...
        min(dx(inb)), max(dx(inb)), h_box, min(abs(xs - pc.x_exit_mm)), rs(1));

% 2. derivatives
for k = 1:2
    if k == 1, m = ms; tag = 'stretched'; else, m = mu; tag = 'uniform  '; end
    u = m.usedInd; x = m.X(u); r = m.Y(u);
    f = sin(3*x) .* cos(2*r); g = sin(3*x) .* sin(2*r);
    far = farfromwall(m, 4);
    e = @(num, ex) max(abs(num(far) - ex(far))) / max(abs(ex(far)));
    fprintf('2. %s: Dx %.1e  D2x %.1e  Dr_s %.1e  D2r_s %.1e  Dr_a %.1e  D2r_a %.1e  Dxr_s %.1e\n', tag, ...
        e(m.Dx*f, 3*cos(3*x).*cos(2*r)), e(m.D2x*f, -9*f), e(m.Dy_symm*f, -2*sin(3*x).*sin(2*r)), ...
        e(m.D2y_symm*f, -4*f), e(m.Dy_asymm*g, 2*sin(3*x).*cos(2*r)), e(m.D2y_asymm*g, -4*g), ...
        e(m.Dxy_symm*f, -6*cos(3*x).*sin(2*r)));
end

% 3. weights
Ws = sum(ms.W(ms.usedInd) .* ms.Y(ms.usedInd)); Wu = sum(mu.W(mu.usedInd) .* mu.Y(mu.usedInd));
fprintf('3. integral of r dx dr over the fluid: stretched %.6f, uniform %.6f (rel diff %.1e)\n', Ws, Wu, ...
        abs(Ws - Wu) / Wu);

% 4. base flow along r = 1.3 mm and r = 5 mm (rows nearest), x 1..20 mm
for rl = [1.3 5]
    [~, is] = min(abs(ms.Y(:,1)/mm - rl)); [~, iu] = min(abs(mu.Y(:,1)/mm - rl));
    Ms = mach(Bs, is, pc); Mu = mach(Bu, iu, pc);
    xq = linspace(1, 20, 2000);
    a = interp1(ms.X(1,:)/mm, Ms, xq); b = interp1(mu.X(1,:)/mm, Mu, xq);
    fprintf('4. Mach on r = %.2f / %.2f mm: mean |dM| %.4f, p99 %.4f (stretched vs uniform)\n', ...
            ms.Y(is,1)/mm, mu.Y(iu,1)/mm, mean(abs(a-b), 'omitnan'), prctile(abs(a-b), 99));
end
fprintf('points: stretched %d fluid (%d DOF), uniform %d fluid\n', ms.ngp, 5*ms.ngp, mu.ngp);
end

function far = farfromwall(m, n)
fl = false(size(m.X)); fl(m.usedInd) = true;
d = bwdist(~fl); d(:, [1:n, end-n+1:end]) = 0; d(end-n+1:end, :) = 0;
far = d(m.usedInd) > n;
end

function M = mach(B, i, pc)
M = hypot(B.U(i,:), B.V(i,:)) ./ sqrt(max(B.T(i,:), eps)) * pc.ref.Ma;
M(B.T(i,:) == 0) = nan;
end
