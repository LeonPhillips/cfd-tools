function [mesh, BF, info] = build_srp_case(p)
% BUILD_SRP_CASE  LinStab2D grid, masked body, boundary sets, sponges and base flow for the srp_M2.5 GSA,
% following the supervisor's method (uniform masked grid, local shock smoothing, effective viscosity).
%
%   [mesh, BF, info] = build_srp_case(p)      p = srp_case_defaults(); edit; call.
%
% Coordinates are nondimensionalised by the body diameter D; velocities by U_inf, density by rho_inf,
% temperature by T_inf, viscosities by mu_inf (freestream, from the SU2 t3 config).
% mesh.X = axial x/D, mesh.Y = radial r/D (LinStab2D 'axy' with y_symmetry: first row at dr/2).
% BF.U axial, BF.V radial, BF.W = 0 (azimuthal), BF.T, BF.RHO, BF.MU = MU_Sutherland + sigma_visc*MU_t.
% Boundary sets in mesh.idx: ri (freestream inflow, x max), li (outflow, x min), ti (top, r max),
%   ji (jet inflow: nozzle exit plane, r < r_exit), wi (no-slip walls: the rest of the masked edge).

t0 = tic;
D = p.D;
mm = 1e-3 / D;                              % 1 mm in units of D

%% --- grid: x grid aligned with the nozzle exit plane; uniform, or stretched (fine box, see srp_case_defaults)
xe = p.x_exit_mm * mm;
Rm = p.r_max_mm * mm;
stretched = isfield(p, 'grid') && isfield(p.grid, 'h_box_mm') && ~isempty(p.grid.h_box_mm) && ...
            abs(p.grid.h_box_mm - p.h_mm) > 1e-12;
if ~stretched
    h  = p.h_mm * mm;
    x0 = xe - floor((xe - p.x_range_mm(1)*mm) / h + 1e-9) * h;
    x1 = xe + floor((p.x_range_mm(2)*mm - xe) / h + 1e-9) * h;
    Nx = round((x1 - x0) / h) + 1;
    Nr = round(Rm / h + 0.5);               % CreateMesh with y_symmetry: dr = R/(Nr-0.5), first row at dr/2
    cmesh = CreateMesh([x0 x1], [0 Rm], Nx, Nr, 4, true, false, 0);
    X = cmesh.X; Y = cmesh.Y;
    fprintf('grid: h = %.4f mm, Nx = %d, Nr = %d, x %.3f..%.3f mm, dr = %.4f mm\n', p.h_mm, Nx, Nr, ...
            x0/mm, x1/mm, (Y(2,1)-Y(1,1))/mm);
else
    % computational grid in index units (dxi = deta = 1, first row at eta = 1/2), mapped to the physical
    % grid by 1-D maps x(xi), r(eta) (applied with LinStab2D's DeformMesh after masking)
    g  = p.grid; hb = g.h_box_mm * mm; ho = p.h_mm * mm;
    xs = stretched_axis(xe, p.x_range_mm*mm, hb, ho, g.box_x_mm*mm, g.max_growth, false);
    rs = stretched_axis(hb/2, [0 Rm], hb, ho, [0 g.box_r_mm*mm], g.max_growth, true);
    Nx = numel(xs); Nr = numel(rs);
    cmesh = CreateMesh([0 Nx-1], [0 Nr-0.5], Nx, Nr, 4, true, false, 0);
    X = repmat(xs(:)', Nr, 1); Y = repmat(rs(:), 1, Nx);
    x0 = xs(1); x1 = xs(end);
    h  = hb;                                % smallest spacing (tolerances and wall tests below)
    fprintf(['grid (stretched): box %.4f mm in x %.2f..%.2f, r 0..%.2f mm; %.4f mm outside; growth <= %.0f%%;\n' ...
             '  Nx = %d, Nr = %d (%d points), x %.3f..%.3f mm, r max %.3f mm\n'], g.h_box_mm, g.box_x_mm, ...
            g.box_r_mm, p.h_mm, 100*g.max_growth, Nx, Nr, Nx*Nr, x0/mm, x1/mm, rs(end)/mm);
end

%% --- mask: body cylinder behind the face, plus the lip/nozzle block up to the exit plane (solid = true)
tol = 1e-6 * h;
solid = (X < p.x_face_mm*mm & Y < p.r_body_mm*mm) | (X < xe - tol & Y < p.r_lip_out_mm*mm);
use_nozzle = isfield(p, 'nozzle_interior') && p.nozzle_interior;
if use_nozzle
    % open the nozzle bore (plenum inlet x_inlet .. exit plane) using the wall contour from the mesh
    rin = nozzle_inner_radius(p.baseflow_prefix, X(1,:) / mm, p.x_inlet_mm, p.x_exit_mm) * mm;   % per column
    bore = X >= p.x_inlet_mm*mm - tol & X < xe - tol & Y < repmat(rin, size(X,1), 1);
    solid = solid & ~bore;
    solid = clean_mask(solid);              % no fluid runs shorter than 3 points (FD stencils)
end
mesh  = MeshMask(cmesh, double(solid));
mesh  = fix_offaxis_symm(mesh);            % LinStab2D bug: see the function below
if stretched
    mesh = DeformMesh(mesh, X, []);         % x map, then r map (1-D branches: exact metric terms)
    mesh = DeformMesh(mesh, [], Y);
    mesh.X = X; mesh.Y = Y;
end
mesh.mask = solid;
Xu = X(mesh.usedInd); Yu = Y(mesh.usedInd);

% boundary sets
mi = mesh.idx.mi;
if use_nozzle
    % plenum inlet face: the first open column of the bore (subsonic total-conditions inlet in SU2)
    xcol = min(Xu(Yu < p.r_lip_out_mm*mm & Xu < p.x_face_mm*mm));
    ji = mi(abs(Xu(mi) - xcol) < tol & Yu(mi) < p.r_plenum_mm*mm - h);
    jname = 'ni';
else
    ji = mi(abs(Xu(mi) - xe) < tol & Yu(mi) < p.r_exit_mm*mm);
    jname = 'ji';
end
wi = setdiff(mi, ji);
idx = struct('li', mesh.idx.li, 'ri', mesh.idx.ri, 'ti', mesh.idx.ti, 'wi', wi(:));
idx.(jname) = ji(:);
idx.li = setdiff(idx.li, [ji; wi]);         % edges that are also walls stay walls
idx.ti = setdiff(idx.ti, idx.ri);           % corner (x max, r max) belongs to the inflow
mesh.idx = idx;
mesh.inflow_set = jname(1);                 % 'n' (plenum inlet) or 'j' (exit plane)
fprintf('mask: %d fluid points (%d masked); %s %d, wi %d, li %d, ri %d, ti %d\n', mesh.ngp, nnz(solid), ...
        jname, numel(ji), numel(idx.wi), numel(idx.li), numel(idx.ri), numel(idx.ti));

% wall normals for Neumann conditions (axis-aligned, pointing into the fluid, from the 4 neighbours)
mesh.wall_normal = wall_normals(solid, mesh.usedInd(idx.wi));

% body-side slip wall (user, 2026-11-10): the cylinder r = r_body behind the face (normal +r, x < x_face)
use_slip = isfield(p, 'bodyside_slip') && p.bodyside_slip;
if use_slip
    wi = mesh.idx.wi; n = mesh.wall_normal;
    s = n(:,1) == 0 & n(:,2) == 1 & Xu(wi) < p.x_face_mm*mm - tol & Yu(wi) > p.r_body_mm*mm - h;
    mesh.idx.si = wi(s);
    mesh.idx.wi = wi(~s);
    mesh.slip_normal = n(s, :);
    mesh.wall_normal = n(~s, :);
    fprintf('body-side slip wall: %d points (no-slip walls left: %d)\n', numel(mesh.idx.si), numel(mesh.idx.wi));
end

use_aux = stretched || (isfield(p, 'bf_h_mm') && ~isempty(p.bf_h_mm) && abs(p.bf_h_mm - p.h_mm) > 1e-12);
if use_aux
    % base flow on a uniform helper grid (spacing bf_h_mm; default the box spacing), then interpolated here
    [BF, raw, sensor, shock, wblend] = baseflow_from_aux(p, mesh, X, Y, solid);
    thr = NaN;
    fluid = ~solid;
else
%% --- base flow: scattered SU2 nodes -> grid (natural neighbour, as the method)
S = read_bin(p.baseflow_prefix);
ref = p.ref;
xs = S.x / D; rs = S.r / D;
fields = struct('U', S.Velocity_x / ref.U, 'V', S.Velocity_r / ref.U, 'T', S.Temperature / ref.T, ...
                'RHO', S.Density / ref.rho, 'P', S.Pressure / (ref.rho * ref.U^2), ...
                'MUt', S.Eddy_Viscosity / ref.mu, 'MUlam', S.Laminar_Viscosity / ref.mu, 'Mach', S.Mach);
F = scatteredInterpolant(xs, rs, fields.U, 'natural', 'none');
BF = struct();
raw = struct();
for f = fieldnames(fields)'
    F.Values = fields.(f{1});
    g = nan(size(X));
    g(mesh.usedInd) = F(Xu, Yu);
    nbad = nnz(isnan(g(mesh.usedInd)));
    if nbad > 0
        Fn = scatteredInterpolant(xs, rs, fields.(f{1}), 'nearest', 'nearest');
        b = mesh.usedInd(isnan(g(mesh.usedInd)));
        g(b) = Fn(X(b), Y(b));
        fprintf('  %s: %d points outside the CFD hull filled by nearest neighbour\n', f{1}, nbad);
    end
    raw.(f{1}) = g;
end
fprintf('interpolated %d fields at %.0f s\n', numel(fieldnames(fields)), toc(t0));

% remove the body-side boundary layer from the base flow (with the slip wall): in each column behind the face,
% every point within bodyside_bl_mm of the cylinder takes the values at that distance
if use_slip
    rb = p.r_body_mm*mm; dbl = p.bodyside_bl_mm*mm;
    cols = find(X(1,:) < p.x_face_mm*mm - tol);
    rows = find(Y(:,1) > rb & Y(:,1) < rb + dbl);
    rref = find(Y(:,1) >= rb + dbl, 1);
    for f = fieldnames(raw)'
        a = raw.(f{1});
        a(rows, cols) = repmat(a(rref, cols), numel(rows), 1);
        raw.(f{1}) = a;
    end
    fprintf('body-side boundary layer replaced: %d rows x %d columns (r < %.2f mm)\n', numel(rows), ...
            numel(cols), (rb + dbl)/mm);
end

%% --- local shock smoothing (method's smooth_shock, with a masked Gaussian and sigmas in mm)
sm = p.smooth;
fluid = ~solid;
[sensor, shock, thr] = shock_sensor(raw.P, mesh, fluid, sm, h, mm);
dist = double(bwdist(shock)) * h / mm;      % mm (bwdist returns single)
wblend = exp(-dist.^2 / (2 * sm.blend_sigma_mm^2));
wblend(solid) = 0;
sig = sm.filter_sigma_mm * mm / h;          % in cells
for f = {'U', 'V', 'T', 'RHO', 'P'}
    a = raw.(f{1});
    s = masked_gauss(a, fluid, sig);
    BF.(f{1}) = a .* (1 - wblend) + s .* wblend;
end
fprintf('shock smoothing: %d sensor points, filter sigma %.3f mm, blend sigma %.3f mm\n', ...
        nnz(shock), sm.filter_sigma_mm, sm.blend_sigma_mm);

%% --- viscosity: Sutherland + sigma_visc * mu_t (frozen)
BF.W = zeros(size(X));
BF.Re = ref.Re; BF.Ma = ref.Ma; BF.Pr = ref.Pr; BF.kappa = ref.kappa; BF.T_0 = ref.T;
BF = sutherland_air(BF);                    % MU, dmudT, d2mudT2 (LinStab2D; also sets cv, c1, c2)
mut = max(raw.MUt, 0);
sc = p.mut_sigma_cells;
if isfield(p, 'mut_sigma_mm') && ~isempty(p.mut_sigma_mm), sc = p.mut_sigma_mm * mm / h; end
if sc > 0
    mut = masked_gauss(mut, fluid, sc);
end
mut(solid) = 0;
BF.MU_sutherland = BF.MU;
BF.MU_t = mut;
BF.MU = BF.MU + p.sigma_visc * mut;
for f = {'U','V','W','T','RHO','P','MU','MU_sutherland','MU_t','dmudT','d2mudT2'}
    a = double(BF.(f{1})); a(solid) = 0; BF.(f{1}) = a;
end
BF.sigma_visc = p.sigma_visc;
end

%% --- sponges (method's GetSponge: tanh ramps at outflow, inflow and top, combined by max)
sp = p.sponge;
right = sp.amp_in/2  * (tanh((X - x1 + sp.width_in_mm*mm) / (sp.trans_in_mm*mm)) + 1);    % freestream inflow, x max
left  = sp.amp_out/2 * (tanh(-(X - x0 - sp.width_out_mm*mm) / (sp.trans_out_mm*mm)) + 1);  % outflow, x min
top   = sp.amp_top/2 * (tanh((Y - Rm + sp.width_top_mm*mm) / (sp.trans_top_mm*mm)) + 1);
sponge = max(left, max(right, top));
sponge(solid) = 0;
mesh.sponge = sponge;

info = struct('p', p, 'raw', raw, 'sensor', sensor, 'shock', shock, 'wblend', wblend, 'sensor_thr', thr, ...
              'x_range', [x0 x1], 'r_max', Rm, 'h', h, 'Nx', Nx, 'Nr', Nr, 'time_s', toc(t0));
fprintf('case built in %.0f s: %d fluid points, %d DOF\n', toc(t0), mesh.ngp, 5*mesh.ngp);
end

% ------------------------------------------------------------------------------------------------
function rin = nozzle_inner_radius(prefix, xcols_mm, x_inlet_mm, x_exit_mm)
% Inner radius of the nozzle/plenum at each grid column (mm), from the wall_noslip segments: the smallest r
% at which a wall segment crosses x = const (the lip outer surface is larger, so min picks the inner wall).
m = jsondecode(fileread([prefix '_wall.json']));
fid = fopen([prefix '_wall.bin'], 'r'); a = fread(fid, [m.N, 4], 'double'); fclose(fid);
a = a * 1e3;                                % mm: x0 r0 x1 r1
rin = zeros(size(xcols_mm));
for j = 1:numel(xcols_mm)
    x = xcols_mm(j);
    if x < x_inlet_mm - 1e-9 || x >= x_exit_mm, continue; end
    lo = min(a(:,1), a(:,3)); hi = max(a(:,1), a(:,3));
    s = find(lo <= x & hi >= x & max(a(:,2), a(:,4)) < 3.2);
    if isempty(s), continue; end
    t = (x - a(s,1)) ./ (a(s,3) - a(s,1)); t(~isfinite(t)) = 0;
    r = a(s,2) + t .* (a(s,4) - a(s,2));
    rin(j) = min(r);
end
end

function solid = clean_mask(solid)
% make fluid runs shorter than 3 points (in x or r) solid, so every FD segment has >= 3 points
for it = 1:20
    old = solid;
    for dim = 1:2
        F = ~solid; if dim == 2, F = F.'; end
        for i = 1:size(F, 1)
            d = diff([false, F(i,:), false]);
            s = find(d == 1); e = find(d == -1) - 1;
            for k = find(e - s + 1 < 3), F(i, s(k):e(k)) = false; end
        end
        if dim == 2, F = F.'; end
        solid = ~F;
    end
    if isequal(old, solid), break; end
end
end

function S = read_bin(prefix)
m = jsondecode(fileread([prefix '_flow.json']));
fid = fopen([prefix '_flow.bin'], 'r');
a = fread(fid, [m.N, numel(m.names)], 'double');
fclose(fid);
S = struct();
for k = 1:numel(m.names), S.(m.names{k}) = a(:, k); end
end

function s = masked_gauss(a, fluid, sigma_cells)
% Gaussian filter of a over the fluid points only (normalised convolution: no leakage from the solid)
a0 = a; a0(~fluid | isnan(a0)) = 0;
w  = double(fluid & ~isnan(a));
num = imgaussfilt(a0, sigma_cells, 'Padding', 'replicate');
den = imgaussfilt(w,  sigma_cells, 'Padding', 'replicate');
s = num ./ max(den, eps);
s(~fluid) = a(~fluid);
end

function [sensor, shock, thr] = shock_sensor(P, mesh, fluid, sm, h, mm)
pv = P(mesh.usedInd);
g  = hypot(mesh.Dx * pv, mesh.Dy * pv);
sensor = zeros(size(P));
sensor(mesh.usedInd) = g;
sensor = masked_gauss(sensor, fluid, sm.sensor_sigma_mm * mm / h);
sensor(~fluid) = 0;
if isfield(sm, 'threshold_abs') && ~isempty(sm.threshold_abs), thr = sm.threshold_abs;   % set by baseflow_from_aux
else, thr = sm.threshold * max(sensor(fluid)); end
shock = fluid & sensor > thr & mesh.X > sm.x_min_mm*mm & mesh.Y > sm.r_min_mm*mm;
end

function n = wall_normals(solid, ind)
% unit normal into the fluid at each wall point (linear indices into the full grid), from which of the
% 4 neighbours (and, at convex corners, the diagonals) are solid
[ny, nx] = size(solid);
[iy, ix] = ind2sub([ny nx], ind(:));
S = @(i, j) solid(sub2ind([ny nx], min(max(i,1),ny), min(max(j,1),nx)));
nxv = double(S(iy, ix-1)) - double(S(iy, ix+1));      % solid on the left -> normal +x
nyv = double(S(iy-1, ix)) - double(S(iy+1, ix));      % solid below -> normal +r
dg = (nxv == 0 & nyv == 0);                            % only a diagonal neighbour is solid (convex corner)
nxv(dg) = double(S(iy(dg)-1, ix(dg)-1) | S(iy(dg)+1, ix(dg)-1)) - double(S(iy(dg)-1, ix(dg)+1) | S(iy(dg)+1, ix(dg)+1));
nyv(dg) = double(S(iy(dg)-1, ix(dg)-1) | S(iy(dg)-1, ix(dg)+1)) - double(S(iy(dg)+1, ix(dg)-1) | S(iy(dg)+1, ix(dg)+1));
nrm = hypot(nxv, nyv); nrm(nrm == 0) = 1;
n = [nxv ./ nrm, nyv ./ nrm];
end

function xs = stretched_axis(xa, lim, hb, ho, box, gmax, from_start)
% Grid coordinates along one axis: spacing hb inside box = [lo hi], growing smoothly (quintic smoothstep in
% position, cell-to-cell growth <= gmax) to ho outside. One node sits at xa (exit plane, or r = hb/2 with
% from_start = true: then nodes run from xa upwards only). Nodes stay within lim (the last one <= lim(2)).
L = 1.875 * (ho - hb) / (gmax * hb);         % transition length (units of hb) for max d(spacing)/dx = gmax
L = max(L, 1) * hb;
S = @(t) (t <= 0) .* 0 + (t >= 1) .* 1 + (t > 0 & t < 1) .* (6*t.^5 - 15*t.^4 + 10*t.^3);
sp = @(x) hb + (ho - hb) * max(S((x - box(2)) / L), S((box(1) - x) / L));
up = xa; x = xa;
while true
    s = sp(x + 0.5 * sp(x));
    if x + s > lim(2) + 1e-9 * hb, break; end
    x = x + s; up(end+1) = x; %#ok<AGROW>
end
if from_start, xs = up; return; end
dn = []; x = xa;
while true
    s = sp(x - 0.5 * sp(x));
    if x - s < lim(1) - 1e-9 * hb, break; end
    x = x - s; dn(end+1) = x; %#ok<AGROW>
end
xs = [fliplr(dn), up];
end

function [BF, raw, sensor, shock, wblend] = baseflow_from_aux(p, mesh, X, Y, solid)
% Build the base flow (interpolation, body-side layer, shock smoothing, mu_t filter) on two uniform helper
% grids, then interpolate it (linear) to this grid's points:
%  - fine:   spacing p.bf_h_mm (default: the box spacing) over the refined box p.grid.box_x_mm / box_r_mm plus
%            a margin (the SU2 refinement region), and
%  - coarse: spacing h_mm (this grid's outer spacing) over the whole domain.
% Points inside the box take the fine values and points outside the coarse ones, blended linearly over the
% last blend_mm inside the box edge. The shock-sensor threshold is set on the fine helper (relative to its
% max, as before) and passed to the coarse one as an absolute value, so both flag the same shocks.
% (One uniform fine helper over the whole domain is ~4M points at 0.0177 mm: it took > 2.5 h to build.)
% Values inside a helper's solid are first filled from its nearest fluid point, so wall-adjacent points never
% see solid values. Sutherland viscosity is recomputed here from the interpolated T.
mm = 1e-3 / p.D;
hf = p.h_mm;
if isfield(p, 'bf_h_mm') && ~isempty(p.bf_h_mm), hf = p.bf_h_mm;
elseif ~isempty(p.grid.h_box_mm),               hf = p.grid.h_box_mm; end
margin_mm = 1.5; blend_mm = 0.5;
bx = p.grid.box_x_mm; br = p.grid.box_r_mm;
pf = p; pf.h_mm = hf; pf.grid.h_box_mm = []; pf.bf_h_mm = [];
pf.x_range_mm = [max(p.x_range_mm(1), bx(1) - margin_mm), min(p.x_range_mm(2), bx(2) + margin_mm)];
pf.r_max_mm   = min(p.r_max_mm, br + margin_mm);
fprintf('--- base flow, fine helper grid: h = %.4f mm over x %.2f..%.2f, r 0..%.2f mm\n', hf, ...
        pf.x_range_mm, pf.r_max_mm);
[mF, BFf, iF] = build_srp_case(pf);
pc = p; pc.grid.h_box_mm = []; pc.bf_h_mm = [];
pc.smooth.threshold_abs = iF.sensor_thr;
fprintf('--- base flow, coarse helper grid: h = %.4f mm over the whole domain (sensor threshold %.4g from the fine grid)\n', ...
        p.h_mm, iF.sensor_thr);
[mC, BFc, iC] = build_srp_case(pc);

% weight of the fine helper at this grid's points (1 inside the box, 0 outside, linear over blend_mm)
u = mesh.usedInd; Xu = X(u); Yu = Y(u);
d = min(Xu/mm - bx(1), br - Yu/mm);
if bx(2) < p.x_range_mm(2) - 1e-9, d = min(d, bx(2) - Xu/mm); end
wf = min(max(d / blend_mm, 0), 1);
kf = wf > 0;
Fi = helper_interp(mF, Xu(kf), Yu(kf)); Ci = helper_interp(mC, Xu, Yu);
mix = @(af, ac) mix_fields(Fi(af, 'linear'), Ci(ac, 'linear'), wf, kf, size(X), u);
BF = struct();
for f = {'U', 'V', 'T', 'RHO', 'P', 'MU_t'}
    BF.(f{1}) = mix(BFf.(f{1}), BFc.(f{1}));
end
raw = struct();
for f = fieldnames(iF.raw)'
    raw.(f{1}) = mix(iF.raw.(f{1}), iC.raw.(f{1}));
end
wblend = mix(iF.wblend, iC.wblend);
sensor = mix(double(iF.sensor), double(iC.sensor));
sf = Fi(double(iF.shock), 'nearest'); sc = Ci(double(iC.shock), 'nearest');
sh = sc; sh(kf) = sf .* (wf(kf) >= 0.5) + sc(kf) .* (wf(kf) < 0.5);
shock = false(size(X)); shock(u) = sh > 0.5;
fprintf('--- base flow interpolated from the helper grids (fine %d, coarse %d points) to %d points (%d in the box)\n', ...
        mF.ngp, mC.ngp, mesh.ngp, nnz(kf));
ref = p.ref;
BF.W = zeros(size(X));
BF.Re = ref.Re; BF.Ma = ref.Ma; BF.Pr = ref.Pr; BF.kappa = ref.kappa; BF.T_0 = ref.T;
BF.T(solid) = 1;                              % keep Sutherland finite in the solid (zeroed below)
BF = sutherland_air(BF);
BF.MU_sutherland = BF.MU;
BF.MU = BF.MU + p.sigma_visc * BF.MU_t;
for f = {'U','V','W','T','RHO','P','MU','MU_sutherland','MU_t','dmudT','d2mudT2'}
    a = double(BF.(f{1})); a(solid) = 0; BF.(f{1}) = a;
end
BF.sigma_visc = p.sigma_visc;
end

function mesh = fix_offaxis_symm(mesh)
% Fix for LinStab2D's AddDW2mesh (eefce30) on masked axisymmetric meshes. It decides whether a radial fluid
% segment touches the axis from the first point of the whole column (yl(1)), not of the segment. So:
%  - in columns whose first fluid point is off the axis (e.g. behind the body, r > r_body), Dy_symm, Dy_asymm,
%    D2y_symm and D2y_asymm are left ZERO: no radial derivatives at all in the operator there;
%  - in columns that do start on the axis, every other segment (e.g. above the body, behind the nozzle bore)
%    gets the axis-reflection closure at its lower end, as if that wall were the axis.
% Off-axis segments need no parity: their rows are replaced by the plain Dy / D2y rows, and the mixed
% derivatives are rebuilt with LinStab2D's own definitions (Dxy_symm = Dy_symm*Dx, Dyx_symm = Dx*Dy_symm).
if ~isfield(mesh, 'Dy_symm'), return; end
[ny, nx] = size(mesh.X);
fl = false(ny, nx); fl(mesh.usedInd) = true;
onaxis = false(ny, nx);
for j = 1:nx
    k = 1;
    while k <= ny && fl(k, j), onaxis(k, j) = true; k = k + 1; end
end
off = ~onaxis(mesh.usedInd);
P = spdiags(double(off), 0, mesh.ngp, mesh.ngp); Q = spdiags(double(~off), 0, mesh.ngp, mesh.ngp);
mesh.Dy_symm   = Q * mesh.Dy_symm   + P * mesh.Dy;
mesh.Dy_asymm  = Q * mesh.Dy_asymm  + P * mesh.Dy;
mesh.D2y_symm  = Q * mesh.D2y_symm  + P * mesh.D2y;
mesh.D2y_asymm = Q * mesh.D2y_asymm + P * mesh.D2y;
mesh.Dxy_asymm = mesh.Dy_asymm * mesh.Dx;  mesh.Dyx_asymm = mesh.Dx * mesh.Dy_asymm;
mesh.Dxy_symm  = mesh.Dy_symm  * mesh.Dx;  mesh.Dyx_symm  = mesh.Dx * mesh.Dy_symm;
fprintf('off-axis radial segments: %d points get plain Dy/D2y in the (a)symmetric matrices\n', nnz(off));
end

function F = helper_interp(mA, Xq, Yq)
% F(a, method): helper-grid field a (solid filled from its nearest fluid point) at the query points
xa = mA.X(1, :); ra = mA.Y(:, 1);
[~, L] = bwdist(~mA.mask);
F = @(a, method) interp_at(griddedInterpolant({ra, xa}, double(a(L)), method, 'nearest'), Yq, Xq);
end

function v = interp_at(G, Yq, Xq)
v = G(Yq, Xq);
end

function g = mix_fields(vf, vc, wf, kf, sz, u)
v = vc; v(kf) = wf(kf) .* vf + (1 - wf(kf)) .* vc(kf);
g = zeros(sz); g(u) = v;
end
