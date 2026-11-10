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

%% --- grid: uniform, x grid aligned with the nozzle exit plane
h  = p.h_mm * mm;
xe = p.x_exit_mm * mm;
x0 = xe - floor((xe - p.x_range_mm(1)*mm) / h + 1e-9) * h;
x1 = xe + floor((p.x_range_mm(2)*mm - xe) / h + 1e-9) * h;
Nx = round((x1 - x0) / h) + 1;
Rm = p.r_max_mm * mm;
Nr = round(Rm / h + 0.5);                   % CreateMesh with y_symmetry: dr = R/(Nr-0.5), first row at dr/2
cmesh = CreateMesh([x0 x1], [0 Rm], Nx, Nr, 4, true, false, 0);
X = cmesh.X; Y = cmesh.Y;
fprintf('grid: h = %.4f mm, Nx = %d, Nr = %d, x %.3f..%.3f mm, dr = %.4f mm\n', p.h_mm, Nx, Nr, ...
        x0/mm, x1/mm, (Y(2,1)-Y(1,1))/mm);

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
[sensor, shock] = shock_sensor(raw.P, mesh, fluid, sm, h, mm);
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
if p.mut_sigma_cells > 0
    mut = masked_gauss(mut, fluid, p.mut_sigma_cells);
end
mut(solid) = 0;
BF.MU_sutherland = BF.MU;
BF.MU_t = mut;
BF.MU = BF.MU + p.sigma_visc * mut;
for f = {'U','V','W','T','RHO','P','MU','MU_sutherland','MU_t','dmudT','d2mudT2'}
    a = double(BF.(f{1})); a(solid) = 0; BF.(f{1}) = a;
end
BF.sigma_visc = p.sigma_visc;

%% --- sponges (method's GetSponge: tanh ramps at outflow, inflow and top, combined by max)
sp = p.sponge;
right = sp.amp_in/2  * (tanh((X - x1 + sp.width_in_mm*mm) / (sp.trans_in_mm*mm)) + 1);    % freestream inflow, x max
left  = sp.amp_out/2 * (tanh(-(X - x0 - sp.width_out_mm*mm) / (sp.trans_out_mm*mm)) + 1);  % outflow, x min
top   = sp.amp_top/2 * (tanh((Y - Rm + sp.width_top_mm*mm) / (sp.trans_top_mm*mm)) + 1);
sponge = max(left, max(right, top));
sponge(solid) = 0;
mesh.sponge = sponge;

info = struct('p', p, 'raw', raw, 'sensor', sensor, 'shock', shock, 'wblend', wblend, ...
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

function [sensor, shock] = shock_sensor(P, mesh, fluid, sm, h, mm)
pv = P(mesh.usedInd);
g  = hypot(mesh.Dx * pv, mesh.Dy * pv);
sensor = zeros(size(P));
sensor(mesh.usedInd) = g;
sensor = masked_gauss(sensor, fluid, sm.sensor_sigma_mm * mm / h);
sensor(~fluid) = 0;
thr = sm.threshold * max(sensor(fluid));
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
