function [L0, idx] = GetLinProblemVarMu(mesh, BF, m, opts)
% GETLINPROBLEMVARMU  LinStab2D's GetLinProblem for the axisymmetric ('axy') model, with a spatially
% varying mean viscosity (frozen eddy viscosity), following the supervisor's GSA method.
%
%   [L0, idx] = GetLinProblemVarMu(mesh, BF, m, opts)
%
% Same inputs/outputs and conventions as GetLinProblem(mesh, BF, 'axy', m): LHS q = RHS dq/dt,
% L0 = (RHS/1i)\LHS - sponge, q = [rho, u(axial), v(radial), w(azimuthal), T] on return.
%
% Base flow fields used (nondimensional, on mesh.X/mesh.Y):
%   BF.RHO, BF.U (axial), BF.V (radial), BF.W (azimuthal), BF.T
%   BF.MU       total mean viscosity  = MU_Sutherland + sigma_visc*MU_t
%   BF.dmudT    d(mu_Sutherland)/dT   (perturbation mu' = dmudT*T'; mu_t is frozen)
%   BF.d2mudT2  d2(mu_Sutherland)/dT2
%   BF.dMUdx, BF.dMUdr (optional) gradient of BF.MU; default mesh.Dx*MU, mesh.Dy*MU
%   BF.Re, BF.Ma, BF.Pr, BF.kappa, BF.cv, BF.c1, BF.c2 as in LinStab2D
% opts.coeffs: coefficient script (default 'getCoeffsAxyVarMu'; 'getCoeffsAxyMuT' or
%              'getCoeffsLaminar' reproduce upstream LinStab2D).
% opts.alpha_rho: density damping coefficient (default 0; method values 3e-5, 1e-4), see below.

if nargin < 4, opts = struct(); end
if ~isfield(opts, 'coeffs'), opts.coeffs = 'getCoeffsAxyVarMu'; end

fprintf('--- Creating axisymmetric operator (%s), m = %g...', opts.coeffs, m); tic();

% mesh (CS2 -> CS1 as in GetLinProblem: x -> z, y -> r)
z   = mesh.X(mesh.usedInd); %#ok<NASGU>
Dz  = mesh.Dx;  D2z = mesh.D2x;
r   = mesh.Y(mesh.usedInd);
Dr  = mesh.Dy;  D2r = mesh.D2y;
nGridPoints = mesh.ngp;

% base flow, CS1: U radial, V azimuthal, W axial
U   = BF.V(mesh.usedInd);
V   = BF.W(mesh.usedInd);
W   = BF.U(mesh.usedInd);
RHO = BF.RHO(mesh.usedInd);
T   = BF.T(mesh.usedInd);
MU  = BF.MU(mesh.usedInd);
P   = RHO.*T./(BF.Ma^2*BF.kappa); %#ok<NASGU>
dmudT   = getf(BF, 'dmudT',   mesh, zeros(size(MU)));
d2mudT2 = getf(BF, 'd2mudT2', mesh, zeros(size(MU)));
if isfield(BF, 'dMUdr'), dMUdr = BF.dMUdr(mesh.usedInd); else, dMUdr = Dr*MU; end
if isfield(BF, 'dMUdx'), dMUdz = BF.dMUdx(mesh.usedInd); else, dMUdz = Dz*MU; end

cv = BF.cv; c1 = BF.c1; c2 = BF.c2; kappa = BF.kappa; Re = BF.Re; %#ok<NASGU>

Z = zeros(nGridPoints, 1);
I = ones(nGridPoints, 1); %#ok<NASGU>
R = reshape(r, nGridPoints, 1);

% base-flow derivatives (as GetLinProblem)
dUdr   = Dr*U;   d2Udr2   = D2r*U;
dVdr   = Dr*V;   d2Vdr2   = D2r*V;
dWdr   = Dr*W;   d2Wdr2   = D2r*W;
dRHOdr = Dr*RHO; d2RHOdr2 = D2r*RHO; %#ok<NASGU>
dTdr   = Dr*T;   d2Tdr2   = D2r*T;   %#ok<NASGU>
dUdz   = Dz*U;   d2Udz2   = D2z*U;
dVdz   = Dz*V;   d2Vdz2   = D2z*V;
dWdz   = Dz*W;   d2Wdz2   = D2z*W;
dRHOdz = Dz*RHO; d2RHOdz2 = D2z*RHO; %#ok<NASGU>
dTdz   = Dz*T;   d2Tdz2   = D2z*T;   %#ok<NASGU>
d2Udrz = Dz*dUdr;
d2Vdrz = Dz*dVdr;
d2Wdrz = Dz*dWdr;
d2RHOdrz = Dz*dRHOdr; d2Tdrz = Dz*dTdr; %#ok<NASGU>

[DR, D2R, DZ, D2Z, D2RZ] = CreateDiffMatrices_Axy(mesh, m);

% coefficients (script sets A0_ij, Ar_ij, Az_ij, Arr_ij, Azz_ij, Arz_ij, B_ij)
eval(opts.coeffs);

% assemble the 5x5 block-diagonal coefficient matrices (eval: names come from the script)
C = struct();
for p = {'A0','Ar','Az','Arr','Azz','Arz','B'}
    rows = cell(5, 1);
    for i = 1:5
        cols = cell(1, 5);
        for j = 1:5
            cols{j} = spdiags(full(eval(sprintf('%s_%d%d', p{1}, i, j))), 0, nGridPoints, nGridPoints);
        end
        rows{i} = [cols{:}];
    end
    C.(p{1}) = vertcat(rows{:});
end
LHS = C.A0 + C.Ar*DR + C.Az*DZ + C.Arr*D2R + C.Azz*D2Z + C.Arz*D2RZ;
RHS = C.B;

idx = mesh.idx;
ngp = mesh.ngp;
boundaries = fields(mesh.idx);
variables  = {'rho','u','v','w','T'};
for ib = 1:length(boundaries)
    b = boundaries{ib};
    idx.([b '_all']) = [];
    for iv = 1:length(variables)
        vv = variables{iv};
        idx.([b '_' vv]) = idx.(b) + ngp*(iv-1);
        idx.([b '_all']) = [idx.([b '_all']); idx.([b '_' vv])];
    end
end
for iv = 1:length(variables)
    idx.([variables{iv} '_j']) = (1:ngp) + ngp*(iv-1);
end

L0 = (RHS/1i)\LHS;

if isfield(mesh, 'sponge')
    sponge  = mesh.sponge(mesh.usedInd);
    L0 = L0 - spdiags(repmat(sponge(:), 5, 1), 0, ngp*5, ngp*5);
end

% density damping (supervisor's alpha_rho term, extended to both directions): added to the rho' equation
%   D_rho = -alpha_rho * ( hz^3 H^-1 D2z' H D2z + hr^3 H^-1 D2r' H D2r ),  H = diag(W r)  (axisymmetric weights)
% so that rho'* H D_rho rho' = -alpha_rho (hz^3 |D2z rho'|^2_H + hr^3 |D2r rho'|^2_H) <= 0 (dissipative).
% D2r is the parity-correct r second derivative of a scalar for this m (block 1 of D2R). The method applies the
% y term only; the scaling h^3 makes the damping vanish as h^-1 at grid scale and as h^3 k^4 for resolved waves.
if isfield(opts, 'alpha_rho') && opts.alpha_rho ~= 0
    a  = opts.alpha_rho;
    hz = median(diff(unique(round(mesh.X(mesh.usedInd), 12))));
    hr = median(diff(unique(round(mesh.Y(mesh.usedInd), 12))));
    Hv = mesh.W(mesh.usedInd) .* mesh.Y(mesh.usedInd); Hv = Hv(:);
    assert(all(isfinite(Hv) & Hv > 0), 'GetLinProblemVarMu: invalid quadrature weights for alpha_rho');
    H  = spdiags(Hv, 0, ngp, ngp); Hi = spdiags(1./Hv, 0, ngp, ngp);
    D2rs = D2R(1:ngp, 1:ngp);
    Drho = -a * (hz^3 * (Hi * (D2z' * H * D2z)) + hr^3 * (Hi * (D2rs' * H * D2rs)));
    L0(1:ngp, 1:ngp) = L0(1:ngp, 1:ngp) + Drho;
    % report the damping rate of a grid-scale checkerboard and of a smooth field (Rayleigh quotients in H)
    [iy, ix] = ind2sub(size(mesh.X), mesh.usedInd(:));
    rq = @(q) -real(q' * H * (Drho * q)) / real(q' * H * q);
    qc = (-1).^(ix + iy);
    qs = cos(2*pi*mesh.X(mesh.usedInd(:)) / 0.11);          % axial wavelength 0.11 D (~2 mm)
    fprintf('\n--- alpha_rho = %.1e (hz %.3e, hr %.3e): damping rate checkerboard %.3f, 2 mm wave %.2e', ...
            a, hz, hr, rq(qc), rq(qs));
end

% back to the mesh CS: q = [rho, u(axial), v(radial), w(azimuthal), T]
II = speye(size(Dz)); ZZ = II*0;
Tm = [II, ZZ, ZZ, ZZ, ZZ;
      ZZ, ZZ, ZZ, II, ZZ;
      ZZ, II, ZZ, ZZ, ZZ;
      ZZ, ZZ, II, ZZ, ZZ;
      ZZ, ZZ, ZZ, ZZ, II];
L0 = Tm*L0*Tm';
fprintf(' Done in %.0f seconds.\n', toc);
end

function v = getf(BF, name, mesh, default)
if isfield(BF, name), v = BF.(name)(mesh.usedInd); else, v = default; end
end
