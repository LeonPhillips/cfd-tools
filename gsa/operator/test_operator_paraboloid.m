function test_operator_paraboloid(repo, outdir)
% TEST_OPERATOR_PARABOLOID  Operator-level checks of GetLinProblemVarMu on LinStab2D's
% example_5_Paraboloid (axisymmetric, Sutherland, curvilinear mesh), for m = 0 and m = 5.
%
%  A  coeffs 'getCoeffsLaminar' (upstream's script through our assembly) == GetLinProblem 'axy'
%     -> checks the wiring (assembly, CS change, idx, sponge): should be identical
%  B  'getCoeffsAxyMuT' == upstream, with dmudT = 0      -> identical (the only difference is mu'*div(tau_bar))
%  C  'getCoeffsAxyVarMu', MU = Sutherland, dMU from chain rule == 'getCoeffsAxyMuT'   -> identical
%  D  'getCoeffsAxyVarMu', MU = Sutherland (dMU by FD of MU) vs upstream: size of the change
%     (full mu' terms + FD-of-MU vs chain rule), reported only
%  E  eddy viscosity: MU + mu_t with a smooth bump; residual-free sanity: L0 changes only where mu_t /= 0

if ~exist(outdir, 'dir'), mkdir(outdir); end
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(genpath(fullfile(repo, 'aux_matlab')));

% build mesh + base flow with the example's own code, up to GetLinProblem
code = fileread(fullfile(repo, 'example_5_Paraboloid.m'));
cut  = regexp(code, '(?m)^.*GetLinProblem\(', 'start', 'once');
assert(~isempty(cut), 'GetLinProblem call not found in example_5');
code = code(1:cut-1);
rx = @(s, p, r) regexprep(s, p, r, 'lineanchors', 'dotexceptnewline');
code = rx(code, '^[ \t]*(clear all|clear variables|close all)[ \t]*;?', '');
code = rx(code, '^verbose[ \t]*=.*$', 'verbose = false;');
fname = 'parab_setup';
fid = fopen(fullfile(outdir, [fname '.m']), 'w');
fprintf(fid, 'function [mesh, baseFlow] = %s()\n%s\nend\n', fname, code);
fclose(fid);
addpath(outdir); rehash;
old = cd(repo); restore = onCleanup(@() cd(old));
[mesh, BF] = feval(fname);
clear restore

rel = @(A, B) full(max(max(abs(A - B)))) / full(max(max(abs(B))));
fid = fopen(fullfile(outdir, 'test_operator_paraboloid.txt'), 'w');
pr = @(varargin) [fprintf(1, varargin{:}), fprintf(fid, varargin{:})];
pr('paraboloid: ngp = %d, Re = %g, Ma = %g\n', mesh.ngp, BF.Re, BF.Ma);

for m = [0 5]
    L_up = GetLinProblem(mesh, BF, 'axy', m);
    L_A  = GetLinProblemVarMu(mesh, BF, m, struct('coeffs', 'getCoeffsLaminar'));
    BF0 = BF; BF0.dmudT = 0*BF.dmudT;
    L_up0 = GetLinProblem(mesh, BF0, 'axy', m);
    L_B   = GetLinProblemVarMu(mesh, BF0, m, struct('coeffs', 'getCoeffsAxyMuT'));
    L_MT  = GetLinProblemVarMu(mesh, BF, m, struct('coeffs', 'getCoeffsAxyMuT'));
    BFc = BF;
    BFc.dMUdr = zeros(size(BF.MU)); BFc.dMUdx = zeros(size(BF.MU));
    BFc.dMUdr(mesh.usedInd) = BF.dmudT(mesh.usedInd) .* (mesh.Dy*BF.T(mesh.usedInd));
    BFc.dMUdx(mesh.usedInd) = BF.dmudT(mesh.usedInd) .* (mesh.Dx*BF.T(mesh.usedInd));
    L_C   = GetLinProblemVarMu(mesh, BFc, m, struct('coeffs', 'getCoeffsAxyVarMu'));
    L_D   = GetLinProblemVarMu(mesh, BF, m);
    pr('m = %d: A wiring |L - L_up|/|L_up| = %.1e\n', m, rel(L_A, L_up));
    pr('       B mu''=0       = %.1e\n', rel(L_B, L_up0));
    pr('       C VarMu(chain rule) vs MuT = %.1e\n', rel(L_C, L_MT));
    pr('       D VarMu(FD of MU) vs upstream = %.1e (full linearisation vs upstream: %.1e)\n', ...
       rel(L_D, L_up), rel(L_MT, L_up));

    % E: eddy viscosity bump
    BFe = BF;
    xc = mean(mesh.X(:)); yc = mean(mesh.Y(:)); s = 0.1*(max(mesh.Y(:)) - min(mesh.Y(:)));
    mut = 50 * exp(-((mesh.X - xc).^2 + (mesh.Y - yc).^2) / s^2);
    BFe.MU = BF.MU + mut;
    L_E = GetLinProblemVarMu(mesh, BFe, m);
    dE  = L_E - L_D;
    rows = mod(find(any(dE ~= 0, 2)) - 1, mesh.ngp) + 1;
    far  = mut(mesh.usedInd) < 1e-12 * max(mut(:));
    stencil_reach = any(far(rows));
    pr('       E eddy bump (peak 50): |dL|/|L| = %.1e, rows changed %d; any changed row far from the bump: %d\n', ...
       rel(L_E, L_D), numel(unique(rows)), stencil_reach);
end
fclose(fid);
end
