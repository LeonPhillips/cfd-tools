function derive_axy_lns(outdir)
% DERIVE_AXY_LNS  Symbolic linearisation of the axisymmetric compressible Navier-Stokes equations
% in LinStab2D's form, and generation of the coefficient scripts.
%
%   derive_axy_lns(outdir)   (needs the Symbolic Math Toolbox; ~minutes)
%
% Writes two scripts with LinStab2D's variable names (drop-in for getCoeffsLaminar.m):
%   getCoeffsAxyMuT.m    mu = mu(T) only, grad(mu_bar) = dmudT*grad(T_bar).
%                        Must reproduce LinStab2D's getCoeffsLaminar.m (checked by test_coeffs_vs_upstream).
%   getCoeffsAxyVarMu.m  mu_bar = MU(r,z) is any field (e.g. Sutherland + sigma*mu_t). Its gradient enters
%                        through the inputs dMUdr, dMUdz. The perturbation is mu' = dmudT*T' (Sutherland only,
%                        i.e. frozen eddy viscosity), with grad(dmudT) = d2mudT2*grad(T_bar).
%                        Heat conduction uses the same total mu: q = -(MU/c1) grad T (Pr_eff = Pr).
%
% Conventions (as LinStab2D 'axy', coordinate system CS1 in GetLinProblem):
%   r radial, theta azimuthal (perturbations ~ exp(1i*m*theta)), z axial
%   q = [rho, U (radial), V (azimuthal), W (axial), T]
%   p = rho*T/c2, c2 = kappa*Ma^2;  E = cv*T + |u|^2/2, cv = 1/(kappa*(kappa-1)*Ma^2)
%   heat flux -(mu/c1) grad T, c1 = (kappa-1)*Re*Pr*Ma^2;  stress tau/Re (Stokes hypothesis)
%   Conservative form dQ/dt = F(q), Q = [rho, rho*u, rho*v, rho*w, rho*E].
%   LinStab2D: LHS q = RHS dq/dt, LHS = dF/dq (coefficients A0, Ar, Az, Arr, Azz, Arz), RHS = 1i*dQ/dq (B).

if nargin < 1, outdir = pwd; end
t0 = tic;

syms r z th m ep real
Re = sym('Re', 'positive'); c1 = sym('c1', 'positive');   % sym(), not syms: 'cv' is a MATLAB function
c2 = sym('c2', 'positive'); cv = sym('cv', 'positive');
syms RHOb(r,z) Ub(r,z) Vb(r,z) Wb(r,z) Tb(r,z) MUb(r,z) dmudTb(r,z)
syms rh(r,z) uh(r,z) vh(r,z) wh(r,z) Th(r,z)

E   = exp(1i*m*th);
rho = RHOb + ep*rh*E;
ur  = Ub   + ep*uh*E;
ut  = Vb   + ep*vh*E;
uz  = Wb   + ep*wh*E;
T   = Tb   + ep*Th*E;
mu  = MUb  + ep*dmudTb*Th*E;          % mu' = dmu/dT * T' (only the T-dependent part is perturbed)
p   = rho*T/c2;

dr  = @(f) diff(f, r);
dt  = @(f) diff(f, th);
dz  = @(f) diff(f, z);
div = @(ar, at, az) (1/r)*dr(r*ar) + (1/r)*dt(at) + dz(az);

% strain rate and stress
divu = div(ur, ut, uz);
Srr = dr(ur);
Stt = (dt(ut) + ur)/r;
Szz = dz(uz);
Srt = (r*dr(ut/r) + dt(ur)/r)/2;
Srz = (dz(ur) + dr(uz))/2;
Stz = (dz(ut) + dt(uz)/r)/2;
trr = mu*(2*Srr - 2/3*divu);
ttt = mu*(2*Stt - 2/3*divu);
tzz = mu*(2*Szz - 2/3*divu);
trt = 2*mu*Srt;
trz = 2*mu*Srz;
ttz = 2*mu*Stz;

% fluxes
F = sym(zeros(5,1));
F(1) = -div(rho*ur, rho*ut, rho*uz);
F(2) = -( (1/r)*dr(r*rho*ur*ur) + (1/r)*dt(rho*ur*ut) + dz(rho*ur*uz) - rho*ut^2/r ) - dr(p) ...
       + ( (1/r)*dr(r*trr) + (1/r)*dt(trt) + dz(trz) - ttt/r )/Re;
F(3) = -( (1/r^2)*dr(r^2*rho*ur*ut) + (1/r)*dt(rho*ut*ut) + dz(rho*ut*uz) ) - dt(p)/r ...
       + ( (1/r^2)*dr(r^2*trt) + (1/r)*dt(ttt) + dz(ttz) )/Re;
F(4) = -( (1/r)*dr(r*rho*uz*ur) + (1/r)*dt(rho*uz*ut) + dz(rho*uz*uz) ) - dz(p) ...
       + ( (1/r)*dr(r*trz) + (1/r)*dt(ttz) + dz(tzz) )/Re;
Et = rho*(cv*T + (ur^2 + ut^2 + uz^2)/2);
ar = trr*ur + trt*ut + trz*uz;
at = trt*ur + ttt*ut + ttz*uz;
az = trz*ur + ttz*ut + tzz*uz;
F(5) = -div((Et + p)*ur, (Et + p)*ut, (Et + p)*uz) + div(ar, at, az)/Re ...
       + div(mu/c1*dr(T), mu/c1*dt(T)/r, mu/c1*dz(T));
Q = [rho; rho*ur; rho*ut; rho*uz; Et];

% linearise: d/d(ep) at ep = 0, then strip exp(1i*m*th) (every term is linear in it)
G  = subs(subs(diff(F, ep), ep, 0), th, 0);
JB = subs(subs(diff(Q, ep), ep, 0), th, 0);
fprintf('linearised in %.0f s\n', toc(t0));

% placeholders for perturbation derivatives
pert  = {rh, uh, vh, wh, Th};
vname = {'1','2','3','4','5'};
kinds = {'rr','zz','rz','r','z','0'};
P = sym(zeros(5, 6));
for j = 1:5
    f = pert{j};
    P(j,:) = sym({sprintf('P%s_rr',vname{j}), sprintf('P%s_zz',vname{j}), sprintf('P%s_rz',vname{j}), ...
                  sprintf('P%s_r',vname{j}),  sprintf('P%s_z',vname{j}),  sprintf('P%s_0',vname{j})});
    G  = subs(G,  [diff(f,r,r), diff(f,z,z), diff(f,r,z), diff(f,z,r)], [P(j,1), P(j,2), P(j,3), P(j,3)]);
    G  = subs(G,  [diff(f,r), diff(f,z)], [P(j,4), P(j,5)]);
    G  = subs(G,  f, P(j,6));
    JB = subs(JB, f, P(j,6));
end

% base-flow names (LinStab2D's)
mk = @(s) sym(s, 'real');
U = mk('U'); V = mk('V'); W = mk('W'); RHO = mk('RHO'); Tn = mk('Tn'); MU = mk('MU');
dmudT = mk('dmudT'); d2mudT2 = mk('d2mudT2'); dMUdr = mk('dMUdr'); dMUdz = mk('dMUdz'); R = mk('R');
base = {Ub,'U'; Vb,'V'; Wb,'W'; RHOb,'RHO'; Tb,'T'};
for k = 1:size(base,1)
    f = base{k,1}; nm = base{k,2};
    G = subs(G, [diff(f,r,r), diff(f,z,z), diff(f,r,z), diff(f,z,r)], ...
                [sym(['d2' nm 'dr2']), sym(['d2' nm 'dz2']), sym(['d2' nm 'drz']), sym(['d2' nm 'drz'])]);
    G = subs(G, [diff(f,r), diff(f,z)], [sym(['d' nm 'dr']), sym(['d' nm 'dz'])]);
end
G  = subs(G,  [diff(dmudTb,r), diff(dmudTb,z)], [d2mudT2*sym('dTdr'), d2mudT2*sym('dTdz')]);
G  = subs(G,  [diff(MUb,r), diff(MUb,z)], [dMUdr, dMUdz]);
G  = subs(G,  [Ub, Vb, Wb, RHOb, Tb, MUb, dmudTb], [U, V, W, RHO, Tn, MU, dmudT]);
JB = subs(JB, [Ub, Vb, Wb, RHOb, Tb], [U, V, W, RHO, Tn]);
G  = subs(G, r, R);
left = setdiff(arrayfun(@char, symvar(G), 'UniformOutput', false), ...
               [arrayfun(@char, P(:).', 'UniformOutput', false), {'U','V','W','RHO','Tn','MU','dmudT','d2mudT2', ...
                'dMUdr','dMUdz','R','m','Re','c1','c2','cv'}, ...
                strcat('d', {'U','V','W','RHO','T'}, 'dr'), strcat('d', {'U','V','W','RHO','T'}, 'dz'), ...
                strcat('d2', {'U','V','W','RHO','T'}, 'dr2'), strcat('d2', {'U','V','W','RHO','T'}, 'dz2'), ...
                strcat('d2', {'U','V','W','RHO','T'}, 'drz')]);
assert(isempty(left), 'unexpected symbols left: %s', strjoin(left, ', '));
assert(isempty(regexp(char(G), '\(r, z\)', 'once')), 'unsubstituted symbolic functions remain');

% coefficients
pre = {'Arr','Azz','Arz','Ar','Az','A0'};    % same order as kinds / columns of P
C = struct();
for i = 1:5
    for j = 1:5
        for k = 1:6
            C.(sprintf('%s_%d%d', pre{k}, i, j)) = diff(G(i), P(j,k));
        end
        C.(sprintf('B_%d%d', i, j)) = 1i*diff(JB(i), P(j,6));
    end
end
fprintf('coefficients extracted at %.0f s\n', toc(t0));

% version 1: mu = mu(T)  ->  grad(mu_bar) = dmudT grad(T_bar)
Cmt = structfun(@(c) subs(c, [dMUdr, dMUdz], [dmudT*sym('dTdr'), dmudT*sym('dTdz')]), C, 'UniformOutput', false);
write_script(fullfile(outdir, 'getCoeffsAxyMuT.m'), Cmt, ...
    'mu = mu(T): grad(mu_bar) = dmudT*grad(T_bar). Must match LinStab2D getCoeffsLaminar.m.');
% version 2: variable mu_bar field, frozen eddy viscosity
write_script(fullfile(outdir, 'getCoeffsAxyVarMu.m'), C, ...
    ['mu_bar = MU(r,z) any field (Sutherland + sigma*mu_t); grad(mu_bar) from inputs dMUdr, dMUdz; ' ...
     'mu'' = dmudT*T'' (frozen eddy viscosity); heat conduction with total MU.']);
fprintf('done in %.0f s\n', toc(t0));
end

function write_script(fname, C, note)
fid = fopen(fname, 'w');
[~, nm] = fileparts(fname);
fprintf(fid, '%% %s  GENERATED by derive_axy_lns.m on %s - do not edit by hand.\n', nm, char(datetime('now','Format','yyyy-MM-dd HH:mm')));
fprintf(fid, '%% %s\n', note);
fprintf(fid, '%% Script with LinStab2D''s variable names (see GetLinProblem.m / getCoeffsLaminar.m), CS1:\n');
fprintf(fid, '%% R radial, theta azimuthal (exp(1i*m*theta)), z axial; q = [rho, U(radial), V(azimuthal), W(axial), T].\n\n');
f = fieldnames(C);
for k = 1:numel(f)
    c = C.(f{k});
    if isAlways(c == 0, 'Unknown', 'false')
        fprintf(fid, '%s = Z;\n', f{k});
    else
        s = char(simplify(c, 'Steps', 10));
        s = regexprep(s, '\<Tn\>', 'T');
        s = regexprep(s, '(?<!\.)([\*/\^])', '.$1');
        fprintf(fid, '%s = %s + Z;\n', f{k}, s);
    end
end
fclose(fid);
fprintf('wrote %s (%d coefficients)\n', fname, numel(f));
end
