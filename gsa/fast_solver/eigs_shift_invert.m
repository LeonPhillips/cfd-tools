function [V, omega, V_adj, omega_adj, info] = eigs_shift_invert(L0, nEig, sigma, opts, solver_opts, W, invW)
% EIGS_SHIFT_INVERT  Shift-invert eigs with one factorisation of (L0 - sigma I), shared by
% the direct and the (discrete, W-norm) adjoint problem.
%
%   [V, omega, V_adj, omega_adj, info] = eigs_shift_invert(L0, nEig, sigma, opts, solver_opts, W, invW)
%
% Replaces, in LinStab2D's example_4_Cavity:
%   [V    ,omega    ] = eigs(L0        ,nEig,sigma      ,opts);
%   [V_adj,omega_adj] = eigs(invW*L0'*W,nEig,conj(sigma),opts);
%
% Direct:  eigs(@(x) (L0 - sigma I)\x, n, nEig, sigma, opts)  (MATLAB's function-handle shift-invert form)
% Adjoint: (invW L0' W - conj(sigma) I)^-1 = invW (L0 - sigma I)^-H W, so no second factorisation.
%
% omega and omega_adj are returned as vectors. Pass W = [] (or omit W, invW) to skip the adjoint.
% solver_opts: see make_shift_invert_solver (default type 'mumps').
% info: solver info, timings, and true residuals ||L v - w v|| / (||L0||_1 ||v||) of every mode,
%       evaluated with the original L0 (independent of the solver).

if nargin < 4 || isempty(opts), opts = struct(); end
if nargin < 5 || isempty(solver_opts), solver_opts = struct('type', 'mumps'); end
if nargin < 6, W = []; end
if nargin < 7, invW = []; end
do_adj = ~isempty(W);
if do_adj && isempty(invW), invW = spdiags(1 ./ diag(W), 0, size(W,1), size(W,2)); end

n = size(L0, 1);
info = struct();
t_all = tic;

S = make_shift_invert_solver(L0 - sigma * speye(n), solver_opts);
info.solver = S.info;

% direct
t = tic;
[V, D] = eigs(@(x) S.solve(x), n, nEig, sigma, opts);
omega = diag(D);
info.t_eigs_direct_s = toc(t);

normL = norm(L0, 1);
info.resid_direct = vecnorm(L0*V - V.*omega.') ./ (normL * vecnorm(V));

V_adj = []; omega_adj = [];
if do_adj
    t = tic;
    [V_adj, Da] = eigs(@(x) invW * S.solve_ctranspose(W * x), n, nEig, conj(sigma), opts);
    omega_adj = diag(Da);
    info.t_eigs_adjoint_s = toc(t);
    LadjV = invW * (L0' * (W * V_adj));
    info.resid_adjoint = vecnorm(LadjV - V_adj.*omega_adj.') ./ (normL * vecnorm(V_adj));
end

S.release();                         % free the factors (MUMPS) before returning
info.t_total_s = toc(t_all);
fprintf('  eigs_shift_invert: %d modes, direct %.1f s%s, total %.1f s, max residual %.1e%s\n', ...
    nEig, info.t_eigs_direct_s, ternary(do_adj, sprintf(', adjoint %.1f s', getfield_or(info, 't_eigs_adjoint_s', 0)), ''), ...
    info.t_total_s, max(info.resid_direct), ...
    ternary(do_adj, sprintf(' (adjoint %.1e)', max(getfield_or(info, 'resid_adjoint', 0))), ''));
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function v = getfield_or(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
