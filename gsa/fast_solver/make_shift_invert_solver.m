function S = make_shift_invert_solver(A, opts)
% MAKE_SHIFT_INVERT_SOLVER  Factorise A once; return handles for A\b and A'\b.
%
%   S = make_shift_invert_solver(A, opts)
%   S.solve(b)             A  \ b
%   S.solve_ctranspose(b)  A' \ b   (same factors)
%   S.info                 timings / memory / factor size
%   S.release()            free the factors
%
% opts.type: 'mumps' (default) or 'umfpack' (MATLAB's built-in sparse LU, the same
%            factorisation LinStab2D uses; useful for checking the wiring).
% opts.allow_fallback: if true, 'mumps' falls back to 'umfpack' with a warning when the
%            MEX is missing (default false: a run sized for MUMPS must not silently fall back).
% Other opts fields are passed to MumpsFactor.

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'type'), opts.type = 'mumps'; end
if ~isfield(opts, 'allow_fallback'), opts.allow_fallback = false; end

type = lower(opts.type);
if strcmp(type, 'mumps') && ~fls_mumps_available()
    if opts.allow_fallback
        warning('make_shift_invert_solver:fallback', 'MUMPS not available, using UMFPACK');
        type = 'umfpack';
    else
        error('make_shift_invert_solver:noMumps', 'MUMPS (zmumpsmex) not on the MATLAB path');
    end
end

switch type
    case 'mumps'
        mopts = rmfield(opts, intersect(fieldnames(opts), {'type', 'allow_fallback'}));
        F = MumpsFactor(A, mopts);
        S.solve            = @(b) F.solve(b);
        S.solve_ctranspose = @(b) F.solve_ctranspose(b);
        S.release          = @() F.release();
        S.factor           = F;
        S.info             = F.info;
        S.info.type        = 'mumps';

    case 'umfpack'
        t0 = tic;
        % P*(R\A)*Q = L*U  (UMFPACK, as in LinStab2D's GetInverseFunction 'lu')
        [L, U, P, Q, R] = lu(A);
        S.solve            = @(b) Q * (U \ (L \ (P * (R \ b))));
        S.solve_ctranspose = @(b) R' \ (P' * (L' \ (U' \ (Q' * b))));
        S.release          = @() [];
        S.info = struct('type', 'umfpack', 'n', size(A,1), 'nnz', nnz(A), ...
                        't_factorise_s', toc(t0), 'factor_entries', nnz(L) + nnz(U));
        fprintf('  UMFPACK: n = %d, factor entries = %.3g, factorise %.1f s\n', ...
                size(A,1), S.info.factor_entries, S.info.t_factorise_s);

    otherwise
        error('make_shift_invert_solver:type', 'unknown solver type ''%s''', opts.type);
end
end
