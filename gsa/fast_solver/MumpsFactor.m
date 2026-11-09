classdef MumpsFactor < handle
% MUMPSFACTOR  Factorise a complex sparse matrix once with MUMPS and solve many times.
%
%   F = MumpsFactor(A, opts)    factorise A (n x n, complex unsymmetric)
%   x = F.solve(b)              A  \ b
%   x = F.solve_ctranspose(b)   A' \ b   (conjugate transpose, same factors)
%   F.release()                 free the MUMPS instance (also done on delete)
%
% Settings (written for LinStab2D shift-invert, following M. Wong's 2026 note):
%   - SYM = 0 (unsymmetric), METIS nested dissection (ICNTL(7) = 5)
%   - threshold partial pivoting CNTL(1) = 0.1, MUMPS default matching/scaling
%   - rows equilibrated by their largest |entry| before factorising, so penalty
%     rows from Dirichlet BCs (1e4..1e8 on the diagonal) don't dominate pivoting
%   - A'\b uses MUMPS's transposed solve (ICNTL(9) ~= 1) with conjugation
%
% opts (all optional):
%   .mem_relax   ICNTL(14), % extra working space over the analysis estimate (default 30)
%   .pivot_tol   CNTL(1) (default 0.1)
%   .verbose     print MUMPS statistics (default true)
%   .check       do one random solve and report the relative residual (default true)

    properties (SetAccess = private)
        n
        info            % struct: timings, memory (MB), factor entries, check residuals
    end
    properties (Access = private)
        id              % MUMPS instance struct
        d               % row scaling: A_scaled = diag(d) * A
        mat_arg         % matrix argument passed on solve calls (empty sparse: no copy of A)
    end

    methods
        function F = MumpsFactor(A, opts)
            if nargin < 2, opts = struct(); end
            o = struct('mem_relax', 30, 'pivot_tol', 0.1, 'verbose', true, 'check', true);
            fn = fieldnames(opts);
            for k = 1:numel(fn), o.(fn{k}) = opts.(fn{k}); end
            if exist('zmumpsmex', 'file') ~= 3
                error('MumpsFactor:noMex', 'zmumpsmex not on the MATLAB path');
            end
            F.n = size(A, 1);
            t0 = tic;

            % --- row equilibration
            rmax = full(max(abs(A), [], 2));
            rmax(rmax == 0) = 1;
            F.d = 1 ./ rmax;
            As = spdiags(F.d, 0, F.n, F.n) * A;

            % --- init
            id = initmumps();
            id.SYM = 0;
            id = zmumps(id);
            id.ICNTL(1:4) = [-1 -1 -1 0];      % silent
            id.ICNTL(7)  = 5;                  % METIS ordering
            id.ICNTL(14) = o.mem_relax;
            id.CNTL(1)   = o.pivot_tol;

            % --- analysis + factorisation (JOB = 4)
            id.JOB = 4;
            t1 = tic;
            id = zmumps(id, As);
            t_fact = toc(t1);
            if id.INFOG(1) < 0
                err = id.INFOG(1:2);
                id.JOB = -2; id = zmumps(id); %#ok<NASGU>   % zmumps needs an output, or it returns without releasing
                error('MumpsFactor:factor', 'MUMPS factorisation failed: INFOG(1:2) = [%d %d]', err);
            end
            F.id = id;
            F.mat_arg = sparse(F.n, F.n);      % solve phase needs no matrix entries

            F.info = struct();
            F.info.n               = F.n;
            F.info.nnz             = nnz(A);
            F.info.t_factorise_s   = t_fact;
            F.info.t_total_s       = toc(t0);
            F.info.factor_entries  = factor_entries(id);
            F.info.mem_est_MB      = id.INFOG(17);   % total memory estimated by the analysis
            F.info.mem_used_MB     = id.INFOG(22);   % total memory effectively used in factorisation
            F.info.n_solves        = 0;

            if o.verbose
                fprintf(['  MUMPS: n = %d, nnz = %d, factor entries = %.3g, factorise %.1f s, ', ...
                         'memory est %d MB / used %d MB\n'], F.n, nnz(A), F.info.factor_entries, ...
                         t_fact, F.info.mem_est_MB, F.info.mem_used_MB);
            end

            if o.check
                b  = complex(randn(F.n,1), randn(F.n,1));
                x  = F.solve(b);
                xh = F.solve_ctranspose(b);
                F.info.check_resid     = norm(A*x  - b) / (norm(A,1)*norm(x)  + norm(b));
                F.info.check_resid_adj = norm(A'*xh - b) / (norm(A,1)*norm(xh) + norm(b));
                if o.verbose
                    fprintf('  MUMPS check: rel. residual %.2e (A\\b), %.2e (A''\\b)\n', ...
                            F.info.check_resid, F.info.check_resid_adj);
                end
            end
        end

        function x = solve(F, b)
            % A x = b  <=>  As x = d .* b
            x = F.run_solve(F.d .* b, 1);
        end

        function x = solve_ctranspose(F, b)
            % A' x = b, A = D^-1 As  =>  As' y = b, x = D y.
            % MUMPS ICNTL(9) ~= 1 solves As.' z = c, so As' y = b  <=>  As.' conj(y) = conj(b).
            y = conj(F.run_solve(conj(b), 2));
            x = F.d .* y;
        end

        function release(F)
            if ~isempty(F.id) && F.id.INST ~= -9999
                F.id.JOB = -2;
                F.id = zmumps(F.id);   % zmumps returns without releasing if nargout = 0
                F.id = [];
            end
        end

        function delete(F)
            F.release();
        end
    end

    methods (Access = private)
        function x = run_solve(F, rhs, icntl9)
            id = F.id;
            id.JOB = 3;
            id.ICNTL(9) = icntl9;
            id.RHS = complex(full(rhs));
            [inform, ~, sol] = zmumpsmex(id.SYM, id.JOB, id.ICNTL, id.CNTL, id.PERM_IN, ...
                id.COLSCA, id.ROWSCA, id.RHS, id.VAR_SCHUR, id.INST, id.REDRHS, id.KEEP, id.DKEEP, F.mat_arg);
            if inform(1) < 0
                error('MumpsFactor:solve', 'MUMPS solve failed: INFOG(1:2) = [%d %d]', inform(1:2));
            end
            x = sol;
            F.info.n_solves = F.info.n_solves + 1;
        end
    end
end

function e = factor_entries(id)
% INFOG(29): entries in the factors (negative = value in millions)
e = id.INFOG(29);
if e < 0, e = -e * 1e6; end
end
