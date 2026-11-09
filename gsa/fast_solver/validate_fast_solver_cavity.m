function validate_fast_solver_cavity(repo, outdir, F, modes)
% VALIDATE_FAST_SOLVER_CAVITY  Check eigs_shift_invert against LinStab2D's own eigs calls on
% example_4_Cavity (50 global modes near omega = 1.5+1i, direct and adjoint).
%
%   validate_fast_solver_cavity(repo, outdir, F, modes)
%     repo   LinStab2D clone (the example is run up to its solve step, unmodified apart from
%            dropping 'clear'/'close' lines, verbose = false and the grid factor F)
%     outdir where the summary (.txt) and eigenvalues (.mat) go
%     F      grid refinement factor on Nx_cavity = 80, Ny_cavity = 50 (default 1 = upstream)
%     modes  cell array of solvers to run, from {'builtin','umfpack','mumps'} (default all three)
%
% 'builtin' = the example's own two eigs calls (two UMFPACK factorisations).
% Matches each eigenvalue to the nearest builtin one; reports max |d omega| and the minimum
% W-norm alignment |v1' W v2| / (|v1|_W |v2|_W) of the matched modes.

if nargin < 3 || isempty(F), F = 1; end
if nargin < 4 || isempty(modes), modes = {'builtin', 'umfpack', 'mumps'}; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
here = fileparts(mfilename('fullpath'));
addpath(here);
fls_mumps_available();   % adds FLS_MUMPS_DIR if set

% --- build the cavity operator with the example's own code, up to its solve step
code = fileread(fullfile(repo, 'example_4_Cavity.m'));
cut  = strfind(code, '%setup eigs options');
assert(~isempty(cut), 'marker "%%setup eigs options" not found in example_4_Cavity.m');
code = code(1:cut(1)-1);
% (MATLAB's '.' matches newlines unless 'dotexceptnewline' is given)
rx = @(s, p, r) regexprep(s, p, r, 'lineanchors', 'dotexceptnewline');
code = rx(code, '^[ \t]*(clear all|clear variables|close all)[ \t]*;?', '');
code = rx(code, '^verbose[ \t]*=.*$', 'verbose = false;');
code = rx(code, '^Nx_cavity[ \t]*=.*$', sprintf('Nx_cavity = round(80*%g);', F));
code = rx(code, '^Ny_cavity[ \t]*=.*$', sprintf('Ny_cavity = round(50*%g);', F));
% Wrap it in a generated function (variables created by eval inside a function are not
% visible to the rest of the function in current MATLAB).
fname = ['cav_setup_F' strrep(sprintf('%g', F), '.', 'p')];
fid = fopen(fullfile(outdir, [fname '.m']), 'w');
fprintf(fid, 'function [L0, W, invW, nEig, omega_target] = %s()\n%s\nend\n', fname, code);
fclose(fid);
addpath(outdir);
rehash;
old = cd(repo);                              % the example's relative paths need the repo as cwd
restore = onCleanup(@() cd(old));
[L0, W, invW, nEig, omega_target] = feval(fname);
clear restore

opts = struct('tol', 1e-8, 'disp', 0, 'p', 300);   % the example's eigs settings
n = size(L0, 1);
fprintf('cavity F = %g: n = %d, nnz = %d, nEig = %d, sigma = %s\n', F, n, nnz(L0), nEig, num2str(omega_target));

R = struct();
for k = 1:numel(modes)
    m = modes{k};
    fprintf('--- %s\n', m);
    t = tic;
    switch m
        case 'builtin'
            [V, D]   = eigs(L0, nEig, omega_target, opts);
            [Va, Da] = eigs(invW*L0'*W, nEig, conj(omega_target), opts);
            R.(m) = struct('omega', diag(D), 'V', V, 'omega_adj', diag(Da), 'V_adj', Va, 'info', struct());
        otherwise
            [V, w, Va, wa, info] = eigs_shift_invert(L0, nEig, omega_target, opts, struct('type', m), W, invW);
            R.(m) = struct('omega', w, 'V', V, 'omega_adj', wa, 'V_adj', Va, 'info', info);
    end
    R.(m).time_s = toc(t);
    fprintf('    %s: %.1f s\n', m, R.(m).time_s);
end

% --- compare against the reference (builtin if run, else the first solver)
ref = modes{1};
fid = fopen(fullfile(outdir, sprintf('validate_cavity_F%g.txt', F)), 'w');
pr = @(varargin) cellfun(@(f) fprintf(f, varargin{:}), {1, fid});
pr('cavity F = %g, n = %d, nnz = %d, nEig = %d, sigma = %s, reference = %s\n', F, n, nnz(L0), nEig, num2str(omega_target), ref);
for k = 1:numel(modes)
    m = modes{k};
    pr('%-8s time %7.1f s', m, R.(m).time_s);
    if isfield(R.(m).info, 'solver')
        s = R.(m).info.solver;
        pr(' | factorise %.1f s, factor entries %.3g', s.t_factorise_s, s.factor_entries);
        if isfield(s, 'mem_used_MB'), pr(', MUMPS mem used %d MB', s.mem_used_MB); end
        pr(' | max residual direct %.1e, adjoint %.1e', max(R.(m).info.resid_direct), max(R.(m).info.resid_adjoint));
    end
    pr('\n');
    if ~strcmp(m, ref)
        [dw, al]   = compare(R.(ref).omega,     R.(ref).V,     R.(m).omega,     R.(m).V,     W);
        [dwa, ala] = compare(R.(ref).omega_adj, R.(ref).V_adj, R.(m).omega_adj, R.(m).V_adj, W);
        % modes away from the window edge: the 40 closest to the shift
        [~, ix] = sort(abs(R.(m).omega - omega_target)); in = ix(1:min(40, end));
        pr('         vs %s: direct max|dw| %.1e (inner 40: %.1e), min alignment 1-%.1e | adjoint max|dw| %.1e, min alignment 1-%.1e\n', ...
           ref, max(dw), max(dw(in)), 1 - min(al), max(dwa), 1 - min(ala));
    end
end
fclose(fid);

S = struct();
for k = 1:numel(modes)
    S.(modes{k}) = rmfield(R.(modes{k}), {'V', 'V_adj'});
end
save(fullfile(outdir, sprintf('validate_cavity_F%g.mat', F)), '-struct', 'S');
end

function [dw, al] = compare(w1, V1, w2, V2, W)
dw = zeros(numel(w2), 1); al = zeros(numel(w2), 1);
for j = 1:numel(w2)
    [dw(j), i] = min(abs(w1 - w2(j)));
    a = V1(:, i); b = V2(:, j);
    al(j) = abs(a' * W * b) / sqrt(real(a' * W * a) * real(b' * W * b));
end
end
