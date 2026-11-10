function test_alpha_rho(h_mm, alpha)
% TEST_ALPHA_RHO  Check the density-damping term of GetLinProblemVarMu on the srp case (operator only, no eigs).
%   test_alpha_rho(0.1, 3e-5)
% Checks: (1) the operators with and without alpha_rho differ only in the rho'-rho' block;
%         (2) that block difference is symmetric in the H = W r inner product and negative semi-definite;
%         (3) prints the damping rate of a grid-scale checkerboard and of a smooth wave (from GetLinProblemVarMu).
if nargin < 1, h_mm = 0.1; end
if nargin < 2, alpha = 3e-5; end
pc = srp_case_defaults(); pc.h_mm = h_mm;
[mesh, BF] = build_srp_case(pc);
ngp = mesh.ngp;
for m = [0 1]
    L_0 = GetLinProblemVarMu(mesh, BF, m);
    L_a = GetLinProblemVarMu(mesh, BF, m, struct('alpha_rho', alpha));
    dL = L_a - L_0;
    [i, j] = find(dL);
    fprintf('\nm = %d: nnz(dL) = %d, all in the rho-rho block: %d\n', m, nnz(dL), all(i <= ngp & j <= ngp));
    D = dL(1:ngp, 1:ngp);
    Hv = mesh.W(mesh.usedInd) .* mesh.Y(mesh.usedInd); H = spdiags(Hv(:), 0, ngp, ngp);
    S = H * D;
    fprintf('  H-symmetry: |HD - (HD)''|/|HD| = %.1e\n', normest(S - S') / normest(S));
    q = randn(ngp, 20);
    e = real(diag(q' * S * q));
    fprintf('  q''HDq for 20 random q: max %.2e (must be <= 0), min %.2e\n', max(e), min(e));
    fprintf('  nnz: operator %d -> %d (+%.1f%%)\n', nnz(L_0), nnz(L_a), 100*(nnz(L_a)/nnz(L_0) - 1));
end
end
