function report = test_coeffs_vs_upstream(n)
% TEST_COEFFS_VS_UPSTREAM  Compare the generated coefficient scripts with LinStab2D's getCoeffsLaminar.m
% at n random points (all base-flow inputs random), m = 3.
%
% Test 1: getCoeffsAxyMuT   vs getCoeffsLaminar           -> should agree to round-off
% Test 2: getCoeffsAxyVarMu with dMUdr = dmudT*dTdr etc.  -> identical to Test 1 (consistency of the eddy version)
% Test 3: getCoeffsAxyVarMu with MU + c (constant mu_t)   -> equals getCoeffsLaminar with MU + c
% Needs LinStab2D's aux_matlab and this folder on the path.

if nargin < 1, n = 200; end
rng(1);
v = random_inputs(n);

up  = run_coeffs('getCoeffsLaminar',  v);
mt  = run_coeffs('getCoeffsAxyMuT',   v);
v2 = v; v2.dMUdr = v.dmudT .* v.dTdr; v2.dMUdz = v.dmudT .* v.dTdz;
vm  = run_coeffs('getCoeffsAxyVarMu', v2);
c = 7.3;
v3 = v2; v3.MU = v.MU + c;            % constant mu_t: grad(mu_bar) unchanged
vm3 = run_coeffs('getCoeffsAxyVarMu', v3);
v4 = v;  v4.MU = v.MU + c;
up3 = run_coeffs('getCoeffsLaminar', v4);

report = struct();
report.test1 = compare(up, mt);
report.test2 = compare(mt, vm);
report.test3 = compare(up3, vm3);
for t = {'test1','test2','test3'}
    R = report.(t{1});
    fprintf('%s: max relative difference %.2e over %d coefficients (worst %s)\n', t{1}, R.max, R.n, R.worst);
    if ~isempty(R.bad)
        fprintf('   coefficients above 1e-10: %s\n', strjoin(R.bad, ', '));
    end
end
end

function v = random_inputs(n)
v = struct();
names = {'U','V','W','RHO','T','MU','dmudT','d2mudT2'};
for k = 1:numel(names), v.(names{k}) = 0.5 + rand(n,1); end
for f = {'U','V','W','RHO','T'}
    for s = {'d%sdr','d%sdz','d2%sdr2','d2%sdz2','d2%sdrz'}
        v.(sprintf(s{1}, f{1})) = randn(n,1);
    end
end
v.dMUdr = randn(n,1); v.dMUdz = randn(n,1);
v.R  = 0.1 + rand(n,1);
v.m  = 3; v.Re = 1234.5; v.cv = 1/(1.4*0.4*2.5^2); v.c2 = 1.4*2.5^2; v.c1 = 0.4*1234.5*0.72*2.5^2;
v.kappa = 1.4;
v.Z = zeros(n,1); v.I = ones(n,1);
end

function S = run_coeffs(script, v) %#ok<INUSD>
f = fieldnames(v);
for k = 1:numel(f), eval([f{k} ' = v.(f{k});']); end
eval(script);
S = struct();
vars = who('-regexp', '^(A0|Ar|Az|Arr|Azz|Arz|B)_\d\d$');
for k = 1:numel(vars), S.(vars{k}) = eval(vars{k}); end
end

function R = compare(a, b)
f = fieldnames(a);
assert(isequal(sort(f), sort(fieldnames(b))), 'coefficient sets differ');
R = struct('max', 0, 'worst', '', 'bad', {{}}, 'n', numel(f));
for k = 1:numel(f)
    x = a.(f{k}); y = b.(f{k});
    if isscalar(x), x = x + zeros(size(y)); end
    if isscalar(y), y = y + zeros(size(x)); end
    d = max(abs(x - y)) / max(1, max(abs(x)));
    if d > R.max, R.max = d; R.worst = f{k}; end
    if d > 1e-10, R.bad{end+1} = sprintf('%s (%.1e)', f{k}, d); end
end
end
