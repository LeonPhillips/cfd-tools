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

% Test 4: what the differences in Test 1 are. With dmudT = 0 (no mu' = dmudT*T' term, only its
% gradient d2mudT2*grad(T) T') they should vanish if upstream only omits mu' * (derivatives of the
% base stress), i.e. the mu'*div(tau_bar) and mu'*(viscous work) terms.
v5 = v; v5.dmudT = zeros(n,1);
up5 = run_coeffs('getCoeffsLaminar', v5);
mt5 = run_coeffs('getCoeffsAxyMuT',  v5);
% Test 5: are the differences independent of m and d2mudT2?
v6 = v; v6.m = 7; v6.d2mudT2 = 3*v.d2mudT2;
up6 = run_coeffs('getCoeffsLaminar', v6);
mt6 = run_coeffs('getCoeffsAxyMuT',  v6);
dd = struct();
for f = {'A0_25','A0_35','A0_45','A0_55'}
    dd.(f{1}) = max(abs((mt6.(f{1}) - up6.(f{1})) - (mt.(f{1}) - up.(f{1})))) / max(abs(mt.(f{1}) - up.(f{1})));
end

report = struct();
report.test1 = compare(up, mt);
report.test2 = compare(mt, vm);
report.test3 = compare(up3, vm3);
report.test4 = compare(up5, mt5);
report.test5 = dd;
fprintf('test5 (change of the difference when m and d2mudT2 change, relative): A0_25 %.1e, A0_35 %.1e, A0_45 %.1e, A0_55 %.1e\n', ...
        dd.A0_25, dd.A0_35, dd.A0_45, dd.A0_55);
for t = {'test1','test2','test3','test4'}
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
