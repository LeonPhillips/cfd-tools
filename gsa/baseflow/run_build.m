function run_build(outdir, h_mm, varargin)
% RUN_BUILD  Build the srp_M2.5 GSA case at grid spacing h_mm, write check plots and save it.
%   run_build(outdir, h_mm, 'field', value, ...)   overrides top-level fields of srp_case_defaults
% Saves <outdir>/case_h<h>.mat with BF, the mesh index data, sponge and parameters (the mesh's
% derivative matrices are rebuilt from these when the operator is made).
p = srp_case_defaults();
p.h_mm = h_mm;
for k = 1:2:numel(varargin)
    f = strsplit(varargin{k}, '.');      % 'smooth.x_min_mm' sets a nested field
    p = setfield(p, f{:}, varargin{k+1});
end
if ~exist(outdir, 'dir'), mkdir(outdir); end
[mesh, BF, info] = build_srp_case(p);
tag = sprintf('h%s', strrep(sprintf('%g', h_mm), '.', 'p'));
if isfield(p, 'tag') && ~isempty(p.tag), tag = [tag '_' char(p.tag)]; end
check_srp_case(mesh, BF, info, outdir, tag);
C = struct('BF', BF, 'p', p, 'idx', mesh.idx, 'usedInd', mesh.usedInd, 'mask', mesh.mask, ...
           'sponge', mesh.sponge, 'wall_normal', mesh.wall_normal, 'X', mesh.X, 'Y', mesh.Y, ...
           'shock', info.shock, 'x_range', info.x_range, 'r_max', info.r_max, 'Nx', info.Nx, 'Nr', info.Nr);
save(fullfile(outdir, ['case_' tag '.mat']), '-struct', 'C', '-v7.3');
fprintf('saved %s\n', fullfile(outdir, ['case_' tag '.mat']));
end
