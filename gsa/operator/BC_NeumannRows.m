function [L0, dofs] = BC_NeumannRows(L0, mesh, idx, borders, variables, penalty)
% BC_NEUMANNROWS  Zero normal-gradient conditions: each constrained row becomes penalty*(n_x Dx + n_r Dr)
% acting on the same variable (supervisor's BC_Neumann_v3 approach, for masked LinStab2D meshes).
%
%   [L0, dofs] = BC_NeumannRows(L0, mesh, idx, borders, variables, penalty)
%     borders    cell array of border names without the 'i': 'l','r' (normal x), 't' (normal r), or 'w'
%                (masked walls; per-node normals from mesh.wall_normal, rows matching idx.wi)
%     variables  cell array from {'rho','u','v','w','T'}
%     penalty    default 1e4 (as the method)
% Uses mesh.Dx / mesh.Dy, which are one-sided (SBP boundary closures) at segment ends, so the normal
% derivative at a wall or edge point uses fluid points only.

if nargin < 6 || isempty(penalty), penalty = 1e4; end
ngp = mesh.ngp;
vars = {'rho','u','v','w','T'};
dofs = [];
for b = borders
    pts = idx.([b{1} 'i'])(:);
    switch b{1}
        % normals point INTO the fluid: with dq/dt = penalty*(n.grad q) the boundary value then relaxes towards
        % the interior (stable); an outward normal makes these rows anti-relax (unstable)
        case 'l', n = repmat([ 1  0], numel(pts), 1);   % x min
        case 'r', n = repmat([-1  0], numel(pts), 1);   % x max
        case 'b', n = repmat([ 0  1], numel(pts), 1);   % r min
        case 't', n = repmat([ 0 -1], numel(pts), 1);   % r max
        case 'w',       n = mesh.wall_normal;
            assert(size(n,1) == numel(pts), 'mesh.wall_normal must match idx.wi');
        otherwise, error('BC_NeumannRows:border', 'unknown border %s', b{1});
    end
    Dn = spdiags(n(:,1), 0, numel(pts), numel(pts)) * mesh.Dx(pts, :) + ...
         spdiags(n(:,2), 0, numel(pts), numel(pts)) * mesh.Dy(pts, :);
    for v = variables
        k = find(strcmp(vars, v{1})) - 1;
        rows = pts + k*ngp;
        L0(rows, :) = 0;
        L0(rows, (1:ngp) + k*ngp) = penalty * Dn;
        dofs = [dofs; rows]; %#ok<AGROW>
    end
end
fprintf('--- Neumann (rows, penalty %.0e): %d dofs on %s for %s\n', penalty, numel(dofs), ...
        strjoin(borders, ','), strjoin(variables, ','));
end
