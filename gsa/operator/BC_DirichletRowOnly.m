function [L0, dofs] = BC_DirichletRowOnly(L0, idx, borders, variables, penalty)
% BC_DIRICHLETROWONLY  Zero-value Dirichlet conditions by zeroing rows only (supervisor's method,
% M. Wong's BC_Dirichlet_RowOnly), with LinStab2D's index names.
%
%   [L0, dofs] = BC_DirichletRowOnly(L0, idx, borders, variables, penalty)
%     borders    cell array of LinStab2D border names without the 'i', e.g. {'l','t','m'}
%                (idx.<b>i_<var> must exist, as made by GetLinProblem / GetLinProblemVarMu)
%     variables  cell array from {'rho','u','v','w','T'}
%     penalty    diagonal value (default 1e8, as the method; LinStab2D's BC_Dirichlet uses 1e4
%                and also zeroes the columns)
% Rows of the constrained dofs become penalty*I; columns are left as they are.

if nargin < 5 || isempty(penalty), penalty = 1e8; end
dofs = [];
for b = borders
    for v = variables
        f = [b{1} 'i_' v{1}];
        if ~isfield(idx, f), error('BC_DirichletRowOnly:noField', 'idx.%s not found', f); end
        dofs = [dofs; idx.(f)(:)]; %#ok<AGROW>
    end
end
dofs = unique(dofs);
L0(dofs, :) = 0;
L0(dofs, dofs) = speye(numel(dofs)) * penalty;
fprintf('--- Dirichlet (row only, penalty %.0e): %d dofs on %s for %s\n', penalty, numel(dofs), ...
        strjoin(borders, ','), strjoin(variables, ','));
end
