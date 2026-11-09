function ok = fls_mumps_available()
% FLS_MUMPS_AVAILABLE  True if the MUMPS MATLAB interface (zmumpsmex + zmumps.m) is on the path.
% Add it with addpath(getenv('FLS_MUMPS_DIR')) or addpath('<dir with zmumpsmex.mexa64>').
d = getenv('FLS_MUMPS_DIR');
if ~isempty(d) && exist(d, 'dir') && exist('zmumpsmex', 'file') ~= 3
    addpath(d);
end
ok = exist('zmumpsmex', 'file') == 3 && exist('zmumps', 'file') == 2 && exist('initmumps', 'file') == 2;
end
