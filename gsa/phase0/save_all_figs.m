function save_all_figs(outdir, prefix)
% Save every open figure as PNG: <outdir>/<prefix>_NN_<figure name>.png
figs = findall(0, 'type', 'figure');
[~, order] = sort([figs.Number]);
figs = figs(order);
for k = 1:numel(figs)
    name = regexprep(figs(k).Name, '[^A-Za-z0-9]+', '_');
    f = fullfile(outdir, sprintf('%s_%02d_%s.png', prefix, k, name));
    try
        exportgraphics(figs(k), f, 'Resolution', 120);
    catch err
        fprintf('could not save figure %d: %s\n', k, err.message);
    end
end
fprintf('saved %d figures with prefix %s\n', numel(figs), prefix);
end
