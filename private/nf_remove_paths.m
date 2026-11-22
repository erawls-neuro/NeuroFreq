function nf_remove_paths(root)
% Remove exactly this directory tree from the active path; never delete files.
if isempty(root), return; end
entries = strsplit(path, pathsep);
keep = ~cellfun(@isempty, entries) & ~strcmp(entries, '.');
drop = entries(keep);
drop = drop(cellfun(@(entry) nf_under(entry, root), drop));
if ~isempty(drop), rmpath(drop{:}); end
end
