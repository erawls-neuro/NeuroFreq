function inside = nf_under(candidate, root)
% True for the root itself or its descendants. The separator boundary keeps
% /data/NeuroFreq_deprecated distinct from /data/NeuroFreq.
inside = false;
if isempty(candidate) || isempty(root), return; end
candidate = nf_canonical(candidate);
root = nf_canonical(root);
if ispc
    candidate = lower(candidate);
    root = lower(root);
end
inside = strcmp(candidate, root);
if ~inside
    if root(end) ~= filesep, root(end+1) = filesep; end
    inside = strncmp(candidate, root, numel(root));
end
end
