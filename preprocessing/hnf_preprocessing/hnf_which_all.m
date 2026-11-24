function paths = hnf_which_all(name)
% Resolve all vendor entry points so shadowed code cannot pass provenance checks.
if exist('OCTAVE_VERSION','builtin')
    suffix = {'.m','.oct',['.' mexext()]};
    paths = {};
    for k = 1:numel(suffix)
        hits = file_in_loadpath([name suffix{k}],'all');
        if ~isempty(hits), paths = [paths; hits(:)]; end %#ok<AGROW>
    end
else
    paths = which(name,'-all');
end
end
