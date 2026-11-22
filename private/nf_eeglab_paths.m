function paths = nf_eeglab_paths(folder)
% These are the core path directories configured by EEGLAB's startup.
% List them explicitly. Recursive genpath would expose private functions,
% optional plugin copies, and shadowing compatibility implementations.
paths = {};
if isempty(folder), return; end
names = {'','functions','functions/adminfunc','functions/sigprocfunc', ...
    'functions/guifunc','functions/studyfunc','functions/popfunc', ...
    'functions/statistics','functions/timefreqfunc','functions/miscfunc', ...
    'functions/supportfiles'};
for k = 1:numel(names)
    member = fullfile(folder,names{k});
    if isfolder(member), paths{end+1} = member; end %#ok<AGROW>
end
end
