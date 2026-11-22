function varargout = nf_run(target, varargin)
% NF_RUN  Run a NeuroFreq function or pipeline with its selected dependencies.
%
%   [TF, settings] = nf_run('nf_tftransform', EEG, 'method', 'stft');
%   results = nf_run(@my_analysis, EEG);
%
% Use the first form when only the NeuroFreq root is on the MATLAB path.
% The second form runs an existing function that calls several NF routines
% and shares the same temporary search path across the whole analysis.
% Inputs and requested outputs are passed through unchanged. On return or
% error, the caller gets exactly the MATLAB path they had on entry.
%
% nf_install records dependency choices in dependencies/_nf_local/runtime.mat.
% Paths inside NeuroFreq can be stored relative to its root so that moving
% the toolbox does not break a local installation. Without that file, NF's
% source folders and the user's current MATLAB dependencies are available.

if isstring(target) && isscalar(target)
    target = char(target);
end
named = ischar(target) && isrow(target);
if named
    if ~isvarname(target) || ~strncmp(target, 'nf_', 3) || ...
            any(strcmp(target, {'nf_install', 'nf_run'}))
        error('NeuroFreq:InvalidTarget', 'Expected a public nf_ function name.');
    end
elseif ~isa(target, 'function_handle')
    error('NeuroFreq:InvalidTarget', 'Pass an NF function name or a pipeline function handle.');
end

root = fileparts(mfilename('fullpath'));
sourceDirs = nf_source_dirs(root);
dependencyRoot = fullfile(root, 'dependencies');
settings = fullfile(dependencyRoot, '_nf_local', 'runtime.mat');
eeglabRoot = '';
providerDirs = {};
if isfile(settings)
    saved = load(settings, 'config');
    if ~isfield(saved, 'config') || ~isfield(saved.config, 'scopePaths') || ...
            ~isfield(saved.config, 'eeglabRoot')
        error('NeuroFreq:InvalidConfiguration', 'Run nf_install to rebuild %s.', settings);
    end
    providerDirs = cellfun(@(d) nf_resolve_dir(root, d), ...
        saved.config.scopePaths(:)', 'UniformOutput', false);
    providerDirs = providerDirs(cellfun(@isfolder, providerDirs));
    if ~isempty(saved.config.eeglabRoot)
        eeglabRoot = nf_resolve_dir(root, saved.config.eeglabRoot);
        if ~isfile(fullfile(eeglabRoot, 'eeglab.m'))
            error('NeuroFreq:MissingEEGLAB', ...
                'The selected EEGLAB has moved from %s. Run nf_install again.', eeglabRoot);
        end
        providerDirs{end + 1} = eeglabRoot;
    end
end

original = path;
% A dormant cached copy must stay dormant even if the user previously
% added the entire NF tree using genpath. Only the selected provider
% directories below are reintroduced.
priorDirs = nf_prior_dirs(original, pwd, dependencyRoot, eeglabRoot);
if nf_under(pwd, dependencyRoot) && ~any(cellfun(@(d) nf_same_file(pwd, d), providerDirs))
    error('NeuroFreq:ConflictingCurrentFolder', ...
        'Change out of the inactive dependency folder %s before running NeuroFreq.', pwd);
end
% The path is restored even if a dependency changes it during execution.
restore = onCleanup(@() path(original)); %#ok<NASGU>
scope = unique([sourceDirs, providerDirs, priorDirs], 'stable');
path(strjoin(scope, pathsep));

if ~isempty(eeglabRoot) && ~nf_same_file(which('eeglab'), fullfile(eeglabRoot, 'eeglab.m'))
    error('NeuroFreq:ShadowedEEGLAB', ...
        'Expected eeglab.m from %s but MATLAB found %s. Check the current folder.', ...
        eeglabRoot, which('eeglab'));
end
if named
    resolved = which(target);
    if isempty(resolved) || ~any(cellfun(@(d) nf_same_file(resolved, ...
            fullfile(d, [target '.m'])), sourceDirs))
        error('NeuroFreq:UnknownFunction', ...
            'Could not resolve %s to an NF source file in %s. MATLAB found "%s".', ...
            target, root, resolved);
    end
end
[varargout{1:nargout}] = feval(target, varargin{:});
end

function folder = nf_resolve_dir(root, folder)
% External choices are absolute paths; installer-owned choices are relative.
if isempty(regexp(folder, '^([A-Za-z]:[/\\]|[/\\])', 'once'))
    folder = fullfile(root, folder);
end
folder = nf_canonical(folder);
end

function same = nf_same_file(first, second)
same = false;
if isempty(first) || isempty(second), return, end
first = nf_canonical(first);
second = nf_canonical(second);
if ispc
    same = strcmpi(first, second);
else
    same = strcmp(first, second);
end
end

function keep = nf_prior_dirs(pathText, cwd, dependencyRoot, eeglabRoot)
% Partial installations can expose popfunc or plugins without eeglab.m on
% the path. Check each enclosing directory, then add back only explicitly
% selected providers. Cache this scan for repeated calls from a pipeline.
persistent previousPath previousCwd previousEEG previousKeep
if isequal(pathText, previousPath) && isequal(cwd, previousCwd) && ...
        isequal(eeglabRoot, previousEEG)
    keep = previousKeep;
    return
end
entries = strsplit(pathText, pathsep, 'CollapseDelimiters', false);
entries = entries(~cellfun('isempty', entries));
cache = containers.Map('KeyType', 'char', 'ValueType', 'char');
allowed = true(size(entries));
for i = 1:numel(entries)
    directory = nf_canonical(entries{i});
    allowed(i) = ~nf_under(directory, dependencyRoot) && ...
        (isempty(eeglabRoot) || isempty(find_owner(directory)));
end
if ~isempty(eeglabRoot)
    owner = find_owner(nf_canonical(cwd));
    if ~isempty(owner) && ~nf_same_file(owner, eeglabRoot)
        error('NeuroFreq:ConflictingCurrentFolder', ...
            ['MATLAB searches its current folder first. Change out of %s ' ...
             'before using the EEGLAB selected for NeuroFreq.'], owner);
    end
end
keep = entries(allowed);
previousPath = pathText;
previousCwd = cwd;
previousEEG = eeglabRoot;
previousKeep = keep;

    function owner = find_owner(directory)
        % Visit each ancestor once, even for an EEGLAB installation with
        % dozens of plugin directories or a large MATLAB toolbox tree.
        visited = {};
        owner = '';
        while ~isempty(directory)
            if isKey(cache, directory)
                owner = cache(directory);
                break
            end
            visited{end + 1} = directory; %#ok<AGROW>
            if isfile(fullfile(directory, 'eeglab.m'))
                owner = directory;
                break
            end
            parent = fileparts(directory);
            if strcmp(parent, directory), break, end
            directory = parent;
        end
        for j = 1:numel(visited)
            cache(visited{j}) = owner;
        end
    end
end
