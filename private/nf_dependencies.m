function [rows, scopePaths] = nf_dependencies(root, eeglabRoot, privateEEGLAB)
% A package has two independent locations: the copy kept with NF, and the
% provider NF will actually use. An external provider wins over the copy.
% The paths returned here are for nf_run; downloading never activates code.
catalog = nf_catalog();
rows = cell(numel(catalog), 2);
scopePaths = {};
% Temporarily hide old NF copies while searching for external providers.
% The installer commits the root-only path after the configuration is saved.
originalPath = path;
restore = onCleanup(@() path(originalPath)); %#ok<NASGU>
nf_remove_paths(fullfile(root, 'dependencies'));
prior = nf_prior_providers(root, eeglabRoot);
fprintf('  Checking MATLAB user folders for off-path installations.\n');
profile = nf_profile_probes(root, catalog);
for k = 1:numel(catalog)
    d = catalog(k);
    fprintf('\n[%d/%d] %s -- %s\n', k, numel(catalog), d.name, d.features);
    fprintf('  Searching MATLAB, the selected EEGLAB plugins, and NF dependencies.\n');
    rows{k,1} = d.name;
    try
        [selected, rows{k,2}] = nf_select_package(root, eeglabRoot, privateEEGLAB, d, prior, profile{k});
        scopePaths = unique([scopePaths selected], 'stable');
    catch problem
        rows{k,2} = ['Unavailable: ' problem.message];
    end
    fprintf('  %s\n', rows{k,2});
end
end


function [selected, note] = nf_select_package(root, eeglabRoot, privateEEGLAB, d, prior, profileFiles)
depot = fullfile(root, 'dependencies', '_nf_local', 'packages', d.key);
cache = nf_cached_package(root, depot, d);
isPlugin = strcmp(d.kind, 'plugin');
inEEGLAB = {};
if isPlugin && ~isempty(eeglabRoot)
    inEEGLAB = nf_find_package_in_children(fullfile(eeglabRoot, 'plugins'), d);
end
onPath = nf_installed_package(root, d);
selected = onPath;
origin = 'MATLAB path';
if isempty(selected)
    external = prior(~cellfun(@(p) nf_under(p, root), prior));
    previousFiles = cellfun(@(p) fullfile(p, [d.probes{1} '.m']), external, 'UniformOutput', false);
    previousFiles = previousFiles(cellfun(@isfile, previousFiles));
    selected = nf_provider_for_files(previousFiles, d);
    if ~isempty(selected)
        origin = 'previous NF selection';
    elseif ~isempty(cache) && any(ismember(cache, prior))
        selected = cache;
        origin = 'previous NF copy';
    end
end
% Older NF installations often placed packages directly in dependencies/.
% Recognize those layouts even after hiding their path entries. Keep the
% managed download separate so the user's existing files are preserved.
if isempty(selected)
    selected = nf_find_package_in_children(fullfile(root, 'dependencies'), d);
    if ~isempty(selected), origin = 'existing NF dependencies folder'; end
end
% A compatible active package remains the default. A private EEGLAB may
% supply a second copy; it is an explicit alternative and also a local cache.
if ~isempty(inEEGLAB) && isempty(onPath)
    selected = inEEGLAB;
    origin = 'selected EEGLAB';
end
if privateEEGLAB && ~isempty(inEEGLAB) && isempty(cache) ...
        && nf_under(inEEGLAB{1}, fullfile(root, 'dependencies'))
    cache = inEEGLAB;
end
if ~isempty(selected)
    fprintf('  Found on %s: %s\n', origin, selected{1});
end

% A plugin on the ordinary MATLAB path may still be missing from the selected
% EEGLAB distribution. Offer to install there explicitly or keep the provider.
missingPlugin = isPlugin && ~isempty(eeglabRoot) && isempty(inEEGLAB);
if isempty(inEEGLAB) || ~isempty(onPath) || ~isequal(selected, inEEGLAB)
    if isempty(selected)
        options = {'Use NF copy', 'Locate local folder', 'Skip'};
        default = 'Skip';
    else
        options = {'Use existing provider', 'Use NF copy', 'Locate local folder', 'Skip'};
        default = 'Use existing provider';
    end
    if ~isempty(inEEGLAB) && ~isequal(inEEGLAB, selected)
        fprintf('  Also available in the selected EEGLAB: %s\n', inEEGLAB{1});
        options = [{'Use selected EEGLAB plugin'}, options];
    end
    located = {};
    if isempty(selected) && ~isempty(profileFiles)
        located = nf_provider_for_files(profileFiles, d);
        if ~isempty(located)
            fprintf('  Found off path in your MATLAB user folder: %s\n', located{1});
            options = [{'Use detected folder'}, options];
        end
    end
    if missingPlugin
        options = [{'Install using EEGLAB'}, options];
    end
    if ~isempty(cache)
        fprintf('  Compatible NF copy: %s\n', cache{1});
    end
    prompt = sprintf('%s provides %s. Choose which copy NF should use.', d.name, d.features);
    if missingPlugin
        prompt = sprintf('%s\nThe selected EEGLAB has no compatible %s plugin in %s.', ...
            prompt, d.name, fullfile(eeglabRoot, 'plugins'));
    end
    action = nf_choose(prompt, d.name, options, default);
    switch action
        case 'Skip'
            selected = {}; note = 'Skipped by user.'; return
        case 'Use selected EEGLAB plugin'
            selected = inEEGLAB;
            origin = 'selected EEGLAB';
        case 'Use detected folder'
            selected = located;
            origin = 'MATLAB user folder';
        case 'Locate local folder'
            folder = nf_pick_folder(['Locate the ' d.name ' package root'], nf_user_folder());
            if isempty(folder)
                note = 'Folder selection canceled.';
                return
            end
            selected = nf_package_paths(folder, d, true);
            if isempty(selected)
                error('NeuroFreq:IncompletePackage', ...
                    'The selected folder does not supply %s.', strjoin(d.probes, ', '));
            end
            origin = 'selected folder';
        case 'Install using EEGLAB'
            fprintf('  Installing %s into %s through EEGLAB.\n', d.name, eeglabRoot);
            selected = nf_install_plugin(eeglabRoot, d);
            origin = 'selected EEGLAB';
        case 'Use NF copy'
            selected = {};
            origin = 'NF copy';
    end
end

% Keep a complete copy even when a different provider was selected.
if isempty(cache)
    try
        cache = nf_acquire_package(depot, d);
    catch problem
        if isempty(selected), rethrow(problem); end
        note = sprintf('Using %s (%s). NF copy unavailable: %s', ...
            selected{1}, origin, problem.message);
        warning('NeuroFreq:MissingDependencyCopy', '%s', note);
        return
    end
end
if isempty(selected)
    selected = cache;
    origin = 'NF copy';
end
if ~isequal(selected, cache)
    fprintf('  NF copy: %s\n', cache{1});
    if isequal(selected, onPath)
        warning('NeuroFreq:ExternalDependency', ...
            ['%s was already installed on the MATLAB path. Its copy in /dependencies ' ...
             'is off the active path. If errors arise, rerun nf_install and choose ' ...
             'Use NF copy for this dependency.'], d.name);
    end
end
note = sprintf('Using %s (%s).', selected{1}, origin);
end


function paths = nf_installed_package(nfRoot, d)
% which -all includes shadowed installations. Require the probes to belong
% to one complete installation; matching unrelated files is insufficient.
paths = {};
if exist('OCTAVE_VERSION','builtin')
    % Octave's which ignores '-all'. The load-path primitive also finds
    % shadowed providers and native Octave implementations.
    probe = d.probes{1};
    found = file_in_loadpath({[probe '.m'],[probe '.oct'],[probe '.' mexext]},'all');
else
    found = which(d.probes{1}, '-all');
end
if ischar(found), found = cellstr(found); end
found = found(~cellfun(@(p) isempty(p) || nf_under(p, nfRoot) || nf_in_backup(p, nfRoot), found));
paths = nf_provider_for_files(found, d);
end


function paths = nf_provider_for_files(found, d)
paths = {};
for k = 1:numel(found)
    candidate = fileparts(found{k});
    % A probe may live below a package root (e.g. HAPPE or PREP). Stop
    % climbing at the maximum depth named by the catalog.
    depth = max(cellfun(@(s) numel(strfind(strrep(s, '\', '/'), '/')) + 1, d.paths));
    for level = 0:depth
        if isempty(candidate), break; end
        paths = nf_package_paths(candidate, d, false);
        if ~isempty(paths) && any(cellfun(@(p) strcmp(nf_canonical(p), ...
                nf_canonical(fileparts(found{k}))), paths))
            return
        end
        paths = {};
        parent = fileparts(candidate);
        if strcmp(parent, candidate), break; end
        candidate = parent;
    end
end
end


function paths = nf_find_package_in_children(parent, d)
% Check direct children of an EEGLAB plugins folder or NF's legacy
% dependencies folder. Only the package's catalogued code paths qualify.
paths = {};
children = dir(parent);
for k = 1:numel(children)
    if ~children(k).isdir || startsWith(children(k).name, '.') || ...
            strcmp(children(k).name, '_nf_local'), continue; end
    paths = nf_package_paths(fullfile(children(k).folder, children(k).name), d, true);
    if ~isempty(paths), return; end
end
end


function paths = nf_cached_package(root, depot, d)
paths = {};
if ~isempty(d.folder)
    paths = nf_package_paths(fullfile(root, d.folder), d, true);
    if ~isempty(paths), return; end
end
copies = dir(fullfile(depot, 'copy-*'));
[~, order] = sort([copies.datenum], 'descend');
for k = reshape(order, 1, [])
    paths = nf_package_paths(fullfile(copies(k).folder, copies(k).name), d, true);
    if ~isempty(paths), return; end
end
end


function paths = nf_package_paths(folder, d, allowWrapper)
% A catalog entry lists only its public code directories. In particular,
% never use genpath on HAPPE or MADE: they carry other toolboxes inside.
paths = {};
if ~isfolder(folder), return; end
for k = 1:numel(d.paths)
    rel = d.paths{k};
    if strcmp(rel, '.')
        paths{end+1} = folder; %#ok<AGROW>
    elseif any(rel == '*')
        hits = dir(fullfile(folder, rel));
        for j = reshape(find([hits.isdir]), 1, [])
            if ~startsWith(hits(j).name, '.')
                paths{end+1} = fullfile(hits(j).folder, hits(j).name); %#ok<AGROW>
            end
        end
    elseif isfolder(fullfile(folder, rel))
        paths{end+1} = fullfile(folder, rel); %#ok<AGROW>
    end
end
if isfield(d, 'requiredFiles') && ...
        ~all(cellfun(@(asset) isfile(fullfile(folder, asset)), d.requiredFiles))
    paths = {};
end
for k = 1:numel(d.probes)
    if ~any(cellfun(@(p) nf_has_function(p, d.probes{k}), paths))
        paths = {};
        break
    end
end
% Both GEDAI 1.7 and 1.8 expose the same MATLAB function signature.
% Read its declared plugin release without launching the plugin or its GUI.
if ~isempty(paths) && strcmpi(d.plugin, 'GEDAI') && ...
        isfield(d, 'pluginVersion') && ~isempty(d.pluginVersion)
    marker = fullfile(folder, 'eegplugin_GEDAI.m');
    version = {};
    if isfile(marker)
        version = regexp(fileread(marker), ...
            '(?m)^\s*vers\s*=\s*''GEDAI v([0-9.]+)[\s'']', 'tokens', 'once');
    end
    if isempty(version) || ~strcmp(version{1}, d.pluginVersion), paths = {}; end
end
if ~isempty(paths) || ~allowWrapper, return; end
% Distribution ZIPs commonly contain one top-level release folder.
items = dir(folder);
items = items([items.isdir] & ~startsWith({items.name}, '.') & ~strcmp({items.name}, '__MACOSX'));
if numel(items) == 1
    paths = nf_package_paths(fullfile(folder, items(1).name), d, false);
end
end


function yes = nf_has_function(folder, name)
base = fullfile(folder, name);
yes = isfile([base '.m']) || isfile([base '.' mexext]);
if exist('OCTAVE_VERSION','builtin')
    yes = yes || isfile([base '.oct']);
else
    yes = yes || isfile([base '.p']);
end
end


function folder = nf_user_folder()
if exist('OCTAVE_VERSION','builtin'), folder = pwd; return; end
folders = strsplit(userpath, pathsep);
folder = pwd;
for k = 1:numel(folders)
    if isfolder(folders{k}), folder = folders{k}; return; end
end
end


function paths = nf_prior_providers(root, eeglabRoot)
% Discover previous external selections without temporarily activating code.
% Cached providers are found separately and can be reused without network.
paths = {};
file = nf_runtime_file(root);
if ~isfile(file), return; end
try
    saved = load(file, 'config');
    previous = saved.config.scopePaths;
    priorEEGLAB = '';
    if isfield(saved.config, 'eeglabRoot') && ~isempty(saved.config.eeglabRoot)
        priorEEGLAB = nf_prior_absolute(saved.config.eeglabRoot, root);
    end
    for k = 1:numel(previous)
        candidate = nf_prior_absolute(previous{k}, root);
        if ~isempty(priorEEGLAB) && ~strcmp(eeglabRoot, priorEEGLAB) && nf_under(candidate, priorEEGLAB)
            continue  % A different EEGLAB was selected on this run.
        end
        if isfolder(candidate) && ~nf_in_backup(candidate, root)
            paths{end+1} = candidate; %#ok<AGROW>
        end
    end
catch problem
    fprintf('  Could not read previous providers: %s\n', problem.message);
end
end


function candidate = nf_prior_absolute(candidate, root)
if isempty(regexp(candidate, '^([A-Za-z]:[/\\]|[/\\])', 'once'))
    candidate = fullfile(root, candidate);
end
candidate = nf_canonical(candidate);
end


function found = nf_profile_probes(root, catalog)
% Search MATLAB's configured userpath once. Store only entry points from
% the catalog; an off-path installation is used only after explicit choice.
% Scan userpath, not the entire user's home or MATLAB's installation tree.
anchors = arrayfun(@(d) [d.probes{1} '.m'], catalog, 'UniformOutput', false);
found = repmat({{}}, size(catalog));
% Octave has no MATLAB userpath. Existing providers on its load path,
% cached selections and explicit folder selection are checked separately.
if exist('OCTAVE_VERSION','builtin'), return; end
folders = strsplit(userpath, pathsep);
for k = 1:numel(folders)
    folder = folders{k};
    if ~isfolder(folder) || nf_under(folder, root), continue; end
    try
        files = dir(fullfile(folder, '**', '*.m'));
    catch problem
        fprintf('  Could not inspect %s: %s\n', folder, problem.message);
        continue
    end
    for j = 1:numel(anchors)
        hits = files(strcmp({files.name}, anchors{j}));
        for h = 1:numel(hits)
            source = fullfile(hits(h).folder, hits(h).name);
            if ~nf_under(source, root) && ~nf_in_backup(source, root)
                found{j}{end+1} = source;
            end
        end
    end
end
end


function inside = nf_in_backup(candidate, root)
% Never reactivate a provider from the renamed NeuroFreq sibling, even
% when that folder still appears on a pre-existing MATLAB search path.
candidate = nf_canonical(candidate);
prefix = [root '_deprecated'];
if ispc
    candidate = lower(candidate);
    prefix = lower(prefix);
end
inside = startsWith(candidate, prefix) && ...
    ~isempty(regexp([candidate(numel(prefix)+1:end) filesep], '^(_[0-9]+)?[/\\]', 'once'));
end


function paths = nf_acquire_package(depot, d)
% Stage beside the destination so a failed network request, canceled
% transfer, or incomplete distribution cannot replace a working copy.
url = d.url;
if isempty(url) && strcmp(d.kind, 'plugin')
    release = nf_plugin_release(d);
    url = release.link;
end
if isempty(url)
    error('NeuroFreq:NoDownload', 'No redistributable download URL is available for %s.', d.name);
end
if ~strcmp(nf_canonical(depot), depot)
    error('NeuroFreq:LinkedDepot', 'The package destination points outside its expected location: %s.', depot);
end
if ~isfolder(depot), mkdir(depot); end
stage = tempname(depot);
mkdir(stage);
cleanup = onCleanup(@() nf_discard_stage(stage)); %#ok<NASGU>
fprintf('  Downloading %s from %s\n', d.name, url);
unpacked = nf_fetch_package(url, stage);
paths = nf_package_paths(unpacked, d, true);
if isempty(paths)
    error('NeuroFreq:IncompleteDownload', ...
        'Downloaded %s does not contain all required functions: %s.', d.name, strjoin(d.probes, ', '));
end
% Each successful download gets its own directory. An incomplete or edited
% older copy remains available for inspection even if this is a retry.
[~, suffix] = fileparts(stage);
destination = fullfile(depot, ['copy-' suffix]);
[ok, message] = movefile(unpacked, destination);
if ~ok, error('NeuroFreq:PackageMove', '%s', message); end
paths = nf_package_paths(destination, d, true);
end


function release = nf_plugin_release(d)
list = nf_eeglab_releases();
index = find(strcmpi({list.plugin}, d.plugin), 1);
if isempty(index)
    error('NeuroFreq:UnknownPlugin', '%s is absent from the EEGLAB release registry.', d.plugin);
end
release = list(index);
end


function paths = nf_install_plugin(eeglabRoot, d)
% EEGLAB's installer chooses its destination via which('eeglab'). Verify
% both entry points before asking it to write anything, then restore the
% caller's path even if the plugin manager or the download fails.
original = path;
restore = onCleanup(@() path(original)); %#ok<NASGU>
dirs = nf_eeglab_paths(eeglabRoot);
addpath(dirs{:}, '-begin');
if ~nf_under(which('eeglab'), eeglabRoot) || ~nf_under(which('plugin_install'), eeglabRoot)
    error('NeuroFreq:WrongEEGLAB', 'Cannot use the plugin manager from %s.', eeglabRoot);
end
release = nf_plugin_release(d);
if isfield(d, 'pluginVersion') && ~isempty(d.pluginVersion) && ...
        ~strcmp(regexprep(release.curversion, '^v', ''), d.pluginVersion)
    error('NeuroFreq:PluginVersion', ...
        'EEGLAB offers %s %s; NF requires %s. Choose Use NF copy for the compatible release.', ...
        d.name, release.curversion, d.pluginVersion);
end
pluginRoot = fullfile(eeglabRoot, 'plugins');
if ~isfolder(pluginRoot), mkdir(pluginRoot); end
if isfile(fullfile(pluginRoot, 'tmp.zip')) || isfolder(fullfile(pluginRoot, [release.plugin release.curversion]))
    error('NeuroFreq:PluginCollision', ...
        'EEGLAB already has a download or a folder for this release. Inspect %s or choose Use NF copy.', pluginRoot);
end
ok = plugin_install(nf_https(release.link), release.plugin, release.curversion, true);
if ok ~= 1
    error('NeuroFreq:PluginManager', 'EEGLAB could not install %s.', d.name);
end
paths = nf_find_package_in_children(fullfile(eeglabRoot, 'plugins'), d);
if isempty(paths)
    error('NeuroFreq:PluginIncomplete', 'EEGLAB did not install the required %s functions.', d.name);
end
end


function nf_discard_stage(folder)
if isfolder(folder)
    [ok, message] = rmdir(folder, 's');
    if ~ok, warning('NeuroFreq:TempCleanup', 'Could not remove %s: %s', folder, message); end
end
end
