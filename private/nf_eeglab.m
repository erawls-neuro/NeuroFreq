function [eeglabRoot, privateEEGLAB, detail] = nf_eeglab(root)
% Find EEGLAB without running its startup script or changing the active path.
% A selected installation is used by nf_run. A private installation is a
% complete SCCN release kept alongside the NF dependency downloads.
eeglabRoot = '';
privateEEGLAB = false;
detail = 'EEGLAB skipped. Functions that call EEGLAB may be unavailable.';
fprintf('\n[EEGLAB] Looking for eeglab.m and checking the current SCCN release.\n');
cache = fullfile(root,'dependencies','_nf_local','eeglab');
current = struct('version','','url','');
try
    entries = nf_eeglab_releases();
    release = entries(find(strcmpi({entries.plugin},'EEGLAB'),1));
    current.version = release.curversion;
    current.url = regexprep(release.link,'^http:','https:');
    fprintf('[EEGLAB] SCCN currently offers %s.\n',current.version);
catch err
    fprintf('[EEGLAB] Could not determine the current release: %s\n',err.message);
end

% Remember the last selected distribution, including one that was never
% put on the user's global path. Store private paths relative to NF so the
% whole tree can move. A newer copy can be selected at the prompt below.
candidate = '';
settings = fullfile(root,'dependencies','_nf_local','runtime.mat');
if isfile(settings)
    try
        prior = load(settings,'config');
        candidate = prior.config.eeglabRoot;
        if ~isempty(candidate)
            selectedDir = javaObject('java.io.File',candidate);
            if ~selectedDir.isAbsolute(), candidate = fullfile(root,candidate); end
        end
    catch err
        fprintf('[EEGLAB] Could not read the last selection: %s\n',err.message);
    end
end
if ~nf_valid_eeglab(candidate), candidate = fileparts(which('eeglab')); end
if ~nf_valid_eeglab(candidate), candidate = nf_cached_eeglab(cache); end
if isempty(candidate)
    choice = nf_choose('No complete EEGLAB installation was found on the MATLAB path.', ...
        'Find EEGLAB',{'Locate EEGLAB','Install NF-only copy','Skip EEGLAB'},'Skip EEGLAB');
else
    [label,~] = nf_eeglab_version(candidate);
    message = sprintf('Found EEGLAB %s at\n%s\nWhich installation should NF use?',label,candidate);
    fprintf('[EEGLAB] Found %s at %s\n',label,candidate);
    choice = nf_choose(message,'Select EEGLAB', ...
        {'Use this copy','Locate EEGLAB','Install NF-only copy'},'Use this copy');
end
if strcmp(choice,'Locate EEGLAB')
    candidate = nf_pick_folder('Select the EEGLAB folder containing eeglab.m',pwd);
    if ~nf_valid_eeglab(candidate)
        fprintf('[EEGLAB] No folder with eeglab.m and its core functions was selected.\n');
        return
    end
elseif strcmp(choice,'Install NF-only copy')
    [installed,downloadDetail] = nf_private_eeglab(cache,current);
    if ~isempty(installed)
        eeglabRoot = installed; privateEEGLAB = true; detail = downloadDetail;
        return
    end
    detail = downloadDetail;
    if isempty(candidate), return; end
    % An existing installation may still be used after a failed download.
    fallback = nf_choose(sprintf('%s\nUse the existing EEGLAB at %s?',detail,candidate), ...
        'EEGLAB download',{'Use existing','Skip EEGLAB'},'Skip EEGLAB');
    if ~strcmp(fallback,'Use existing'), return; end
elseif ~strcmp(choice,'Use this copy')
    return
end

[label,year] = nf_eeglab_version(candidate);
today = datevec(now);
isOld = isfinite(year) && today(1)-year > 3;
if isOld || ~isfinite(year)
    if isOld
        ageText = sprintf('EEGLAB %s was released in %d or earlier, more than three calendar years ago.',label,year);
    else
        ageText = sprintf('The release year of EEGLAB %s could not be verified.',label);
    end
    message = sprintf(['%s\nNF uses recent EEGLAB methods. Some NF functions may fail with this copy.\n' ...
        'An NF-only release leaves your regular EEGLAB installation in place.'],ageText);
    warning('NeuroFreq:EEGLABAge','%s',message);
    choices = {'Install NF-only copy','EEGLAB updater','Keep this copy'};
    if nf_under(candidate,cache), choices(2) = []; end
    answer = nf_choose(message,'EEGLAB compatibility',choices,'Keep this copy');
    if strcmp(answer,'Install NF-only copy')
        [installed,detail] = nf_private_eeglab(cache,current);
        if ~isempty(installed)
            eeglabRoot = installed; privateEEGLAB = true;
            return
        end
        fallback = nf_choose(sprintf('%s\nKeep the selected EEGLAB for NF?',detail), ...
            'EEGLAB download',{'Keep this copy','Skip EEGLAB'},'Keep this copy');
        if ~strcmp(fallback,'Keep this copy'), return; end
    elseif strcmp(answer,'EEGLAB updater')
        candidate = nf_update_eeglab(candidate);
        if ~nf_valid_eeglab(candidate)
            detail = 'EEGLAB updater finished without a usable installation. Run nf_install to select a copy.';
            return
        end
        [label,year] = nf_eeglab_version(candidate);
    end
end
eeglabRoot = candidate;
privateEEGLAB = nf_under(candidate,cache);
detail = sprintf('Selected EEGLAB %s at %s.',label,candidate);
if ~isfinite(year) || today(1)-year > 3
    detail = [detail ' Some recent EEGLAB methods may be unavailable.'];
end
end


function [candidate,detail] = nf_private_eeglab(cache,current)
% Install side-by-side. An interrupted download or invalid ZIP never
% replaces an existing distribution (including a user-modified copy).
candidate = '';
try
    if isempty(current.version)
        candidate = nf_cached_eeglab(cache);
        if isempty(candidate), error('No SCCN release metadata or usable cached release.'); end
        detail = sprintf('Using the cached NF-only EEGLAB at %s. The latest release could not be verified.',candidate);
        return
    end
    if isempty(regexp(current.version,'^\d+[A-Za-z0-9._-]*$','once')) || ...
            ~strncmpi(current.url,'https://',8)
        error('The SCCN catalog did not provide a valid release version and HTTPS download.');
    end
    dest = fullfile(cache,current.version);
    if ~nf_complete_private_eeglab(dest,current.version)
        if exist(dest,'dir') || exist(dest,'file')
            error('The existing folder %s is incomplete. It was preserved.',dest);
        end
        if ~isfolder(cache), mkdir(cache); end
        stage = tempname(cache);
        mkdir(stage);
        cleanup = onCleanup(@()rmdir(stage,'s')); %#ok<NASGU>
        fprintf('[EEGLAB] Downloading the complete SCCN distribution %s.\n',current.version);
        unpacked = nf_fetch_package(current.url,stage);
        source = unpacked;
        if ~nf_valid_eeglab(source)
            % SCCN packages its distribution inside a single top-level folder.
            files = dir(unpacked);
            options = files([files.isdir] & ~ismember({files.name},{'.','..'}));
            valid = arrayfun(@(folder) nf_valid_eeglab(fullfile(unpacked,folder.name)),options);
            if nnz(valid) ~= 1, error('The download contains no unique complete EEGLAB installation.'); end
            source = fullfile(unpacked,options(valid).name);
        end
        if ~nf_complete_private_eeglab(source,current.version)
            error('The extracted EEGLAB version or bundled plugins do not match the release.');
        end
        % Use an exclusive rename, so two installers cannot nest a second
        % copy inside the destination or overwrite it.
        sourceDir = javaObject('java.io.File',source);
        if ~sourceDir.renameTo(javaObject('java.io.File',dest))
            if ~nf_complete_private_eeglab(dest,current.version), error('Could not install %s.',dest); end
        end
    end
    candidate = dest;
    detail = sprintf('Selected NF-only EEGLAB %s at %s.',current.version,dest);
catch err
    detail = sprintf('NF-only EEGLAB installation failed: %s',err.message);
    warning('NeuroFreq:EEGLABDownload','%s',detail);
end
end


function candidate = nf_cached_eeglab(cache)
% Sort version components numerically (2024.10 should follow 2024.2).
candidate = '';
if ~isfolder(cache), return; end
folders = dir(cache);
folders = folders([folders.isdir]);
keys = zeros(0,4); roots = {};
for k = 1:numel(folders)
    folder = fullfile(cache,folders(k).name);
    if isempty(regexp(folders(k).name,'^\d','once')) || ~nf_complete_private_eeglab(folder,folders(k).name), continue; end
    parts = regexp(folders(k).name,'\d+','match');
    key = zeros(1,4);
    n = min(4,numel(parts));
    key(1:n) = str2double(parts(1:n));
    keys(end+1,:) = key; roots{end+1} = folder; %#ok<AGROW>
end
if ~isempty(roots)
    [~,order] = sortrows(keys);
    candidate = roots{order(end)};
end
end


function yes = nf_valid_eeglab(folder)
% An EEGLAB root needs both its entry point and the core functions that NF
% calls. Checking just the folder name would accept incomplete source trees.
yes = ischar(folder) && ~isempty(folder) && isfile(fullfile(folder,'eeglab.m')) && ...
    isfile(fullfile(folder,'functions','adminfunc','eeg_checkset.m')) && ...
    isfile(fullfile(folder,'functions','popfunc','pop_loadset.m'));
end


function yes = nf_complete_private_eeglab(folder,version)
% An SCCN release ZIP includes plugins. Its core version should also match
% the catalog label (or its retained version folder during an offline run).
yes = nf_valid_eeglab(folder);
if yes
    [installed,~] = nf_eeglab_version(folder);
    yes = strcmpi(installed,version) && ...
        ~isempty(dir(fullfile(folder,'plugins','*','eegplugin_*.m')));
end
end


function [label,year] = nf_eeglab_version(folder)
% Read metadata from the *selected* installation. Calling eeg_getversion
% on the current path can instead report another copy. Modification times
% (including the date the ZIP was unpacked) do not establish release age.
label = 'unknown'; year = NaN;
versionFile = fullfile(folder,'functions','adminfunc','eeg_getversion.m');
if isfile(versionFile)
    source = fileread(versionFile);
    value = regexp(source,'(?m)^\s*vers\s*=\s*''([^'']+)''','tokens','once');
    if ~isempty(value), label = value{1}; end
    value = regexp(source,'(?m)^\s*releaseDate\s*=\s*''\d+-[a-zA-Z]+-(\d{2,4})[ '']','tokens','once');
    if ~isempty(value)
        year = str2double(value{1});
        if year >= 90 && year < 100, year = year + 1900; end
        if year < 90, year = year + 2000; end
    end
end
contents = fullfile(folder,'Contents.m');
if strcmp(label,'unknown') && isfile(contents)
    value = regexp(fileread(contents),'(?im)^%\s*Version\s+([0-9][.\w-]*)','tokens','once');
    if ~isempty(value), label = value{1}; end
end
% Calendar-versioned releases encode their own year more reliably than
% copied source dates. Numbered releases 4 through 14 predate EEGLAB 2019.
value = regexp(label,'^((?:19|20)\d{2})(?:\D|$)','tokens','once');
if ~isempty(value)
    year = str2double(value{1});
elseif ~isfinite(year) && ~isempty(regexp(label,'^(?:[4-9]|1[0-4])(?:\.|$)','once'))
    year = 2018; % latest possible year for purposes of the age warning
end
end


function selected = nf_update_eeglab(original)
% This option delegates the update to the same command used by EEGLAB's
% Help menu. Its own dialog controls the destination, path and old version.
selected = original;
startingPath = path; startingFolder = pwd;
restoreFolder = onCleanup(@()nf_restore_eeglab_cwd(startingFolder)); %#ok<NASGU>
try
    cd(original);
    needed = nf_eeglab_paths(original);
    addpath(needed{:},'-begin');
    updater = which('eeglab_update');
    if isempty(updater) || ~nf_under(updater,original)
        error('This EEGLAB installation has no built-in eeglab_update. Select an NF-only release or locate a newer distribution.');
    end
    if ~usejava('desktop')
        error('EEGLAB uses graphical dialogs for this update. Run nf_install in the MATLAB desktop.');
    end
    fprintf('[EEGLAB] Opening the selected distribution''s built-in updater.\n');
    eeglab_update;
    % Upstream changes to the new directory even if the user declines a
    % permanent path change. On cancellation the working directory is old.
    if nf_valid_eeglab(pwd) && ~strcmp(pwd,original)
        selected = pwd;
    else
        active = fileparts(which('eeglab'));
        if nf_valid_eeglab(active), selected = active; end
    end
catch err
    warning('NeuroFreq:EEGLABUpdater','EEGLAB updater: %s',err.message);
end
if strcmp(selected,original)
    path(startingPath);
    fprintf('[EEGLAB] No newly selected installation was detected.\n');
else
    % Remove temporary additions to the old release but respect actual path
    % changes the user chose in the EEGLAB updater (e.g., its new folder).
    before = strsplit(startingPath,pathsep);
    after = strsplit(path,pathsep);
    extra = intersect(setdiff(after,before),nf_eeglab_paths(original));
    if ~isempty(extra), rmpath(extra{:}); end
end
end


function nf_restore_eeglab_cwd(folder)
% The updater can optionally rename the user's original working directory.
while ~isfolder(folder)
    parent = fileparts(folder);
    if isempty(parent) || strcmp(parent,folder), return; end
    folder = parent;
end
cd(folder);
end
