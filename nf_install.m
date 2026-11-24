function report = nf_install(varargin)
% NF_INSTALL  Check this NeuroFreq installation and configure its dependencies.
%
%   nf_install
%   report = nf_install
%   report = nf_install('Native','off')  % use MATLAB implementations only
%
% Run from the NeuroFreq folder, or call nf_install with that folder on the
% MATLAB path. This function uses its own location to identify the installation.
% It compares the full distribution with the current GitHub default branch.
% Replacing a distribution always requires permission. Your existing files,
% including uncommitted work and local data, stay together in a sibling backup.
%
% At the end, only the NeuroFreq root needs to be on the MATLAB path. Use
%   EEG = nf_run('nf_filter', EEG);
% or wrap a whole analysis with nf_run(@my_analysis, ...). A private EEGLAB
% and the selected dependencies are available for the duration of that call.
%
% The package list, EEGLAB logic and update transaction live in private/.
% To add a dependency, start with private/nf_catalog.m.

options = inputParser;
addParameter(options,'Native','auto',@(s) ischar(s) || (isstring(s) && isscalar(s)));
parse(options,varargin{:});
nativeMode = validatestring(options.Results.Native,{'auto','off'});
root = fileparts(mfilename('fullpath'));
fprintf('\nNeuroFreq installer\n  %s\n\n', root);
lock = nf_install_lock(root); %#ok<NASGU> % released on success, error or Ctrl+C

% Stage 1. Check GitHub before changing paths or configuring dependencies.
[status, detail] = nf_check_repository(root);
report = struct('root', root, 'repository', status, 'detail', detail, ...
    'eeglab', '', 'dependencies', {cell(0, 2)}, 'features', [], ...
    'configured', false, 'pathSaved', false, 'nativeMode',nativeMode,'native',{cell(0,3)});
fprintf('%s\n', detail);
if strcmp(status, 'updated')
    % MATLAB may still have this version of nf_install in memory. Finish with
    % the configuration function shipped in the newly verified distribution.
    rehash;
    if ~isfile(fullfile(root, 'nf_install.m'))
        fprintf('The installed GitHub revision predates nf_install. The review installer is in the backup above.\n');
        fprintf('The new NeuroFreq root is on the path for this session.\n');
        return
    end
    if exist(fullfile(root, 'private', 'nf_configure.m'), 'file') ~= 2
        fprintf('Run nf_install again to continue with the updated installer.\n');
        return
    end
end

% Stage 2. Inspect existing installations before removing NF subdirectories
% from the path, so copies added by an earlier setup can still be detected.
if ~isfile(fullfile(root, 'private', 'nf_configure.m'))
    fprintf('Dependency setup requires the private folder shipped with NeuroFreq.\n');
    return
end
root = char(javaObject('java.io.File', root).getCanonicalPath());
report = nf_configure(root, report);
nf_repo_clean_path(root);
addpath(root, '-begin');
report.pathSaved = nf_offer_savepath(root);
end

function release = nf_install_lock(root)
% Separate MATLAB sessions can run the installer, but only one may change a
% particular NF folder at a time. The OS releases the lock on process exit.
% Keep its small, reusable file in tempdir, outside the compared distribution.
if ~usejava('jvm')
    error('NeuroFreq:InstallerJVM', ...
        'The installer requires a Java-enabled MATLAB or Octave with a JRE. Core functions can also be called directly with nf_run.');
end
key = char(javaObject('java.io.File',root).getCanonicalPath());
if ispc, key = lower(key); end
digest = javaMethod('getInstance', 'java.security.MessageDigest', 'SHA-256');
bytes = digest.digest(typecast(unicode2native(key, 'UTF-8'), 'int8'));
key = lower(reshape(dec2hex(typecast(bytes, 'uint8'), 2).', 1, []));
file = javaObject('java.io.RandomAccessFile', ...
    fullfile(tempdir, ['neurofreq-install-' key '.lock']), 'rw');
channel = file.getChannel();
try
    handle = channel.tryLock();
    if isempty(handle)
        error('NeuroFreq:InstallerBusy', ...
            'Another installer is using %s. Wait for it to finish, then retry.', root);
    end
catch failure
    file.close();
    rethrow(failure)
end
release = onCleanup(@() file.close()); % closing the channel releases its lock
end

function saved = nf_offer_savepath(root)
saved = false;
% savepath writes the user's entire current pathdef. Make that operation
% explicit. Choosing the session-only option still saves the NF configuration.
answer = nf_choose(['Keep this path for future sessions? Saving records ' ...
    'your current search path, including the NeuroFreq root.'], ...
    'Save path?', {'Save path', 'This session only'}, 'This session only');
if strcmp(answer, 'Save path')
    saved = savepath == 0;
    if ~saved
        warning('NeuroFreq:SavePath', ...
            'MATLAB could not save its path. Add this folder at startup: %s', root);
    end
else
    fprintf('For future sessions, add this folder to your MATLAB path: %s\n', root);
end
end

%% GitHub verification and non-destructive replacement
% Keep the repair code in this file. It must work even if private/ is missing.
function [status, detail] = nf_check_repository(root)
% GitHub's default branch is the reference distribution. Pin its commit once
% so that a push during installation cannot change which files we compare.
api = 'https://api.github.com/repos/erawls-neuro/NeuroFreq';
status = 'unverified';
if ~nf_repo_recognized(root)
    error('NeuroFreq:Location',['%s is not an active NeuroFreq installation folder. ' ...
        'Run nf_install from the intended distribution.'],root);
end
root = char(javaObject('java.io.File',root).getCanonicalPath());
try
    options = weboptions('ContentType', 'json', 'Timeout', 45, ...
        'UserAgent', 'NeuroFreq MATLAB installer');
    head = webread([api '/commits?per_page=1'], options);
    if ischar(head), head = jsondecode(head); end
    if iscell(head), head = head{1}; end
    revision = head(1).sha;
    tree = webread([api '/git/trees/' head(1).commit.tree.sha '?recursive=1'], options);
    if ischar(tree), tree = jsondecode(tree); end
    if tree.truncated || ~strcmp(tree.sha, head(1).commit.tree.sha)
        error('NeuroFreq:IncompleteTree', 'GitHub did not return the complete requested tree.');
    end
    fprintf('Comparing %s with GitHub commit %s ...\n', root, revision);
    [changes, hasInstaller] = nf_repo_compare(root, tree.tree);
catch err
    detail = sprintf('Could not verify NeuroFreq against GitHub. %s', err.message);
    return
end

if isempty(changes)
    status = 'current';
    detail = sprintf('The local distribution matches GitHub commit %s.', revision);
    return
end
fprintf('Differences relative to %s:\n', revision(1:12));
fprintf('  %s\n', changes{:});
question = sprintf(['Found %d differences from the latest NeuroFreq. Replace this folder?\n\n%s\n\n' ...
    'Your entire current folder (including personal files and Git history) will be kept beside it ' ...
    'with a _deprecated suffix. See the Command Window for the differences.'], numel(changes),root);
if ~hasInstaller
    % A review/development checkout can contain a newer installer than the
    % published branch. Make the consequence explicit before replacing it.
    question = sprintf(['%s\n\nThe published version predates nf_install. Updating will ' ...
        'archive the installer too. Choose Keep this copy to configure this checkout.'],question);
end
% Confirm locally even if the private installer helpers are missing.
choice = 'Keep this copy';
if (usejava('desktop') && usejava('awt')) || ...
        (exist('OCTAVE_VERSION','builtin') && isguirunning())
    choice = questdlg(question,'Update NeuroFreq?','Update','Keep this copy','Keep this copy');
end
if ~strcmp(choice,'Update')
    status = 'kept';
    detail = 'Kept the current NeuroFreq distribution.';
    if ~hasInstaller
        status = 'preview';
        detail = 'Kept this checkout. The published distribution predates nf_install.m.';
    end
    return
end

try
    % Keep staging on the same volume as the installation. Both moves below
    % are then directory renames, even when the system temp folder is elsewhere.
    scratch = tempname(fileparts(root));
    if ~javaObject('java.io.File',scratch).mkdir()
        error('NeuroFreq:Stage', 'Could not create a staging directory beside %s.', root);
    end
    cleanup = onCleanup(@() nf_repo_remove_stage(scratch)); %#ok<NASGU>
    fprintf('Downloading and verifying NeuroFreq %s ...\n', revision(1:12));
    archive = fullfile(scratch,'NeuroFreq.zip');
    url = [api '/zipball/' revision];
    if exist('OCTAVE_VERSION','builtin')
        urlwrite(url,archive);
    else
        websave(archive,url,weboptions('Timeout',120,'UserAgent','NeuroFreq MATLAB installer'));
    end
    nf_repo_unpack(archive,scratch);
    folders = dir(scratch);
    folders = folders([folders.isdir] & ~ismember({folders.name},{'.','..'}));
    if numel(folders) ~= 1
        error('NeuroFreq:Archive', 'Expected one NeuroFreq directory in the GitHub ZIP.');
    end
    staged = fullfile(scratch,folders.name);
    differences = nf_repo_compare(staged,tree.tree);
    if ~isempty(differences)
        error('NeuroFreq:Archive', 'The downloaded release failed verification (%s).', differences{1});
    end
    % Recheck just before swapping. A newly created or removed user file merits
    % another review. All bytes still go into the backup, including local edits.
    if ~isequal(changes,nf_repo_compare(root,tree.tree))
        error('NeuroFreq:Changed', 'The local file inventory changed during the download. Run nf_install again.');
    end
    backup = nf_repo_swap(root,staged);
    status = 'updated';
    detail = sprintf('Installed NeuroFreq %s. Previous distribution: %s',revision(1:12),backup);
catch err
    if strcmp(err.identifier,'NeuroFreq:Rollback'), rethrow(err); end
    detail = sprintf('Could not finish the NeuroFreq update. %s',err.message);
end
end


function nf_repo_unpack(archive, folder)
% The extraction destination stays within the dedicated staging directory.
% Keep this check local so an incomplete installation can still repair itself.
zipFile = javaObject('java.util.zip.ZipFile', archive);
closeZip = onCleanup(@() zipFile.close());
entries = zipFile.entries();
while entries.hasMoreElements()
    entry = entries.nextElement();
    name = strrep(char(entry.getName()), '\', '/');
    if ~isempty(regexp(name,'(^/|^[a-zA-Z]:|(^|/)\.\.(/|$))','once'))
        error('NeuroFreq:Archive','Invalid release archive member: %s.',name);
    end
end
clear closeZip
unzip(archive,folder);
end


function recognized = nf_repo_recognized(root)
% A copied installer in Downloads or the user's home must never rename that
% directory. Accept a partial distribution with a couple of distinctive NF
% files, or an explicitly named NeuroFreq directory being bootstrapped.
item = javaObject('java.io.File',root);
canonical = char(item.getCanonicalPath());
home = char(javaMethod('getProperty','java.lang.System','user.home'));
[~,name] = fileparts(canonical);
recognized = false;
if ~isempty(regexpi(name,'_deprecated(_[0-9]+)?$','once')) || ...
        isempty(item.getCanonicalFile().getParent()) || ...
        javaMethod('isSymbolicLink','java.nio.file.Files',item.toPath()) || ...
        nf_repo_inside(home,canonical) || nf_repo_inside(matlabroot,canonical)
    return
end
markers = {'nf_viewerapp.mlapp','nf_logo.png','timefrequency/nf_cwt.m', ...
    'dependencies/readme_dependencies.txt','setutilities/nf_epoch.m'};
markerCount = sum(cellfun(@(p) isfile(fullfile(root,p)),markers));
namedNF = ~isempty(regexpi(name,'^NeuroFreq([-_].*)?$','once'));
recognized = markerCount >= 2 || namedNF;
end


function [changes, hasInstaller] = nf_repo_compare(root,entries)
% Compare actual bytes. Neither a Git checkout nor Git on PATH is required.
% GitHub sometimes omits fields on directory entries; MATLAB's JSON decoder
% may therefore return a cell array instead of a homogeneous struct array.
if ~iscell(entries), entries = num2cell(entries); end
expected = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(entries)
    e = entries{i};
    if strcmp(e.type,'tree'), continue; end
    if ~strcmp(e.type,'blob') || ~ismember(e.mode,{'100644','100755'})
        error('NeuroFreq:TreeType','Unsupported Git entry at %s (mode %s).',e.path,e.mode);
    end
    expected(e.path) = e.sha;
    % Register ancestor directories too. Empty string means "directory".
    for slash = reshape(find(e.path == '/'),1,[])
        expected(e.path(1:slash-1)) = '';
    end
end
if expected.Count == 0 || isKey(expected,'dependencies/_nf_local')
    error('NeuroFreq:Tree','Invalid or conflicting upstream file inventory.');
end
hasInstaller = isKey(expected,'nf_install.m');
changes = nf_repo_walk(root,'',expected);
missing = keys(expected);
for i = 1:numel(missing)
    if ~isempty(expected(missing{i}))
        changes{end+1,1} = ['Missing: ' missing{i}]; %#ok<AGROW>
    end
end
changes = sort(changes); % Stable output and a reliable second inventory check.
end


function changes = nf_repo_walk(folder,relative,expected)
% Walk known directories only. An unexpected directory is reported once,
% regardless of how many personal datasets or linked folders it contains.
changes = cell(0,1);
directory = javaObject('java.io.File',folder);
if ~directory.isDirectory() || ~directory.canRead()
    error('NeuroFreq:Read','Cannot read directory %s.',folder);
end
files = directory.listFiles();
for i = 1:numel(files)
    item = files(i);
    key = [relative char(item.getName())];
    if ismember(key,{'.git','dependencies/_nf_local'}), continue; end
    if ~isKey(expected,key)
        changes{end+1,1} = ['Additional: ' key]; %#ok<AGROW>
        continue
    end
    sha = expected(key);
    remove(expected,key);
    if javaMethod('isSymbolicLink','java.nio.file.Files',item.toPath())
        changes{end+1,1} = ['Local symbolic link: ' key]; %#ok<AGROW>
    elseif isempty(sha) && item.isDirectory()
        nested = nf_repo_walk(char(item.getPath()),[key '/'],expected);
        changes = [changes; nested]; %#ok<AGROW>
    elseif isempty(sha) || ~item.isFile()
        changes{end+1,1} = ['Wrong file type: ' key]; %#ok<AGROW>
    elseif ~strcmp(nf_repo_sha(char(item.getPath())),sha)
        changes{end+1,1} = ['Modified: ' key]; %#ok<AGROW>
    end
end
end


function sha = nf_repo_sha(filename)
% A Git blob hashes "blob <byte count>", a NUL, then the original bytes.
% Chunk reads keep memory use independent of the largest example dataset.
[fid,message] = fopen(filename,'rb');
if fid < 0, error('NeuroFreq:Read','%s: %s',filename,message); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-1');
header = [sprintf('blob %d',javaObject('java.io.File',filename).length()) char(0)];
digest.update(typecast(uint8(header),'int8'));
while true
    bytes = fread(fid,1024*1024,'*uint8');
    if isempty(bytes), break; end
    digest.update(typecast(bytes,'int8'));
end
[message,number] = ferror(fid);
if number ~= 0, error('NeuroFreq:Read','%s: %s',filename,message); end
sha = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end


function backup = nf_repo_swap(root,staged)
% Once the new tree has been verified, preserve the *whole* old tree.
% Never merge directories or overwrite a backup. Files.move without
% REPLACE_EXISTING fails if the destination has been claimed by another run.
backup = [root '_deprecated'];
suffix = 1;
while javaObject('java.io.File',backup).exists()
    suffix = suffix + 1;
    backup = sprintf('%s_deprecated_%d',root,suffix);
end
originalPath = path;
originalFolder = pwd;
rootFile = javaObject('java.io.File',root);
inside = nf_repo_inside(originalFolder,root);
archived = false;
installed = false;
try
    nf_repo_clean_path(root);
    if inside, cd(char(rootFile.getParent())); end
    nf_repo_move(root,backup);
    archived = true;
    nf_repo_move(staged,root);
    installed = true;
    addpath(root,'-begin');
    if inside, cd(root); end
    rehash;
catch failure
    try
        if installed
            if nf_repo_inside(pwd,root), cd(char(rootFile.getParent())); end
            nf_repo_move(root,staged);
        end
        if archived, nf_repo_move(backup,root); end
        cd(originalFolder);
        path(originalPath);
        rehash;
    catch recovery
        error('NeuroFreq:Rollback',['Update failed: %s\nRecovery also failed: %s\n' ...
            'Inspect %s and %s before retrying.'],failure.message,recovery.message,root,backup);
    end
    rethrow(failure);
end
end


function nf_repo_move(source,target)
opts = javaArray('java.nio.file.CopyOption',0);
javaMethod('move','java.nio.file.Files',javaObject('java.io.File',source).toPath(),javaObject('java.io.File',target).toPath(),opts);
end


function nf_repo_clean_path(root)
% Remove paths into this NF tree and its numbered sibling backups only.
[parent,name] = fileparts(root);
pattern = ['^' regexptranslate('escape',name) '_deprecated(_[0-9]+)?$'];
entries = strsplit(path,pathsep);
removeEntries = false(size(entries));
for i = 1:numel(entries)
    entry = entries{i};
    if isempty(entry) || strcmp(entry,'.'), continue; end
    item = javaObject('java.io.File',entry);
    if ~item.isAbsolute(), item = javaObject('java.io.File',fullfile(pwd,entry)); end
    full = char(item.getCanonicalPath());
    removeEntry = nf_repo_inside(full,root);
    if ~removeEntry && nf_repo_inside(full,parent) && numel(full) > numel(parent)
        firstChild = fullfile(parent,'x');
        relative = strsplit(full(numel(firstChild):end),filesep);
        if ispc
            removeEntry = ~isempty(regexpi(relative{1},pattern,'once'));
        else
            removeEntry = ~isempty(regexp(relative{1},pattern,'once'));
        end
    end
    removeEntries(i) = removeEntry;
end
if any(removeEntries), path(strjoin(entries(~removeEntries),pathsep)); end
end


function nf_repo_remove_stage(scratch)
% Only delete the private staging directory created by this invocation.
if isfolder(scratch)
    [ok,message] = rmdir(scratch,'s');
    if ~ok, warning('NeuroFreq:Cleanup','Could not remove %s: %s',scratch,message); end
end
end

function inside = nf_repo_inside(name, root)
% Compare full filesystem path components (also handles a Windows drive root).
item = javaObject('java.io.File', name);
base = javaObject('java.io.File', root);
inside = item.getCanonicalFile().toPath().startsWith(base.getCanonicalFile().toPath());
end
