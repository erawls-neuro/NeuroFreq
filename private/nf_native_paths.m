function [dirs,allowed] = nf_native_paths(root,native)
% Verify the installed ABI, current source (including shared headers), and
% compiled bytes before exposing a cache directory on the nf_run path.
dirs = {};
allowed = struct();
if ~strcmp(native.abi,nf_native_abi()) || ~strcmp(native.mode,'auto') || isempty(native.targets)
    return
end
cache = fullfile(root,'dependencies','_nf_local','native');
if ~strcmp(nf_canonical(cache),cache), return; end
persistent recent
if isempty(recent) || ~strcmp(recent.root,root) || ~isequal(recent.native,native)
    recent = struct('root',root,'native',native,'watch',[],'stamp',[]);
end

% Even a warm call checks actual bytes. This detects equal-size edits with
% preserved modification times (common when switching Git checkouts).
% Reuse the previous decision only if the same manifest, paths, headers,
% and content are present. Do the more expensive per-target diagnosis only
% when something changed; an unreadable or stale target then falls back.
if ~isempty(recent.watch)
    stamp = content_stamp(recent.watch);
    if ~isempty(stamp) && isequal(stamp,recent.stamp)
        dirs = recent.dirs;
        allowed = recent.allowed;
        return
    end
end

% Snapshot the candidate files before checking each independent receipt.
% A second snapshot below prevents a concurrent edit from being memorized
% as if it were the content we just verified.
[watch,firstStamp] = snapshot(root,native.targets);
for i = 1:numel(native.targets)
    item = native.targets(i);
    source = fullfile(root,item.source);
    binary = fullfile(root,item.binary);
    try
        if ~isfile(source) || ~isfile(binary) || ~nf_under(binary,cache), continue; end
        % MATLAB can resolve any executable sibling in a path directory.
        % Only publish an isolated generation containing this MEX + receipt.
        parent = javaObject('java.io.File',fileparts(binary));
        children = cellstr(parent.list());
        if ~isequal(sort(children(:)),sort({[item.name '.' mexext];'sha256.txt'})), continue; end
        inputs = nf_native_inputs(source);
        if ~strcmp(nf_native_digest(inputs,'NF source v1'),item.sourceHash) || ...
                ~strcmp(nf_native_digest({binary},'NF binary v1'),item.binaryHash)
            continue
        end
        allowed.(item.name) = nf_canonical(binary);
        dirs{end+1} = fileparts(allowed.(item.name)); %#ok<AGROW>
    catch
        % One broken cached target must not disable other kernels or MATLAB.
        continue
    end
end

% Create an aggregate fingerprint only after every file has independently
% matched its install receipt. If any target is invalid, continue checking
% it individually on later calls so repaired files can recover immediately.
recent.watch = [];
if numel(dirs) == numel(native.targets) && ~isempty(firstStamp) && ...
        isequal(firstStamp,content_stamp(watch))
    recent.watch = watch;
    recent.stamp = firstStamp;
    recent.dirs = dirs;
    recent.allowed = allowed;
end
end

function [watch,stamp] = snapshot(root,targets)
watch = [];
stamp = [];
try
    sources = cellfun(@(f) fullfile(root,f),{targets.source},'UniformOutput',false);
    binaries = cellfun(@(f) fullfile(root,f),{targets.binary},'UniformOutput',false);
    parents = cellfun(@fileparts,sources,'UniformOutput',false);
    [folders,first] = unique(parents,'stable');
    files = [sources,binaries];
    for k = 1:numel(first)
        inputs = nf_native_inputs(sources{first(k)});
        files = [files,inputs(2:end)]; %#ok<AGROW>
    end
    toFile = @(p) javaObject('java.io.File',p);
    watch.files = cellfun(toFile,files,'UniformOutput',false);
    watch.paths = cellfun(@(f) f.toPath(),watch.files,'UniformOutput',false);
    watch.folders = cellfun(toFile,folders,'UniformOutput',false);
    watch.binFolders = cellfun(@(b) toFile(fileparts(b)),binaries,'UniformOutput',false);
    stamp = content_stamp(watch);
catch
    % Unavailable inputs are diagnosed independently by the caller.
end
end

function stamp = content_stamp(watch)
% Collect the names of shared headers and hash every tracked file. Reuse
% the Java file objects and one digest to avoid repeated directory walks,
% canonical-path construction and separate digest-to-hex conversions.
try
    md = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
    headerNames = cell(size(watch.folders));
    for k = 1:numel(watch.folders)
        entries = cellstr(watch.folders{k}.list());
        matches = ~cellfun('isempty',regexpi(entries,'\.(h|hpp|cuh)$','once'));
        headerNames{k} = sort(entries(matches));
    end
    locations = cell(size(watch.files));
    for k = 1:numel(watch.paths)
        locations{k} = char(watch.files{k}.getCanonicalPath());
        bytes = javaMethod('readAllBytes','java.nio.file.Files',watch.paths{k});
        md.update(typecast(int64(numel(bytes)),'int8'));
        if ~isempty(bytes), md.update(bytes); end
    end
    bins = cell(size(watch.binFolders));
    for k = 1:numel(bins)
        bins{k} = sort(cellstr(watch.binFolders{k}.list()));
    end
    stamp = {headerNames,locations,bins,md.digest()};
catch
    stamp = []; % Missing, unreadable or concurrently replaced.
end
end
