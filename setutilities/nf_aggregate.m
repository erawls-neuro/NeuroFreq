function tfagg = nf_aggregate(varargin)
%NF_AGGREGATE  Aggregate NeuroFreq TF sets across .mat files.
%
% Supported calling patterns:
%
%   1) GUI selection mode (no inputs):
%       tfagg = nf_aggregate()
%
%   2) DIR-struct mode (one required struct input):
%       files = dir('mydir/*subject*_tf.mat');
%       tfagg  = nf_aggregate(files)
%       tfagg  = nf_aggregate(files, 'aggtype', 'singletrial')
%
% The first input (if present) MUST be the struct output of dir().
%
% Aggregation modes:
%   - 'avg'        : stacks condition-averaged TF fields across subjects
%   - 'singletrial': concatenates single-trial TF fields across subjects
%   - 'stat'       : stacks statistical/regression outputs across subjects
%
% Auto-detection:
%   1) If tf.scale == 'stat' -> 'stat'
%   2) Else, if isfield(tf,'conds') -> 'avg'  (NF: .conds exists => averaged)
%   3) Else -> 'singletrial'                  (NF: no .conds => singletrial, even if ntrial==1)
%
% NF shape rules:
%   - singletrial TF fields: chan x freq x time x trial   (4D only)
%   - avg TF fields:         chan x freq x time x cond    (4D only)
%   - stat measure fields:   top-level structs (e.g., power, power_induced_ap, phase_itc)
%                            containing numeric metrics (e.g., corrcoef) and metadata
%   - stat numeric fields:   any other numeric arrays (stacked across subjects)
%
% Output:
%   tfagg has common header fields: freqs, times, chanlocs, conds, aggtype,
%   plus stacked data in top-level measure fields (power*/phase*) and/or tfagg.stat.
%
% Notes:
%   - Field names in outputs are sanitized with matlab.lang.makeValidName.
%     Raw names are stored in *Fields_raw.
%   - This function is intentionally strict about layout consistency across files.
%   - Canonical output "scale" is set to the aggregation type (avg/singletrial/stat).
%     The original source scale (if present) is stored in tfagg.source_scale.
%
% IMPORTANT NF CONVENTION (enforced here):
%   - Power and phase measures are TOP-LEVEL siblings (no tf.power.(...) containers).
%   - This function will ignore and NOT create any nested group containers.

filesIn = [];

if nargin >= 1
    filesIn = varargin{1};
end

if nargin >= 2
    opts = local_parse_args(varargin{2:end});
else
    opts = local_parse_args();
end

if isempty(filesIn) == 0
    if isstruct(filesIn) == 0
        error('nf_aggregate:BadInputs', 'If provided, first input must be a dir() struct.');
    end
    [files, source_dir] = local_files_from_dirstruct(filesIn);
else
    [files, source_dir] = local_select_matfiles_gui();
end

nFiles = numel(files);

if nFiles == 0
    error('nf_aggregate:NoFiles', 'No usable .mat files were provided/selected.');
end

if opts.verbose == 1
    fprintf('[nf_aggregate]: aggregating %d file(s).\n', nFiles);
end

nm = cell(nFiles, 1);
fullPaths = cell(nFiles, 1);

for iFile = 1:nFiles
    nm{iFile} = files(iFile).name;
    fullPaths{iFile} = local_direntry_to_fullpath(files(iFile), source_dir);

    if exist(fullPaths{iFile}, 'file') ~= 2
        error('nf_aggregate:BadPath', 'File does not exist: %s', fullPaths{iFile});
    end
end

source_dir = local_common_dir(fullPaths);

% Load all TF sets
tfs = cell(nFiles, 1);

for iFile = 1:nFiles
    tf = local_load_tfset(fullPaths{iFile});
    tf = local_validate_tfset_minimum(tf, fullPaths{iFile}, opts);
    tfs{iFile} = tf;
end

% Determine aggregation type
if isempty(opts.aggtype) == 0
    aggType = lower(opts.aggtype);
else
    aggType = local_detect_aggtype(tfs{1});
end

if strcmpi(aggType, 'avg') == 0 && strcmpi(aggType, 'singletrial') == 0 && strcmpi(aggType, 'stat') == 0
    error('nf_aggregate:BadAggType', 'Unknown aggtype: %s', aggType);
end

if opts.verbose == 1
    fprintf('[nf_aggregate]: aggtype = %s\n', aggType);
end

% Validate shared header consistency
ref = tfs{1};

checkConds = false;

if strcmpi(aggType, 'avg') == 1
    checkConds = true;
end

valid = NaN(1,nFiles);
valid(1) = 1;
for iFile = 2:nFiles
    valid(iFile) = local_validate_layout_match(ref, tfs{iFile}, fullPaths{1}, fullPaths{iFile}, checkConds, opts.strict);
end

% remove bad sets
tfs = tfs(find(valid)); %#ok
fullPaths = fullPaths(find(valid)); %#ok
nFiles = sum(valid);

% Build payloads per file (includes NF rule enforcement per mode)
payloads = cell(nFiles, 1);

for iFile = 1:nFiles

    tf = tfs{iFile};

    local_validate_tf_rules_per_mode(tf, aggType, fullPaths{iFile}, opts);

    if strcmpi(aggType, 'avg') == 1
        payloads{iFile} = local_prepare_avg_payload(tf, fullPaths{iFile}, opts);
    elseif strcmpi(aggType, 'singletrial') == 1
        payloads{iFile} = local_prepare_singletrial_payload(tf, fullPaths{iFile}, opts);
    elseif strcmpi(aggType, 'stat') == 1
        payloads{iFile} = local_prepare_stat_payload(tf, fullPaths{iFile}, opts);
    else
        error('nf_aggregate:InternalError', 'Unexpected aggType branch.');
    end
end

% Guarantee subid_value exists
for iFile = 1:nFiles
    if isfield(payloads{iFile}, 'subid_value') == 0
        payloads{iFile}.subid_value = '';
    end

    if isempty(payloads{iFile}.subid_value) == 1
        payloads{iFile}.subid_value = local_make_subid_value(payloads{iFile}.filepath, iFile);
        payloads{iFile}.subid_value_is_dummy = true;
    end
end

% Build aggregate
if strcmpi(aggType, 'avg') || strcmpi(aggType, 'singletrial')
    tfagg = local_build_singletrial_aggregate(payloads, opts);
elseif strcmpi(aggType, 'stat') == 1
    tfagg = local_build_stat_aggregate(payloads, opts);
else
    error('nf_aggregate:InternalError', 'Unknown aggType in builder.');
end

% Attach top-level metadata
tfagg.aggtype = aggType;
tfagg.nsubjects = nFiles;

subid = cell(nFiles, 1);

for iFile = 1:nFiles
    subid{iFile} = payloads{iFile}.subid_value;
end

tfagg.subid = subid;
tfagg.source_dir = source_dir;
tfagg.source_paths = fullPaths;

tfagg.freqs = ref.freqs;
tfagg.times = ref.times;
tfagg.chanlocs = ref.chanlocs;

if isfield(ref, 'conds') == 1
    tfagg.conds = ref.conds;
else
    tfagg.conds = {};
end

tfagg.source_scale = '';

if isfield(ref, 'scale') == 1
    if ischar(ref.scale) == 1
        tfagg.source_scale = ref.scale;
    end
end

tfagg.scale = aggType;

if opts.verbose == 1
    fprintf('[nf_aggregate]: aggregation complete.\n');
end

end

% Local: parse args
function opts = local_parse_args(varargin)

opts = struct();
opts.verbose = false;
opts.strict = true;
opts.aggtype = '';

if mod(numel(varargin), 2) ~= 0
    error('nf_aggregate:BadArgs', 'Arguments must be name/value pairs.');
end

for i = 1:2:numel(varargin)
    key = varargin{i};
    val = varargin{i + 1};

    if isstring(key) == 1
        key = char(key);
    end

    if ischar(key) == 0
        error('nf_aggregate:BadArgs', 'Parameter names must be strings.');
    end

    key = lower(key);

    if strcmp(key, 'verbose') == 1
        opts.verbose = logical(val);
    elseif strcmp(key, 'strict') == 1
        opts.strict = logical(val);
    elseif strcmp(key, 'aggtype') == 1
        if isstring(val) == 1
            val = char(val);
        end
        opts.aggtype = val;
    else
        error('nf_aggregate:BadArgs', 'Unknown parameter: %s', key);
    end
end

if isempty(opts.aggtype) == 0
    if ischar(opts.aggtype) == 0
        error('nf_aggregate:BadArgs', 'aggtype must be a string.');
    end
end

end

% Local: resolve files from dir() struct
function [files, source_dir] = local_files_from_dirstruct(d)

if isempty(d) == 1
    error('nf_aggregate:BadInputs', 'dir() struct input is empty.');
end

if isfield(d, 'name') == 0
    error('nf_aggregate:BadInputs', 'dir() struct input must contain a .name field.');
end

files = d;
files = local_filter_dirstruct(files);

if isempty(files) == 1
    error('nf_aggregate:NoFiles', 'dir() struct contained no usable .mat files.');
end

if isfield(files, 'folder') == 1
    if isempty(files(1).folder) == 0
        source_dir = files(1).folder;
    else
        source_dir = pwd;
    end
else
    source_dir = pwd;
end

end

% Local: GUI selection of .mat files
function [files, source_dir] = local_select_matfiles_gui()

[sel, pth] = uigetfile('*.mat', 'Select NeuroFreq TF file(s)', 'MultiSelect', 'on');

if isequal(sel, 0)
    error('nf_aggregate:NoFiles', 'File selection cancelled.');
end

if ischar(sel) == 1
    sel = {sel};
end

if iscell(sel) == 0
    error('nf_aggregate:NoFiles', 'Unexpected selection output from uigetfile.');
end

n = numel(sel);

files = repmat(struct('name', '', 'folder', '', 'isdir', false), n, 1);

for i = 1:n
    files(i).name = sel{i};
    files(i).folder = pth;
    files(i).isdir = false;
end

files = local_filter_dirstruct(files);
source_dir = pth;

if isempty(files) == 1
    error('nf_aggregate:NoFiles', 'No usable .mat files selected.');
end

end

% Local: filter dir-struct entries to usable .mat files
function files = local_filter_dirstruct(files)

if isempty(files) == 1
    return
end

keep = false(numel(files), 1);

for i = 1:numel(files)

    if isfield(files(i), 'isdir') == 1
        if files(i).isdir == 1
            keep(i) = false;
            continue
        end
    end

    nm = files(i).name;

    if isempty(nm) == 1
        keep(i) = false;
        continue
    end

    if strcmp(nm, '.') == 1 || strcmp(nm, '..') == 1
        keep(i) = false;
        continue
    end

    [~, ~, ext] = fileparts(nm);

    if strcmpi(ext, '.mat') == 1
        keep(i) = true;
    else
        keep(i) = false;
    end
end

files = files(keep);

end

% Local: build a full path from a dir entry
function fullp = local_direntry_to_fullpath(d, fallback_dir)

nm = d.name;

hasSep = false;

if contains(nm, '/') == 1
    hasSep = true;
end

if contains(nm, '\') == 1
    hasSep = true;
end

if hasSep == true
    fullp = nm;
    return
end

if isfield(d, 'folder') == 1
    if isempty(d.folder) == 0
        fullp = fullfile(d.folder, nm);
        return
    end
end

if isempty(fallback_dir) == 0
    fullp = fullfile(fallback_dir, nm);
    return
end

fullp = fullfile(pwd, nm);

end

% Local: common directory across all full paths
function common = local_common_dir(fullPaths)

common = '';

if isempty(fullPaths) == 1
    return
end

dirs = cell(numel(fullPaths), 1);

for i = 1:numel(fullPaths)
    dirs{i} = fileparts(fullPaths{i});
end

first = dirs{1};
allSame = true;

for i = 2:numel(dirs)
    if strcmp(first, dirs{i}) == 0
        allSame = false;
    end
end

if allSame == 1
    common = first;
else
    common = '';
end

end

% Local: load tf set
function tf = local_load_tfset(filepath)

S = load(filepath);

fns = fieldnames(S);

if isempty(fns) == 1
    error('nf_aggregate:LoadError', 'No variables in file: %s', filepath);
end

if numel(fns) == 1 %#ok
    tf = S.(fns{1});
    return
end

if isfield(S, 'tf') == 1
    tf = S.tf;
    return
end

if isfield(S, 'tfset') == 1
    tf = S.tfset;
    return
end

error('nf_aggregate:LoadError', 'Could not identify TF struct in file: %s', filepath);

end

% Local: minimal validation (do NOT auto-insert .conds)
function tf = local_validate_tfset_minimum(tf, filepath, opts)

if isstruct(tf) == 0
    error('nf_aggregate:BadTF', 'TF set is not a struct. File: %s', filepath);
end

req = {'freqs', 'times', 'chanlocs'};

for i = 1:numel(req)
    if isfield(tf, req{i}) == 0
        error('nf_aggregate:BadTF', 'Missing field %s. File: %s', req{i}, filepath);
    end
end

if isempty(tf.freqs) == 1 || isnumeric(tf.freqs) == 0
    error('nf_aggregate:BadTF', 'freqs must be a numeric vector. File: %s', filepath);
end

if isempty(tf.times) == 1 || isnumeric(tf.times) == 0
    error('nf_aggregate:BadTF', 'times must be a numeric vector. File: %s', filepath);
end

tf.freqs = tf.freqs(:)';
tf.times = tf.times(:)';

if isfield(tf, 'conds') == 1

    if isstring(tf.conds) == 1
        tf.conds = cellstr(tf.conds(:));
    elseif iscategorical(tf.conds) == 1
        tf.conds = cellstr(string(tf.conds(:)));
    elseif iscell(tf.conds) == 1

        for k = 1:numel(tf.conds)

            if isstring(tf.conds{k}) == 1
                tf.conds{k} = char(tf.conds{k});
            elseif iscategorical(tf.conds{k}) == 1
                tf.conds{k} = char(string(tf.conds{k}));
            end

        end

    elseif isnumeric(tf.conds) == 1
        % leave numeric (stat payloads may use numeric codes)
    else
        error('nf_aggregate:BadTF', 'conds has unsupported type. File: %s', filepath);
    end

end

if isfield(tf, 'ntrial') == 1
    if isempty(tf.ntrial) == 1
        if opts.strict == 1
            error('nf_aggregate:BadTF', 'ntrial is present but empty. File: %s', filepath);
        end
    else
        if ~isnumeric(tf.ntrial)
            error('nf_aggregate:BadTF', 'ntrial must be numeric. File: %s', filepath);
        end
        if ~isscalar(tf.ntrial)
            error('nf_aggregate:BadTF', 'ntrial must be a scalar. File: %s', filepath);
        end
        if tf.ntrial < 1
            error('nf_aggregate:BadTF', 'ntrial must be >= 1. File: %s', filepath);
        end
    end
end

end

% Local: detect agg type
function aggType = local_detect_aggtype(tf)

aggType = '';

if isfield(tf, 'scale') == 1
    if ischar(tf.scale) == 1
        if strcmpi(tf.scale, 'stat') == 1
            aggType = 'stat';
        end
    end
end

if isempty(aggType) == 1
    if isfield(tf, 'conds') == 1
        aggType = 'avg';
    else
        aggType = 'singletrial';
    end
end

end

% Local: enforce NF structural rules per mode (STRICT)
function local_validate_tf_rules_per_mode(tf, aggType, filepath, opts)

if opts.strict == 0
    return
end

if strcmpi(aggType, 'avg') == 1

    if isfield(tf, 'conds') == 0
        error('nf_aggregate:NFRuleViolation', ...
            'AVG mode requires .conds to exist (NF rule). File: %s', filepath);
    end
    %
    % if iscell(tf.conds) == 0
    %     error('nf_aggregate:NFRuleViolation', ...
    %         'AVG mode requires .conds to be a cell array. File: %s', filepath);
    % end

    if isempty(tf.conds) == 1
        error('nf_aggregate:NFRuleViolation', ...
            'AVG mode requires non-empty .conds. File: %s', filepath);
    end

elseif strcmpi(aggType, 'singletrial') == 1

    if isfield(tf, 'conds') == 1
        error('nf_aggregate:NFRuleViolation', ...
            'SINGLETRIAL mode forbids .conds (NF rule). File: %s', filepath);
    end

    if isfield(tf, 'ntrial') == 0
        error('nf_aggregate:NFRuleViolation', ...
            'SINGLETRIAL mode requires .ntrial to exist (NF rule). File: %s', filepath);
    end

    if isempty(tf.ntrial) == 1
        error('nf_aggregate:NFRuleViolation', ...
            'SINGLETRIAL mode requires non-empty .ntrial. File: %s', filepath);
    end

elseif strcmpi(aggType, 'stat') == 1

    % no extra structural rule here

else
    error('nf_aggregate:InternalError', 'Unknown aggType in rule enforcement.');
end

end

% Local: validate layout match
function valid = local_validate_layout_match(tfRef, tfCur, pathRef, pathCur, checkConds, strict)

tol = 1e-10;

if numel(tfRef.freqs) ~= numel(tfCur.freqs)
    error('nf_aggregate:LayoutMismatch', 'freq length mismatch: %s vs %s', pathRef, pathCur);
end

if max(abs(tfRef.freqs(:) - tfCur.freqs(:))) > tol
    error('nf_aggregate:LayoutMismatch', 'freq values mismatch: %s vs %s', pathRef, pathCur);
end

if numel(tfRef.times) ~= numel(tfCur.times)
    error('nf_aggregate:LayoutMismatch', 'time length mismatch: %s vs %s', pathRef, pathCur);
end

if max(abs(tfRef.times(:) - tfCur.times(:))) > tol
    error('nf_aggregate:LayoutMismatch', 'time values mismatch: %s vs %s', pathRef, pathCur);
end

labelsRef = local_chan_labels(tfRef.chanlocs);
labelsCur = local_chan_labels(tfCur.chanlocs);

if numel(labelsRef) ~= numel(labelsCur)
    error('nf_aggregate:LayoutMismatch', 'chan count mismatch: %s vs %s', pathRef, pathCur);
end

if strict == 1
    for i = 1:numel(labelsRef)
        if strcmpi(labelsRef{i}, labelsCur{i}) == 0
            error('nf_aggregate:LayoutMismatch', 'chan labels mismatch: %s vs %s', pathRef, pathCur);
        end
    end
end

if checkConds == 1
    if isfield(tfRef, 'conds') == 0 || isfield(tfCur, 'conds') == 0
        error('nf_aggregate:LayoutMismatch', 'conds missing for avg layout check: %s vs %s', pathRef, pathCur);
    end

    % if tfRef.conds ~= tfCur.conds
    %     warning('nf_aggregate:LayoutMismatch', 'conds length mismatch: %s vs %s', pathRef, pathCur);
    %     warning('continuing WITHOUT this file');
    %     valid = 0;
    %     return
    % else
    valid = 1;
    % end

    % for i = 1:numel(tfRef.conds)
    %     if strcmpi(tfRef.conds{i}, tfCur.conds{i}) == 0
    %         error('nf_aggregate:LayoutMismatch', 'conds mismatch: %s vs %s', pathRef, pathCur);
    %     end
    % end
end

end

% Local: extract channel labels robustly
function labels = local_chan_labels(chanlocs)

nChan = numel(chanlocs);
labels = cell(nChan, 1);

for i = 1:nChan
    if isstruct(chanlocs(i)) == 1
        if isfield(chanlocs(i), 'labels') == 1
            labels{i} = chanlocs(i).labels;
        elseif isfield(chanlocs(i), 'label') == 1
            labels{i} = chanlocs(i).label;
        else
            labels{i} = sprintf('Ch%d', i);
        end
    else
        labels{i} = sprintf('Ch%d', i);
    end
end

end

% Local: prepare avg payload (TOP-LEVEL power*/phase* only)
function payload = local_prepare_avg_payload(tf, filepath, opts)

payload = struct();
payload.filepath = filepath;

payload.freqs = tf.freqs;
payload.times = tf.times;
payload.chanlocs = tf.chanlocs;
payload.conds = tf.conds;

payload.subid_value = local_get_subid_value(tf, filepath, 1);

nChan = numel(tf.chanlocs);
nFreq = numel(tf.freqs);
nTime = numel(tf.times);
nConds = tf.conds;

[powAll, powRaw] = local_collect_top_level_tf_numeric(tf, {'power', 'pow_'}, filepath, opts);
[phAll, phRaw] = local_collect_top_level_tf_numeric(tf, {'phase', 'ph_'}, filepath, opts);

if opts.strict == 1
    if isempty(fieldnames(powAll)) == 1 && isempty(fieldnames(phAll)) == 1
        error('nf_aggregate:NoDataFields', 'No 4D top-level power/phase fields found for avg aggregation. File: %s', filepath);
    end
end

payload.powerFields_raw = powRaw;
payload.phaseFields_raw = phRaw;

pf = fieldnames(powAll);

for i = 1:numel(pf)
    fn = pf{i};
    X = powAll.(fn);
    X = local_standardize_avg_field(X, nConds, nChan, nFreq, nTime, filepath, fn);
    payload.(fn) = X;
end

phf = fieldnames(phAll);

for i = 1:numel(phf)
    fn = phf{i};
    X = phAll.(fn);
    X = local_standardize_avg_field(X, nConds, nChan, nFreq, nTime, filepath, fn);
    payload.(fn) = X;
end

payload.powerFields = fieldnames(powAll);
payload.phaseFields = fieldnames(phAll);

if opts.strict == 1
    overlap = intersect(payload.powerFields, payload.phaseFields);

    if isempty(overlap) == 0
        error('nf_aggregate:FieldNameCollision', ...
            'Power/phase field-name collision after makeValidName in %s: %s', filepath, overlap{1});
    end
end

if isfield(tf,'behavior')
    payload.behavior = tf.behavior;
end

end

% Local: prepare singletrial payload (TOP-LEVEL power*/phase* only)
function payload = local_prepare_singletrial_payload(tf, filepath, opts)

payload = struct();
payload.filepath = filepath;

payload.freqs = tf.freqs;
payload.times = tf.times;
payload.chanlocs = tf.chanlocs;

payload.subid_value = local_get_subid_value(tf, filepath, 1);

nChan = numel(tf.chanlocs);
nFreq = numel(tf.freqs);
nTime = numel(tf.times);

[powAll, powRaw] = local_collect_top_level_tf_numeric(tf, {'power', 'pow_'}, filepath, opts);
[phAll, phRaw] = local_collect_top_level_tf_numeric(tf, {'phase', 'ph_'}, filepath, opts);

if opts.strict == 1
    % if isempty(fieldnames(powAll)) == 1 && isempty(fieldnames(phAll)) == 1
    %     error('nf_aggregate:NoDataFields', 'No 4D top-level power/phase fields found for singletrial aggregation. File: %s', filepath);
    % end
end

payload.powerFields_raw = powRaw;
payload.phaseFields_raw = phRaw;

pf = fieldnames(powAll);

for i = 1:numel(pf)
    fn = pf{i};
    X = powAll.(fn);
    X = local_standardize_singletrial_field(X, tf.ntrial, nChan, nFreq, nTime, filepath, fn);
    payload.(fn) = X;
end

phf = fieldnames(phAll);

for i = 1:numel(phf)
    fn = phf{i};
    X = phAll.(fn);
    X = local_standardize_singletrial_field(X, tf.ntrial, nChan, nFreq, nTime, filepath, fn);
    payload.(fn) = X;
end

payload.powerFields = fieldnames(powAll);
payload.phaseFields = fieldnames(phAll);

if opts.strict == 1
    overlap = intersect(payload.powerFields, payload.phaseFields);

    if isempty(overlap) == 0
        error('nf_aggregate:FieldNameCollision', ...
            'Power/phase field-name collision after makeValidName in %s: %s', filepath, overlap{1});
    end
end

payload.ntrials = tf.ntrial;

payload.trial_condindex = [];
payload.trial_condlabel = {};

if isfield(tf, 'trial_condindex') == 1
    if numel(tf.trial_condindex) == payload.ntrials
        payload.trial_condindex = tf.trial_condindex(:);
    end
end

if isfield(tf, 'trial_condlabel') == 1
    if numel(tf.trial_condlabel) == payload.ntrials
        payload.trial_condlabel = tf.trial_condlabel(:);
    end
end

end

% Local: prepare stat payload (TOP-LEVEL measure structs + numeric stat arrays)
function payload = local_prepare_stat_payload(tf, filepath, opts)

payload = struct();
payload.filepath = filepath;

payload.freqs = tf.freqs;
payload.times = tf.times;
payload.chanlocs = tf.chanlocs;

payload.subid_value = local_get_subid_value(tf, filepath, 1);

payload.conds = {};

if isfield(tf, 'conds') == 1
    payload.conds = tf.conds;
end

% (A) numeric "stat" fields (misc outputs) stacked into payload.stat
payload.stat = struct();

fns = fieldnames(tf);

skip = {'freqs', 'times', 'chanlocs', 'conds', 'scale', 'trial_condindex', 'trial_condlabel', 'ntrial', ...
    'ntrls_in', 'ntrls_used', 'regressors', 'df', 'Fs', 'permutations', 'tfce_par', ...
    'prednames', 'model', 'powerFields', 'phaseFields', 'statFields', ...
    'powerFields_raw', 'phaseFields_raw', 'statFields_raw'};

rawKept = {};

for i = 1:numel(fns)
    fn = fns{i};

    if any(strcmp(fn, skip))
        continue
    end

    v = tf.(fn);

    if isnumeric(v) == 0
        continue
    end

    if isempty(v) == 1
        continue
    end

    key = matlab.lang.makeValidName(fn);

    if isfield(payload.stat, key) == 1
        error('nf_aggregate:NameCollision', ...
            'makeValidName collision (stat numeric) in %s: "%s"', filepath, key);
    end

    payload.stat.(key) = v;
    rawKept{end + 1, 1} = fn; %#ok<AGROW>
end

payload.statFields_raw = rawKept;
payload.statFields = fieldnames(payload.stat);

% (B) top-level measure structs (power*/phase*) stacked separately
measures = local_collect_top_level_measure_structs(tf);
payload.measures = measures;
payload.measureFields = fieldnames(measures);

payload.measureFields_raw = payload.measureFields;

if isfield(tf, 'prednames') == 1
    payload.prednames = tf.prednames;
else
    payload.prednames = {};
end

if isfield(tf, 'model') == 1
    payload.model = tf.model;
else
    payload.model = '';
end

if opts.strict == 1
    if isempty(payload.statFields) == 1 && isempty(payload.measureFields) == 1
        error('nf_aggregate:NoStatFields', 'No numeric stat fields or measure structs found in %s', filepath);
    end
end

end

% Local: standardize avg fields to chan x freq x time x cond (4D only)
function X4 = local_standardize_avg_field(X, nConds, nChan, nFreq, nTime, filename, fname)

if ndims(X) ~= 4
    error('nf_aggregate:BadAvgShape', ...
        'Avg field %s must be 4D (chan x freq x time x cond). File: %s', fname, filename);
end

if size(X, 1) ~= nChan
    error('nf_aggregate:BadAvgShape', 'Field %s channel dim mismatch. File: %s', fname, filename);
end

if size(X, 2) ~= nFreq
    error('nf_aggregate:BadAvgShape', 'Field %s freq dim mismatch. File: %s', fname, filename);
end

if size(X, 3) ~= nTime
    error('nf_aggregate:BadAvgShape', 'Field %s time dim mismatch. File: %s', fname, filename);
end

if size(X, 4) ~= nConds
    error('nf_aggregate:BadAvgShape', 'Field %s cond dim mismatch. File: %s', fname, filename);
end

X4 = X;

end

% Local: standardize singletrial fields to chan x freq x time x trial (4D only)
function X4 = local_standardize_singletrial_field(X, nTrial, nChan, nFreq, nTime, filename, fname)

if ndims(X) ~= 4
    error('nf_aggregate:BadSingleTrialShape', ...
        'Single-trial field %s must be 4D (chan x freq x time x trial). File: %s', fname, filename);
end

if size(X, 1) ~= nChan
    error('nf_aggregate:BadSingleTrialShape', 'Field %s channel dim mismatch. File: %s', fname, filename);
end

if size(X, 2) ~= nFreq
    error('nf_aggregate:BadSingleTrialShape', 'Field %s freq dim mismatch. File: %s', fname, filename);
end

if size(X, 3) ~= nTime
    error('nf_aggregate:BadSingleTrialShape', 'Field %s time dim mismatch. File: %s', fname, filename);
end

if size(X, 4) ~= nTrial
    error('nf_aggregate:BadSingleTrialShape', 'Field %s trial dim mismatch. File: %s', fname, filename);
end

X4 = X;

end


% Local: build singletrial aggregate (concat trials across subjects)
function tfagg = local_build_singletrial_aggregate(payloads, opts)

tfagg = struct();

p0 = payloads{1};

pf = p0.powerFields;
phf = p0.phaseFields;

nSub = numel(payloads);

tfagg.powerFields_raw = p0.powerFields_raw;
tfagg.phaseFields_raw = p0.phaseFields_raw;

tfagg.powerFields = pf;
tfagg.phaseFields = phf;

if opts.strict == 1
    for s = 2:nSub
        if ~isequal(pf, payloads{s}.powerFields)
            error('nf_aggregate:FieldMismatch', ...
                'Power field set mismatch between %s and %s', payloads{1}.filepath, payloads{s}.filepath);
        end

        if ~isequal(phf, payloads{s}.phaseFields)
            error('nf_aggregate:FieldMismatch', ...
                'Phase field set mismatch between %s and %s', payloads{1}.filepath, payloads{s}.filepath);
        end
    end
end

nChan = numel(p0.chanlocs);
nFreq = numel(p0.freqs);
nTime = numel(p0.times);

subid = cell(nSub, 1);

for s = 1:nSub
    subid{s} = payloads{s}.subid_value;
end

tfagg.subid = subid;

nTrialsAll = 0;

for s = 1:nSub
    if isfield(payloads{s},'conds')
        nTrialsAll = nTrialsAll + payloads{s}.conds;
    else
        nTrialsAll = nTrialsAll + payloads{s}.ntrials;
    end
end

for i = 1:numel(pf)
    fn = pf{i};
    tfagg.(fn) = nan(nChan, nFreq, nTime, nTrialsAll);
end

for i = 1:numel(phf)
    fn = phf{i};
    tfagg.(fn) = nan(nChan, nFreq, nTime, nTrialsAll);
end

tfagg.trial_subjectindex = nan(nTrialsAll, 1);
tfagg.trial_subid = cell(nTrialsAll, 1);

tfagg.trial_fileindex = nan(nTrialsAll, 1);

cursor = 0;

for s = 1:nSub
    p = payloads{s};
    if isfield(payloads{s},'conds')
        nT = p.conds;
    else
        nT = p.ntrials;
    end

    idx = (cursor + 1):(cursor + nT);

    for i = 1:numel(pf)
        fn = pf{i};
        tfagg.(fn)(:, :, :, idx) = p.(fn);
    end

    for i = 1:numel(phf)
        fn = phf{i};
        tfagg.(fn)(:, :, :, idx) = p.(fn);
    end

    if (s==1) && isfield(p,'behavior')
        tfagg.behavior = p.behavior;
    elseif s>1 && isfield(p,'behavior')
        tfagg.behavior = [tfagg.behavior p.behavior];
    end

    tfagg.trial_subjectindex(idx) = s;

    for k = 1:nT
        tfagg.trial_subid{idx(k)} = subid{s};
    end

    tfagg.trial_fileindex(idx) = s;

    tfagg.trial_condindex(idx) = 1;

    cursor = cursor + nT;
end

if isfield(payloads{s},'conds')
    tfagg.conds = nTrialsAll;
else
    tfagg.ntrials = nTrialsAll;
end

end

% Local: build stat aggregate
% - tfagg.stat.(numericStatField) stacked as subjects-first
% - tfagg.(measureName).(metricName) stacked as subjects-first for numeric metrics
function tfagg = local_build_stat_aggregate(payloads, opts)

tfagg = struct();

p0 = payloads{1};

nSub = numel(payloads);

% ---- (A) numeric stat fields -> tfagg.stat ----
tfagg.stat = struct();

sf = p0.statFields;

tfagg.statFields_raw = p0.statFields_raw;
tfagg.statFields = sf;

if isfield(p0, 'prednames') == 1
    tfagg.prednames = p0.prednames;
else
    tfagg.prednames = {};
end

if isfield(p0, 'model') == 1
    tfagg.model = p0.model;
else
    tfagg.model = '';
end

if opts.strict == 1
    for s = 2:nSub

        if ~isequal(sf, payloads{s}.statFields)
            error('nf_aggregate:StatFieldMismatch', ...
                'Stat numeric fields mismatch between %s and %s', payloads{1}.filepath, payloads{s}.filepath);
        end

        if isfield(p0, 'prednames') == 1 && isfield(payloads{s}, 'prednames') == 1
            if ~isequal(p0.prednames, payloads{s}.prednames)
                error('nf_aggregate:PredMismatch', ...
                    'Predictor names mismatch between %s and %s', payloads{1}.filepath, payloads{s}.filepath);
            end
        end

    end
end

for i = 1:numel(sf)

    fn = sf{i};

    X0 = double(p0.stat.(fn));

    baseSz = size(X0);

    X = nan([baseSz nSub]);

    for s = 1:nSub

        p = payloads{s};

        if isfield(p.stat, fn) == 0
            error('nf_aggregate:MissingField', 'Missing stat numeric field %s in %s', fn, p.filepath);
        end

        Xs = double(p.stat.(fn));

        if ~isequal(size(Xs), size(X0))
            error('nf_aggregate:StatShapeMismatch', ...
                'Stat numeric field %s size mismatch between %s and %s', fn, p0.filepath, p.filepath);
        end

        idx = repmat({':'}, 1, ndims(X0));

        X(idx{:}, s) = Xs;

    end

    X = local_move_subject_last_to_first(X);

    tfagg.stat.(fn) = X;

end

% ---- (B) top-level measure structs -> tfagg.(measure).(metric) ----
mFields = p0.measureFields;

tfagg.measureFields_raw = p0.measureFields_raw;
tfagg.measureFields = mFields;

if opts.strict == 1
    for s = 2:nSub
        if ~isequal(mFields, payloads{s}.measureFields)
            error('nf_aggregate:MeasureFieldMismatch', ...
                'Measure struct field set mismatch between %s and %s', payloads{1}.filepath, payloads{s}.filepath);
        end
    end
end

for iM = 1:numel(mFields)

    mName = mFields{iM};

    mAgg = local_stack_measure_struct_metrics(payloads, mName, opts);

    tfagg.(mName) = mAgg;

end

end

function Y = local_move_subject_last_to_first(X)

nd = ndims(X);

if nd < 2
    Y = X;
    return
end

perm = 1:nd;
perm = [nd perm(1:(nd - 1))];

Y = permute(X, perm);

end

function mOut = local_stack_measure_struct_metrics(payloads, mName, opts)

nSub = numel(payloads);

m0 = payloads{1}.measures.(mName);

metricNames = fieldnames(m0);

mOut = struct();

for iF = 1:numel(metricNames)

    fn = metricNames{iF};

    v0 = m0.(fn);

    if isnumeric(v0) == 0 || isempty(v0) == 1
        mOut.(fn) = v0;

        if opts.strict == 1
            for s = 2:nSub
                mS = payloads{s}.measures.(mName);

                if isfield(mS, fn) == 0
                    error('nf_aggregate:MissingMeasureMetric', ...
                        'Missing metric %s.%s in %s', mName, fn, payloads{s}.filepath);
                end

                vS = mS.(fn);

                if ~isequal(vS, v0)
                    error('nf_aggregate:MeasureMetaMismatch', ...
                        '%s.%s mismatch between %s and %s', mName, fn, payloads{1}.filepath, payloads{s}.filepath);
                end
            end
        end

        continue
    end

    X0 = double(v0);

    baseSz = size(X0);

    X = nan([baseSz nSub]);

    for s = 1:nSub

        mS = payloads{s}.measures.(mName);

        if isfield(mS, fn) == 0
            error('nf_aggregate:MissingMeasureMetric', ...
                'Missing metric %s.%s in %s', mName, fn, payloads{s}.filepath);
        end

        Xs = double(mS.(fn));

        if ~isequal(size(Xs), baseSz)
            error('nf_aggregate:MeasureMetricShapeMismatch', ...
                '%s.%s size mismatch between %s and %s', mName, fn, payloads{1}.filepath, payloads{s}.filepath);
        end

        idx = repmat({':'}, 1, ndims(X0));

        X(idx{:}, s) = Xs;

    end

    X = local_move_subject_last_to_first(X);

    mOut.(fn) = X;

end

end

function subid = local_get_subid_value(tf, filepath, fallbackIndex)

subid = '';

candidates = {'subid', 'subID', 'subject', 'subjectID', 'participant', 'participantID', 'id'};

for iC = 1:numel(candidates)
    fn = candidates{iC};

    if isfield(tf, fn) == 1
        v = tf.(fn);

        if ischar(v) == 1
            subid = v;
            break
        end

        if isstring(v) == 1
            subid = char(v);
            break
        end

        if isnumeric(v) == 1
            if isscalar(v)
                subid = num2str(v);
                break
            end
        end

        if iscategorical(v) == 1
            subid = char(string(v(1)));
            break
        end
    end
end

if isempty(subid) == 1
    [~, nm, ~] = fileparts(filepath);

    if isempty(nm) == 0
        subid = nm;
    else
        subid = ['sub' num2str(fallbackIndex)];
    end
end

subid = strtrim(subid);

if isempty(subid) == 1
    subid = ['sub' num2str(fallbackIndex)];
end

end

function sub = local_make_subid_value(filename, fallbackIdx)

subParsed = local_parse_subid_from_filename(filename);

if isempty(subParsed) == 0
    sub = subParsed;
else
    [~, base, ~] = fileparts(filename);

    if isempty(base) == 0
        sub = base;
    else
        sub = sprintf('sub%d', fallbackIdx);
    end
end

sub = char(string(sub));
sub = strtrim(sub);

if isempty(sub) == 1
    sub = sprintf('sub%d', fallbackIdx);
end

end

function sub = local_parse_subid_from_filename(filename)

sub = [];

[~, base, ~] = fileparts(filename);

tok = regexp(base, 'sub-([A-Za-z0-9]+)', 'tokens', 'once');

if isempty(tok) == 0
    sub = tok{1};
end

end

% Collect top-level numeric TF fields whose names match power/phase prefixes.
% Enforces "TOP-LEVEL SIBLINGS" convention and ignores any nested tf.power / tf.phase containers.
function [outStruct, rawNames] = local_collect_top_level_tf_numeric(tf, prefixes, filepath, opts)

outStruct = struct();
rawNames = {};

fns = fieldnames(tf);

for i = 1:numel(fns)

    fn = fns{i};

    % Never treat nested containers as data sources
    if strcmpi(fn, 'power') == 1 || strcmpi(fn, 'phase') == 1
        v = tf.(fn);

        if isstruct(v) == 1
            continue
        end
    end

    isMatch = false;

    for j = 1:numel(prefixes)
        pfx = prefixes{j};

        if strncmpi(fn, pfx, numel(pfx)) == 1
            isMatch = true;
        end
    end

    if isMatch == 0
        continue
    end

    X = tf.(fn);

    if isnumeric(X) == 0
        continue
    end

    if isempty(X) == 1
        continue
    end

    if ndims(X) ~= 4
        continue
    end

    key = matlab.lang.makeValidName(fn);

    if isfield(outStruct, key) == 1
        if opts.strict == 1
            error('nf_aggregate:NameCollision', ...
                'Name collision for %s in %s', key, filepath);
        else
            continue
        end
    end

    outStruct.(key) = X;
    rawNames{end + 1, 1} = fn; %#ok<AGROW>

end

end

% Collect top-level measure structs (power*/phase*) for STAT aggregation.
% Example expected layout:
%   tf.power.corrcoef (numeric)
%   tf.power.tstat (numeric)
%   tf.power.df (numeric or scalar)
%   tf.power.type (string/char)
% and similarly for tf.power_induced_ap, tf.phase_itc, etc.
function measures = local_collect_top_level_measure_structs(tf)

measures = struct();

fn = fieldnames(tf);

for i = 1:numel(fn)

    name = fn{i};

    if strcmpi(name, 'stat') == 1
        continue
    end

    if contains(lower(name), 'power') == 0 && contains(lower(name), 'phase') == 0
        continue
    end

    v = tf.(name);

    if isstruct(v) == 0
        continue
    end

    if isempty(fieldnames(v)) == 1
        continue
    end

    key = matlab.lang.makeValidName(name);

    if isfield(measures, key) == 1
        error('nf_aggregate:NameCollision', ...
            'makeValidName collision (measure struct) in %s: "%s"', 'tf', key);
    end

    measures.(key) = v;

end

end






