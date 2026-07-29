function [NF, report] = nf_concat(files, varargin)
%NF_CONCAT Concatenate compatible NeuroFreq files or structures across trials.
%
%   [NF, REPORT] = NF_CONCAT(FILES)
%
%   FILES can be one of the following:
%       1. A cell array of NeuroFreq structures.
%       2. A cell array of .mat files containing one NeuroFreq structure each.
%       3. A string array or cell array of filenames.
%       4. A dir() structure array.
%       5. A wildcard pattern resolving to files.
%       6. A cell array of EEGLAB .set files, if setTransformArgs is supplied.
%
%   The function concatenates trial-resolved NeuroFreq fields along dimension 4.
%   By default, it keeps the modal time vector across runs and drops runs with
%   a nonmatching time vector. If two or more runs have the usual long epoch and
%   one run has a short or otherwise strange epoch, the strange run is dropped.
%
%   Name/value options:
%
%       timeMode
%           'majority'   default. Keep runs matching the modal time vector.
%                        Ties are resolved by choosing the longest time vector.
%           'intersect'  keep all runs, but keep only time points present in
%                        every run.
%           'strict'     require all time vectors to match.
%           'first'      keep runs matching the first file's time vector.
%           'longest'    keep runs matching the longest time vector.
%
%       timeTolerance
%           Tolerance used when comparing time vectors. Default is 1e-6.
%
%       freqTolerance
%           Tolerance used when comparing frequency vectors and numeric scalar
%           metadata. Default is 1e-6.
%
%       chanMode
%           'strict' requires matching channel labels when labels are present.
%           'count' only requires the same number of channels.
%           'none' skips channel checks. Default is 'strict'.
%
%       dataFields
%           Cell array of trial-resolved fields to concatenate. Default is
%           automatic detection, preferring power and phase.
%
%       setTransformArgs
%           Cell array passed to nf_tftransform when FILES are .set files.
%           Example: {'method', 'stft', 'freqs', 1:30}
%
%       matVariable
%           Variable name to load from each .mat file. Default is automatic
%           detection of one NeuroFreq-like scalar structure.
%
%       runLabels
%           Numeric, string, or cell vector giving run labels. Used in REPORT
%           and by behaviorTrialMode = 'fixed' when numeric.
%
%       behaviorTrialMode
%           'none'       default. Preserve behavior trial fields.
%           'cumulative' add the cumulative number of retained trials to
%                        behaviorTrialFields.
%           'fixed'      add behaviorTrialOffset * (runLabel - 1) when numeric
%                        runLabels are supplied; otherwise add
%                        behaviorTrialOffset * (inputIndex - 1).
%
%       behaviorTrialOffset
%           Scalar offset used when behaviorTrialMode is 'fixed'.
%
%       behaviorTrialFields
%           Behavior fields to offset. Default is {'trial','prevtrial','nexttrial'}.
%
%       offsetEventEpoch
%           Offset event.epoch by cumulative retained trials. Default is true.
%
%       offsetEpochEvent
%           Offset epoch.event by cumulative event count. Default is true.
%
%       addMergeInfo
%           Add REPORT to NF.merge. Default is true.
%
%       verbose
%           Show warnings about dropped runs and nonfatal decisions. Default is true.
%
%   Example, from EEGLAB .set files:
%
%       filestmp = dir('/path/to/sub-*_task-ContrastChangeDetection*_stim*.set');
%
%       TF = nf_concat(filestmp, ...
%           'setTransformArgs', {'method', 'stft', 'freqs', 1:1:30}, ...
%           'timeMode', 'majority');
%
%   Example, keep all runs but only full-overlap time points:
%
%       TF = nf_concat(filestmp, ...
%           'setTransformArgs', {'method', 'stft', 'freqs', 1:1:30}, ...
%           'timeMode', 'intersect');

parser = inputParser;
parser.FunctionName = mfilename;

addRequired(parser, 'files');
addParameter(parser, 'timeMode', 'majority');
addParameter(parser, 'timeTolerance', 1e-6);
addParameter(parser, 'freqTolerance', 1e-6);
addParameter(parser, 'chanMode', 'strict');
addParameter(parser, 'dataFields', {});
addParameter(parser, 'trialDim', 4);
addParameter(parser, 'timeDim', 3);
addParameter(parser, 'setTransformArgs', {});
addParameter(parser, 'matVariable', '');
addParameter(parser, 'runLabels', []);
addParameter(parser, 'behaviorTrialMode', 'none');
addParameter(parser, 'behaviorTrialOffset', []);
addParameter(parser, 'behaviorTrialFields', {'trial', 'prevtrial', 'nexttrial'});
addParameter(parser, 'offsetEventEpoch', true);
addParameter(parser, 'offsetEpochEvent', true);
addParameter(parser, 'addMergeInfo', true);
addParameter(parser, 'verbose', true);

parse(parser, files, varargin{:});

opts = parser.Results;
opts.timeMode = normalize_mode_string(opts.timeMode);
opts.chanMode = normalize_mode_string(opts.chanMode);
opts.behaviorTrialMode = normalize_mode_string(opts.behaviorTrialMode);
opts.dataFields = normalize_field_list(opts.dataFields, 'dataFields');
opts.behaviorTrialFields = normalize_field_list(opts.behaviorTrialFields, 'behaviorTrialFields');
opts.matVariable = char(string(opts.matVariable));

if strcmp(opts.timeMode, 'intersection')
    opts.timeMode = 'intersect';
end

if strcmp(opts.chanMode, 'size')
    opts.chanMode = 'count';
end

if strcmp(opts.behaviorTrialMode, 'preserve')
    opts.behaviorTrialMode = 'none';
end

if strcmp(opts.behaviorTrialMode, 'cumulativentrials')
    opts.behaviorTrialMode = 'cumulative';
end

validate_options(opts);

inputItems = normalize_nf_inputs(files);

if isempty(inputItems)
    error('nf_concat:NoInputs', 'No input files or structures were supplied.');
end

nInput = numel(inputItems);

if isempty(opts.runLabels)
    opts.runLabels = 1:nInput;
else
    if numel(opts.runLabels) ~= nInput
        error('nf_concat:BadRunLabels', 'runLabels must have one element per input.');
    end
end

TF = cell(1, nInput);
sources = cell(1, nInput);

for inputIndex = 1:nInput
    [TF{inputIndex}, sources{inputIndex}] = load_nf_input(inputItems{inputIndex}, opts, inputIndex);
    validate_nf_core(TF{inputIndex}, sources{inputIndex});
end

timeInfo = choose_time_policy(TF, opts, sources);

if ~any(timeInfo.keep)
    error('nf_concat:NoKeptInputs', 'No inputs survived the time compatibility policy.');
end

templateIndex = timeInfo.templateIndex;

if isempty(opts.dataFields)
    dataFields = resolve_data_fields(TF{templateIndex}, opts);
else
    dataFields = opts.dataFields;
end

if isempty(dataFields)
    error('nf_concat:NoDataFields', 'No trial-resolved NeuroFreq data fields were found.');
end

keptInputIndices = find(timeInfo.keep);
nKept = numel(keptInputIndices);
TFkept = cell(1, nKept);
keptSources = cell(1, nKept);
keptRunLabels = get_indexed_labels(opts.runLabels, keptInputIndices);

for keptIndex = 1:nKept
    inputIndex = keptInputIndices(keptIndex);
    TFtmp = TF{inputIndex};
    TFtmp = apply_time_selection(TFtmp, dataFields, timeInfo.timeIndices{inputIndex}, timeInfo.times, opts, sources{inputIndex});
    TFkept{keptIndex} = TFtmp;
    keptSources{keptIndex} = sources{inputIndex};
end

check_compatibility(TFkept, dataFields, opts, keptSources);

nTrialsByRun = zeros(1, nKept);
eventCountsByRun = zeros(1, nKept);

for keptIndex = 1:nKept
    nTrialsByRun(keptIndex) = get_ntrls_from_data(TFkept{keptIndex}, dataFields{1}, opts.trialDim, keptSources{keptIndex});

    if isfield(TFkept{keptIndex}, 'event')
        eventCountsByRun(keptIndex) = numel(TFkept{keptIndex}.event);
    else
        eventCountsByRun(keptIndex) = 0;
    end
end

NF = TFkept{1};

for fieldIndex = 1:numel(dataFields)
    fieldName = dataFields{fieldIndex};
    pieces = cell(1, nKept);

    for keptIndex = 1:nKept
        pieces{keptIndex} = TFkept{keptIndex}.(fieldName);
    end

    NF.(fieldName) = cat(opts.trialDim, pieces{:});
end

NF.times = cast_vector_like(timeInfo.times, TFkept{1}.times);
NF.ntrls = sum(nTrialsByRun);

if isfield(NF, 'nsensor')
    NF.nsensor = size(NF.(dataFields{1}), 1);
end

if any_field_present(TFkept, 'event')
    NF.event = merge_event_fields(TFkept, nTrialsByRun, opts);
end

if any_field_present(TFkept, 'epoch')
    NF.epoch = merge_epoch_fields(TFkept, eventCountsByRun, opts);
end

if any_field_present(TFkept, 'behavior')
    NF.behavior = merge_behavior_fields(TFkept, nTrialsByRun, keptInputIndices, opts);
end

report = make_merge_report(sources, keptSources, keptInputIndices, keptRunLabels, timeInfo, dataFields, nTrialsByRun, eventCountsByRun, opts);

if opts.addMergeInfo
    NF.merge = report;
end

end

function validate_options(opts)

validTimeModes = {'majority', 'intersect', 'strict', 'first', 'longest'};
validChanModes = {'strict', 'count', 'none'};
validBehaviorModes = {'none', 'cumulative', 'fixed'};

if ~any(strcmp(opts.timeMode, validTimeModes))
    error('nf_concat:BadTimeMode', 'Unknown timeMode: %s.', opts.timeMode);
end

if ~any(strcmp(opts.chanMode, validChanModes))
    error('nf_concat:BadChanMode', 'Unknown chanMode: %s.', opts.chanMode);
end

if ~any(strcmp(opts.behaviorTrialMode, validBehaviorModes))
    error('nf_concat:BadBehaviorTrialMode', 'Unknown behaviorTrialMode: %s.', opts.behaviorTrialMode);
end

if ~isnumeric(opts.timeTolerance)
    error('nf_concat:BadTimeTolerance', 'timeTolerance must be numeric.');
end

if ~isscalar(opts.timeTolerance)
    error('nf_concat:BadTimeTolerance', 'timeTolerance must be scalar.');
end

if opts.timeTolerance < 0
    error('nf_concat:BadTimeTolerance', 'timeTolerance must be nonnegative.');
end

if ~isnumeric(opts.freqTolerance)
    error('nf_concat:BadFreqTolerance', 'freqTolerance must be numeric.');
end

if ~isscalar(opts.freqTolerance)
    error('nf_concat:BadFreqTolerance', 'freqTolerance must be scalar.');
end

if opts.freqTolerance < 0
    error('nf_concat:BadFreqTolerance', 'freqTolerance must be nonnegative.');
end

if ~isnumeric(opts.trialDim)
    error('nf_concat:BadTrialDim', 'trialDim must be numeric.');
end

if ~isscalar(opts.trialDim)
    error('nf_concat:BadTrialDim', 'trialDim must be scalar.');
end

if opts.trialDim < 1
    error('nf_concat:BadTrialDim', 'trialDim must be positive.');
end

if ~isnumeric(opts.timeDim)
    error('nf_concat:BadTimeDim', 'timeDim must be numeric.');
end

if ~isscalar(opts.timeDim)
    error('nf_concat:BadTimeDim', 'timeDim must be scalar.');
end

if opts.timeDim < 1
    error('nf_concat:BadTimeDim', 'timeDim must be positive.');
end

if opts.timeDim == opts.trialDim
    error('nf_concat:BadDimensions', 'timeDim and trialDim must be different.');
end

if ~iscell(opts.setTransformArgs)
    error('nf_concat:BadSetTransformArgs', 'setTransformArgs must be a cell array.');
end

if strcmp(opts.behaviorTrialMode, 'fixed')
    if isempty(opts.behaviorTrialOffset)
        error('nf_concat:MissingBehaviorTrialOffset', 'behaviorTrialOffset is required when behaviorTrialMode is fixed.');
    end

    if ~isnumeric(opts.behaviorTrialOffset)
        error('nf_concat:BadBehaviorTrialOffset', 'behaviorTrialOffset must be numeric.');
    end

    if ~isscalar(opts.behaviorTrialOffset)
        error('nf_concat:BadBehaviorTrialOffset', 'behaviorTrialOffset must be scalar.');
    end
end

end

function inputItems = normalize_nf_inputs(files)

inputItems = {};

if iscell(files)
    for fileIndex = 1:numel(files)
        inputItems = append_nf_input(inputItems, files{fileIndex});
    end
elseif isstring(files)
    fileCell = cellstr(files(:));

    for fileIndex = 1:numel(fileCell)
        inputItems = append_nf_input(inputItems, fileCell{fileIndex});
    end
elseif ischar(files)
    inputItems = append_nf_input(inputItems, files);
elseif isstruct(files)
    if is_dir_listing(files)
        for fileIndex = 1:numel(files)
            inputItems{end + 1} = fullfile(files(fileIndex).folder, files(fileIndex).name);
        end
    else
        for fileIndex = 1:numel(files)
            inputItems{end + 1} = files(fileIndex);
        end
    end
else
    error('nf_concat:BadInputType', 'FILES must be paths, a dir listing, NeuroFreq structures, or a cell array containing those items.');
end

inputItems = inputItems(:)';

end

function inputItems = append_nf_input(inputItems, item)

if isstring(item)
    item = char(item);
end

if ischar(item)
    if has_wildcard(item)
        listing = dir(item);

        if isempty(listing)
            error('nf_concat:NoWildcardMatches', 'No files matched pattern: %s.', item);
        end

        listing = sort_dir_listing(listing);

        for listingIndex = 1:numel(listing)
            inputItems{end + 1} = fullfile(listing(listingIndex).folder, listing(listingIndex).name);
        end
    else
        inputItems{end + 1} = item;
    end
elseif isstruct(item)
    if is_dir_listing(item)
        for listingIndex = 1:numel(item)
            inputItems{end + 1} = fullfile(item(listingIndex).folder, item(listingIndex).name);
        end
    else
        for structIndex = 1:numel(item)
            inputItems{end + 1} = item(structIndex);
        end
    end
else
    error('nf_concat:BadInputCell', 'Each cell element must be a path, a dir listing, or a NeuroFreq structure.');
end

end

function tf = has_wildcard(textValue)

starHit = strfind(textValue, '*');
questionHit = strfind(textValue, '?');

if isempty(starHit)
    hasStar = false;
else
    hasStar = true;
end

if isempty(questionHit)
    hasQuestion = false;
else
    hasQuestion = true;
end

tf = hasStar || hasQuestion;

end

function listing = sort_dir_listing(listing)

names = cell(1, numel(listing));

for listingIndex = 1:numel(listing)
    names{listingIndex} = listing(listingIndex).name;
end

[~, order] = sort(names);
listing = listing(order);

end

function tf = is_dir_listing(value)

tf = false;

if ~isstruct(value)
    return;
end

if ~isfield(value, 'folder')
    return;
end

if ~isfield(value, 'name')
    return;
end

tf = true;

end

function [TF, source] = load_nf_input(inputItem, opts, inputIndex)

if isstruct(inputItem)
    TF = inputItem;
    source = sprintf('struct input %d', inputIndex);
    return;
end

if isstring(inputItem)
    inputItem = char(inputItem);
end

if ~ischar(inputItem)
    error('nf_concat:BadInputItem', 'Input %d is neither a file path nor a structure.', inputIndex);
end

fileName = inputItem;
source = fileName;

if exist(fileName, 'file') ~= 2
    error('nf_concat:MissingFile', 'File does not exist: %s.', fileName);
end

[filePath, baseName, extension] = fileparts(fileName);
extension = lower(extension);

if strcmp(extension, '.mat')
    TF = load_nf_mat(fileName, opts.matVariable);
elseif strcmp(extension, '.set')
    if isempty(opts.setTransformArgs)
        error('nf_concat:SetNeedsTransformArgs', 'Input is a .set file, so setTransformArgs must be supplied: %s.', fileName);
    end

    if exist('pop_loadset', 'file') ~= 2
        error('nf_concat:MissingPopLoadset', 'pop_loadset is not on the MATLAB path.');
    end

    if exist('nf_tftransform', 'file') ~= 2
        error('nf_concat:MissingTfTransform', 'nf_tftransform is not on the MATLAB path.');
    end

    EEG = pop_loadset('filename', [baseName extension], 'filepath', filePath);
    TF = nf_tftransform(EEG, opts.setTransformArgs{:});
else
    error('nf_concat:UnsupportedFileType', 'Unsupported file type for %s. Use .mat, .set, or pass structures directly.', fileName);
end

end

function TF = load_nf_mat(fileName, matVariable)

if ~isempty(matVariable)
    loaded = load(fileName, matVariable);

    if ~isfield(loaded, matVariable)
        error('nf_concat:MissingMatVariable', 'Variable %s was not found in %s.', matVariable, fileName);
    end

    TF = loaded.(matVariable);

    if ~is_neurofreq_like(TF)
        error('nf_concat:BadMatVariable', 'Variable %s in %s does not look like a NeuroFreq structure.', matVariable, fileName);
    end

    return;
end

loaded = load(fileName);
variableNames = fieldnames(loaded);
candidateNames = {};

for variableIndex = 1:numel(variableNames)
    candidate = loaded.(variableNames{variableIndex});

    if is_neurofreq_like(candidate)
        candidateNames{end + 1} = variableNames{variableIndex};
    end
end

if isempty(candidateNames)
    error('nf_concat:NoNeuroFreqVariable', 'No NeuroFreq-like scalar structure was found in %s.', fileName);
end

if numel(candidateNames) == 1
    TF = loaded.(candidateNames{1});
    return;
end

preferredNames = {'NF', 'TF', 'tf', 'neurofreq', 'neuroFreq', 'contrastTF', 'contrastTFF'};

for preferredIndex = 1:numel(preferredNames)
    thisName = preferredNames{preferredIndex};

    if isfield(loaded, thisName)
        if is_neurofreq_like(loaded.(thisName))
            TF = loaded.(thisName);
            return;
        end
    end
end

error('nf_concat:AmbiguousMatFile', 'More than one NeuroFreq-like structure was found in %s. Supply matVariable.', fileName);

end

function tf = is_neurofreq_like(value)

tf = false;

if ~isstruct(value)
    return;
end

if ~isscalar(value)
    return;
end

if ~isfield(value, 'freqs')
    return;
end

if ~isfield(value, 'times')
    return;
end

tf = true;

end

function validate_nf_core(TF, source)

if ~isstruct(TF)
    error('nf_concat:BadNeuroFreqStruct', 'Input %s is not a structure.', source);
end

if ~isscalar(TF)
    error('nf_concat:BadNeuroFreqStruct', 'Input %s is not a scalar NeuroFreq structure.', source);
end

if ~isfield(TF, 'freqs')
    error('nf_concat:MissingFreqs', 'Input %s is missing field freqs.', source);
end

if ~isfield(TF, 'times')
    error('nf_concat:MissingTimes', 'Input %s is missing field times.', source);
end

if ~isnumeric(TF.freqs)
    error('nf_concat:BadFreqs', 'Input %s has nonnumeric freqs.', source);
end

if ~isnumeric(TF.times)
    error('nf_concat:BadTimes', 'Input %s has nonnumeric times.', source);
end

if ~isvector(TF.freqs)
    error('nf_concat:BadFreqs', 'Input %s has nonvector freqs.', source);
end

if ~isvector(TF.times)
    error('nf_concat:BadTimes', 'Input %s has nonvector times.', source);
end

if isempty(TF.freqs)
    error('nf_concat:EmptyFreqs', 'Input %s has empty freqs.', source);
end

if isempty(TF.times)
    error('nf_concat:EmptyTimes', 'Input %s has empty times.', source);
end

if any(isnan(double(TF.freqs(:))))
    error('nf_concat:BadFreqs', 'Input %s has NaN freqs.', source);
end

if any(isnan(double(TF.times(:))))
    error('nf_concat:BadTimes', 'Input %s has NaN times.', source);
end

if isfield(TF, 'ntrls')
    if ~isnumeric(TF.ntrls)
        error('nf_concat:BadNtrls', 'Input %s has nonnumeric ntrls.', source);
    end

    if ~isscalar(TF.ntrls)
        error('nf_concat:BadNtrls', 'Input %s has nonscalar ntrls.', source);
    end

    if TF.ntrls < 0
        error('nf_concat:BadNtrls', 'Input %s has negative ntrls.', source);
    end
end

end

function dataFields = resolve_data_fields(TF, opts)

fieldNames = fieldnames(TF);
dataFields = {};
knownTrialFields = {'power', 'phase'};

for knownIndex = 1:numel(knownTrialFields)
    fieldName = knownTrialFields{knownIndex};

    if isfield(TF, fieldName)
        if is_candidate_data_field(TF, fieldName, opts, true)
            dataFields{end + 1} = fieldName;
        end
    end
end

for fieldIndex = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIndex};

    if any(strcmp(fieldName, dataFields))
        continue;
    end

    if is_candidate_data_field(TF, fieldName, opts, false)
        dataFields{end + 1} = fieldName;
    end
end

end

function tf = is_candidate_data_field(TF, fieldName, opts, allowCollapsedSingletonTrial)

tf = false;

value = TF.(fieldName);

if ~isnumeric(value)
    return;
end

nFreq = numel(TF.freqs);
nTime = numel(TF.times);

if isfield(TF, 'nsensor')
    nChan = TF.nsensor;
else
    nChan = size(value, 1);
end

if size(value, 1) ~= nChan
    return;
end

if size(value, 2) ~= nFreq
    return;
end

if size(value, opts.timeDim) ~= nTime
    return;
end

nTrials = get_ntrls_from_struct(TF, fieldName, opts.trialDim);

if isempty(nTrials)
    return;
end

if size(value, opts.trialDim) ~= nTrials
    return;
end

if ndims(value) >= opts.trialDim
    tf = true;
    return;
end

if allowCollapsedSingletonTrial
    if nTrials == 1
        tf = true;
        return;
    end
end

end

function nTrials = get_ntrls_from_struct(TF, fieldName, trialDim)

nTrials = [];

if isfield(TF, 'ntrls')
    if isnumeric(TF.ntrls)
        if isscalar(TF.ntrls)
            nTrials = double(TF.ntrls);
            return;
        end
    end
end

if nargin > 1
    if ~isempty(fieldName)
        if isfield(TF, fieldName)
            nTrials = size(TF.(fieldName), trialDim);
            return;
        end
    end
end

end

function timeInfo = choose_time_policy(TF, opts, sources)

nInput = numel(TF);
timeVectors = cell(1, nInput);

for inputIndex = 1:nInput
    timeVectors{inputIndex} = double(TF{inputIndex}.times(:)');
end

if strcmp(opts.timeMode, 'intersect')
    [commonTimes, commonIndices] = intersect_time_vectors(timeVectors, opts.timeTolerance);

    if isempty(commonTimes)
        error('nf_concat:EmptyTimeIntersection', 'No time points are common to all inputs.');
    end

    timeInfo.keep = true(1, nInput);
    timeInfo.timeIndices = commonIndices;
    timeInfo.times = commonTimes;
    timeInfo.templateIndex = 1;
    timeInfo.templateCount = nInput;
    timeInfo.nTimeOriginal = get_time_lengths(timeVectors);
    timeInfo.nTimeOutput = numel(commonTimes);
    timeInfo.dropReason = repmat({''}, 1, nInput);
    timeInfo.policyNote = 'All inputs kept; time dimension reduced to the across-input intersection.';
    return;
end

[groupIds, groupFirstIndices, groupCounts, groupLengths] = group_time_vectors(timeVectors, opts.timeTolerance);

if strcmp(opts.timeMode, 'strict')
    if numel(groupCounts) > 1
        message = make_time_mismatch_message(sources, timeVectors);
        error('nf_concat:TimeMismatch', 'timeMode is strict, but time vectors differ.\n%s', message);
    end

    selectedGroup = 1;
elseif strcmp(opts.timeMode, 'first')
    selectedGroup = groupIds(1);
elseif strcmp(opts.timeMode, 'longest')
    selectedGroup = select_longest_group(groupCounts, groupLengths);
else
    selectedGroup = select_majority_group(groupCounts, groupLengths);
end

templateIndex = groupFirstIndices(selectedGroup);
keep = groupIds == selectedGroup;
timeIndices = cell(1, nInput);

templateTimes = timeVectors{templateIndex};

for inputIndex = 1:nInput
    if keep(inputIndex)
        timeIndices{inputIndex} = 1:numel(timeVectors{inputIndex});
    else
        timeIndices{inputIndex} = [];
    end
end

dropReason = repmat({''}, 1, nInput);

for inputIndex = 1:nInput
    if ~keep(inputIndex)
        dropReason{inputIndex} = 'time vector did not match the selected time template';
    end
end

if opts.verbose
    if any(~keep)
        warning('nf_concat:TimeTemplateDrop', 'Dropping %d of %d inputs because their time vectors do not match the selected template. Kept template has %d time points.', sum(~keep), nInput, numel(templateTimes));
    end

    if strcmp(opts.timeMode, 'majority')
        if groupCounts(selectedGroup) == 1
            if nInput > 1
                warning('nf_concat:SingletonTimeTemplate', 'No time vector was repeated across inputs. The longest time vector was selected because timeMode is majority. Use timeMode = intersect to keep all inputs.');
            end
        end
    end
end

timeInfo.keep = keep;
timeInfo.timeIndices = timeIndices;
timeInfo.times = templateTimes;
timeInfo.templateIndex = templateIndex;
timeInfo.templateCount = groupCounts(selectedGroup);
timeInfo.nTimeOriginal = get_time_lengths(timeVectors);
timeInfo.nTimeOutput = numel(templateTimes);
timeInfo.dropReason = dropReason;
timeInfo.policyNote = sprintf('Selected time template from input %d with %d matching input(s).', templateIndex, groupCounts(selectedGroup));

end

function lengths = get_time_lengths(timeVectors)

lengths = zeros(1, numel(timeVectors));

for inputIndex = 1:numel(timeVectors)
    lengths(inputIndex) = numel(timeVectors{inputIndex});
end

end

function [groupIds, groupFirstIndices, groupCounts, groupLengths] = group_time_vectors(timeVectors, tolerance)

nInput = numel(timeVectors);
groupIds = zeros(1, nInput);
groupFirstIndices = [];
groupCounts = [];
groupLengths = [];

for inputIndex = 1:nInput
    assigned = false;

    for groupIndex = 1:numel(groupFirstIndices)
        representativeIndex = groupFirstIndices(groupIndex);

        if vectors_equal_with_tolerance(timeVectors{inputIndex}, timeVectors{representativeIndex}, tolerance)
            groupIds(inputIndex) = groupIndex;
            groupCounts(groupIndex) = groupCounts(groupIndex) + 1;
            assigned = true;
            break;
        end
    end

    if ~assigned
        groupFirstIndices(end + 1) = inputIndex;
        groupCounts(end + 1) = 1;
        groupLengths(end + 1) = numel(timeVectors{inputIndex});
        groupIds(inputIndex) = numel(groupFirstIndices);
    end
end

end

function selectedGroup = select_majority_group(groupCounts, groupLengths)

selectedGroup = 1;

for groupIndex = 2:numel(groupCounts)
    if groupCounts(groupIndex) > groupCounts(selectedGroup)
        selectedGroup = groupIndex;
    elseif groupCounts(groupIndex) == groupCounts(selectedGroup)
        if groupLengths(groupIndex) > groupLengths(selectedGroup)
            selectedGroup = groupIndex;
        end
    end
end

end

function selectedGroup = select_longest_group(groupCounts, groupLengths)

selectedGroup = 1;

for groupIndex = 2:numel(groupLengths)
    if groupLengths(groupIndex) > groupLengths(selectedGroup)
        selectedGroup = groupIndex;
    elseif groupLengths(groupIndex) == groupLengths(selectedGroup)
        if groupCounts(groupIndex) > groupCounts(selectedGroup)
            selectedGroup = groupIndex;
        end
    end
end

end

function message = make_time_mismatch_message(sources, timeVectors)

lines = cell(1, numel(sources));

for inputIndex = 1:numel(sources)
    thisTimes = timeVectors{inputIndex};
    lines{inputIndex} = sprintf('  %d: %s | nTimes = %d | first = %.12g | last = %.12g', inputIndex, sources{inputIndex}, numel(thisTimes), thisTimes(1), thisTimes(end));
end

message = strjoin(lines, newline);

end

function [commonTimes, commonIndices] = intersect_time_vectors(timeVectors, tolerance)

nInput = numel(timeVectors);
commonTimes = timeVectors{1};
commonIndices = cell(1, nInput);
commonIndices{1} = 1:numel(commonTimes);

for inputIndex = 2:nInput
    thisTimes = timeVectors{inputIndex};
    keepPosition = false(1, numel(commonTimes));
    thisIndices = nan(1, numel(commonTimes));

    for timeIndex = 1:numel(commonTimes)
        difference = abs(thisTimes - commonTimes(timeIndex));
        matchedIndex = find(difference <= tolerance, 1, 'first');

        if ~isempty(matchedIndex)
            keepPosition(timeIndex) = true;
            thisIndices(timeIndex) = matchedIndex;
        end
    end

    commonTimes = commonTimes(keepPosition);

    for previousIndex = 1:(inputIndex - 1)
        commonIndices{previousIndex} = commonIndices{previousIndex}(keepPosition);
    end

    commonIndices{inputIndex} = thisIndices(keepPosition);
end

end

function TF = apply_time_selection(TF, dataFields, timeIndices, timesOut, opts, source)

if isempty(timeIndices)
    error('nf_concat:InternalTimeSelection', 'Attempted to time-subset a dropped input: %s.', source);
end

for fieldIndex = 1:numel(dataFields)
    fieldName = dataFields{fieldIndex};

    if ~isfield(TF, fieldName)
        error('nf_concat:MissingDataField', 'Input %s is missing data field %s.', source, fieldName);
    end

    nTime = size(TF.(fieldName), opts.timeDim);

    if max(timeIndices) > nTime
        error('nf_concat:BadTimeIndex', 'Time index exceeds field %s time dimension in %s.', fieldName, source);
    end

    TF.(fieldName) = subset_along_dimension(TF.(fieldName), opts.timeDim, timeIndices);
end

TF.times = cast_vector_like(timesOut, TF.times);

end

function data = subset_along_dimension(data, dimension, indices)

nDimensions = max(ndims(data), dimension);
subscripts = repmat({':'}, 1, nDimensions);
subscripts{dimension} = indices;
data = data(subscripts{:});

end

function out = cast_vector_like(values, template)

if isa(template, 'single')
    castValues = single(values);
else
    castValues = double(values);
end

if isrow(template)
    out = reshape(castValues, 1, []);
else
    out = reshape(castValues, [], 1);
end

end

function check_compatibility(TFlist, dataFields, opts, sources)

if isempty(TFlist)
    error('nf_concat:NoCompatibilityInputs', 'No inputs were supplied for compatibility checking.');
end

check_data_field_compatibility(TFlist, dataFields, opts, sources);
check_frequency_compatibility(TFlist, opts, sources);
check_time_compatibility_after_selection(TFlist, opts, sources);
check_channel_compatibility(TFlist, dataFields{1}, opts, sources);
check_scalar_metadata_compatibility(TFlist, opts, sources);
check_ntrls_compatibility(TFlist, dataFields, opts, sources);

end

function check_data_field_compatibility(TFlist, dataFields, opts, sources)

nInputs = numel(TFlist);

template = TFlist{1};

for fieldIndex = 1:numel(dataFields)
    fieldName = dataFields{fieldIndex};

    if ~isfield(template, fieldName)
        error('nf_concat:MissingTemplateField', 'The template input is missing data field %s.', fieldName);
    end

    templateData = template.(fieldName);
    templateSize = full_size(templateData, max([ndims(templateData), opts.trialDim, opts.timeDim]));

    for inputIndex = 1:nInputs
        thisTF = TFlist{inputIndex};

        if ~isfield(thisTF, fieldName)
            error('nf_concat:MissingDataField', 'Input %s is missing data field %s.', sources{inputIndex}, fieldName);
        end

        thisData = thisTF.(fieldName);

        if ~isnumeric(thisData)
            error('nf_concat:BadDataField', 'Field %s in %s is not numeric.', fieldName, sources{inputIndex});
        end

        thisSize = full_size(thisData, numel(templateSize));

        if numel(thisSize) ~= numel(templateSize)
            error('nf_concat:BadDataFieldSize', 'Field %s in %s has incompatible dimensionality.', fieldName, sources{inputIndex});
        end

        for dimension = 1:numel(templateSize)
            if dimension == opts.trialDim
                continue;
            end

            if thisSize(dimension) ~= templateSize(dimension)
                error('nf_concat:BadDataFieldSize', 'Field %s in %s has size %d in dimension %d; expected %d.', fieldName, sources{inputIndex}, thisSize(dimension), dimension, templateSize(dimension));
            end
        end
    end
end

end

function dims = full_size(data, nDimensions)

baseSize = size(data);
nDimensions = max(nDimensions, numel(baseSize));
dims = ones(1, nDimensions);
dims(1:numel(baseSize)) = baseSize;

end

function check_frequency_compatibility(TFlist, opts, sources)

templateFreqs = double(TFlist{1}.freqs(:)');

for inputIndex = 2:numel(TFlist)
    thisFreqs = double(TFlist{inputIndex}.freqs(:)');

    if ~vectors_equal_with_tolerance(templateFreqs, thisFreqs, opts.freqTolerance)
        error('nf_concat:FrequencyMismatch', 'Frequency vectors do not match between %s and %s.', sources{1}, sources{inputIndex});
    end
end

end

function check_time_compatibility_after_selection(TFlist, opts, sources)

templateTimes = double(TFlist{1}.times(:)');

for inputIndex = 2:numel(TFlist)
    thisTimes = double(TFlist{inputIndex}.times(:)');

    if ~vectors_equal_with_tolerance(templateTimes, thisTimes, opts.timeTolerance)
        error('nf_concat:InternalTimeMismatch', 'Time vectors still differ after time selection between %s and %s.', sources{1}, sources{inputIndex});
    end
end

end

function check_channel_compatibility(TFlist, referenceField, opts, sources)

if strcmp(opts.chanMode, 'none')
    return;
end

nInputs = numel(TFlist);
templateNChan = size(TFlist{1}.(referenceField), 1);
templateLabels = get_channel_labels(TFlist{1});

for inputIndex = 1:nInputs
    thisNChan = size(TFlist{inputIndex}.(referenceField), 1);

    if thisNChan ~= templateNChan
        error('nf_concat:ChannelCountMismatch', 'Channel count differs in %s. Expected %d, found %d.', sources{inputIndex}, templateNChan, thisNChan);
    end

    if isfield(TFlist{inputIndex}, 'nsensor')
        if TFlist{inputIndex}.nsensor ~= thisNChan
            error('nf_concat:NsensorMismatch', 'nsensor does not match data size in %s.', sources{inputIndex});
        end
    end
end

if strcmp(opts.chanMode, 'count')
    return;
end

if isempty(templateLabels)
    return;
end

for inputIndex = 2:nInputs
    thisLabels = get_channel_labels(TFlist{inputIndex});

    if isempty(thisLabels)
        continue;
    end

    if numel(thisLabels) ~= numel(templateLabels)
        error('nf_concat:ChannelLabelMismatch', 'Channel label count differs in %s.', sources{inputIndex});
    end

    for channelIndex = 1:numel(templateLabels)
        if ~strcmp(templateLabels{channelIndex}, thisLabels{channelIndex})
            error('nf_concat:ChannelLabelMismatch', 'Channel label mismatch at channel %d between %s and %s.', channelIndex, sources{1}, sources{inputIndex});
        end
    end
end

end

function labels = get_channel_labels(TF)

labels = {};

if ~isfield(TF, 'chanlocs')
    return;
end

chanlocs = TF.chanlocs;

if isempty(chanlocs)
    return;
end

if ~isstruct(chanlocs)
    return;
end

if ~isfield(chanlocs, 'labels')
    return;
end

labels = cell(1, numel(chanlocs));

for channelIndex = 1:numel(chanlocs)
    labelValue = chanlocs(channelIndex).labels;

    if isstring(labelValue)
        labelValue = char(labelValue);
    end

    if ischar(labelValue)
        labels{channelIndex} = labelValue;
    else
        labels{channelIndex} = '';
    end
end

end

function check_scalar_metadata_compatibility(TFlist, opts, sources)

fieldsToCheck = {'Fs', 'window', 'overlap', 'method', 'scale'};
template = TFlist{1};

for fieldIndex = 1:numel(fieldsToCheck)
    fieldName = fieldsToCheck{fieldIndex};

    if ~isfield(template, fieldName)
        continue;
    end

    for inputIndex = 2:numel(TFlist)
        if ~isfield(TFlist{inputIndex}, fieldName)
            continue;
        end

        if ~values_equal(template.(fieldName), TFlist{inputIndex}.(fieldName), opts.freqTolerance)
            error('nf_concat:MetadataMismatch', 'Metadata field %s differs between %s and %s.', fieldName, sources{1}, sources{inputIndex});
        end
    end
end

end

function check_ntrls_compatibility(TFlist, dataFields, opts, sources)

for inputIndex = 1:numel(TFlist)
    referenceCount = get_ntrls_from_data(TFlist{inputIndex}, dataFields{1}, opts.trialDim, sources{inputIndex});

    for fieldIndex = 2:numel(dataFields)
        thisCount = get_ntrls_from_data(TFlist{inputIndex}, dataFields{fieldIndex}, opts.trialDim, sources{inputIndex});

        if thisCount ~= referenceCount
            error('nf_concat:TrialCountMismatch', 'Data fields have different trial counts in %s.', sources{inputIndex});
        end
    end

    if isfield(TFlist{inputIndex}, 'ntrls')
        if double(TFlist{inputIndex}.ntrls) ~= referenceCount
            error('nf_concat:NtrlsMismatch', 'ntrls is %d but data field %s has %d trials in %s.', TFlist{inputIndex}.ntrls, dataFields{1}, referenceCount, sources{inputIndex});
        end
    end
end

end

function nTrials = get_ntrls_from_data(TF, fieldName, trialDim, source)

if ~isfield(TF, fieldName)
    error('nf_concat:MissingDataField', 'Input %s is missing data field %s.', source, fieldName);
end

nTrials = size(TF.(fieldName), trialDim);

end

function tf = vectors_equal_with_tolerance(a, b, tolerance)

a = double(a(:)');
b = double(b(:)');

if numel(a) ~= numel(b)
    tf = false;
    return;
end

if isempty(a)
    tf = true;
    return;
end

difference = abs(a - b);

if max(difference) <= tolerance
    tf = true;
else
    tf = false;
end

end

function tf = values_equal(a, b, tolerance)

tf = false;

if isnumeric(a) && isnumeric(b)
    if ~isequal(size(a), size(b))
        return;
    end

    if isempty(a) && isempty(b)
        tf = true;
        return;
    end

    difference = abs(double(a(:)) - double(b(:)));

    if max(difference) <= tolerance
        tf = true;
    end

    return;
end

if ischar(a) || isstring(a)
    if ischar(b) || isstring(b)
        if strcmp(char(string(a)), char(string(b)))
            tf = true;
        end
    end

    return;
end

if islogical(a) && islogical(b)
    tf = isequal(a, b);
    return;
end

tf = isequaln(a, b);

end

function tf = any_field_present(TFlist, fieldName)

tf = false;

for inputIndex = 1:numel(TFlist)
    if isfield(TFlist{inputIndex}, fieldName)
        tf = true;
        return;
    end
end

end

function eventOut = merge_event_fields(TFlist, nTrialsByRun, opts)

parts = cell(1, numel(TFlist));
cumulativeTrials = 0;

for inputIndex = 1:numel(TFlist)
    if isfield(TFlist{inputIndex}, 'event')
        part = TFlist{inputIndex}.event;
    else
        part = struct([]);
    end

    part = reshape_struct_row(part);

    if opts.offsetEventEpoch
        part = offset_struct_field(part, 'epoch', cumulativeTrials);
    end

    parts{inputIndex} = part;
    cumulativeTrials = cumulativeTrials + nTrialsByRun(inputIndex);
end

eventOut = concatenate_struct_parts(parts);

end

function epochOut = merge_epoch_fields(TFlist, eventCountsByRun, opts)

parts = cell(1, numel(TFlist));
cumulativeEvents = 0;

for inputIndex = 1:numel(TFlist)
    if isfield(TFlist{inputIndex}, 'epoch')
        part = TFlist{inputIndex}.epoch;
    else
        part = struct([]);
    end

    part = reshape_struct_row(part);

    if opts.offsetEpochEvent
        part = offset_struct_field(part, 'event', cumulativeEvents);
    end

    parts{inputIndex} = part;
    cumulativeEvents = cumulativeEvents + eventCountsByRun(inputIndex);
end

epochOut = concatenate_struct_parts(parts);

end

function behaviorOut = merge_behavior_fields(TFlist, nTrialsByRun, keptInputIndices, opts)

parts = cell(1, numel(TFlist));
cumulativeTrials = 0;

for keptIndex = 1:numel(TFlist)
    if isfield(TFlist{keptIndex}, 'behavior')
        part = TFlist{keptIndex}.behavior;
    else
        part = struct([]);
    end

    part = reshape_struct_row(part);
    offsetAmount = get_behavior_trial_offset(cumulativeTrials, keptInputIndices(keptIndex), opts);

    if offsetAmount ~= 0
        for fieldIndex = 1:numel(opts.behaviorTrialFields)
            part = offset_struct_field(part, opts.behaviorTrialFields{fieldIndex}, offsetAmount);
        end
    end

    parts{keptIndex} = part;
    cumulativeTrials = cumulativeTrials + nTrialsByRun(keptIndex);
end

behaviorOut = concatenate_struct_parts(parts);

end

function offsetAmount = get_behavior_trial_offset(cumulativeTrials, originalInputIndex, opts)

offsetAmount = 0;

if strcmp(opts.behaviorTrialMode, 'none')
    return;
end

if strcmp(opts.behaviorTrialMode, 'cumulative')
    offsetAmount = cumulativeTrials;
    return;
end

if strcmp(opts.behaviorTrialMode, 'fixed')
    if isnumeric(opts.runLabels)
        runLabel = opts.runLabels(originalInputIndex);
        offsetAmount = opts.behaviorTrialOffset * (double(runLabel) - 1);
    else
        offsetAmount = opts.behaviorTrialOffset * (originalInputIndex - 1);
    end

    return;
end

end

function part = reshape_struct_row(part)

if isempty(part)
    part = struct([]);
    return;
end

if ~isstruct(part)
    error('nf_concat:BadMetadataField', 'Metadata fields event, epoch, and behavior must be structures when present.');
end

part = reshape(part, 1, []);

end

function part = offset_struct_field(part, fieldName, offsetAmount)

if isempty(part)
    return;
end

if ~isfield(part, fieldName)
    return;
end

for itemIndex = 1:numel(part)
    value = part(itemIndex).(fieldName);
    value = offset_numeric_or_cell(value, offsetAmount);
    part(itemIndex).(fieldName) = value;
end

end

function value = offset_numeric_or_cell(value, offsetAmount)

if isnumeric(value)
    value = value + offsetAmount;
    return;
end

if iscell(value)
    for cellIndex = 1:numel(value)
        if isnumeric(value{cellIndex})
            value{cellIndex} = value{cellIndex} + offsetAmount;
        end
    end

    return;
end

end

function out = concatenate_struct_parts(parts)

allFields = {};

for partIndex = 1:numel(parts)
    if isempty(parts{partIndex})
        continue;
    end

    thisFields = fieldnames(parts{partIndex});

    for fieldIndex = 1:numel(thisFields)
        if ~any(strcmp(allFields, thisFields{fieldIndex}))
            allFields{end + 1} = thisFields{fieldIndex};
        end
    end
end

if isempty(allFields)
    out = struct([]);
    return;
end

for partIndex = 1:numel(parts)
    parts{partIndex} = align_struct_fields(parts{partIndex}, allFields);
end

out = [parts{:}];

end

function part = align_struct_fields(part, allFields)

if isempty(part)
    values = cell(numel(allFields), 1);
    part = cell2struct(values, allFields, 1);
    part = part([]);
    return;
end

for fieldIndex = 1:numel(allFields)
    fieldName = allFields{fieldIndex};

    if ~isfield(part, fieldName)
        for itemIndex = 1:numel(part)
            part(itemIndex).(fieldName) = [];
        end
    end
end

part = orderfields(part, allFields);

end

function labels = get_indexed_labels(labelsIn, indices)

if isnumeric(labelsIn)
    labels = labelsIn(indices);
elseif isstring(labelsIn)
    labels = labelsIn(indices);
elseif iscell(labelsIn)
    labels = labelsIn(indices);
else
    labels = indices;
end

end

function report = make_merge_report(sources, keptSources, keptInputIndices, keptRunLabels, timeInfo, dataFields, nTrialsByRun, eventCountsByRun, opts)

report = struct();
report.function = mfilename;
report.created = datestr(now, 30);
report.nInput = numel(sources);
report.nKept = numel(keptInputIndices);
report.nDropped = report.nInput - report.nKept;
report.sources = sources;
report.keptSources = keptSources;
report.keptInputIndices = keptInputIndices;
report.droppedInputIndices = find(~timeInfo.keep);
report.dropReason = timeInfo.dropReason;
report.runLabels = opts.runLabels;
report.keptRunLabels = keptRunLabels;
report.dataFields = dataFields;
report.nTrialsByRun = nTrialsByRun;
report.nTrialsTotal = sum(nTrialsByRun);
report.eventCountsByRun = eventCountsByRun;
report.timeMode = opts.timeMode;
report.timePolicyNote = timeInfo.policyNote;
report.timeTemplateInputIndex = timeInfo.templateIndex;
report.timeTemplateCount = timeInfo.templateCount;
report.nTimeOriginal = timeInfo.nTimeOriginal;
report.nTimeOutput = timeInfo.nTimeOutput;
report.timeStart = timeInfo.times(1);
report.timeEnd = timeInfo.times(end);
report.settings = struct();
report.settings.timeTolerance = opts.timeTolerance;
report.settings.freqTolerance = opts.freqTolerance;
report.settings.chanMode = opts.chanMode;
report.settings.trialDim = opts.trialDim;
report.settings.timeDim = opts.timeDim;
report.settings.behaviorTrialMode = opts.behaviorTrialMode;
report.settings.behaviorTrialOffset = opts.behaviorTrialOffset;
report.settings.behaviorTrialFields = opts.behaviorTrialFields;
report.settings.offsetEventEpoch = opts.offsetEventEpoch;
report.settings.offsetEpochEvent = opts.offsetEpochEvent;

end

function fields = normalize_field_list(fields, optionName)

if isempty(fields)
    fields = {};
    return;
end

if ischar(fields)
    fields = {fields};
    return;
end

if isstring(fields)
    fields = cellstr(fields(:));
    fields = fields(:)';
    return;
end

if iscell(fields)
    for fieldIndex = 1:numel(fields)
        if isstring(fields{fieldIndex})
            fields{fieldIndex} = char(fields{fieldIndex});
        end

        if ~ischar(fields{fieldIndex})
            error('nf_concat:BadFieldList', '%s must contain only character vectors or strings.', optionName);
        end
    end

    fields = fields(:)';
    return;
end

error('nf_concat:BadFieldList', '%s must be a character vector, string, or cell array.', optionName);

end

function mode = normalize_mode_string(mode)

if isstring(mode)
    mode = char(mode);
end

if ~ischar(mode)
    error('nf_concat:BadModeString', 'Mode values must be character vectors or strings.');
end

mode = lower(strtrim(mode));

end
