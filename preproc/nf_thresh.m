function [EEG, rejected, info] = nf_thresh(EEG, voltageThreshold, ...
    powerThreshold, frequencyRange, times, interpolate, maxBadChannels, ...
    frontalChannels, varargin)
% NF_THRESH  Mark, locally repair, and reject artifactual EEG epochs.
%
% [EEG, REJECTED, INFO] = NF_THRESH(EEG, VOLTAGE, POWER, FREQUENCIES,
% TIMES, INTERPOLATE, MAXBAD, FRONTAL, ...)
%
% VOLTAGE is in microvolts. The defaults reproduce the NeuroFreq/MADE
% criteria: +/-125 uV and FFT power outside [-100 30] dB at 20-40 Hz.
% Epochs containing a flagged frontal channel or more than MAXBAD flagged
% channels are rejected. Sparse nonfrontal flags are locally interpolated
% when INTERPOLATE is true. This helper never performs global interpolation
% or rereferencing; those operations belong once at the pipeline exit.

nf_validate_epoched_eeg(EEG);

if nargin < 2 || isempty(voltageThreshold)
    voltageThreshold = 125;
end
if nargin < 3 || isempty(powerThreshold)
    powerThreshold = [-100 30];
end
if nargin < 4 || isempty(frequencyRange)
    frequencyRange = [20 40];
end
if nargin < 5 || isempty(times)
    times = [EEG.xmin EEG.xmax];
end
if nargin < 6 || isempty(interpolate)
    interpolate = true;
end
if nargin < 7 || isempty(maxBadChannels)
    maxBadChannels = floor(EEG.nbchan / 10);
end
if nargin < 8
    frontalChannels = {};
end

parser = inputParser;
parser.FunctionName = 'nf_thresh';
addParameter(parser, 'interpolationMethod', 'sphericalKang', @nf_is_text);
parse(parser, varargin{:});
options = parser.Results;
options.interpolationMethod = nf_normalize_interpolation_method( ...
    options.interpolationMethod);

nf_validate_settings(EEG, voltageThreshold, powerThreshold, ...
    frequencyRange, times, interpolate, maxBadChannels);
powerThreshold = reshape(powerThreshold, 1, 2);
frequencyRange = reshape(frequencyRange, 1, 2);
times = reshape(times, 1, 2);
if interpolate
    EEG = nf_normalize_locations(EEG);
    [frontalIndices, frontalLabels, frontalSource] = ...
        nf_resolve_frontal_channels(EEG, frontalChannels);
else
    frontalIndices = [];
    frontalLabels = {};
    frontalSource = 'not-used-without-local-interpolation';
end

originalTrials = EEG.trials;
originalEpochIds = 1:originalTrials;
behavior = [];
hasBehavior = isfield(EEG, 'etc') && isfield(EEG.etc, 'behavior');
if hasBehavior
    nf_validate_trial_metadata(EEG.etc.behavior, originalTrials, ...
        'EEG.etc.behavior');
end
if isfield(EEG, 'etc') && isfield(EEG.etc, 'nf_epoch_ids')
    nf_validate_trial_metadata(EEG.etc.nf_epoch_ids, originalTrials, ...
        'EEG.etc.nf_epoch_ids');
    originalEpochIds = EEG.etc.nf_epoch_ids;
end

EEG = nf_clear_threshold_marks(EEG);
EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -voltageThreshold, ...
    voltageThreshold, times(1), times(2), 0, 0);
EEG = eeg_checkset(EEG);
voltageFlags = nf_get_rejection_matrix(EEG, 'rejthreshE', ...
    EEG.nbchan, originalTrials);

EEG.specdata = [];
EEG = pop_rejspec(EEG, 1, 'elecrange', 1:EEG.nbchan, ...
    'method', 'fft', 'threshold', powerThreshold, ...
    'freqlimits', frequencyRange, 'specdata', [], ...
    'eegplotplotallrej', 0, 'eegplotreject', 0);
EEG = eeg_checkset(EEG);
spectralFlags = nf_get_rejection_matrix(EEG, 'rejfreqE', ...
    EEG.nbchan, originalTrials);
EEG.specdata = [];

badChannelEpoch = voltageFlags | spectralFlags;
frontalVoltageRejected = any(voltageFlags(frontalIndices, :), 1);
frontalSpectralRejected = any(spectralFlags(frontalIndices, :), 1);
tooManyChannelsRejected = sum(badChannelEpoch, 1) > maxBadChannels;

if interpolate
    rejected = frontalVoltageRejected | frontalSpectralRejected | ...
        tooManyChannelsRejected;
else
    rejected = any(badChannelEpoch, 1);
end

if all(rejected)
    error('nf_thresh:AllRejected', ...
        'All epochs failed the requested artifact criteria.');
end

localInterpolationMask = false(EEG.nbchan, originalTrials);
icaInvalidated = false;
for trialIndex = 1:originalTrials
    badChannels = find(badChannelEpoch(:, trialIndex));
    if isempty(badChannels) || rejected(trialIndex) || ~interpolate
        continue
    end
    singleEpoch = pop_select(EEG, 'trial', trialIndex, 'sorttrial', 'off');
    singleEpoch = eeg_interp(singleEpoch, badChannels, ...
        char(options.interpolationMethod));
    EEG.data(:, :, trialIndex) = singleEpoch.data;
    localInterpolationMask(badChannels, trialIndex) = true;
end
if any(localInterpolationMask(:))
    EEG = nf_clear_ica(EEG);
    icaInvalidated = true;
    EEG = eeg_checkset(EEG);
end

retainedIndices = find(~rejected);
if hasBehavior
    behavior = nf_subset_trial_metadata(EEG.etc.behavior, ...
        retainedIndices, originalTrials, 'EEG.etc.behavior');
end
retainedEpochIds = nf_subset_trial_metadata(originalEpochIds, ...
    retainedIndices, originalTrials, 'EEG.etc.nf_epoch_ids');

EEG = pop_rejepoch(EEG, rejected, 0);
EEG = eeg_checkset(EEG);
EEG = nf_clear_threshold_marks(EEG);
if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
if hasBehavior
    EEG.etc.behavior = behavior;
end
EEG.etc.nf_epoch_ids = retainedEpochIds;

rejectionReason = repmat({''}, 1, originalTrials);
for trialIndex = 1:originalTrials
    if ~rejected(trialIndex) && any(localInterpolationMask(:, trialIndex))
        rejectionReason{trialIndex} = 'retained-local-interpolation';
    elseif ~rejected(trialIndex)
        rejectionReason{trialIndex} = 'retained-clean';
    elseif frontalVoltageRejected(trialIndex) && ...
            frontalSpectralRejected(trialIndex)
        rejectionReason{trialIndex} = 'frontal-voltage-and-spectral';
    elseif frontalVoltageRejected(trialIndex)
        rejectionReason{trialIndex} = 'frontal-voltage';
    elseif frontalSpectralRejected(trialIndex)
        rejectionReason{trialIndex} = 'frontal-spectral';
    elseif tooManyChannelsRejected(trialIndex)
        rejectionReason{trialIndex} = 'too-many-channels';
    else
        rejectionReason{trialIndex} = 'artifact-no-local-interpolation';
    end
end

info = struct();
info.schemaVersion = '2.0.0';
info.units.voltage = 'microvolts';
info.units.spectralPower = 'dB';
info.criteria.voltageThreshold = voltageThreshold;
info.criteria.powerThreshold = powerThreshold;
info.criteria.frequencyRangeHz = frequencyRange;
info.criteria.timesSeconds = times;
info.criteria.spectralMethod = 'fft';
info.criteria.maxBadChannels = maxBadChannels;
info.criteria.interpolateSparseNonfrontal = logical(interpolate);
info.frontal.indices = frontalIndices;
info.frontal.labels = frontalLabels;
info.frontal.source = frontalSource;
info.masks.voltage = voltageFlags;
info.masks.spectral = spectralFlags;
info.masks.anyArtifact = badChannelEpoch;
info.masks.localInterpolation = localInterpolationMask;
info.masks.rejected = rejected;
info.masks.frontalVoltageRejected = frontalVoltageRejected;
info.masks.frontalSpectralRejected = frontalSpectralRejected;
info.masks.tooManyChannelsRejected = tooManyChannelsRejected;
info.rejectionReason = rejectionReason;
info.localBadCount = sum(badChannelEpoch, 1);
info.originalEpochIds = originalEpochIds;
info.retainedIndices = retainedIndices;
info.retainedEpochIds = retainedEpochIds;
info.nOriginal = originalTrials;
info.nRejected = sum(rejected);
info.nRetained = numel(retainedIndices);
info.nLocallyRepairedEpochs = sum(any(localInterpolationMask, 1));
info.nLocallyRepairedChannelEpochs = sum(localInterpolationMask(:));
info.interpolationMethod = char(options.interpolationMethod);
info.icaDecompositionInvalidated = icaInvalidated;
EEG.etc.nf_thresh = info;

end

function EEG = nf_clear_threshold_marks(EEG)
if ~isfield(EEG, 'reject') || isempty(EEG.reject)
    EEG.reject = struct();
end
fields = {'rejthresh', 'rejthreshE', 'rejfreq', 'rejfreqE', ...
    'rejglobal', 'rejglobalE'};
for index = 1:numel(fields)
    EEG.reject.(fields{index}) = [];
end
EEG.specdata = [];
end

function EEG = nf_clear_ica(EEG)
fields = {'icaact', 'icachansind', 'icasphere', 'icaweights', 'icawinv', ...
    'stats', 'specicaact', 'dipfit'};
for index = 1:numel(fields)
    if isfield(EEG, fields{index})
        EEG.(fields{index}) = [];
    end
end
if isfield(EEG, 'reject') && isstruct(EEG.reject)
    rejectionFields = fieldnames(EEG.reject);
    for index = 1:numel(rejectionFields)
        fieldName = rejectionFields{index};
        if strcmpi(fieldName, 'gcompreject') || ...
                strncmpi(fieldName, 'icarej', 6)
            EEG.reject.(fieldName) = [];
        end
    end
end
if isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification')
    EEG.etc = rmfield(EEG.etc, 'ic_classification');
end
if isfield(EEG, 'etc') && isfield(EEG.etc, 'amica')
    EEG.etc = rmfield(EEG.etc, 'amica');
end
if isfield(EEG, 'mods')
    EEG = rmfield(EEG, 'mods');
end
end

function matrix = nf_get_rejection_matrix(EEG, fieldName, channels, trials)
if ~isfield(EEG, 'reject') || ~isfield(EEG.reject, fieldName)
    error('nf_thresh:MissingRejectionFlags', ...
        'EEGLAB did not create EEG.reject.%s.', fieldName);
end
matrix = logical(EEG.reject.(fieldName));
if trials ~= channels && isequal(size(matrix), [trials channels])
    matrix = matrix';
end
if ~isequal(size(matrix), [channels trials])
    error('nf_thresh:InvalidRejectionFlags', ...
        'EEG.reject.%s has unexpected dimensions.', fieldName);
end
end

function [indices, labels, source] = nf_resolve_frontal_channels(EEG, requested)
allLabels = {EEG.chanlocs.labels};
for index = 1:numel(allLabels)
    if ~(ischar(allLabels{index}) || ...
            (isstring(allLabels{index}) && isscalar(allLabels{index}))) || ...
            isempty(strtrim(char(allLabels{index})))
        error('nf_thresh:InvalidChannelLabels', ...
            'Every channel must have a nonempty scalar-text label.');
    end
    allLabels{index} = char(allLabels{index});
end
normalizedLabels = cellfun(@lower, allLabels, 'UniformOutput', false);
if numel(unique(normalizedLabels)) ~= EEG.nbchan
    error('nf_thresh:DuplicateChannelLabels', ...
        'Channel labels must be unique, ignoring case.');
end

if ischar(requested)
    requested = {requested};
elseif isstring(requested)
    requested = cellstr(requested);
end

if ~isempty(requested)
    requested = requested(:)';
    if ~( iscellstr(requested) || istring(requested) ) %#ok<ISCLSTR>
        error('nf_thresh:InvalidFrontalChannels', ...
            'frontalChannels must contain scalar-text channel labels.');
    end
    requestedNormalized = cellfun(@lower, requested, 'UniformOutput', false);
    [present, indices] = ismember(requestedNormalized, normalizedLabels);
    if any(~present)
        missing = strjoin(requested(~present), ', ');
        error('nf_thresh:MissingFrontalChannels', ...
            'Requested frontal channels were not found: %s', missing);
    end
    labels = allLabels(indices);
    source = 'requested-labels';
    return
end

standard = {'Fp1', 'Fpz', 'Fp2', 'AF7', 'AF3', 'AFz', 'AF4', 'AF8'};
standardNormalized = cellfun(@lower, standard, 'UniformOutput', false);
[present, standardIndices] = ismember(standardNormalized, normalizedLabels);
indices = standardIndices(present);
if ~isempty(indices)
    labels = allLabels(indices);
    source = 'standard-frontopolar-labels';
    return
end

x = nan(1, EEG.nbchan);
for channelIndex = 1:EEG.nbchan
    if isfield(EEG.chanlocs, 'X') && ...
            isnumeric(EEG.chanlocs(channelIndex).X) && ...
            isscalar(EEG.chanlocs(channelIndex).X) && ...
            isfinite(EEG.chanlocs(channelIndex).X)
        x(channelIndex) = EEG.chanlocs(channelIndex).X;
    end
end
validIndices = find(isfinite(x));
minimumCoordinates = max(4, ceil(0.8 * EEG.nbchan));
if numel(validIndices) < minimumCoordinates
    error('nf_thresh:UnknownFrontalChannels', ...
        ['No standard frontopolar labels were found and too few valid X ' ...
        'coordinates are available. Supply frontalChannels explicitly.']);
end
count = max(1, ceil(0.1 * numel(validIndices)));
[~, order] = sort(x(validIndices), 'descend');
indices = validIndices(order(1:count));
labels = allLabels(indices);
source = 'largest-positive-X-coordinates';
end

function nf_validate_settings(EEG, voltageThreshold, powerThreshold, ...
    frequencyRange, times, interpolate, maxBadChannels)
if ~nf_is_positive_scalar(voltageThreshold)
    error('nf_thresh:InvalidVoltageThreshold', ...
        'voltageThreshold must be a positive scalar in microvolts.');
end
if ~nf_is_increasing_pair(powerThreshold)
    error('nf_thresh:InvalidPowerThreshold', ...
        'powerThreshold must be a finite increasing pair in dB.');
end
if ~nf_is_increasing_pair(frequencyRange) || frequencyRange(1) < 0 || ...
        frequencyRange(2) >= EEG.srate / 2
    error('nf_thresh:InvalidFrequencyRange', ...
        'frequencyRange must be nonnegative and below Nyquist.');
end
if ~nf_is_increasing_pair(times) || times(1) < EEG.xmin || ...
        times(2) > EEG.xmax
    error('nf_thresh:InvalidTimes', ...
        'times must be a finite increasing pair inside the epoch.');
end
if ~nf_is_logical_scalar(interpolate)
    error('nf_thresh:InvalidInterpolationFlag', ...
        'interpolate must be a logical scalar.');
end
if ~nf_is_nonnegative_integer(maxBadChannels) || ...
        maxBadChannels >= EEG.nbchan
    error('nf_thresh:InvalidMaximum', ...
        'maxBadChannels must be a nonnegative integer below EEG.nbchan.');
end
if interpolate && maxBadChannels > EEG.nbchan - 3
    error('nf_thresh:InvalidMaximum', ...
        ['maxBadChannels must leave at least three donor channels for ' ...
        'spherical interpolation.']);
end
end

function nf_validate_trial_metadata(value, trials, fieldName)
if istable(value)
    count = height(value);
elseif isvector(value)
    count = numel(value);
else
    count = size(value, 1);
end
if count ~= trials
    error('nf_thresh:MetadataMismatch', ...
        '%s has %d trial entries but EEG has %d trials.', ...
        fieldName, count, trials);
end
end

function value = nf_subset_trial_metadata(value, indices, trials, fieldName)
nf_validate_trial_metadata(value, trials, fieldName);
if istable(value)
    value = value(indices, :);
elseif isvector(value)
    value = value(indices);
else
    dimensions = repmat({':'}, 1, ndims(value));
    dimensions{1} = indices;
    value = value(dimensions{:});
end
end

function nf_validate_epoched_eeg(EEG)
if ~isstruct(EEG) || numel(EEG) ~= 1 || ~isfield(EEG, 'data') || ...
        ~isfield(EEG, 'srate') || ~isfield(EEG, 'nbchan') || ...
        ~isfield(EEG, 'pnts') || ~isfield(EEG, 'trials') || ...
        ~isfield(EEG, 'xmin') || ~isfield(EEG, 'xmax') || ...
        ~isfield(EEG, 'chanlocs')
    error('nf_thresh:InvalidEEG', ...
        'EEG must be one valid EEGLAB dataset structure.');
end
if ~isnumeric(EEG.nbchan) || ~isscalar(EEG.nbchan) || ...
        ~isfinite(EEG.nbchan) || EEG.nbchan < 3 || ...
        EEG.nbchan ~= round(EEG.nbchan) || ...
        ~isnumeric(EEG.pnts) || ~isscalar(EEG.pnts) || ...
        ~isfinite(EEG.pnts) || EEG.pnts < 2 || EEG.pnts ~= round(EEG.pnts) || ...
        ~isnumeric(EEG.trials) || ~isscalar(EEG.trials) || ...
        ~isfinite(EEG.trials) || EEG.trials < 2 || ...
        EEG.trials ~= round(EEG.trials) || ...
        ~isnumeric(EEG.srate) || ~isscalar(EEG.srate) || ...
        ~isfinite(EEG.srate) || EEG.srate <= 0 || ...
        ~isnumeric(EEG.xmin) || ~isscalar(EEG.xmin) || ...
        ~isfinite(EEG.xmin) || ~isnumeric(EEG.xmax) || ...
        ~isscalar(EEG.xmax) || ~isfinite(EEG.xmax) || ...
        EEG.xmin >= EEG.xmax || ...
        ~isnumeric(EEG.data) || ~isreal(EEG.data) || isempty(EEG.data) || ...
        size(EEG.data, 3) ~= EEG.trials || ...
        size(EEG.data, 1) ~= EEG.nbchan || size(EEG.data, 2) ~= EEG.pnts || ...
        numel(EEG.chanlocs) ~= EEG.nbchan || any(~isfinite(EEG.data(:)))
    error('nf_thresh:ContinuousData', ...
        'nf_thresh requires a consistent epoched EEGLAB dataset.');
end
end

function EEG = nf_normalize_locations(EEG)
if exist('convertlocs', 'file') ~= 2
    error('nf_thresh:MissingLocationConverter', ...
        'EEGLAB convertlocs.m is required for local interpolation.');
end
hasCartesian = true;
hasTopographic = true;
for index = 1:numel(EEG.chanlocs)
    hasCartesian = hasCartesian && nf_has_xyz(EEG.chanlocs(index));
    hasTopographic = hasTopographic && ...
        nf_has_topography(EEG.chanlocs(index));
end
if hasCartesian
    EEG.chanlocs = convertlocs(EEG.chanlocs, 'cart2all');
elseif hasTopographic
    EEG.chanlocs = convertlocs(EEG.chanlocs, 'topo2all');
else
    error('nf_thresh:MissingChannelCoordinates', ...
        ['Every interpolated channel needs finite nonzero X/Y/Z or finite ' ...
        'theta/radius coordinates.']);
end
nf_validate_locations(EEG.chanlocs);
end

function valid = nf_has_xyz(chanloc)
valid = isfield(chanloc, 'X') && isfield(chanloc, 'Y') && ...
    isfield(chanloc, 'Z');
if ~valid
    return
end
coordinates = [chanloc.X chanloc.Y chanloc.Z];
valid = isnumeric(coordinates) && isreal(coordinates) && ...
    numel(coordinates) == 3 && all(isfinite(coordinates)) && ...
    norm(double(coordinates)) > 0;
end

function valid = nf_has_topography(chanloc)
valid = isfield(chanloc, 'theta') && isfield(chanloc, 'radius') && ...
    isnumeric(chanloc.theta) && isreal(chanloc.theta) && ...
    isscalar(chanloc.theta) && isfinite(chanloc.theta) && ...
    isnumeric(chanloc.radius) && isreal(chanloc.radius) && ...
    isscalar(chanloc.radius) && isfinite(chanloc.radius) && ...
    chanloc.radius >= 0;
end

function nf_validate_locations(chanlocs)
coordinates = zeros(numel(chanlocs), 3);
for index = 1:numel(chanlocs)
    if ~nf_has_xyz(chanlocs(index)) || ...
            ~nf_has_topography(chanlocs(index))
        error('nf_thresh:MissingChannelCoordinates', ...
            'Channel %s lacks complete interpolation geometry.', ...
            char(chanlocs(index).labels));
    end
    coordinates(index, :) = double([chanlocs(index).X ...
        chanlocs(index).Y chanlocs(index).Z]);
    coordinates(index, :) = coordinates(index, :) ./ ...
        norm(coordinates(index, :));
end
for first = 1:size(coordinates, 1) - 1
    distances = sqrt(sum((coordinates(first + 1:end, :) - ...
        coordinates(first, :)) .^ 2, 2));
    if any(distances < 1e-10)
        second = first + find(distances < 1e-10, 1);
        error('nf_thresh:DuplicateChannelCoordinates', ...
            'Channels %s and %s occupy the same scalp coordinate.', ...
            char(chanlocs(first).labels), char(chanlocs(second).labels));
    end
end
end

function method = nf_normalize_interpolation_method(value)
normalized = lower(regexprep(strtrim(char(value)), '[ _-]', ''));
switch normalized
    case 'spherical'
        method = 'spherical';
    case 'sphericalkang'
        method = 'sphericalKang';
    case 'sphericalcrd'
        method = 'sphericalCRD';
    case 'invdist'
        method = 'invdist';
    case 'v4'
        method = 'v4';
    otherwise
        error('nf_thresh:UnknownInterpolationMethod', ...
            ['interpolationMethod must be spherical, sphericalKang, ' ...
            'sphericalCRD, invdist, or v4.']);
end
end

function valid = nf_is_positive_scalar(value)
valid = isnumeric(value) && isscalar(value) && isfinite(value) && value > 0;
end

function valid = nf_is_increasing_pair(value)
valid = isnumeric(value) && numel(value) == 2 && all(isfinite(value)) && ...
    value(1) < value(2);
end

function valid = nf_is_nonnegative_integer(value)
valid = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value >= 0 && value == round(value);
end

function valid = nf_is_logical_scalar(value)
valid = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
    ismember(value, [0 1]);
end

function valid = nf_is_text(value)
valid = ischar(value) || (isstring(value) && isscalar(value));
end
