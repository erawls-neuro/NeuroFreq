function [EEG, info] = nf_badchans(EEG, maxBad, interpolate, method, varargin)
% NF_BADCHANS  Detect and remove bad scalp channels before spatial cleaning.
%
% [EEG, INFO] = NF_BADCHANS(EEG, MAXBAD, INTERPOLATE, METHOD, ...)
%
% METHOD is 'faster' (default) or 'cleanrawdata'. The FASTER path calls
% the published channel_properties and min_z functions. MAXBAD is a QC
% ceiling: the function errors when a detector exceeds it; the detector is
% never relaxed to force a dataset through the ceiling.
% Standard mean/SD z-scoring is enforced across compatible min_z variants.
% For average-referenced data without a recoverable recording reference,
% the intended no-reference channel properties are calculated locally to
% avoid the undefined dist_inds branch in official FASTER 1.2.4.
%
% Name/value inputs:
%   reference             [] for automatic, a channel label, or a channel index
%   fasterOptions         overrides for the standard FASTER min_z fields:
%                             measure=[1 1 1] and z=[3 3 3]
%   cleanCorrelation      clean_channels correlation threshold (default 0.8)
%   cleanHighpass         diagnostic clean_channels high-pass (default 1 Hz)
%   interpolationMethod   'sphericalKang' (default)
%
% In a preprocessing pipeline, pass INTERPOLATE=false. Globally removed
% channels should remain absent through GEDAI and ICA, then be interpolated
% once at the end. INTERPOLATE=true is retained for direct-call backwards
% compatibility.

nf_validate_eeg(EEG);

if nargin < 2 || isempty(maxBad)
    maxBad = floor(EEG.nbchan / 5);
end
if nargin < 3 || isempty(interpolate)
    interpolate = true;
end
if nargin < 4 || isempty(method)
    method = 'faster';
end

parser = inputParser;
parser.FunctionName = 'nf_badchans';
addParameter(parser, 'reference', [], @nf_is_reference);
addParameter(parser, 'fasterOptions', struct(), @nf_is_faster_options);
addParameter(parser, 'cleanCorrelation', 0.8, @nf_is_correlation);
addParameter(parser, 'cleanHighpass', 1, @nf_is_positive_scalar);
addParameter(parser, 'interpolationMethod', 'sphericalKang', @nf_is_text);
parse(parser, varargin{:});
options = parser.Results;
options.fasterOptions = nf_normalize_faster_options(options.fasterOptions);
options.interpolationMethod = nf_normalize_interpolation_method( ...
    options.interpolationMethod);

if ~nf_is_nonnegative_integer(maxBad)
    error('nf_badchans:InvalidMaximum', 'maxBad must be a nonnegative integer.');
end
if ~nf_is_logical_scalar(interpolate)
    error('nf_badchans:InvalidInterpolationFlag', 'interpolate must be a logical scalar.');
end

method = lower(char(method));
if ismember(method, {'faster', 'cleanrawdata'}) || logical(interpolate)
    EEG = nf_normalize_locations(EEG);
end
originalLocations = EEG.chanlocs;
originalLabels = nf_channel_labels(originalLocations);
nf_require_unique_labels(originalLabels);

if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
if ~isfield(EEG.etc, 'ogchan') || isempty(EEG.etc.ogchan)
    EEG.etc.ogchan = originalLocations;
end

[flatMask, flatDetails] = nf_flat_channel_mask(EEG);
nonfiniteMask = nf_nonfinite_channel_mask(EEG);
[referenceIndex, referenceMode, referenceIsRemoved, referenceSource] = ...
    nf_reference_channel( ...
    EEG, options.reference, flatMask, nonfiniteMask);

artifactFlatMask = flatMask;
if referenceIsRemoved
    artifactFlatMask(referenceIndex) = false;
end
mandatoryMask = artifactFlatMask | nonfiniteMask;

if strcmp(method, 'faster')
    [detectorMask, detector] = nf_run_faster(EEG, mandatoryMask, ...
        referenceIndex, referenceMode, options.fasterOptions);
elseif strcmp(method, 'cleanrawdata')
    [detectorMask, detector] = nf_run_cleanraw(EEG, mandatoryMask, ...
        referenceIndex, referenceIsRemoved, options.cleanCorrelation, ...
        options.cleanHighpass);
else
    error('nf_badchans:UnknownMethod', ...
        'method must be ''faster'' or ''cleanrawdata''.');
end

detectedArtifactMask = mandatoryMask | detectorMask;
detectedArtifactIndices = find(detectedArtifactMask);
if numel(detectedArtifactIndices) > maxBad
    badList = strjoin(originalLabels(detectedArtifactIndices), ', ');
    error('nf_badchans:TooManyBadChannels', ...
        ['The detector identified %d artifact channels, above maxBad=%d. ' ...
        'The dataset was not modified. Channels: %s'], ...
        numel(detectedArtifactIndices), maxBad, badList);
end

referenceRemovalMask = false(1, EEG.nbchan);
if referenceIsRemoved
    referenceRemovalMask(referenceIndex) = true;
end
removedMask = detectedArtifactMask | referenceRemovalMask;
removedIndices = find(removedMask);

if EEG.nbchan - numel(removedIndices) < 3
    error('nf_badchans:TooFewChannelsRemain', ...
        'Fewer than three channels would remain after channel removal.');
end

artifactLabels = originalLabels(detectedArtifactIndices);
removedLabels = originalLabels(removedIndices);
if ~isempty(removedIndices)
    EEG = pop_select(EEG, 'nochannel', removedIndices);
end
EEG = eeg_checkset(EEG);

interpolated = false;
referenceRestored = false;
orderRestored = false;
if interpolate && ~isempty(removedIndices)
    if ~isempty(detectedArtifactIndices)
        interpolationLocations = originalLocations(~referenceRemovalMask);
        EEG = eeg_interp(EEG, interpolationLocations, ...
            options.interpolationMethod);
        interpolated = true;
    end
    if referenceIsRemoved
        EEG = pop_reref(EEG, [], ...
            'refloc', originalLocations(referenceIndex), 'refica', 'remove');
        referenceRestored = true;
    else
        EEG = pop_reref(EEG, [], 'refica', 'remove');
    end
    EEG = eeg_checkset(EEG);
    [EEG, orderRestored] = nf_restore_channel_order(EEG, originalLocations);
end

info = struct();
info.schemaVersion = '2.0.0';
info.method = method;
info.nOriginal = numel(originalLabels);
info.maxBad = maxBad;
info.reference.mode = referenceMode;
info.reference.index = referenceIndex;
info.reference.label = nf_label_at(originalLabels, referenceIndex);
info.reference.removedZeroReference = referenceIsRemoved;
info.reference.restored = referenceRestored;
info.reference.source = referenceSource;
info.flat.indices = find(flatMask);
info.flat.labels = originalLabels(flatMask);
info.flat.standardDeviation = flatDetails.standardDeviation;
info.flat.range = flatDetails.range;
info.flat.tolerance = flatDetails.tolerance;
info.nonfinite.indices = find(nonfiniteMask);
info.nonfinite.labels = originalLabels(nonfiniteMask);
info.detector = detector;
info.artifact.indices = detectedArtifactIndices;
info.artifact.labels = artifactLabels;
info.artifact.nDetected = numel(detectedArtifactIndices);
info.removed.indices = removedIndices;
info.removed.labels = removedLabels;
info.removed.nRemoved = numel(removedIndices);
info.interpolation.requested = logical(interpolate);
info.interpolation.performed = interpolated;
info.interpolation.method = char(options.interpolationMethod);
info.interpolation.originalOrderRestored = orderRestored;
info.labelsBefore = originalLabels;
info.labelsAfter = {EEG.chanlocs.labels};

EEG.etc.badchans = numel(detectedArtifactIndices);
EEG.etc.badchanindices = detectedArtifactIndices;
EEG.etc.badchanlabels = artifactLabels;
EEG.etc.nf_removed_reference = struct('removed', referenceIsRemoved, ...
    'index', referenceIndex, 'label', nf_label_at(originalLabels, referenceIndex));
EEG.etc.nf_badchans = info;

end

function [EEG, restored] = nf_restore_channel_order(EEG, originalLocations)
expected = lower(string({originalLocations.labels}));
current = lower(string({EEG.chanlocs.labels}));
[present, order] = ismember(expected, current);
if any(~present) || numel(current) ~= numel(expected)
    error('nf_badchans:MontageRestorationFailed', ...
        'The interpolated channel set does not match the original montage.');
end
restored = ~isequal(order, 1:numel(order));
if ~restored
    return
end
EEG.data = EEG.data(order, :, :);
EEG.chanlocs = EEG.chanlocs(order);
if isfield(EEG, 'reject') && isstruct(EEG.reject)
    rejectionFields = {'rejmanualE', 'rejthreshE', 'rejconstE', ...
        'rejfreqE', 'rejjpE', 'rejkurtE', 'rejglobalE'};
    for index = 1:numel(rejectionFields)
        fieldName = rejectionFields{index};
        if ~isfield(EEG.reject, fieldName)
            continue
        end
        value = EEG.reject.(fieldName);
        if isnumeric(value) || islogical(value)
            if size(value, 1) == numel(order)
                EEG.reject.(fieldName) = value(order, :);
            elseif size(value, 2) == numel(order)
                EEG.reject.(fieldName) = value(:, order);
            end
        end
    end
end
EEG = eeg_checkset(EEG);
end

function [detectorMask, detector] = nf_run_faster(EEG, mandatoryMask, ...
    referenceIndex, referenceMode, fasterOptions)
dependencies = {'channel_properties', 'min_z', 'distancematrix', ...
    'hurst_exponent', 'nanmean'};
for dependencyIndex = 1:numel(dependencies)
    if ~nf_function_available(dependencies{dependencyIndex})
        error('nf_badchans:MissingFASTER', ...
            ['FASTER dependency %s was not found. Install the complete ' ...
            'FASTER distribution and its MATLAB dependencies.'], ...
            dependencies{dependencyIndex});
    end
end

usableIndices = find(~mandatoryMask);
if numel(usableIndices) < 3
    error('nf_badchans:TooFewUsableChannels', ...
        'FASTER requires at least three finite, non-flat channels.');
end

if strcmp(referenceMode, 'zero-reference')
    if ~ismember(referenceIndex, usableIndices)
        usableIndices = sort([usableIndices referenceIndex]);
    end
    diagnostic = pop_select(EEG, 'channel', usableIndices);
    diagnosticReference = find(usableIndices == referenceIndex, 1);
elseif strcmp(referenceMode, 'named-reference')
    diagnostic = pop_select(EEG, 'channel', usableIndices);
    diagnosticReference = find(usableIndices == referenceIndex, 1);
    diagnostic = pop_reref(diagnostic, diagnosticReference, 'keepref', 'on');
elseif strcmp(referenceMode, 'average-fallback')
    diagnostic = pop_select(EEG, 'channel', usableIndices);
    diagnostic = pop_reref(diagnostic, []);
    diagnosticReference = [];
else
    error('nf_badchans:InternalReferenceError', ...
        'Unrecognized FASTER reference mode.');
end

diagnostic = eeg_checkset(diagnostic);
if isempty(diagnosticReference)
    listProperties = ...
        nf_faster_no_reference_properties(diagnostic);
    propertyImplementation = ...
        'nf_badchans intended FASTER no-reference implementation';
else
    listProperties = channel_properties( ...
        diagnostic, ...
        1:diagnostic.nbchan, ...
        diagnosticReference);
    propertyImplementation = which('channel_properties');
end
minZOptions = fasterOptions;
minZOptions.stat = 'z';
fasterMaskLocal = logical(min_z(listProperties, minZOptions));
fasterMaskLocal = fasterMaskLocal(:)';
if numel(fasterMaskLocal) ~= diagnostic.nbchan
    error('nf_badchans:InvalidFASTEROutput', ...
        'min_z returned an unexpected number of channel decisions.');
end
if ~isempty(diagnosticReference)
    fasterMaskLocal(diagnosticReference) = false;
end

detectorMask = false(1, EEG.nbchan);
detectorMask(usableIndices(fasterMaskLocal)) = true;

detector = struct();
detector.name = 'FASTER';
detector.channelProperties = listProperties;
detector.propertyNames = {'meanCorrelation', 'variance', 'hurstExponent'};
detector.analyzedOriginalIndices = usableIndices;
detector.referenceIndexInDiagnostic = diagnosticReference;
detector.fasterMaskInDiagnostic = fasterMaskLocal;
detector.options = fasterOptions;
detector.backendOptions = minZOptions;
detector.implementation.channelProperties = propertyImplementation;
detector.implementation.minZ = which('min_z');
end

function listProperties = nf_faster_no_reference_properties(EEG)
channelData = reshape(double(EEG.data), EEG.nbchan, []);
zeroMask = false(EEG.nbchan, 1);
for channelIndex = 1:EEG.nbchan
    channel = channelData(channelIndex, :);
    zeroMask(channelIndex) = max(channel) == 0 && min(channel) == 0;
end
usableIndices = find(~zeroMask);
if numel(usableIndices) < 2
    error('nf_badchans:InsufficientFASTERCorrelationChannels', ...
        ['FASTER requires at least two nonzero channels to calculate ' ...
        'interchannel correlations.']);
end

correlationMatrix = abs(corrcoef(channelData(usableIndices, :)'));
meanCorrelation = nan(EEG.nbchan, 1);
meanCorrelation(usableIndices) = mean(correlationMatrix, 2);
if any(zeroMask)
    meanCorrelation(zeroMask) = mean(meanCorrelation(usableIndices));
end

channelVariance = var(channelData, 0, 2);
finiteVariance = isfinite(channelVariance);
if any(~finiteVariance)
    if ~any(finiteVariance)
        error('nf_badchans:InvalidFASTERVariance', ...
            'FASTER could not calculate a finite channel variance.');
    end
    channelVariance(~finiteVariance) = ...
        mean(channelVariance(finiteVariance));
end
hurstExponent = nan(EEG.nbchan, 1);
for channelIndex = 1:EEG.nbchan
    hurstExponent(channelIndex) = ...
        hurst_exponent(channelData(channelIndex, :));
end

listProperties = [meanCorrelation channelVariance hurstExponent];
for propertyIndex = 1:size(listProperties, 2)
    property = listProperties(:, propertyIndex);
    missing = isnan(property);
    if any(missing)
        available = property(~missing);
        if isempty(available)
            error('nf_badchans:InvalidFASTERProperty', ...
                'FASTER property %d contains no valid values.', ...
                propertyIndex);
        end
        property(missing) = mean(available);
    end
    property = property - median(property);
    listProperties(:, propertyIndex) = property;
end
end

function [detectorMask, detector] = nf_run_cleanraw(EEG, mandatoryMask, ...
    referenceIndex, referenceIsRemoved, correlationThreshold, highpass)
if exist('clean_channels', 'file') ~= 2
    error('nf_badchans:MissingCleanRawData', ...
        'clean_channels.m was not found. Install clean_rawdata.');
end
if exist('pop_eegfiltnew', 'file') ~= 2
    error('nf_badchans:MissingFilter', ...
        'pop_eegfiltnew.m is required for clean_channels diagnostic preparation.');
end
if EEG.trials ~= 1
    error('nf_badchans:EpochedCleanRawData', ...
        'cleanrawdata channel detection requires continuous EEG.');
end
if EEG.pnts <= round(5 * EEG.srate)
    error('nf_badchans:InsufficientCleanRawData', ...
        ['cleanrawdata channel detection requires more than five seconds ' ...
        'of continuous data.']);
end
if highpass >= EEG.srate / 2
    error('nf_badchans:InvalidCleanHighpass', ...
        'cleanHighpass must be below the data Nyquist frequency.');
end

excludedMask = mandatoryMask;
if referenceIsRemoved
    excludedMask(referenceIndex) = true;
end
usableIndices = find(~excludedMask);
if numel(usableIndices) < 3
    error('nf_badchans:TooFewUsableChannels', ...
        'clean_channels requires at least three usable channels.');
end

diagnostic = pop_select(EEG, 'channel', usableIndices);
diagnostic = pop_eegfiltnew(diagnostic, 'locutoff', highpass);
diagnostic = eeg_checkset(diagnostic);
previousRandomState = rng;
randomCleanup = onCleanup(@() rng(previousRandomState));
[~, removedLocal] = clean_channels(diagnostic, correlationThreshold);
clear randomCleanup
removedLocal = nf_cleanraw_removed_mask(removedLocal, numel(usableIndices));

detectorMask = false(1, EEG.nbchan);
detectorMask(usableIndices(removedLocal)) = true;

detector = struct();
detector.name = 'clean_rawdata clean_channels';
detector.correlationThreshold = correlationThreshold;
detector.diagnosticHighpassHz = highpass;
detector.minimumDurationSeconds = 5;
detector.analyzedOriginalIndices = usableIndices;
detector.removedMaskInDiagnostic = removedLocal;
end

function mask = nf_cleanraw_removed_mask(value, channelCount)
if isempty(value)
    mask = false(1, channelCount);
    return
end
if islogical(value)
    mask = value(:)';
    if numel(mask) ~= channelCount
        error('nf_badchans:InvalidCleanRawOutput', ...
            'clean_channels returned a logical mask with unexpected length.');
    end
    return
end

if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
        any(~isfinite(value(:))) || any(value(:) ~= round(value(:))) || ...
        any(value(:) < 1) || any(value(:) > channelCount) || ...
        numel(unique(value(:))) ~= numel(value)
    error('nf_badchans:InvalidCleanRawOutput', ...
        ['clean_channels must return either a channel-length logical mask ' ...
        'or unique valid removed-channel indices.']);
end
mask = false(1, channelCount);
mask(value(:)') = true;
end

function [referenceIndex, referenceMode, referenceIsRemoved, source] = ...
    nf_reference_channel(EEG, requestedReference, flatMask, nonfiniteMask)
referenceIndex = [];
referenceMode = '';
referenceIsRemoved = false;
source = '';

if ~isempty(requestedReference)
    referenceIndex = nf_resolve_channel(EEG, requestedReference);
    if nonfiniteMask(referenceIndex)
        error('nf_badchans:InvalidReference', ...
            'The requested reference channel contains nonfinite samples.');
    end
    if flatMask(referenceIndex)
        referenceMode = 'zero-reference';
        referenceIsRemoved = true;
        source = 'explicit-reference-option';
    else
        usableIndices = find(~flatMask & ~nonfiniteMask);
        if ~nf_all_have_xyz(EEG.chanlocs(usableIndices))
            error('nf_badchans:MissingReferenceCoordinates', ...
                ['A named FASTER reference requires finite X, Y, and Z ' ...
                'coordinates for every analyzed channel.']);
        end
        referenceMode = 'named-reference';
        source = 'explicit-reference-option';
    end
    return
end

[metadataIndex, metadataSource] = nf_metadata_reference(EEG);
if ~isempty(metadataIndex) && ~nonfiniteMask(metadataIndex) && ...
        flatMask(metadataIndex)
    referenceIndex = metadataIndex;
    referenceMode = 'zero-reference';
    referenceIsRemoved = true;
    source = metadataSource;
    return
end

labels = {EEG.chanlocs.labels};
reliableLabelMask = false(1, EEG.nbchan);
for index = 1:EEG.nbchan
    reliableLabelMask(index) = nf_is_reference_label(labels{index});
end
reliableZeroIndices = find(flatMask & ~nonfiniteMask & reliableLabelMask);
if numel(reliableZeroIndices) == 1
    referenceIndex = reliableZeroIndices;
    referenceMode = 'zero-reference';
    referenceIsRemoved = true;
    source = 'reference-channel-label';
    return
end

if ~isempty(metadataIndex) && ~flatMask(metadataIndex) && ...
        ~nonfiniteMask(metadataIndex)
    referenceIndex = metadataIndex;
    usableIndices = find(~flatMask & ~nonfiniteMask);
else
    usableIndices = [];
end
if ~isempty(referenceIndex) && nf_all_have_xyz(EEG.chanlocs(usableIndices))
    referenceMode = 'named-reference';
    source = metadataSource;
else
    referenceIndex = [];
    referenceMode = 'average-fallback';
    source = 'diagnostic-average-reference';
end
end

function [index, source] = nf_metadata_reference(EEG)
index = [];
source = '';
if ~isfield(EEG, 'ref') || ...
        ~(ischar(EEG.ref) || (isstring(EEG.ref) && isscalar(EEG.ref)))
    return
end
referenceLabel = strtrim(char(EEG.ref));
if isempty(referenceLabel) || ismember(lower(referenceLabel), ...
        {'common', 'average', 'averef', 'common average', 'unknown', 'none'})
    return
end
labels = {EEG.chanlocs.labels};
matches = find(strcmpi(labels, referenceLabel));
if numel(matches) == 1
    index = matches;
    source = 'EEG.ref-metadata';
end
end

function valid = nf_is_reference_label(label)
normalized = lower(regexprep(strtrim(char(label)), '[^a-z0-9]', ''));
valid = ismember(normalized, {'ref', 'reference', 'onlineref', ...
    'onlinereference', 'recordingref', 'recordingreference'});
end

function index = nf_resolve_channel(EEG, reference)
if isnumeric(reference)
    index = reference;
else
    label = char(reference);
    labels = {EEG.chanlocs.labels};
    index = find(strcmpi(labels, label));
end

if numel(index) ~= 1 || index < 1 || index > EEG.nbchan || index ~= round(index)
    error('nf_badchans:InvalidReference', ...
        'reference must resolve to exactly one valid channel.');
end
end

function [mask, details] = nf_flat_channel_mask(EEG)
data = reshape(EEG.data, EEG.nbchan, []);
standardDeviation = nan(1, EEG.nbchan);
channelRange = nan(1, EEG.nbchan);
for index = 1:EEG.nbchan
    channelData = double(data(index, :));
    if any(~isfinite(channelData))
        continue
    end
    standardDeviation(index) = std(channelData, 0, 2);
    channelRange(index) = max(channelData) - min(channelData);
end
usableScale = standardDeviation(isfinite(standardDeviation) & ...
    standardDeviation > 0);
if isempty(usableScale)
    typicalScale = 1;
else
    typicalScale = median(usableScale);
end
tolerance = max(1e-12, 1e-6 * typicalScale);
mask = isfinite(standardDeviation) & ...
    (standardDeviation <= tolerance | channelRange <= tolerance);
details = struct();
details.standardDeviation = standardDeviation;
details.range = channelRange;
details.tolerance = tolerance;
end

function mask = nf_nonfinite_channel_mask(EEG)
data = reshape(EEG.data, EEG.nbchan, []);
mask = any(~isfinite(data), 2)';
end

function labels = nf_channel_labels(chanlocs)
labels = cell(1, numel(chanlocs));
for index = 1:numel(chanlocs)
    if ~isfield(chanlocs, 'labels') || ...
            ~(ischar(chanlocs(index).labels) || ...
            (isstring(chanlocs(index).labels) && ...
            isscalar(chanlocs(index).labels))) || ...
            isempty(strtrim(char(chanlocs(index).labels)))
        error('nf_badchans:MissingChannelLabel', ...
            'Every channel must have a nonempty label.');
    end
    labels{index} = char(chanlocs(index).labels);
end
end

function nf_require_unique_labels(labels)
normalized = cellfun(@lower, labels, 'UniformOutput', false);
if numel(unique(normalized)) ~= numel(normalized)
    error('nf_badchans:DuplicateChannelLabels', ...
        'Channel labels must be unique, ignoring case.');
end
end

function value = nf_label_at(labels, index)
if isempty(index)
    value = '';
else
    value = labels{index};
end
end

function valid = nf_has_xyz(chanloc)
valid = false;
if ~isfield(chanloc, 'X') || ~isfield(chanloc, 'Y') || ~isfield(chanloc, 'Z')
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

function valid = nf_all_have_xyz(chanlocs)
valid = ~isempty(chanlocs);
for index = 1:numel(chanlocs)
    if ~nf_has_xyz(chanlocs(index))
        valid = false;
        return
    end
end
end

function EEG = nf_normalize_locations(EEG)
if ~nf_function_available('convertlocs')
    error('nf_badchans:MissingLocationConverter', ...
        'EEGLAB convertlocs.m is required for FASTER and interpolation.');
end
hasCartesian = nf_all_have_xyz(EEG.chanlocs);
hasTopographic = true;
for index = 1:numel(EEG.chanlocs)
    hasTopographic = hasTopographic && ...
        nf_has_topography(EEG.chanlocs(index));
end
if hasCartesian
    EEG.chanlocs = convertlocs(EEG.chanlocs, 'cart2all');
elseif hasTopographic
    EEG.chanlocs = convertlocs(EEG.chanlocs, 'topo2all');
else
    error('nf_badchans:MissingChannelCoordinates', ...
        ['Every analyzed channel needs finite nonzero X/Y/Z or finite ' ...
        'theta/radius coordinates.']);
end
nf_validate_locations(EEG.chanlocs);
end

function nf_validate_locations(chanlocs)
coordinates = zeros(numel(chanlocs), 3);
for index = 1:numel(chanlocs)
    if ~nf_has_xyz(chanlocs(index)) || ...
            ~nf_has_topography(chanlocs(index))
        error('nf_badchans:MissingChannelCoordinates', ...
            ['Channel %s lacks geometry required by FASTER and EEGLAB ' ...
            'interpolation.'], char(chanlocs(index).labels));
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
        error('nf_badchans:DuplicateChannelCoordinates', ...
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
        error('nf_badchans:UnknownInterpolationMethod', ...
            ['interpolationMethod must be spherical, sphericalKang, ' ...
            'sphericalCRD, invdist, or v4.']);
end
end

function available = nf_function_available(name)
available = exist(name, 'file') == 2 || exist(name, 'builtin') == 5;
end

function nf_validate_eeg(EEG)
if ~isstruct(EEG) || numel(EEG) ~= 1
    error('nf_badchans:InvalidEEG', 'EEG must be one EEGLAB dataset structure.');
end
required = {'data', 'srate', 'nbchan', 'pnts', 'trials', 'chanlocs'};
for index = 1:numel(required)
    if ~isfield(EEG, required{index})
        error('nf_badchans:InvalidEEG', 'EEG.%s is required.', required{index});
    end
end
if ~isnumeric(EEG.nbchan) || ~isscalar(EEG.nbchan) || ...
        ~isfinite(EEG.nbchan) || EEG.nbchan < 3 || ...
        EEG.nbchan ~= round(EEG.nbchan) || ...
        ~isnumeric(EEG.pnts) || ~isscalar(EEG.pnts) || ...
        ~isfinite(EEG.pnts) || EEG.pnts < 2 || EEG.pnts ~= round(EEG.pnts) || ...
        ~isnumeric(EEG.trials) || ~isscalar(EEG.trials) || ...
        ~isfinite(EEG.trials) || EEG.trials < 1 || ...
        EEG.trials ~= round(EEG.trials) || ...
        ~isnumeric(EEG.data) || ~isreal(EEG.data) || isempty(EEG.data) || ...
        EEG.nbchan ~= size(EEG.data, 1) || EEG.pnts ~= size(EEG.data, 2) || ...
        EEG.trials ~= size(EEG.data, 3) || ...
        ~isnumeric(EEG.srate) || ~isscalar(EEG.srate) || ...
        ~isfinite(EEG.srate) || EEG.srate <= 0
    error('nf_badchans:InvalidEEG', 'EEG.data is empty or inconsistent with EEG.nbchan.');
end
if numel(EEG.chanlocs) ~= EEG.nbchan
    error('nf_badchans:InvalidEEG', 'EEG.chanlocs must contain one entry per channel.');
end
end

function valid = nf_is_reference(value)
valid = isempty(value) || ischar(value) || ...
    (isstring(value) && isscalar(value)) || ...
    (isnumeric(value) && isscalar(value) && isfinite(value));
end

function valid = nf_is_faster_options(value)
valid = isstruct(value) && isscalar(value);
end

function options = nf_normalize_faster_options(options)
propertyCount = 3;
allowed = {'measure', 'z'};
names = fieldnames(options);
for index = 1:numel(names)
    if ~ismember(names{index}, allowed)
        error('nf_badchans:UnknownFASTEROption', ...
            'Unknown fasterOptions field: %s.', names{index});
    end
end
if ~isfield(options, 'measure')
    options.measure = true(1, propertyCount);
else
    value = options.measure;
    if ~(isnumeric(value) || islogical(value)) || ...
            numel(value) ~= propertyCount || ...
            any(~isfinite(value(:))) || any(~ismember(value(:), [0 1]))
        error('nf_badchans:InvalidFASTEROptions', ...
            'fasterOptions.measure must contain three binary values.');
    end
    options.measure = logical(reshape(value, 1, propertyCount));
end
if ~isfield(options, 'z')
    options.z = 3 * ones(1, propertyCount);
else
    value = options.z;
    if ~isnumeric(value) || ~isreal(value) || ...
            numel(value) ~= propertyCount || ...
            any(~isfinite(value(:))) || any(value(:) <= 0)
        error('nf_badchans:InvalidFASTEROptions', ...
            'fasterOptions.z must contain three positive finite values.');
    end
    options.z = reshape(value, 1, propertyCount);
end
end

function valid = nf_is_correlation(value)
valid = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value > 0 && value <= 1;
end

function valid = nf_is_positive_scalar(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value > 0;
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
