function [EEG, info] = nf_badchans(EEG, maxBad, interpolate, method, varargin)
% NF_BADCHANS  Detect and remove bad scalp channels before precleaning.
%
% [EEG, INFO] = NF_BADCHANS(EEG, MAXBAD, INTERPOLATE, METHOD, ...)
%
% METHOD is 'faster' (default), 'cleanrawdata', 'prep', 'happeer',
% 'eeglab', or 'none'. HAPPE+ER may also be written as 'HAPPE+ER'.
% Named vendor methods execute their resolved MATLAB entry points and
% record the resolved files, observed signatures, versions, commits,
% and SHA-256 hashes without assuming an upstream release identity.
% METHOD='none' reports flat/nonfinite diagnostics without removing them.
% MAXBAD is a QC ceiling: the function errors when a detector exceeds it;
% the detector is never relaxed to force a dataset through the ceiling.
% Standard mean/SD z-scoring is enforced across compatible min_z variants.
% For average-referenced data without a recoverable recording reference,
% the intended no-reference channel properties are calculated locally to
% avoid the undefined dist_inds branch in some installed FASTER releases.
%
% Name/value inputs:
%   reference             [] for automatic, a channel label, or a channel index
%   fasterOptions         overrides for the standard FASTER min_z fields:
%                             measure=[1 1 1] and z=[3 3 3]
%   cleanCorrelation      clean_channels correlation threshold (default 0.8)
%   cleanHighpass         diagnostic clean_channels high-pass (default 1 Hz)
%   cleanrawdataOptions   clean_flatlines and clean_channels parameters
%   prepOptions           PREP removeTrend and findNoisyChannels fields
%   happeerOptions        Current HAPPE bad-channel stage parameters:
%                             lowDensity=[] (must be supplied explicitly)
%                             allowLowDensityInference=false
%                             stagePoint=1 (pre-wavelet stage only)
%                             runInitialCleanRawData=true
%                             initialFlatlineCriterion=3
%                             initialChannelCriterion=0.1
%                             initialLineNoiseCriterion=20
%                             distance='Euclidian'
%                             params=struct()
%   eeglabOptions         pop_rejchan measures and parameters
%   interpolationMethod   'sphericalKang' (default)
%
% In a preprocessing pipeline, pass INTERPOLATE=false. Globally removed
% channels should remain absent through GEDAI and ICA, then be interpolated
% once at the end. INTERPOLATE=true is retained for direct-call backwards
% compatibility.

nf_validate_eeg(EEG);

if nargin < 2 || isempty(maxBad)
    maxBad = floor(EEG.nbchan / 10);
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
addParameter(parser, 'cleanrawdataOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'prepOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'happeerOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'eeglabOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'interpolationMethod', 'sphericalKang', @nf_is_text);
parse(parser, varargin{:});
options = parser.Results;
method = nf_normalize_channel_method(method);
if strcmp(method, 'faster')
    options.fasterOptions = ...
        nf_normalize_faster_options(options.fasterOptions);
elseif strcmp(method, 'cleanrawdata')
    options.cleanrawdataOptions = nf_normalize_cleanrawdata_options( ...
        options.cleanrawdataOptions, options.cleanCorrelation, ...
        options.cleanHighpass);
elseif strcmp(method, 'happeer')
    options.happeerOptions = nf_normalize_happeer_options( ...
        options.happeerOptions);
elseif strcmp(method, 'eeglab')
    options.eeglabOptions = nf_normalize_eeglab_options( ...
        options.eeglabOptions, EEG.srate);
end
options.interpolationMethod = nf_normalize_interpolation_method( ...
    options.interpolationMethod);

if ~nf_is_nonnegative_integer(maxBad)
    error('nf_badchans:InvalidMaximum', 'maxBad must be a nonnegative integer.');
end
if ~nf_is_logical_scalar(interpolate)
    error('nf_badchans:InvalidInterpolationFlag', 'interpolate must be a logical scalar.');
end

if ismember(method, {'faster', 'cleanrawdata', 'prep', 'happeer'}) || ...
        logical(interpolate)
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
        referenceIndex, referenceIsRemoved, ...
        options.cleanrawdataOptions);
elseif strcmp(method, 'prep')
    [detectorMask, detector] = nf_run_prep(EEG, referenceIndex, ...
        referenceIsRemoved, options.prepOptions);
elseif strcmp(method, 'happeer')
    [detectorMask, detector] = nf_run_happeer(EEG, mandatoryMask, ...
        referenceIndex, referenceIsRemoved, options.happeerOptions);
elseif strcmp(method, 'eeglab')
    [detectorMask, detector] = nf_run_eeglab(EEG, mandatoryMask, ...
        referenceIndex, referenceIsRemoved, options.eeglabOptions);
elseif strcmp(method, 'none')
    [detectorMask, detector] = nf_run_none(EEG);
else
    error('nf_badchans:UnknownMethod', ...
        ['method must be faster, cleanrawdata, prep, happeer, ' ...
        'eeglab, or none.']);
end

if strcmp(method, 'none')
    detectedArtifactMask = false(1, EEG.nbchan);
else
    detectedArtifactMask = mandatoryMask | detectorMask;
end
detectedArtifactIndices = find(detectedArtifactMask);
if numel(detectedArtifactIndices) > maxBad
    badList = strjoin(originalLabels(detectedArtifactIndices), ', ');
    error('nf_badchans:TooManyBadChannels', ...
        ['The detector identified %d artifact channels, above maxBad=%d. ' ...
        'The dataset was not modified. Channels: %s'], ...
        numel(detectedArtifactIndices), maxBad, badList);
end

referenceRemovalMask = false(1, EEG.nbchan);
effectiveReferenceRemoval = referenceIsRemoved && ~strcmp(method, 'none');
if effectiveReferenceRemoval
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
info.schemaVersion = '3.0.0';
info.method = method;
info.nOriginal = numel(originalLabels);
info.maxBad = maxBad;
info.reference.mode = referenceMode;
info.reference.index = referenceIndex;
info.reference.label = nf_label_at(originalLabels, referenceIndex);
info.reference.identifiedZeroReference = referenceIsRemoved;
info.reference.removedZeroReference = effectiveReferenceRemoval;
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
info.provenance = detector.provenance;
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
EEG.etc.nf_removed_reference = struct('removed', effectiveReferenceRemoval, ...
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
vendorContract = nf_vendor_contract( ...
    'FASTER', ...
    'vendor-exact-stage', ...
    {'channel_properties', 'min_z'}, ...
    {{'faster'}, {'faster'}});

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
    vendorContract.contractLevel = 'vendor-primitive-compatible';
    vendorContract.notes = ...
        ['Official FASTER min_z was used after NeuroFreq calculated ' ...
        'the intended no-reference channel properties.'];
else
    listProperties = channel_properties( ...
        diagnostic, ...
        1:diagnostic.nbchan, ...
        diagnosticReference);
    propertyImplementation = which('channel_properties');
end
minZOptions = fasterOptions;
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
detector.provenance = vendorContract;
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
    referenceIndex, referenceIsRemoved, cleanOptions)
vendorContract = nf_vendor_contract( ...
    'clean_rawdata', ...
    'vendor-exact-stage', ...
    {'clean_flatlines', 'clean_channels'}, ...
    {{'clean_rawdata', 'cleanrawdata'}, ...
    {'clean_rawdata', 'cleanrawdata'}});
nf_require_vendor_signature('clean_flatlines', 3, 1);
nf_require_vendor_signature('clean_channels', 7, 2);
cleanRawRoots = cellfun(@fileparts, ...
    {vendorContract.functions.path}, 'UniformOutput', false);
if numel(unique(lower(string(cleanRawRoots)))) ~= 1
    error('nf_badchans:MixedCleanRawDataInstall', ...
        ['clean_flatlines and clean_channels resolved from different ' ...
        'clean_rawdata installations. Remove the mixed or shadowed copy.']);
end
vendorContract.packageRoot = cleanRawRoots{1};
vendorContract.coLocatedStageFunctionsVerified = true;
nf_require_continuous(EEG, 'cleanrawdata');
minimumSeconds = max(cleanOptions.flatlineDuration, ...
    cleanOptions.windowLength);
nf_require_duration(EEG, minimumSeconds, 'cleanrawdata');
if cleanOptions.diagnosticHighpassHz >= EEG.srate / 2
    error('nf_badchans:InvalidCleanHighpass', ...
        'cleanrawdataOptions.diagnosticHighpassHz must be below Nyquist.');
end
if cleanOptions.diagnosticHighpassHz > 0
    nf_require_function('pop_eegfiltnew', 'EEGLAB FIR filtering');
end

excludedMask = nf_nonfinite_channel_mask(EEG);
if referenceIsRemoved
    excludedMask(referenceIndex) = true;
end
usableIndices = find(~excludedMask);
if numel(usableIndices) < 3
    error('nf_badchans:TooFewUsableChannels', ...
        'clean_rawdata requires at least three usable channels.');
end

diagnostic = pop_select(EEG, 'channel', usableIndices);
diagnostic = eeg_checkset(diagnostic);

beforeFlat = diagnostic;
diagnostic = clean_flatlines( ...
    diagnostic, ...
    cleanOptions.flatlineDuration, ...
    cleanOptions.flatlineJitter);
[currentOriginalIndices, flatRemovedOriginalIndices] = ...
    nf_track_vendor_channel_selection( ...
    usableIndices, beforeFlat, diagnostic, 'clean_flatlines');
if numel(currentOriginalIndices) < 3
    error('nf_badchans:TooFewUsableChannels', ...
        'clean_flatlines left fewer than three channels.');
end
if cleanOptions.diagnosticHighpassHz > 0
    diagnostic = pop_eegfiltnew(diagnostic, ...
        'locutoff', cleanOptions.diagnosticHighpassHz);
    diagnostic = eeg_checkset(diagnostic);
end

previousRandomState = rng;
randomCleanup = onCleanup(@() rng(previousRandomState));
[~, removedLocal] = clean_channels( ...
    diagnostic, ...
    cleanOptions.correlationThreshold, ...
    cleanOptions.lineNoiseThreshold, ...
    cleanOptions.windowLength, ...
    cleanOptions.maxBrokenTime, ...
    cleanOptions.ransacSamples, ...
    cleanOptions.subsetSize);
clear randomCleanup
removedLocal = nf_cleanraw_removed_mask( ...
    removedLocal, numel(currentOriginalIndices));
channelRemovedOriginalIndices = ...
    currentOriginalIndices(removedLocal);

detectorMask = false(1, EEG.nbchan);
detectorMask(flatRemovedOriginalIndices) = true;
detectorMask(channelRemovedOriginalIndices) = true;

detector = struct();
detector.name = ...
    'clean_rawdata clean_flatlines + clean_channels';
detector.options = cleanOptions;
detector.analyzedOriginalIndices = usableIndices;
detector.flatlineRemovedOriginalIndices = ...
    flatRemovedOriginalIndices;
detector.cleanChannelsRemovedOriginalIndices = ...
    channelRemovedOriginalIndices;
detector.nativeSafetyOriginalIndices = ...
    find(mandatoryMask & ~detectorMask);
detector.provenance = vendorContract;
end

function [detectorMask, detector] = nf_run_prep(EEG, referenceIndex, ...
    referenceIsRemoved, prepOptions)
vendorContract = nf_vendor_contract( ...
    'PREP', ...
    'vendor-exact-stage-composition', ...
    {'removeTrend', 'findNoisyChannels'}, ...
    {{'eeg-clean-tools', 'preppipeline'}, ...
    {'eeg-clean-tools', 'preppipeline'}});
nf_require_vendor_signature('removeTrend', 2, 2);
nf_require_vendor_signature('findNoisyChannels', 2, 1);
prepRoots = cellfun(@fileparts, ...
    {vendorContract.functions.path}, 'UniformOutput', false);
if numel(unique(lower(string(prepRoots)))) ~= 1
    error('nf_badchans:MixedPREPInstall', ...
        ['removeTrend and findNoisyChannels resolved from different PREP ' ...
        'installations. Remove the mixed or shadowed installation.']);
end
vendorContract.version = '';
vendorContract.upstreamReleaseVerified = false;
vendorContract.packageStageRoot = prepRoots{1};
vendorContract.coLocatedPrepHelpersVerified = true;
vendorContract.verificationScope = ...
    ['Co-located PREP stage files with recorded SHA-256 identities; no ' ...
    'upstream release label is asserted by this standalone adapter.'];
nf_require_continuous(EEG, 'PREP findNoisyChannels');
resolvedOptions = prepOptions;
if ~isfield(resolvedOptions, 'evaluationChannels') || ...
        isempty(resolvedOptions.evaluationChannels)
    resolvedOptions.evaluationChannels = 1:EEG.nbchan;
end
resolvedOptions.evaluationChannels = nf_validate_channel_indices( ...
    resolvedOptions.evaluationChannels, EEG.nbchan, ...
    'prepOptions.evaluationChannels');
if referenceIsRemoved
    resolvedOptions.evaluationChannels = setdiff( ...
        resolvedOptions.evaluationChannels, referenceIndex, 'stable');
end
if numel(resolvedOptions.evaluationChannels) < 3
    error('nf_badchans:TooFewUsableChannels', ...
        'PREP requires at least three evaluation channels.');
end

previousRandomState = rng;
randomCleanup = onCleanup(@() rng(previousRandomState));
[diagnostic, detrendOutput] = removeTrend(EEG, resolvedOptions);
prepOutput = findNoisyChannels(diagnostic, resolvedOptions);
clear randomCleanup
if ~isstruct(prepOutput) || ...
        ~isfield(prepOutput, 'noisyChannels') || ...
        ~isstruct(prepOutput.noisyChannels) || ...
        ~isfield(prepOutput.noisyChannels, 'all')
    error('nf_badchans:InvalidPREPOutput', ...
        'findNoisyChannels did not return noisyChannels.all.');
end
badIndices = nf_validate_vendor_indices( ...
    prepOutput.noisyChannels.all, EEG.nbchan, ...
    'findNoisyChannels');
if referenceIsRemoved
    badIndices = setdiff(badIndices, referenceIndex, 'stable');
end

detectorMask = false(1, EEG.nbchan);
detectorMask(badIndices) = true;

detector = struct();
detector.name = 'PREP findNoisyChannels';
detector.optionsRequested = prepOptions;
detector.optionsResolved = nf_prep_resolved_options(prepOutput);
detector.detrend = detrendOutput;
detector.analyzedOriginalIndices = ...
    resolvedOptions.evaluationChannels;
detector.vendorOutput = prepOutput;
detector.provenance = vendorContract;
detector.provenance.notes = ...
    ['The installed PREP removeTrend and findNoisyChannels functions were ' ...
    'called in their released order for standalone channel screening. ' ...
    'This is a PREP stage composition, not the complete prepPipeline.'];
end

function [detectorMask, detector] = nf_run_happeer(EEG, mandatoryMask, ...
    referenceIndex, referenceIsRemoved, happeOptions)
vendorContract = nf_vendor_contract( ...
    'HAPPE+ER bad-channel stage', ...
    'vendor-exact-stage', ...
    {'happe_detectBadChans', 'pop_clean_rawdata', 'pop_rejchan'}, ...
    {{'happe'}, ...
    {'clean_rawdata', 'cleanrawdata'}, ...
    {'eeglab'}});
nf_require_vendor_signature('happe_detectBadChans', 3, 1);
vendorContract.notes = ...
    ['NeuroFreq calls the uniquely resolved pop_clean_rawdata and ' ...
    'happe_detectBadChans files. The recorded paths, observed signatures, ' ...
    'and SHA-256 hashes identify the installed code; they do not assert ' ...
    'an upstream release when release metadata is unavailable.'];
nf_require_continuous(EEG, 'HAPPE+ER bad-channel detection');
if EEG.srate / 2 < 100
    error('nf_badchans:HAPPEERSamplingRate', ...
        ['This happe_detectBadChans adapter evaluates through 100 Hz. ' ...
        'The sampling rate must therefore be at least 200 Hz.']);
end

excludedMask = nf_nonfinite_channel_mask(EEG);
if referenceIsRemoved
    excludedMask(referenceIndex) = true;
end
currentOriginalIndices = find(~excludedMask);
if numel(currentOriginalIndices) < 3
    error('nf_badchans:TooFewUsableChannels', ...
        'HAPPE+ER requires at least three finite channels.');
end
diagnostic = pop_select(EEG, ...
    'channel', currentOriginalIndices);
diagnostic = eeg_checkset(diagnostic);

hasTopLevelLowDensity = ~isempty(happeOptions.lowDensity);
hasParameterLowDensity = ...
    isfield(happeOptions.params, 'lowDensity') && ...
    ~isempty(happeOptions.params.lowDensity);
if hasParameterLowDensity && ...
        ~nf_is_logical_scalar(happeOptions.params.lowDensity)
    error('nf_badchans:InvalidHAPPEEROptions', ...
        ['happeerOptions.params.lowDensity must be a logical ' ...
        'scalar when supplied.']);
end
if hasTopLevelLowDensity && hasParameterLowDensity && ...
        logical(happeOptions.lowDensity) ~= ...
        logical(happeOptions.params.lowDensity)
    error('nf_badchans:ConflictingHAPPEEROptions', ...
        ['happeerOptions.lowDensity and params.lowDensity disagree. ' ...
        'Supply one value or make them identical.']);
end

lowDensityResolution = struct();
lowDensityResolution.inferred = false;
lowDensityResolution.channelCount = [];
if hasTopLevelLowDensity
    resolvedLowDensity = happeOptions.lowDensity;
    lowDensityResolution.source = 'happeerOptions.lowDensity';
elseif hasParameterLowDensity
    resolvedLowDensity = happeOptions.params.lowDensity;
    lowDensityResolution.source = ...
        'happeerOptions.params.lowDensity';
elseif happeOptions.allowLowDensityInference
    [inferenceChannelCount, inferenceSource] = ...
        nf_happeer_inference_channel_count(EEG);
    resolvedLowDensity = inferenceChannelCount <= 32;
    lowDensityResolution.inferred = true;
    lowDensityResolution.channelCount = inferenceChannelCount;
    lowDensityResolution.source = inferenceSource;
    vendorContract.contractLevel = ...
        'vendor-exact-stage-with-neurofreq-mode-inference';
    vendorContract.notes = [vendorContract.notes ...
        ' lowDensity was inferred only because ' ...
        'allowLowDensityInference=true; the inference source and ' ...
        'channel count are recorded in optionsResolved.'];
else
    error('nf_badchans:HAPPEERLowDensityRequired', ...
        ['HAPPE+ER requires an explicit acquisition-density mode. ' ...
        'Set happeerOptions.lowDensity (recommended), set ' ...
        'happeerOptions.params.lowDensity, or explicitly opt in to ' ...
        'the non-strict allowLowDensityInference fallback.']);
end
resolvedLowDensity = logical(resolvedLowDensity);
resolvedParameters = happeOptions.params;
resolvedParameters.lowDensity = resolvedLowDensity;

initialRemovedOriginalIndices = [];
if happeOptions.runInitialCleanRawData
    beforeInitial = diagnostic;
    previousRandomState = rng;
    randomCleanup = onCleanup(@() rng(previousRandomState));
    diagnostic = pop_clean_rawdata( ...
        diagnostic, ...
        'FlatlineCriterion', ...
        happeOptions.initialFlatlineCriterion, ...
        'ChannelCriterion', ...
        happeOptions.initialChannelCriterion, ...
        'LineNoiseCriterion', ...
        happeOptions.initialLineNoiseCriterion, ...
        'Highpass', 'off', ...
        'BurstCriterion', 'off', ...
        'WindowCriterion', 'off', ...
        'BurstRejection', 'off', ...
        'Distance', happeOptions.distance);
    clear randomCleanup
    [currentOriginalIndices, initialRemovedOriginalIndices] = ...
        nf_track_vendor_channel_selection( ...
        currentOriginalIndices, beforeInitial, diagnostic, ...
        'HAPPE initial pop_clean_rawdata');
end

nf_require_minimum_diagnostic_channels( ...
    diagnostic, 'HAPPE happe_detectBadChans');
beforeHappe = diagnostic;
previousRandomState = rng;
randomCleanup = onCleanup(@() rng(previousRandomState));
diagnostic = happe_detectBadChans( ...
    diagnostic, resolvedParameters, happeOptions.stagePoint);
clear randomCleanup
[currentOriginalIndices, happeRemovedOriginalIndices] = ...
    nf_track_vendor_channel_selection( ...
    currentOriginalIndices, beforeHappe, diagnostic, ...
    'HAPPE happe_detectBadChans');

detectorMask = false(1, EEG.nbchan);
detectorMask(initialRemovedOriginalIndices) = true;
detectorMask(happeRemovedOriginalIndices) = true;
detectorMask(mandatoryMask) = true;
if referenceIsRemoved
    detectorMask(referenceIndex) = false;
end

resolvedOptions = happeOptions;
resolvedOptions.lowDensity = resolvedLowDensity;
resolvedOptions.params = resolvedParameters;
resolvedOptions.lowDensityResolution = lowDensityResolution;
resolvedOptions.vendorSettings = nf_happeer_vendor_settings( ...
    resolvedLowDensity, happeOptions.stagePoint);
detector = struct();
detector.name = 'HAPPE happe_detectBadChans';
detector.optionsRequested = happeOptions;
detector.optionsResolved = resolvedOptions;
detector.analyzedOriginalIndices = find(~excludedMask);
detector.initialCleanRawDataRemovedOriginalIndices = ...
    initialRemovedOriginalIndices;
detector.happeStageRemovedOriginalIndices = ...
    happeRemovedOriginalIndices;
detector.keptOriginalIndices = currentOriginalIndices;
detector.provenance = vendorContract;
end

function settings = nf_happeer_vendor_settings(lowDensity, stagePoint)
settings = struct();
settings.helper = 'happe_detectBadChans';
settings.lowDensity = lowDensity;
settings.stagePoint = stagePoint;
settings.spectrumRange = [1 100];
settings.spectrumNorm = 'on';
settings.distance = 'Euclidian';
if lowDensity
    settings.order = {'pop_clean_rawdata', 'pop_rejchan'};
    settings.flatlineCriterion = 'off';
    settings.channelCriterion = 0.7;
    settings.lineNoiseCriterion = 2.5;
    settings.spectrumThreshold = [-2.75 2.75];
    return
end
settings.order = {'pop_rejchan', 'pop_clean_rawdata'};
settings.flatlineCriterion = 'off';
settings.spectrumThreshold = [-5 1.8935];
settings.channelCriterion = 0.485;
settings.lineNoiseCriterion = 7.1;
end

function [channelCount, source] = ...
    nf_happeer_inference_channel_count(EEG)
channelCount = EEG.nbchan;
source = 'EEG.nbchan-current-selected-montage';
if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc) || ...
        ~isscalar(EEG.etc) || ~isfield(EEG.etc, 'ogchan') || ...
        ~isstruct(EEG.etc.ogchan) || isempty(EEG.etc.ogchan)
    return
end
originalCount = numel(EEG.etc.ogchan);
if originalCount > channelCount
    channelCount = originalCount;
    source = 'EEG.etc.ogchan-original-montage-count';
end
end

function [detectorMask, detector] = nf_run_eeglab(EEG, mandatoryMask, ...
    referenceIndex, referenceIsRemoved, eeglabOptions)
vendorContract = nf_vendor_contract( ...
    'EEGLAB', ...
    'vendor-exact-stage', ...
    {'pop_rejchan'}, ...
    {{'eeglab'}});
excludedMask = mandatoryMask;
if referenceIsRemoved
    excludedMask(referenceIndex) = true;
end
analyzedOriginalIndices = find(~excludedMask);
if numel(analyzedOriginalIndices) < 3
    error('nf_badchans:TooFewUsableChannels', ...
        'EEGLAB pop_rejchan requires at least three usable channels.');
end
diagnostic = pop_select(EEG, ...
    'channel', analyzedOriginalIndices);
diagnostic = eeg_checkset(diagnostic);

runs = repmat(struct( ...
    'measure', '', ...
    'threshold', [], ...
    'removedOriginalIndices', [], ...
    'values', []), 1, numel(eeglabOptions.measures));
detectorMask = false(1, EEG.nbchan);
for measureIndex = 1:numel(eeglabOptions.measures)
    measureName = eeglabOptions.measures{measureIndex};
    threshold = nf_eeglab_threshold( ...
        eeglabOptions, measureName);
    arguments = { ...
        'elec', 1:diagnostic.nbchan, ...
        'threshold', threshold, ...
        'norm', eeglabOptions.norm, ...
        'measure', measureName, ...
        'indexonly', 'on', ...
        'verbose', eeglabOptions.verbose};
    if strcmp(measureName, 'spec')
        arguments = [arguments ...
            {'freqrange', eeglabOptions.frequencyRange}];
    end
    [~, removedLocal, values] = pop_rejchan( ...
        diagnostic, arguments{:});
    removedLocal = nf_validate_vendor_indices( ...
        removedLocal, diagnostic.nbchan, 'EEGLAB pop_rejchan');
    removedOriginal = analyzedOriginalIndices(removedLocal);
    detectorMask(removedOriginal) = true;
    runs(measureIndex).measure = measureName;
    runs(measureIndex).threshold = threshold;
    runs(measureIndex).removedOriginalIndices = ...
        removedOriginal;
    runs(measureIndex).values = values;
end

detector = struct();
detector.name = 'EEGLAB pop_rejchan';
detector.options = eeglabOptions;
detector.analyzedOriginalIndices = analyzedOriginalIndices;
detector.runs = runs;
detector.provenance = vendorContract;
end

function [detectorMask, detector] = nf_run_none(EEG)
detectorMask = false(1, EEG.nbchan);
detector = struct();
detector.name = 'none';
detector.options = struct();
detector.provenance = struct( ...
    'provider', 'NeuroFreq', ...
    'contractLevel', 'native-safety-only', ...
    'version', '', ...
    'commit', '', ...
    'functions', struct([]), ...
    'notes', ...
    ['No automated channel detector ran. Flat and nonfinite channel ' ...
    'diagnostics are reported but are not removed.']);
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

function [keptOriginalIndices, removedOriginalIndices] = ...
    nf_track_vendor_channel_selection( ...
    beforeOriginalIndices, beforeEEG, afterEEG, stage)
if numel(beforeOriginalIndices) ~= beforeEEG.nbchan
    error('nf_badchans:InternalChannelMappingError', ...
        'The channel map entering %s was inconsistent.', stage);
end
beforeLabels = lower(string(nf_channel_labels(beforeEEG.chanlocs)));
afterLabels = lower(string(nf_channel_labels(afterEEG.chanlocs)));
[present, locations] = ismember(afterLabels, beforeLabels);
if any(~present) || numel(unique(locations)) ~= numel(locations)
    error('nf_badchans:VendorChannelMappingError', ...
        '%s returned unknown or duplicate channel labels.', stage);
end
keptOriginalIndices = beforeOriginalIndices(locations);
removedMask = ~ismember(beforeLabels, afterLabels);
removedOriginalIndices = beforeOriginalIndices(removedMask);
end

function indices = nf_validate_channel_indices(value, channelCount, name)
if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
        isempty(value) || any(~isfinite(value(:))) || ...
        any(value(:) ~= round(value(:))) || any(value(:) < 1) || ...
        any(value(:) > channelCount) || ...
        numel(unique(value(:))) ~= numel(value)
    error('nf_badchans:InvalidChannelIndices', ...
        '%s must contain unique valid channel indices.', name);
end
indices = reshape(value, 1, []);
end

function indices = nf_validate_vendor_indices(value, channelCount, stage)
if isempty(value)
    indices = [];
    return
end
if islogical(value)
    value = value(:)';
    if numel(value) ~= channelCount
        error('nf_badchans:InvalidVendorOutput', ...
            '%s returned a logical mask with unexpected length.', stage);
    end
    indices = find(value);
    return
end
if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
        any(~isfinite(value(:))) || any(value(:) ~= round(value(:))) || ...
        any(value(:) < 1) || any(value(:) > channelCount) || ...
        numel(unique(value(:))) ~= numel(value)
    error('nf_badchans:InvalidVendorOutput', ...
        '%s returned invalid channel indices.', stage);
end
indices = reshape(value, 1, []);
end

function resolved = nf_prep_resolved_options(prepOutput)
parameterNames = { ...
    'evaluationChannels', ...
    'robustDeviationThreshold', ...
    'highFrequencyNoiseThreshold', ...
    'correlationWindowSeconds', ...
    'correlationThreshold', ...
    'badTimeThreshold', ...
    'ransacSampleSize', ...
    'ransacChannelFraction', ...
    'ransacCorrelationThreshold', ...
    'ransacUnbrokenTime', ...
    'ransacWindowSeconds', ...
    'ransacOff'};
resolved = struct();
for index = 1:numel(parameterNames)
    fieldName = parameterNames{index};
    if isfield(prepOutput, fieldName)
        resolved.(fieldName) = prepOutput.(fieldName);
    end
end
end

function threshold = nf_eeglab_threshold(options, measure)
if strcmp(measure, 'prob')
    threshold = options.probabilityThreshold;
elseif strcmp(measure, 'kurt')
    threshold = options.kurtosisThreshold;
elseif strcmp(measure, 'spec')
    threshold = options.spectrumThreshold;
elseif strcmp(measure, 'std')
    threshold = options.standardDeviationThreshold;
else
    error('nf_badchans:UnknownEEGLABMeasure', ...
        'Unknown EEGLAB channel measure: %s.', measure);
end
end

function nf_require_continuous(EEG, stage)
if EEG.trials ~= 1
    error('nf_badchans:ContinuousDataRequired', ...
        '%s requires continuous EEG.', stage);
end
end

function nf_require_duration(EEG, minimumSeconds, stage)
durationSeconds = EEG.pnts / EEG.srate;
if durationSeconds <= minimumSeconds
    error('nf_badchans:InsufficientData', ...
        '%s requires more than %.3g seconds of data.', ...
        stage, minimumSeconds);
end
end

function nf_require_minimum_diagnostic_channels(EEG, stage)
if EEG.nbchan < 3
    error('nf_badchans:TooFewUsableChannels', ...
        '%s requires at least three channels.', stage);
end
end

function nf_require_function(functionName, provider)
if ~nf_function_available(functionName)
    error('nf_badchans:MissingDependency', ...
        '%s requires %s.m.', provider, functionName);
end
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
        ['EEGLAB convertlocs.m is required by the selected channel ' ...
        'detector or interpolation.']);
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
            ['Channel %s lacks geometry required by spatial channel ' ...
            'detection and EEGLAB interpolation.'], ...
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

function method = nf_normalize_channel_method(value)
if ~nf_is_text(value)
    error('nf_badchans:InvalidMethod', ...
        'method must be one text scalar.');
end
normalized = lower(regexprep(strtrim(char(value)), '[ _+-]', ''));
if strcmp(normalized, 'faster')
    method = 'faster';
elseif ismember(normalized, {'cleanrawdata', 'cleanraw'})
    method = 'cleanrawdata';
elseif strcmp(normalized, 'prep')
    method = 'prep';
elseif ismember(normalized, {'happeer', 'happier'})
    method = 'happeer';
elseif ismember(normalized, {'eeglab', 'legacy', 'legacyeeglab'})
    method = 'eeglab';
elseif strcmp(normalized, 'none')
    method = 'none';
else
    error('nf_badchans:UnknownMethod', ...
        ['method must be faster, cleanrawdata, prep, happeer, ' ...
        'eeglab, or none.']);
end
end

function options = nf_normalize_cleanrawdata_options( ...
    options, legacyCorrelation, legacyHighpass)
allowed = { ...
    'flatlineDuration', ...
    'flatlineJitter', ...
    'correlationThreshold', ...
    'lineNoiseThreshold', ...
    'windowLength', ...
    'maxBrokenTime', ...
    'ransacSamples', ...
    'subsetSize', ...
    'diagnosticHighpassHz'};
nf_reject_unknown_fields(options, allowed, 'cleanrawdataOptions');
options = nf_set_default( ...
    options, 'flatlineDuration', 5);
options = nf_set_default( ...
    options, 'flatlineJitter', 20);
options = nf_set_default( ...
    options, 'correlationThreshold', legacyCorrelation);
options = nf_set_default( ...
    options, 'lineNoiseThreshold', 4);
options = nf_set_default( ...
    options, 'windowLength', 5);
options = nf_set_default( ...
    options, 'maxBrokenTime', 0.4);
options = nf_set_default( ...
    options, 'ransacSamples', 50);
options = nf_set_default( ...
    options, 'subsetSize', 0.25);
options = nf_set_default( ...
    options, 'diagnosticHighpassHz', legacyHighpass);
nf_validate_clean_channel_options(options, 'cleanrawdataOptions');
if ~nf_is_nonnegative_scalar(options.diagnosticHighpassHz)
    error('nf_badchans:InvalidCleanRawDataOptions', ...
        ['cleanrawdataOptions.diagnosticHighpassHz must be a ' ...
        'nonnegative finite scalar.']);
end
end

function options = nf_normalize_happeer_options(options)
allowed = { ...
    'lowDensity', ...
    'allowLowDensityInference', ...
    'stagePoint', ...
    'runInitialCleanRawData', ...
    'initialFlatlineCriterion', ...
    'initialChannelCriterion', ...
    'initialLineNoiseCriterion', ...
    'distance', ...
    'params'};
nf_reject_unknown_fields(options, allowed, 'happeerOptions');
options = nf_set_default(options, 'lowDensity', []);
options = nf_set_default(options, 'allowLowDensityInference', false);
options = nf_set_default(options, 'stagePoint', 1);
options = nf_set_default(options, 'runInitialCleanRawData', true);
options = nf_set_default(options, 'initialFlatlineCriterion', 3);
options = nf_set_default(options, 'initialChannelCriterion', 0.1);
options = nf_set_default(options, 'initialLineNoiseCriterion', 20);
options = nf_set_default(options, 'distance', 'Euclidian');
options = nf_set_default(options, 'params', struct());
if ~isempty(options.lowDensity) && ...
        ~nf_is_logical_scalar(options.lowDensity)
    error('nf_badchans:InvalidHAPPEEROptions', ...
        ['happeerOptions.lowDensity must be empty or a logical ' ...
        'scalar. Empty requires explicit opt-in to inference.']);
end
if ~isempty(options.lowDensity)
    options.lowDensity = logical(options.lowDensity);
end
if ~nf_is_logical_scalar(options.allowLowDensityInference)
    error('nf_badchans:InvalidHAPPEEROptions', ...
        ['happeerOptions.allowLowDensityInference must be a ' ...
        'logical scalar.']);
end
options.allowLowDensityInference = ...
    logical(options.allowLowDensityInference);
if ~isnumeric(options.stagePoint) || ...
        ~isscalar(options.stagePoint) || ...
        ~isfinite(options.stagePoint) || ...
        options.stagePoint ~= 1
    error('nf_badchans:InvalidHAPPEEROptions', ...
        ['happeerOptions.stagePoint must be 1. HAPPE stage 2 is a ' ...
        'post-wavelet stage and cannot be executed by this ' ...
        'pre-wavelet bad-channel API.']);
end
if ~nf_is_logical_scalar(options.runInitialCleanRawData)
    error('nf_badchans:InvalidHAPPEEROptions', ...
        ['happeerOptions.runInitialCleanRawData must be a ' ...
        'logical scalar.']);
end
options.runInitialCleanRawData = ...
    logical(options.runInitialCleanRawData);
if ~nf_is_positive_or_off(options.initialFlatlineCriterion)
    error('nf_badchans:InvalidHAPPEEROptions', ...
        ['happeerOptions.initialFlatlineCriterion must be a positive ' ...
        'finite scalar or off.']);
end
if nf_is_text(options.initialFlatlineCriterion)
    options.initialFlatlineCriterion = 'off';
end
if ~nf_is_correlation_or_off(options.initialChannelCriterion)
    error('nf_badchans:InvalidHAPPEEROptions', ...
        ['happeerOptions.initialChannelCriterion must be in (0, 1] ' ...
        'or off.']);
end
if nf_is_text(options.initialChannelCriterion)
    options.initialChannelCriterion = 'off';
end
if ~nf_is_positive_or_off(options.initialLineNoiseCriterion)
    error('nf_badchans:InvalidHAPPEEROptions', ...
        ['happeerOptions.initialLineNoiseCriterion must be a positive ' ...
        'finite scalar or off.']);
end
if nf_is_text(options.initialLineNoiseCriterion)
    options.initialLineNoiseCriterion = 'off';
end
if ~nf_is_text(options.distance) || ...
        isempty(strtrim(char(options.distance)))
    error('nf_badchans:InvalidHAPPEEROptions', ...
        'happeerOptions.distance must be nonempty text.');
end
options.distance = char(options.distance);
if ~nf_is_scalar_struct(options.params)
    error('nf_badchans:InvalidHAPPEEROptions', ...
        'happeerOptions.params must be one scalar struct.');
end
end

function options = nf_normalize_eeglab_options(options, samplingRate)
allowed = { ...
    'measures', ...
    'norm', ...
    'verbose', ...
    'probabilityThreshold', ...
    'kurtosisThreshold', ...
    'spectrumThreshold', ...
    'standardDeviationThreshold', ...
    'frequencyRange'};
nf_reject_unknown_fields(options, allowed, 'eeglabOptions');
options = nf_set_default(options, 'measures', {'kurt'});
options = nf_set_default(options, 'norm', 'on');
options = nf_set_default(options, 'verbose', 'on');
options = nf_set_default(options, 'probabilityThreshold', 5);
options = nf_set_default(options, 'kurtosisThreshold', 5);
options = nf_set_default(options, 'spectrumThreshold', 5);
options = nf_set_default(options, 'standardDeviationThreshold', 100);
defaultUpperFrequency = min(45, samplingRate / 2);
if defaultUpperFrequency <= 1
    defaultLowerFrequency = 0;
else
    defaultLowerFrequency = 1;
end
options = nf_set_default( ...
    options, 'frequencyRange', ...
    [defaultLowerFrequency defaultUpperFrequency]);
options.measures = nf_normalize_eeglab_measures(options.measures);
options.norm = nf_normalize_on_off( ...
    options.norm, 'eeglabOptions.norm');
options.verbose = nf_normalize_on_off( ...
    options.verbose, 'eeglabOptions.verbose');
thresholdFields = { ...
    'probabilityThreshold', ...
    'kurtosisThreshold', ...
    'spectrumThreshold', ...
    'standardDeviationThreshold'};
for index = 1:numel(thresholdFields)
    fieldName = thresholdFields{index};
    if ~nf_is_threshold(options.(fieldName))
        error('nf_badchans:InvalidEEGLABOptions', ...
            ['eeglabOptions.%s must contain one finite threshold or ' ...
            'an increasing pair.'], fieldName);
    end
end
if ~nf_is_nonnegative_ordered_pair(options.frequencyRange) || ...
        options.frequencyRange(2) > samplingRate / 2
    error('nf_badchans:InvalidEEGLABOptions', ...
        ['eeglabOptions.frequencyRange must be an increasing ' ...
        'nonnegative pair at or below Nyquist.']);
end
options.frequencyRange = reshape(options.frequencyRange, 1, 2);
end

function measures = nf_normalize_eeglab_measures(value)
if nf_is_text(value)
    value = {char(value)};
elseif isstring(value)
    value = cellstr(value(:));
end
if ~iscell(value) || isempty(value)
    error('nf_badchans:InvalidEEGLABOptions', ...
        'eeglabOptions.measures must contain at least one measure.');
end
measures = cell(1, numel(value));
for index = 1:numel(value)
    if ~nf_is_text(value{index})
        error('nf_badchans:InvalidEEGLABOptions', ...
            'Every eeglabOptions.measures entry must be text.');
    end
    normalized = lower(regexprep( ...
        strtrim(char(value{index})), '[ _-]', ''));
    if ismember(normalized, {'prob', 'probability', 'jointprobability'})
        measures{index} = 'prob';
    elseif ismember(normalized, {'kurt', 'kurtosis'})
        measures{index} = 'kurt';
    elseif ismember(normalized, {'spec', 'spectrum', 'spectral'})
        measures{index} = 'spec';
    elseif ismember(normalized, {'std', 'standarddeviation', 'sd'})
        measures{index} = 'std';
    else
        error('nf_badchans:InvalidEEGLABOptions', ...
            'Unknown eeglabOptions measure: %s.', char(value{index}));
    end
end
if numel(unique(measures)) ~= numel(measures)
    error('nf_badchans:InvalidEEGLABOptions', ...
        'eeglabOptions.measures cannot contain duplicates.');
end
end

function nf_validate_clean_channel_options(options, prefix)
positiveFields = { ...
    'flatlineDuration', ...
    'flatlineJitter', ...
    'lineNoiseThreshold', ...
    'windowLength', ...
    'maxBrokenTime'};
for index = 1:numel(positiveFields)
    fieldName = positiveFields{index};
    if ~nf_is_positive_scalar(options.(fieldName))
        error('nf_badchans:InvalidChannelOptions', ...
            '%s.%s must be a positive finite scalar.', ...
            prefix, fieldName);
    end
end
if ~nf_is_correlation(options.correlationThreshold)
    error('nf_badchans:InvalidChannelOptions', ...
        '%s.correlationThreshold must be in (0, 1].', prefix);
end
if ~nf_is_positive_integer(options.ransacSamples)
    error('nf_badchans:InvalidChannelOptions', ...
        '%s.ransacSamples must be a positive integer.', prefix);
end
if ~nf_is_fraction(options.subsetSize)
    error('nf_badchans:InvalidChannelOptions', ...
        '%s.subsetSize must be in (0, 1].', prefix);
end
end

function options = nf_set_default(options, fieldName, value)
if ~isfield(options, fieldName) || isempty(options.(fieldName))
    options.(fieldName) = value;
end
end

function nf_reject_unknown_fields(options, allowed, prefix)
names = fieldnames(options);
for index = 1:numel(names)
    if ~ismember(names{index}, allowed)
        error('nf_badchans:UnknownOption', ...
            'Unknown %s field: %s.', prefix, names{index});
    end
end
end

function value = nf_normalize_on_off(value, name)
if nf_is_logical_scalar(value)
    if logical(value)
        value = 'on';
    else
        value = 'off';
    end
    return
end
if ~nf_is_text(value)
    error('nf_badchans:InvalidOnOffOption', ...
        '%s must be on, off, true, or false.', name);
end
value = lower(strtrim(char(value)));
if ~ismember(value, {'on', 'off'})
    error('nf_badchans:InvalidOnOffOption', ...
        '%s must be on, off, true, or false.', name);
end
end

function contract = nf_vendor_contract( ...
    provider, contractLevel, functionNames, acceptedTokens)
if numel(functionNames) ~= numel(acceptedTokens)
    error('nf_badchans:InternalContractError', ...
        'The dependency contract is internally inconsistent.');
end
records = repmat(struct( ...
    'name', '', ...
    'path', '', ...
    'sha256', '', ...
    'release', '', ...
    'commit', '', ...
    'inputArity', NaN, ...
    'outputArity', NaN), 1, numel(functionNames));
for index = 1:numel(functionNames)
    functionName = functionNames{index};
    functionPath = nf_resolve_vendor_function( ...
        functionName, acceptedTokens{index}, provider);
    records(index).name = functionName;
    records(index).path = functionPath;
    records(index).sha256 = nf_file_sha256(functionPath);
    if isempty(records(index).sha256)
        error('nf_badchans:VendorHashFailed', ...
            ['Could not calculate a SHA-256 hash for %s at %s. ' ...
            'The strict vendor contract was not satisfied.'], ...
            functionName, functionPath);
    end
    records(index).release = nf_infer_release(functionPath);
    records(index).commit = nf_nearest_git_commit(functionPath);
    [records(index).inputArity, records(index).outputArity] = ...
        nf_vendor_signature(functionName);
end
contract = struct();
contract.provider = provider;
contract.contractLevel = contractLevel;
contract.identitySource = ...
    'unique MATLAB path resolution plus SHA-256 file identity';
contract.version = nf_provider_version(provider, records);
contract.commit = nf_contract_commits(records);
contract.functions = records;
contract.verified = true;
contract.codeIdentityVerified = true;
contract.upstreamReleaseVerified = false;
contract.verificationScope = ...
    ['Installed-file identity only; upstream release authenticity is ' ...
    'not inferred when release or commit metadata is unavailable.'];
contract.notes = '';
end

function functionPath = nf_resolve_vendor_function( ...
    functionName, acceptedTokens, provider)
resolved = which(functionName, '-all');
if isempty(resolved)
    error('nf_badchans:MissingVendorDependency', ...
        '%s requires an installed %s.m vendor entry point.', ...
        provider, functionName);
end
if ischar(resolved)
    resolved = cellstr(resolved);
elseif isstring(resolved)
    resolved = cellstr(resolved(:));
elseif ~iscell(resolved)
    error('nf_badchans:DependencyResolutionFailed', ...
        'MATLAB returned an unsupported path result for %s.', ...
        functionName);
end
resolved = resolved(~cellfun(@isempty, resolved));
resolved = unique(resolved, 'stable');
if numel(resolved) ~= 1
    error('nf_badchans:ShadowedVendorDependency', ...
        ['%s resolved to %d files. Remove shadowed copies before ' ...
        'running a strict vendor method. Paths: %s'], ...
        functionName, numel(resolved), strjoin(resolved, ' | '));
end
functionPath = resolved{1};
normalizedPath = lower(strrep(functionPath, '\', '/'));
matchesProvider = false;
for index = 1:numel(acceptedTokens)
    matchesProvider = matchesProvider || ...
        contains(normalizedPath, lower(acceptedTokens{index}));
end
if ~matchesProvider
    error('nf_badchans:UnverifiedVendorDependency', ...
        ['%s resolved outside the expected %s installation: %s. ' ...
        'NeuroFreq will not silently treat this file as vendor code.'], ...
        functionName, provider, functionPath);
end
end

function [inputArity, outputArity] = nf_vendor_signature(functionName)
inputArity = NaN;
outputArity = NaN;
try
    inputArity = nargin(functionName);
catch
    inputArity = NaN;
end
try
    outputArity = nargout(functionName);
catch
    outputArity = NaN;
end
end

function nf_require_vendor_signature( ...
    functionName, expectedInputArity, expectedOutputArity)
[inputArity, outputArity] = nf_vendor_signature(functionName);
if ~isfinite(inputArity) || ~isfinite(outputArity) || ...
        inputArity ~= expectedInputArity || ...
        outputArity ~= expectedOutputArity
    error('nf_badchans:UnexpectedVendorSignature', ...
        ['%s has signature nargin=%g, nargout=%g; this adapter ' ...
        'requires nargin=%d, nargout=%d. Use the matching installed ' ...
        'MATLAB release or update the adapter explicitly.'], ...
        functionName, inputArity, outputArity, ...
        expectedInputArity, expectedOutputArity);
end
end

function hash = nf_file_sha256(filePath)
hash = '';
fileId = fopen(filePath, 'rb');
if fileId < 0
    return
end
fileCleanup = onCleanup(@() fclose(fileId));
bytes = fread(fileId, Inf, '*uint8');
clear fileCleanup
try
    digestEngine = java.security.MessageDigest.getInstance('SHA-256');
    digestEngine.update(typecast(bytes, 'int8'));
    digest = typecast(digestEngine.digest(), 'uint8');
    hash = lower(reshape(dec2hex(digest, 2).', 1, []));
catch
    hash = '';
end
end

function release = nf_infer_release(filePath)
normalizedPath = strrep(filePath, '\', '/');
expression = [ ...
    ['(?i)(?:prep(?:pipeline)?|faster|clean[_-]?rawdata|eeglab|' ...
    'happe(?:\+?er)?)'] ...
    '[^/]*?([0-9]+(?:\.[0-9]+)+)'];
tokens = regexp(normalizedPath, expression, 'tokens', 'once');
if isempty(tokens)
    release = '';
else
    release = tokens{1};
end
end

function version = nf_provider_version(provider, records)
version = '';
if contains(lower(provider), 'prep') && ...
        nf_function_available('getPrepVersion')
    try
        version = char(string(getPrepVersion()));
    catch
        version = '';
    end
elseif strcmpi(provider, 'EEGLAB') && ...
        nf_function_available('eeg_getversion')
    try
        version = char(string(eeg_getversion()));
    catch
        version = '';
    end
end
if ~isempty(version)
    return
end
releases = {records.release};
releases = releases(~cellfun(@isempty, releases));
releases = unique(releases, 'stable');
if numel(releases) == 1
    version = releases{1};
elseif numel(releases) > 1
    version = strjoin(releases, ' + ');
end
end

function commits = nf_contract_commits(records)
values = {records.commit};
values = values(~cellfun(@isempty, values));
values = unique(values, 'stable');
if isempty(values)
    commits = '';
elseif numel(values) == 1
    commits = values{1};
else
    commits = values;
end
end

function commit = nf_nearest_git_commit(filePath)
commit = '';
directory = fileparts(filePath);
for level = 1:8
    gitDirectory = fullfile(directory, '.git');
    if exist(gitDirectory, 'dir') == 7
        headPath = fullfile(gitDirectory, 'HEAD');
        if exist(headPath, 'file') ~= 2
            return
        end
        head = strtrim(fileread(headPath));
        if startsWith(head, 'ref:')
            referenceName = strtrim(head(5:end));
            referencePath = fullfile( ...
                gitDirectory, strrep(referenceName, '/', filesep));
            if exist(referencePath, 'file') == 2
                commit = strtrim(fileread(referencePath));
            end
        elseif ~isempty(regexp(head, '^[0-9a-fA-F]{40}$', 'once'))
            commit = lower(head);
        end
        return
    end
    parent = fileparts(directory);
    if strcmp(parent, directory)
        return
    end
    directory = parent;
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

function valid = nf_is_scalar_struct(value)
valid = isstruct(value) && isscalar(value);
end

function options = nf_normalize_faster_options(options)
propertyCount = 3;
allowed = {'measure', 'z'};
if isfield(options, 'stat')
    error('nf_badchans:UnsupportedFASTERStatOption', ...
        ['fasterOptions.stat is not accepted. nf_badchans enforces ' ...
        'the installed min_z measure/z contract; supply only ' ...
        'fasterOptions.measure and fasterOptions.z.']);
end
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

function valid = nf_is_correlation_or_off(value)
valid = nf_is_correlation(value);
if valid
    return
end
valid = nf_is_text(value) && ...
    strcmpi(strtrim(char(value)), 'off');
end

function valid = nf_is_positive_scalar(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value > 0;
end

function valid = nf_is_positive_or_off(value)
valid = nf_is_positive_scalar(value);
if valid
    return
end
valid = nf_is_text(value) && ...
    strcmpi(strtrim(char(value)), 'off');
end

function valid = nf_is_nonnegative_scalar(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value >= 0;
end

function valid = nf_is_positive_integer(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value > 0 && value == round(value);
end

function valid = nf_is_fraction(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value > 0 && value <= 1;
end

function valid = nf_is_threshold(value)
valid = isnumeric(value) && isreal(value) && isvector(value) && ...
    ~isempty(value) && numel(value) <= 2 && ...
    all(isfinite(value(:)));
if valid && numel(value) == 2
    valid = value(1) < value(2);
end
end

function valid = nf_is_ordered_pair(value)
valid = isnumeric(value) && isreal(value) && numel(value) == 2 && ...
    all(isfinite(value(:))) && value(1) < value(2);
end

function valid = nf_is_positive_ordered_pair(value)
valid = nf_is_ordered_pair(value) && all(value(:) > 0);
end

function valid = nf_is_nonnegative_ordered_pair(value)
valid = nf_is_ordered_pair(value) && all(value(:) >= 0);
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
