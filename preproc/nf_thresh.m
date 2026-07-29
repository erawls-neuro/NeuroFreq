function [EEG, rejected, info] = nf_thresh(EEG, voltageThreshold, ...
    powerThreshold, frequencyRange, times, interpolate, maxBadChannels, ...
    frontalChannels, varargin)
% NF_THRESH  Detect, mark, locally repair, and reject artifactual EEG epochs.
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
%
% Name/value options:
%   detectors               Any composable subset of amplitude, fft,
%                           peak2peak, step, gradient, flatline, clipping,
%                           faster, jointprobability, and eeglabstats.
%   amplitudeOptions        thresholdUV, timesSeconds, channelIndices.
%   fftBands                N-by-2 frequencies, N-by-4 frequencies plus
%                           dB limits, or a structure array.
%   fftOptions              bands and channelIndices.
%   peak2peakOptions        thresholdUV, windowSeconds, stepSeconds,
%                           timesSeconds, channelIndices.
%   stepOptions             Same fields as peak2peakOptions.
%   gradientOptions         thresholdUVPerMs, minimumConsecutiveSamples,
%                           timesSeconds, channelIndices.
%   flatlineOptions         maximumRangeUV, windowSeconds, stepSeconds,
%                           timesSeconds, channelIndices.
%   clippingOptions         limitsUV, toleranceUV,
%                           minimumConsecutiveSamples, timesSeconds,
%                           channelIndices.
%   fasterOptions           Global epoch channelIndices and min_z
%                           measure/z; optional FASTER v1.2.4 local
%                           channel detection uses localChannelDetection,
%                           localMeasure, localZ, localChannelIndices,
%                           ignoredChannels, and
%                           exactVendorInterpolation.
%                           Local detection and vendor interpolation default
%                           off; all four local measures default to z=3.
%   jointProbabilityOptions EEGLAB local/global SD thresholds, channels,
%                           superpose, and visualizationType.
%   eeglabStatsOptions      methods (joint probability and/or kurtosis)
%                           plus the EEGLAB statistic settings above.
%   eeglabOptions           Compatibility alias for eeglabStatsOptions.
%   repairOptions           sparseChannelAction, frontalAction, and
%                           excessChannelAction accept interpolate,
%                           rejectepoch, or mark. detectedEpochAction
%                           accepts rejectepoch or mark. epochAction
%                           accepts reject or mark.
%
% Explicit vendor detectors are strict contracts: the installed EEGLAB or
% FASTER functions are called directly and missing dependencies are errors.

nf_validate_epoched_eeg(EEG, 1);

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
addParameter(parser, 'detectors', {'amplitude', 'fft'}, ...
    @nf_is_detector_list);
addParameter(parser, 'amplitudeOptions', struct(), @isstruct);
addParameter(parser, 'fftOptions', struct(), @isstruct);
addParameter(parser, 'fftBands', [], @nf_is_fft_band_input);
addParameter(parser, 'peak2peakOptions', struct(), @isstruct);
addParameter(parser, 'stepOptions', struct(), @isstruct);
addParameter(parser, 'gradientOptions', struct(), @isstruct);
addParameter(parser, 'flatlineOptions', struct(), @isstruct);
addParameter(parser, 'clippingOptions', struct(), @isstruct);
addParameter(parser, 'fasterOptions', struct(), @isstruct);
addParameter(parser, 'jointProbabilityOptions', struct(), @isstruct);
addParameter(parser, 'eeglabStatsOptions', struct(), @isstruct);
addParameter(parser, 'eeglabOptions', struct(), @isstruct);
addParameter(parser, 'repairOptions', struct(), @isstruct);
parse(parser, varargin{:});
options = parser.Results;
options.interpolationMethod = nf_normalize_interpolation_method( ...
    options.interpolationMethod);
options.detectors = nf_normalize_detectors(options.detectors);

if ~isempty(options.detectors)
    minimumTrials = 1;
    populationDetectors = {'faster', 'jointprobability', ...
        'eeglabstats'};
    if any(ismember(options.detectors, populationDetectors))
        minimumTrials = 2;
    end
    nf_validate_epoched_eeg(EEG, minimumTrials);
    nf_validate_settings(EEG, voltageThreshold, powerThreshold, ...
        frequencyRange, times, options.detectors);
    if any(strcmp(options.detectors, 'fft'))
        powerThreshold = reshape(powerThreshold, 1, 2);
        frequencyRange = reshape(frequencyRange, 1, 2);
    end
    timedDetectors = {'amplitude', 'peak2peak', 'step', 'gradient', ...
        'flatline', 'clipping'};
    if any(ismember(options.detectors, timedDetectors))
        times = reshape(times, 1, 2);
    end
end
if isempty(options.detectors)
    detectorSettings = struct();
else
    detectorSettings = nf_resolve_detector_settings(EEG, options, ...
        voltageThreshold, powerThreshold, frequencyRange, times);
end
hasRepairableDetector = nf_has_repairable_detector( ...
    options.detectors, detectorSettings);
if hasRepairableDetector
    if ~nf_is_logical_scalar(interpolate)
        error('nf_thresh:InvalidInterpolationFlag', ...
            'interpolate must be a logical scalar.');
    end
    effectiveInterpolate = logical(interpolate);
    effectiveMaxBadChannels = maxBadChannels;
else
    effectiveInterpolate = false;
    effectiveMaxBadChannels = 0;
end
repairSettings = nf_resolve_repair_settings(options.repairOptions, ...
    effectiveInterpolate);
requiresInterpolation = nf_repair_requires_interpolation( ...
    options.detectors, detectorSettings, repairSettings);
nf_validate_repair_limits(EEG, effectiveMaxBadChannels, ...
    hasRepairableDetector, requiresInterpolation);

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

if ~isempty(options.detectors)
    EEG = nf_clear_threshold_marks(EEG);
end
channelMasks = struct();
epochMasks = struct();
detectorDetails = struct();
provenance = repmat(nf_empty_provenance(), 1, 0);
badChannelEpoch = false(EEG.nbchan, originalTrials);
repairableChannelEpoch = false(EEG.nbchan, originalTrials);
directEpochArtifact = false(1, originalTrials);
for detectorIndex = 1:numel(options.detectors)
    detector = options.detectors{detectorIndex};
    [EEG, channelMask, epochMask, detail, contract] = ...
        nf_run_detector(EEG, detector, detectorSettings);
    if ~isequal(size(channelMask), [EEG.nbchan originalTrials])
        error('nf_thresh:InvalidDetectorMask', ...
            'Detector %s returned an invalid channel-by-epoch mask.', ...
            detector);
    end
    if ~isequal(size(epochMask), [1 originalTrials])
        error('nf_thresh:InvalidDetectorMask', ...
            'Detector %s returned an invalid epoch mask.', detector);
    end
    fieldName = nf_detector_field(detector);
    channelMasks.(fieldName) = logical(channelMask);
    epochMasks.(fieldName) = logical(epochMask);
    detectorDetails.(fieldName) = detail;
    provenance(end + 1) = contract;
    badChannelEpoch = badChannelEpoch | logical(channelMask);
    if nf_detector_channels_are_repairable(detector)
        repairableChannelEpoch = ...
            repairableChannelEpoch | logical(channelMask);
    end
    directEpochArtifact = directEpochArtifact | logical(epochMask);
end
if ~isempty(options.detectors)
    EEG.specdata = [];
end

hasFasterDetector = any(strcmp(options.detectors, 'faster'));
fasterLocalChannelMask = false(EEG.nbchan, originalTrials);
fasterGlobalEpochMask = false(1, originalTrials);
fasterVendorInterpolationRequested = false;
fasterIgnoredChannels = [];
if hasFasterDetector
    fasterField = nf_detector_field('faster');
    fasterLocalChannelMask = channelMasks.(fasterField);
    fasterGlobalEpochMask = epochMasks.(fasterField);
    fasterVendorInterpolationRequested = ...
        detectorSettings.faster.exactVendorInterpolation;
    fasterIgnoredChannels = detectorSettings.faster.ignoredChannels;
end

hasChannelArtifacts = any(repairableChannelEpoch(:));
usesFrontalRule = hasChannelArtifacts && ...
    ~strcmp(repairSettings.frontalAction, ...
    repairSettings.sparseChannelAction);
usesVendorInterpolation = fasterVendorInterpolationRequested && ...
    any(fasterLocalChannelMask(:));
usesInterpolation = hasChannelArtifacts && ...
    (strcmp(repairSettings.sparseChannelAction, 'interpolate') || ...
    strcmp(repairSettings.frontalAction, 'interpolate') || ...
    strcmp(repairSettings.excessChannelAction, 'interpolate'));
usesInterpolation = usesInterpolation || usesVendorInterpolation;
if usesFrontalRule || usesInterpolation
    EEG = nf_normalize_locations(EEG);
    [frontalIndices, frontalLabels, frontalSource] = ...
        nf_resolve_frontal_channels(EEG, frontalChannels);
else
    frontalIndices = [];
    frontalLabels = {};
    frontalSource = 'not-required-by-selected-actions';
end

frontalArtifactRejected = any( ...
    repairableChannelEpoch(frontalIndices, :), 1);
tooManyChannelsRejected = ...
    sum(repairableChannelEpoch, 1) > effectiveMaxBadChannels;
[candidateRejected, selectedActions] = nf_select_epoch_actions( ...
    repairableChannelEpoch, directEpochArtifact, frontalArtifactRejected, ...
    tooManyChannelsRejected, repairSettings);
if strcmp(repairSettings.epochAction, 'reject')
    rejected = candidateRejected;
else
    rejected = false(1, originalTrials);
end

if all(rejected)
    error('nf_thresh:AllRejected', ...
        'All epochs failed the requested artifact criteria.');
end

vendorFasterInterpolationMask = false(EEG.nbchan, originalTrials);
vendorFasterInterpolationCalled = false;
if fasterVendorInterpolationRequested
    badChannelsByEpoch = cell(1, originalTrials);
    for trialIndex = 1:originalTrials
        if rejected(trialIndex)
            badChannelsByEpoch{trialIndex} = [];
            continue
        end
        badChannels = find(fasterLocalChannelMask(:, trialIndex));
        badChannelsByEpoch{trialIndex} = reshape(badChannels, 1, []);
        if isempty(badChannels)
            continue
        end
        donorChannels = setdiff(1:EEG.nbchan, ...
            [badChannelsByEpoch{trialIndex} ...
            reshape(fasterIgnoredChannels, 1, [])]);
        if numel(donorChannels) < 3
            error('nf_thresh:InsufficientFasterInterpolationDonors', ...
                ['Epoch %d has fewer than three FASTER interpolation ' ...
                'donor channels after applying ignoredChannels.'], ...
                trialIndex);
        end
        vendorFasterInterpolationMask(badChannels, trialIndex) = true;
    end
    warningState = warning;
    warningCleanup = onCleanup(@() warning(warningState));
    EEG = h_epoch_interp_spl(EEG, badChannelsByEpoch, ...
        fasterIgnoredChannels);
    clear warningCleanup
    vendorFasterInterpolationCalled = true;
    EEG = nf_eeg_checkset_if_available(EEG);
end

neuroFreqInterpolationMask = false(EEG.nbchan, originalTrials);
for trialIndex = 1:originalTrials
    nativeRepairMask = repairableChannelEpoch(:, trialIndex) & ...
        ~vendorFasterInterpolationMask(:, trialIndex);
    badChannels = find(nativeRepairMask);
    shouldInterpolate = strcmp(selectedActions{trialIndex}, ...
        'interpolate');
    if isempty(badChannels) || rejected(trialIndex) || ...
            ~shouldInterpolate
        continue
    end
    if EEG.nbchan - numel(badChannels) < 3
        error('nf_thresh:InsufficientInterpolationDonors', ...
            ['Epoch %d has fewer than three clean donor channels. ' ...
            'Choose rejectepoch or mark for this artifact class.'], ...
            trialIndex);
    end
    singleEpoch = pop_select(EEG, 'trial', trialIndex, 'sorttrial', 'off');
    singleEpoch = eeg_interp(singleEpoch, badChannels, ...
        char(options.interpolationMethod));
    EEG.data(:, :, trialIndex) = singleEpoch.data;
    neuroFreqInterpolationMask(badChannels, trialIndex) = true;
end
localInterpolationMask = vendorFasterInterpolationMask | ...
    neuroFreqInterpolationMask;
icaInvalidated = false;
if any(localInterpolationMask(:))
    EEG = nf_clear_ica(EEG);
    icaInvalidated = true;
    EEG = eeg_checkset(EEG);
end

if hasFasterDetector
    detectorDetails.faster.vendorInterpolationRequested = ...
        fasterVendorInterpolationRequested;
    detectorDetails.faster.vendorInterpolationCalled = ...
        vendorFasterInterpolationCalled;
    detectorDetails.faster.vendorInterpolationFunction = ...
        'h_epoch_interp_spl';
    detectorDetails.faster.vendorInterpolationMask = ...
        vendorFasterInterpolationMask;
    detectorDetails.faster.vendorInterpolationApplied = ...
        any(vendorFasterInterpolationMask(:));
    detectorDetails.faster.vendorInterpolationScope = ...
        ['Released FASTER interpolation applied only to retained epochs ' ...
        'flagged by its local per-epoch channel detector.'];
end

retainedIndices = find(~rejected);
if hasBehavior
    behavior = nf_subset_trial_metadata(EEG.etc.behavior, ...
        retainedIndices, originalTrials, 'EEG.etc.behavior');
end
retainedEpochIds = nf_subset_trial_metadata(originalEpochIds, ...
    retainedIndices, originalTrials, 'EEG.etc.nf_epoch_ids');

if any(rejected)
    EEG = pop_rejepoch(EEG, rejected, 0);
    EEG = eeg_checkset(EEG);
end
if ~isempty(options.detectors)
    EEG = nf_clear_threshold_marks(EEG);
end
if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
if hasBehavior
    EEG.etc.behavior = behavior;
end
EEG.etc.nf_epoch_ids = retainedEpochIds;

rejectionReason = repmat({''}, 1, originalTrials);
for trialIndex = 1:originalTrials
    if rejected(trialIndex) && directEpochArtifact(trialIndex) && ...
            strcmp(repairSettings.detectedEpochAction, 'rejectepoch')
        rejectionReason{trialIndex} = 'detector-epoch-artifact';
    elseif rejected(trialIndex) && tooManyChannelsRejected(trialIndex) && ...
            strcmp(repairSettings.excessChannelAction, 'rejectepoch')
        rejectionReason{trialIndex} = 'too-many-channels';
    elseif rejected(trialIndex) && frontalArtifactRejected(trialIndex) && ...
            strcmp(repairSettings.frontalAction, 'rejectepoch')
        rejectionReason{trialIndex} = 'frontal-channel-artifact';
    elseif rejected(trialIndex)
        rejectionReason{trialIndex} = 'channel-artifact';
    elseif any(localInterpolationMask(:, trialIndex))
        rejectionReason{trialIndex} = 'retained-local-interpolation';
    elseif candidateRejected(trialIndex)
        rejectionReason{trialIndex} = 'retained-marked-epoch';
    elseif any(badChannelEpoch(:, trialIndex)) || ...
            directEpochArtifact(trialIndex)
        rejectionReason{trialIndex} = 'retained-marked-artifact';
    else
        rejectionReason{trialIndex} = 'retained-clean';
    end
end

info = struct();
info.schemaVersion = '3.0.0';
info.units.voltage = 'microvolts';
info.units.spectralPower = 'dB';
info.criteria.voltageThreshold = voltageThreshold;
info.criteria.powerThreshold = powerThreshold;
info.criteria.frequencyRangeHz = frequencyRange;
info.criteria.timesSeconds = times;
info.criteria.spectralMethod = 'fft';
info.criteria.maxBadChannels = effectiveMaxBadChannels;
info.criteria.requestedMaxBadChannels = maxBadChannels;
info.criteria.channelRepairCapable = hasRepairableDetector;
info.criteria.interpolateSparseNonfrontal = effectiveInterpolate;
info.criteria.interpolationRequiredByResolvedActions = ...
    requiresInterpolation;
if isfield(detectorSettings, 'amplitude')
    info.criteria.resolvedAmplitude = detectorSettings.amplitude;
else
    info.criteria.resolvedAmplitude = struct();
end
if isfield(detectorSettings, 'fft')
    info.criteria.resolvedFFT = detectorSettings.fft;
else
    info.criteria.resolvedFFT = struct();
end
if isempty(options.detectors)
    info.detectors = {'none'};
else
    info.detectors = options.detectors;
end
info.activeDetectors = options.detectors;
info.parameters = detectorSettings;
info.repair = repairSettings;
info.provenance = provenance;
info.detectorDetails = detectorDetails;
info.frontal.indices = frontalIndices;
info.frontal.labels = frontalLabels;
info.frontal.source = frontalSource;
info.masks.byDetectorChannel = channelMasks;
info.masks.byDetectorEpoch = epochMasks;
info.masks.byDetectorAnyEpoch = nf_detector_any_epoch_masks( ...
    channelMasks, epochMasks, originalTrials);
info.masks.voltage = nf_legacy_detector_mask(channelMasks, ...
    'amplitude', EEG.nbchan, originalTrials);
info.masks.spectral = nf_legacy_detector_mask(channelMasks, ...
    'fft', EEG.nbchan, originalTrials);
info.masks.fasterLocalChannel = fasterLocalChannelMask;
info.masks.fasterGlobalEpoch = fasterGlobalEpochMask;
info.masks.frontalVoltageRejected = any( ...
    info.masks.voltage(frontalIndices, :), 1);
info.masks.frontalSpectralRejected = any( ...
    info.masks.spectral(frontalIndices, :), 1);
info.masks.anyArtifact = badChannelEpoch;
info.masks.repairableChannelArtifact = repairableChannelEpoch;
info.masks.directEpochArtifact = directEpochArtifact;
info.masks.localInterpolation = localInterpolationMask;
info.masks.vendorFasterInterpolation = vendorFasterInterpolationMask;
info.masks.neuroFreqInterpolation = neuroFreqInterpolationMask;
info.masks.rejected = rejected;
info.masks.candidateRejected = candidateRejected;
info.masks.frontalArtifactRejected = frontalArtifactRejected;
info.masks.tooManyChannelsRejected = tooManyChannelsRejected;
info.rejectionReason = rejectionReason;
info.selectedActions = selectedActions;
info.localBadCount = sum(badChannelEpoch, 1);
info.repairableLocalBadCount = sum(repairableChannelEpoch, 1);
info.originalEpochIds = originalEpochIds;
info.retainedIndices = retainedIndices;
info.retainedEpochIds = retainedEpochIds;
info.nOriginal = originalTrials;
info.nRejected = sum(rejected);
info.nRetained = numel(retainedIndices);
info.nCandidateRejected = sum(candidateRejected);
info.nArtifactEpochs = sum(any(badChannelEpoch, 1) | ...
    directEpochArtifact);
info.nMarkedArtifactEpochs = sum((any(badChannelEpoch, 1) | ...
    directEpochArtifact) & ~rejected & ...
    ~any(localInterpolationMask, 1));
info.nLocallyRepairedEpochs = sum(any(localInterpolationMask, 1));
info.nLocallyRepairedChannelEpochs = sum(localInterpolationMask(:));
info.interpolationMethod = char(options.interpolationMethod);
info.vendorFasterInterpolation.requested = ...
    fasterVendorInterpolationRequested;
info.vendorFasterInterpolation.called = vendorFasterInterpolationCalled;
info.vendorFasterInterpolation.applied = ...
    any(vendorFasterInterpolationMask(:));
info.vendorFasterInterpolation.function = 'h_epoch_interp_spl';
info.vendorFasterInterpolation.ignoredChannels = fasterIgnoredChannels;
info.vendorFasterInterpolation.fullPipelineEquivalent = false;
info.vendorFasterInterpolation.scope = ...
    ['Exact released FASTER local interpolation component inside the ' ...
    'NeuroFreq composition; not a claim of full FASTER pipeline order.'];
info.icaDecompositionInvalidated = icaInvalidated;
EEG.etc.nf_thresh = info;

end

function settings = nf_resolve_detector_settings(EEG, options, ...
    voltageThreshold, powerThreshold, frequencyRange, times)
allChannels = 1:EEG.nbchan;
settings = struct();

if any(strcmp(options.detectors, 'amplitude'))
defaults = struct();
defaults.thresholdUV = voltageThreshold;
defaults.timesSeconds = times;
defaults.channelIndices = allChannels;
settings.amplitude = nf_merge_known_options(defaults, ...
    options.amplitudeOptions, 'amplitudeOptions');
nf_validate_amplitude_options(EEG, settings.amplitude);
end

if any(strcmp(options.detectors, 'fft'))
defaults = struct();
defaults.bands = nf_normalize_fft_bands([], frequencyRange, ...
    powerThreshold, EEG.srate);
defaults.channelIndices = allChannels;
if ~isempty(options.fftBands) && isfield(options.fftOptions, 'bands')
    error('nf_thresh:ConflictingFFTBands', ...
        ['Specify FFT bands with either fftBands or fftOptions.bands, ' ...
        'not both.']);
end
settings.fft = nf_merge_known_options(defaults, options.fftOptions, ...
    'fftOptions');
if ~isempty(options.fftBands)
    settings.fft.bands = options.fftBands;
end
settings.fft.bands = nf_normalize_fft_bands(settings.fft.bands, ...
    frequencyRange, powerThreshold, EEG.srate);
nf_validate_channel_indices(settings.fft.channelIndices, EEG.nbchan, ...
    'fftOptions.channelIndices');
end

if any(strcmp(options.detectors, 'peak2peak'))
defaults = struct();
defaults.thresholdUV = 200;
defaults.windowSeconds = min(0.200, diff(times));
defaults.stepSeconds = min(0.100, defaults.windowSeconds);
defaults.timesSeconds = times;
defaults.channelIndices = allChannels;
settings.peak2peak = nf_merge_known_options(defaults, ...
    options.peak2peakOptions, 'peak2peakOptions');
nf_validate_window_detector_options(EEG, settings.peak2peak, ...
    'peak2peakOptions');
if ~nf_is_positive_scalar(settings.peak2peak.thresholdUV)
    error('nf_thresh:InvalidPeakToPeakOptions', ...
        'peak2peakOptions.thresholdUV must be positive.');
end
end

if any(strcmp(options.detectors, 'step'))
defaults = struct();
defaults.thresholdUV = 100;
defaults.windowSeconds = min(0.100, diff(times) / 2);
defaults.stepSeconds = min(0.050, defaults.windowSeconds);
defaults.timesSeconds = times;
defaults.channelIndices = allChannels;
settings.step = nf_merge_known_options(defaults, options.stepOptions, ...
    'stepOptions');
nf_validate_window_detector_options(EEG, settings.step, 'stepOptions');
if ~nf_is_positive_scalar(settings.step.thresholdUV)
    error('nf_thresh:InvalidStepOptions', ...
        'stepOptions.thresholdUV must be positive.');
end
end

if any(strcmp(options.detectors, 'gradient'))
defaults = struct();
defaults.thresholdUVPerMs = 50;
defaults.minimumConsecutiveSamples = 1;
defaults.timesSeconds = times;
defaults.channelIndices = allChannels;
settings.gradient = nf_merge_known_options(defaults, ...
    options.gradientOptions, 'gradientOptions');
nf_validate_point_detector_options(EEG, settings.gradient, ...
    'gradientOptions');
if ~nf_is_positive_scalar(settings.gradient.thresholdUVPerMs)
    error('nf_thresh:InvalidGradientOptions', ...
        'gradientOptions.thresholdUVPerMs must be positive.');
end
if ~nf_is_positive_integer(settings.gradient.minimumConsecutiveSamples)
    error('nf_thresh:InvalidGradientOptions', ...
        ['gradientOptions.minimumConsecutiveSamples must be a positive ' ...
        'integer.']);
end
end

if any(strcmp(options.detectors, 'flatline'))
defaults = struct();
defaults.maximumRangeUV = 1;
defaults.windowSeconds = min(0.200, diff(times));
defaults.stepSeconds = min(0.050, defaults.windowSeconds);
defaults.timesSeconds = times;
defaults.channelIndices = allChannels;
settings.flatline = nf_merge_known_options(defaults, ...
    options.flatlineOptions, 'flatlineOptions');
nf_validate_window_detector_options(EEG, settings.flatline, ...
    'flatlineOptions');
if ~nf_is_nonnegative_scalar(settings.flatline.maximumRangeUV)
    error('nf_thresh:InvalidFlatlineOptions', ...
        'flatlineOptions.maximumRangeUV must be nonnegative.');
end
end

if any(strcmp(options.detectors, 'clipping'))
defaults = struct();
defaults.limitsUV = [-500 500];
defaults.toleranceUV = 0.5;
defaults.minimumConsecutiveSamples = 2;
defaults.timesSeconds = times;
defaults.channelIndices = allChannels;
settings.clipping = nf_merge_known_options(defaults, ...
    options.clippingOptions, 'clippingOptions');
nf_validate_point_detector_options(EEG, settings.clipping, ...
    'clippingOptions');
if ~nf_is_increasing_pair(settings.clipping.limitsUV)
    error('nf_thresh:InvalidClippingOptions', ...
        'clippingOptions.limitsUV must be a finite increasing pair.');
end
if ~nf_is_nonnegative_scalar(settings.clipping.toleranceUV)
    error('nf_thresh:InvalidClippingOptions', ...
        'clippingOptions.toleranceUV must be nonnegative.');
end
if settings.clipping.toleranceUV >= ...
        diff(settings.clipping.limitsUV) / 2
    error('nf_thresh:InvalidClippingOptions', ...
        'clippingOptions.toleranceUV is too large for limitsUV.');
end
if ~nf_is_positive_integer(settings.clipping.minimumConsecutiveSamples)
    error('nf_thresh:InvalidClippingOptions', ...
        ['clippingOptions.minimumConsecutiveSamples must be a positive ' ...
        'integer.']);
end
end

if any(strcmp(options.detectors, 'faster'))
defaults = struct();
defaults.channelIndices = allChannels;
defaults.measure = [];
defaults.z = 3;
defaults.localChannelDetection = false;
defaults.localMeasure = [1 1 1 1];
defaults.localZ = [3 3 3 3];
defaults.localChannelIndices = allChannels;
defaults.ignoredChannels = [];
defaults.exactVendorInterpolation = false;
settings.faster = nf_merge_known_options(defaults, ...
    options.fasterOptions, 'fasterOptions');
nf_validate_channel_indices(settings.faster.channelIndices, EEG.nbchan, ...
    'fasterOptions.channelIndices');
if ~(isempty(settings.faster.measure) || ...
        ((isnumeric(settings.faster.measure) || ...
        islogical(settings.faster.measure)) && ...
        isvector(settings.faster.measure) && ...
        all(isfinite(double(settings.faster.measure(:)))) && ...
        all(ismember(double(settings.faster.measure(:)), [0 1]))))
    error('nf_thresh:InvalidFasterOptions', ...
        'fasterOptions.measure must be empty, numeric, or logical.');
end
if ~(isnumeric(settings.faster.z) && isreal(settings.faster.z) && ...
        ~isempty(settings.faster.z) && ...
        all(isfinite(settings.faster.z(:))) && ...
        all(settings.faster.z(:) > 0))
    error('nf_thresh:InvalidFasterOptions', ...
        'fasterOptions.z must contain positive finite values.');
end
if ~nf_is_logical_scalar(settings.faster.localChannelDetection)
    error('nf_thresh:InvalidFasterOptions', ...
        'fasterOptions.localChannelDetection must be binary.');
end
settings.faster.localChannelDetection = ...
    logical(settings.faster.localChannelDetection);
nf_validate_channel_indices(settings.faster.localChannelIndices, ...
    EEG.nbchan, 'fasterOptions.localChannelIndices');
settings.faster.localChannelIndices = ...
    reshape(settings.faster.localChannelIndices, 1, []);
if ~isempty(settings.faster.ignoredChannels)
    nf_validate_channel_indices(settings.faster.ignoredChannels, ...
        EEG.nbchan, 'fasterOptions.ignoredChannels');
end
settings.faster.ignoredChannels = ...
    reshape(settings.faster.ignoredChannels, 1, []);
if ~((isnumeric(settings.faster.localMeasure) || ...
        islogical(settings.faster.localMeasure)) && ...
        isvector(settings.faster.localMeasure) && ...
        ~isempty(settings.faster.localMeasure) && ...
        all(isfinite(double(settings.faster.localMeasure(:)))))
    error('nf_thresh:InvalidFasterOptions', ...
        ['fasterOptions.localMeasure must be a finite numeric or logical ' ...
        'scalar or four-element vector.']);
end
settings.faster.localMeasure = nf_expand_setting_vector( ...
    settings.faster.localMeasure, 4, 'fasterOptions.localMeasure');
if ~all(ismember(double(settings.faster.localMeasure), [0 1]))
    error('nf_thresh:InvalidFasterOptions', ...
        ['fasterOptions.localMeasure must be binary and have one value ' ...
        'for each of the four FASTER local channel measures.']);
end
settings.faster.localMeasure = logical(settings.faster.localMeasure);
if ~(isnumeric(settings.faster.localZ) && ...
        isreal(settings.faster.localZ) && ...
        isvector(settings.faster.localZ) && ...
        ~isempty(settings.faster.localZ) && ...
        all(isfinite(settings.faster.localZ(:))))
    error('nf_thresh:InvalidFasterOptions', ...
        ['fasterOptions.localZ must be a finite numeric scalar or ' ...
        'four-element vector.']);
end
settings.faster.localZ = nf_expand_setting_vector( ...
    settings.faster.localZ, 4, 'fasterOptions.localZ');
if ~(isnumeric(settings.faster.localZ) && ...
        isreal(settings.faster.localZ) && ...
        all(isfinite(settings.faster.localZ)) && ...
        all(settings.faster.localZ > 0))
    error('nf_thresh:InvalidFasterOptions', ...
        ['fasterOptions.localZ must contain four positive finite ' ...
        'thresholds.']);
end
if ~nf_is_logical_scalar(settings.faster.exactVendorInterpolation)
    error('nf_thresh:InvalidFasterOptions', ...
        'fasterOptions.exactVendorInterpolation must be binary.');
end
settings.faster.exactVendorInterpolation = ...
    logical(settings.faster.exactVendorInterpolation);
if settings.faster.exactVendorInterpolation && ...
        ~settings.faster.localChannelDetection
    error('nf_thresh:InvalidFasterOptions', ...
        ['fasterOptions.exactVendorInterpolation requires ' ...
        'localChannelDetection=true.']);
end
end

if any(strcmp(options.detectors, 'jointprobability'))
defaults = struct();
defaults.channelIndices = allChannels;
defaults.localThresholdSD = 5;
defaults.globalThresholdSD = 5;
defaults.superpose = false;
defaults.visualizationType = 0;
settings.jointprobability = nf_merge_known_options(defaults, ...
    options.jointProbabilityOptions, 'jointProbabilityOptions');
nf_validate_eeglab_stat_options(EEG, settings.jointprobability, ...
    'jointProbabilityOptions');
end

if any(strcmp(options.detectors, 'eeglabstats'))
defaults = struct();
defaults.channelIndices = allChannels;
defaults.methods = {'jointprobability', 'kurtosis'};
defaults.localThresholdSD = 5;
defaults.globalThresholdSD = 5;
defaults.superpose = false;
defaults.visualizationType = 0;
eeglabStatsOverrides = nf_merge_override_aliases( ...
    options.eeglabStatsOptions, options.eeglabOptions, ...
    'eeglabStatsOptions', 'eeglabOptions');
settings.eeglabstats = nf_merge_known_options(defaults, ...
    eeglabStatsOverrides, 'eeglabStatsOptions');
settings.eeglabstats.methods = nf_normalize_eeglab_stat_methods( ...
    settings.eeglabstats.methods);
nf_validate_eeglab_stat_options(EEG, settings.eeglabstats, ...
    'eeglabStatsOptions');
end
end

function repair = nf_resolve_repair_settings(overrides, interpolate)
defaults = struct();
if interpolate
    defaults.sparseChannelAction = 'interpolate';
else
    defaults.sparseChannelAction = 'rejectepoch';
end
defaults.frontalAction = 'rejectepoch';
defaults.excessChannelAction = 'rejectepoch';
defaults.detectedEpochAction = 'rejectepoch';
defaults.epochAction = 'reject';
repair = nf_merge_known_options(defaults, overrides, 'repairOptions');
repair.sparseChannelAction = nf_normalize_channel_action( ...
    repair.sparseChannelAction, 'repairOptions.sparseChannelAction');
repair.frontalAction = nf_normalize_channel_action( ...
    repair.frontalAction, 'repairOptions.frontalAction');
repair.excessChannelAction = nf_normalize_channel_action( ...
    repair.excessChannelAction, 'repairOptions.excessChannelAction');
repair.detectedEpochAction = nf_normalize_detected_epoch_action( ...
    repair.detectedEpochAction);
repair.epochAction = nf_normalize_epoch_action(repair.epochAction);
end

function [EEG, channelMask, epochMask, detail, contract] = ...
        nf_run_detector(EEG, detector, settings)
switch detector
    case 'amplitude'
        contract = nf_vendor_contract(detector, 'EEGLAB', ...
            {'pop_eegthresh'}, 'vendor-exact');
        [EEG, channelMask, epochMask, detail] = ...
            nf_detect_amplitude(EEG, settings.amplitude);
    case 'fft'
        contract = nf_vendor_contract(detector, 'EEGLAB', ...
            {'pop_rejspec'}, 'vendor-exact');
        [EEG, channelMask, epochMask, detail] = ...
            nf_detect_fft(EEG, settings.fft);
    case 'peak2peak'
        [channelMask, detail] = nf_detect_peak2peak(EEG, ...
            settings.peak2peak);
        epochMask = false(1, EEG.trials);
        contract = nf_native_contract(detector, ...
            'moving-window peak-to-peak');
    case 'step'
        [channelMask, detail] = nf_detect_step(EEG, settings.step);
        epochMask = false(1, EEG.trials);
        contract = nf_native_contract(detector, ...
            'adjacent-window mean step');
    case 'gradient'
        [channelMask, detail] = nf_detect_gradient(EEG, ...
            settings.gradient);
        epochMask = false(1, EEG.trials);
        contract = nf_native_contract(detector, ...
            'sample gradient in microvolts per millisecond');
    case 'flatline'
        [channelMask, detail] = nf_detect_flatline(EEG, ...
            settings.flatline);
        epochMask = false(1, EEG.trials);
        contract = nf_native_contract(detector, ...
            'moving-window signal range');
    case 'clipping'
        [channelMask, detail] = nf_detect_clipping(EEG, ...
            settings.clipping);
        epochMask = false(1, EEG.trials);
        contract = nf_native_contract(detector, ...
            'user-defined acquisition-rail occupancy');
    case 'faster'
        fasterFunctions = nf_faster_function_names(settings.faster);
        contract = nf_vendor_contract(detector, 'FASTER', ...
            fasterFunctions, 'vendor-exact-components');
        [channelMask, epochMask, detail] = nf_detect_faster(EEG, ...
            settings.faster);
        contract.definition = detail.contractScope;
    case 'jointprobability'
        contract = nf_vendor_contract(detector, 'EEGLAB', ...
            {'pop_jointprob'}, 'vendor-exact');
        [EEG, channelMask, epochMask, detail] = ...
            nf_detect_joint_probability(EEG, ...
            settings.jointprobability);
    case 'eeglabstats'
        functionNames = nf_eeglab_stat_function_names( ...
            settings.eeglabstats.methods);
        contract = nf_vendor_contract(detector, 'EEGLAB', ...
            functionNames, 'vendor-exact');
        [EEG, channelMask, epochMask, detail, functionNames] = ...
            nf_detect_eeglab_stats(EEG, settings.eeglabstats);
    otherwise
        error('nf_thresh:UnknownDetector', ...
            'Unknown artifact detector: %s.', detector);
end
channelMask = logical(channelMask);
epochMask = reshape(logical(epochMask), 1, EEG.trials);
end

function repairable = nf_detector_channels_are_repairable(detector)
repairable = ismember(detector, {'amplitude', 'fft', 'peak2peak', ...
    'step', 'gradient', 'flatline', 'clipping', 'faster'});
end

function capable = nf_has_repairable_detector(detectors, settings)
nativeDetectors = {'amplitude', 'fft', 'peak2peak', 'step', ...
    'gradient', 'flatline', 'clipping'};
capable = any(ismember(detectors, nativeDetectors));
if capable || ~any(strcmp(detectors, 'faster')) || ...
        ~isfield(settings, 'faster')
    return
end
capable = logical(settings.faster.localChannelDetection);
end

function required = nf_repair_requires_interpolation( ...
        detectors, settings, repair)
if ~nf_has_repairable_detector(detectors, settings)
    required = false;
    return
end
actions = {repair.sparseChannelAction, repair.frontalAction, ...
    repair.excessChannelAction};
required = any(strcmp(actions, 'interpolate'));
if any(strcmp(detectors, 'faster')) && isfield(settings, 'faster')
    required = required || ...
        logical(settings.faster.exactVendorInterpolation);
end
end

function [EEG, mask, epochMask, detail] = nf_detect_amplitude(EEG, ...
    settings)
nf_require_functions({'pop_eegthresh'}, 'amplitude', 'EEGLAB');
EEG = nf_clear_rejection_fields(EEG, {'rejthresh', 'rejthreshE'});
EEG = pop_eegthresh(EEG, 1, settings.channelIndices, ...
    -settings.thresholdUV, settings.thresholdUV, ...
    settings.timesSeconds(1), settings.timesSeconds(2), 0, 0);
EEG = nf_eeg_checkset_if_available(EEG);
mask = nf_get_channel_rejection_matrix(EEG, 'rejthreshE', ...
    settings.channelIndices);
epochMask = false(1, EEG.trials);
detail = struct();
detail.thresholdUV = settings.thresholdUV;
detail.timesSeconds = settings.timesSeconds;
detail.channelIndices = settings.channelIndices;
detail.nFlaggedChannelEpochs = sum(mask(:));
end

function [EEG, mask, epochMask, detail] = nf_detect_fft(EEG, settings)
nf_require_functions({'pop_rejspec'}, 'fft', 'EEGLAB');
numberBands = numel(settings.bands);
bandMasks = false(EEG.nbchan, EEG.trials, numberBands);
for bandIndex = 1:numberBands
    band = settings.bands(bandIndex);
    EEG = nf_clear_rejection_fields(EEG, {'rejfreq', 'rejfreqE'});
    EEG.specdata = [];
    EEG = pop_rejspec(EEG, 1, 'elecrange', ...
        settings.channelIndices, 'method', 'fft', 'threshold', ...
        band.powerThresholdDb, 'freqlimits', ...
        band.frequencyRangeHz, 'specdata', [], ...
        'eegplotplotallrej', 0, 'eegplotreject', 0);
    EEG = nf_eeg_checkset_if_available(EEG);
    bandMasks(:, :, bandIndex) = ...
        nf_get_channel_rejection_matrix(EEG, 'rejfreqE', ...
        settings.channelIndices);
end
mask = any(bandMasks, 3);
epochMask = false(1, EEG.trials);
EEG.specdata = [];
detail = struct();
detail.bands = settings.bands;
detail.channelIndices = settings.channelIndices;
detail.bandMasks = bandMasks;
detail.nFlaggedChannelEpochs = sum(mask(:));
end

function [mask, detail] = nf_detect_peak2peak(EEG, settings)
[data, selectedTimes] = nf_selected_data(EEG, ...
    settings.channelIndices, settings.timesSeconds);
windowSamples = nf_seconds_to_samples(settings.windowSeconds, ...
    EEG.srate, 2);
stepSamples = nf_seconds_to_samples(settings.stepSeconds, EEG.srate, 1);
starts = nf_window_starts(size(data, 2), windowSamples, stepSamples);
localMask = false(numel(settings.channelIndices), EEG.trials);
maximumPeakToPeak = zeros(numel(settings.channelIndices), EEG.trials);
for startIndex = 1:numel(starts)
    indices = starts(startIndex):(starts(startIndex) + ...
        windowSamples - 1);
    windowData = data(:, indices, :);
    ranges = squeeze(max(windowData, [], 2) - min(windowData, [], 2));
    ranges = nf_force_channel_trial_matrix(ranges, ...
        numel(settings.channelIndices), EEG.trials);
    maximumPeakToPeak = max(maximumPeakToPeak, ranges);
    localMask = localMask | ranges > settings.thresholdUV;
end
mask = nf_expand_local_mask(localMask, settings.channelIndices, ...
    EEG.nbchan, EEG.trials);
detail = struct();
detail.definition = ...
    'Maximum signal range within any moving window.';
detail.thresholdUV = settings.thresholdUV;
detail.windowSeconds = settings.windowSeconds;
detail.stepSeconds = settings.stepSeconds;
detail.timesSeconds = selectedTimes;
detail.maximumPeakToPeakUV = maximumPeakToPeak;
detail.nFlaggedChannelEpochs = sum(mask(:));
end

function [mask, detail] = nf_detect_step(EEG, settings)
[data, selectedTimes] = nf_selected_data(EEG, ...
    settings.channelIndices, settings.timesSeconds);
windowSamples = nf_seconds_to_samples(settings.windowSeconds, ...
    EEG.srate, 1);
stepSamples = nf_seconds_to_samples(settings.stepSeconds, EEG.srate, 1);
firstBoundary = windowSamples + 1;
lastBoundary = size(data, 2) - windowSamples + 1;
boundaries = firstBoundary:stepSamples:lastBoundary;
if isempty(boundaries)
    error('nf_thresh:DetectorWindowTooLong', ...
        ['stepOptions.windowSeconds requires two adjacent windows ' ...
        'inside the selected time range.']);
end
localMask = false(numel(settings.channelIndices), EEG.trials);
maximumStep = zeros(numel(settings.channelIndices), EEG.trials);
for boundaryIndex = 1:numel(boundaries)
    boundary = boundaries(boundaryIndex);
    before = data(:, boundary - windowSamples:boundary - 1, :);
    after = data(:, boundary:boundary + windowSamples - 1, :);
    difference = squeeze(abs(mean(after, 2) - mean(before, 2)));
    difference = nf_force_channel_trial_matrix(difference, ...
        numel(settings.channelIndices), EEG.trials);
    maximumStep = max(maximumStep, difference);
    localMask = localMask | difference > settings.thresholdUV;
end
mask = nf_expand_local_mask(localMask, settings.channelIndices, ...
    EEG.nbchan, EEG.trials);
detail = struct();
detail.definition = ...
    'Maximum absolute difference between adjacent window means.';
detail.thresholdUV = settings.thresholdUV;
detail.windowSeconds = settings.windowSeconds;
detail.stepSeconds = settings.stepSeconds;
detail.timesSeconds = selectedTimes;
detail.maximumStepUV = maximumStep;
detail.nFlaggedChannelEpochs = sum(mask(:));
end

function [mask, detail] = nf_detect_gradient(EEG, settings)
[data, selectedTimes] = nf_selected_data(EEG, ...
    settings.channelIndices, settings.timesSeconds);
sampleIntervalMs = 1000 / EEG.srate;
gradient = abs(diff(data, 1, 2)) / sampleIntervalMs;
above = gradient > settings.thresholdUVPerMs;
localMask = false(numel(settings.channelIndices), EEG.trials);
for channelIndex = 1:size(above, 1)
    for trialIndex = 1:EEG.trials
        values = reshape(above(channelIndex, :, trialIndex), 1, []);
        localMask(channelIndex, trialIndex) = nf_has_run(values, ...
            settings.minimumConsecutiveSamples);
    end
end
maximumGradient = squeeze(max(gradient, [], 2));
maximumGradient = nf_force_channel_trial_matrix(maximumGradient, ...
    numel(settings.channelIndices), EEG.trials);
mask = nf_expand_local_mask(localMask, settings.channelIndices, ...
    EEG.nbchan, EEG.trials);
detail = struct();
detail.definition = ...
    'Absolute first temporal derivative in microvolts per millisecond.';
detail.thresholdUVPerMs = settings.thresholdUVPerMs;
detail.minimumConsecutiveSamples = ...
    settings.minimumConsecutiveSamples;
detail.timesSeconds = selectedTimes;
detail.maximumGradientUVPerMs = maximumGradient;
detail.nFlaggedChannelEpochs = sum(mask(:));
end

function [mask, detail] = nf_detect_flatline(EEG, settings)
[data, selectedTimes] = nf_selected_data(EEG, ...
    settings.channelIndices, settings.timesSeconds);
windowSamples = nf_seconds_to_samples(settings.windowSeconds, ...
    EEG.srate, 2);
stepSamples = nf_seconds_to_samples(settings.stepSeconds, EEG.srate, 1);
starts = nf_window_starts(size(data, 2), windowSamples, stepSamples);
localMask = false(numel(settings.channelIndices), EEG.trials);
minimumRange = inf(numel(settings.channelIndices), EEG.trials);
for startIndex = 1:numel(starts)
    indices = starts(startIndex):(starts(startIndex) + ...
        windowSamples - 1);
    windowData = data(:, indices, :);
    ranges = squeeze(max(windowData, [], 2) - min(windowData, [], 2));
    ranges = nf_force_channel_trial_matrix(ranges, ...
        numel(settings.channelIndices), EEG.trials);
    minimumRange = min(minimumRange, ranges);
    localMask = localMask | ranges <= settings.maximumRangeUV;
end
mask = nf_expand_local_mask(localMask, settings.channelIndices, ...
    EEG.nbchan, EEG.trials);
detail = struct();
detail.definition = ...
    'Signal range at or below the limit within any moving window.';
detail.maximumRangeUV = settings.maximumRangeUV;
detail.windowSeconds = settings.windowSeconds;
detail.stepSeconds = settings.stepSeconds;
detail.timesSeconds = selectedTimes;
detail.minimumWindowRangeUV = minimumRange;
detail.nFlaggedChannelEpochs = sum(mask(:));
end

function [mask, detail] = nf_detect_clipping(EEG, settings)
[data, selectedTimes] = nf_selected_data(EEG, ...
    settings.channelIndices, settings.timesSeconds);
lower = settings.limitsUV(1) + settings.toleranceUV;
upper = settings.limitsUV(2) - settings.toleranceUV;
atRail = data <= lower | data >= upper;
localMask = false(numel(settings.channelIndices), EEG.trials);
for channelIndex = 1:size(atRail, 1)
    for trialIndex = 1:EEG.trials
        values = reshape(atRail(channelIndex, :, trialIndex), 1, []);
        localMask(channelIndex, trialIndex) = nf_has_run(values, ...
            settings.minimumConsecutiveSamples);
    end
end
mask = nf_expand_local_mask(localMask, settings.channelIndices, ...
    EEG.nbchan, EEG.trials);
detail = struct();
detail.definition = ...
    'Consecutive samples at or beyond user-defined acquisition rails.';
detail.limitsUV = settings.limitsUV;
detail.toleranceUV = settings.toleranceUV;
detail.minimumConsecutiveSamples = ...
    settings.minimumConsecutiveSamples;
detail.timesSeconds = selectedTimes;
detail.nFlaggedChannelEpochs = sum(mask(:));
end

function [channelMask, epochMask, detail] = nf_detect_faster(EEG, ...
    settings)
functionNames = nf_faster_function_names(settings);
nf_require_functions(functionNames, 'faster', 'FASTER');
properties = epoch_properties(EEG, settings.channelIndices);
if size(properties, 1) ~= EEG.trials && ...
        size(properties, 2) == EEG.trials
    properties = properties';
end
if size(properties, 1) ~= EEG.trials
    error('nf_thresh:InvalidFasterOutput', ...
        ['FASTER epoch_properties returned %d rows for %d epochs. ' ...
        'The installed FASTER release is not API-compatible.'], ...
        size(properties, 1), EEG.trials);
end
numberMeasures = size(properties, 2);
if isempty(settings.measure)
    measure = true(1, numberMeasures);
else
    measure = nf_expand_setting_vector(settings.measure, ...
        numberMeasures, 'fasterOptions.measure');
    measure = logical(measure);
end
z = nf_expand_setting_vector(settings.z, numberMeasures, ...
    'fasterOptions.z');
minZOptions = struct();
minZOptions.measure = measure;
minZOptions.z = z;
decisions = min_z(properties, minZOptions);
if numel(decisions) ~= EEG.trials
    error('nf_thresh:InvalidFasterOutput', ...
        'FASTER min_z returned %d decisions for %d epochs.', ...
        numel(decisions), EEG.trials);
end
epochMask = reshape(logical(decisions), 1, EEG.trials);
channelMask = false(EEG.nbchan, EEG.trials);
localProperties = cell(1, EEG.trials);
localMinZOptions = struct();
localMinZOptions.measure = settings.localMeasure;
localMinZOptions.z = settings.localZ;
if settings.localChannelDetection
    localDecisions = false(numel(settings.localChannelIndices), ...
        EEG.trials);
    for trialIndex = 1:EEG.trials
        trialProperties = single_epoch_channel_properties(EEG, ...
            trialIndex, settings.localChannelIndices);
        expectedSize = [numel(settings.localChannelIndices) 4];
        if ~isequal(size(trialProperties), expectedSize)
            error('nf_thresh:InvalidFasterLocalOutput', ...
                ['FASTER single_epoch_channel_properties returned ' ...
                '%d-by-%d for %d selected channels; v1.2.4 requires ' ...
                'one row per channel and four property columns.'], ...
                size(trialProperties, 1), size(trialProperties, 2), ...
                numel(settings.localChannelIndices));
        end
        trialDecisions = min_z(trialProperties, localMinZOptions);
        if numel(trialDecisions) ~= ...
                numel(settings.localChannelIndices)
            error('nf_thresh:InvalidFasterLocalOutput', ...
                ['FASTER min_z returned %d local decisions for %d ' ...
                'selected channels in epoch %d.'], ...
                numel(trialDecisions), ...
                numel(settings.localChannelIndices), trialIndex);
        end
        localDecisions(:, trialIndex) = ...
            reshape(logical(trialDecisions), [], 1);
        localProperties{trialIndex} = trialProperties;
    end
    channelMask(settings.localChannelIndices, :) = localDecisions;
end
detail = struct();
detail.channelIndices = settings.channelIndices;
detail.properties = properties;
detail.minZOptions = minZOptions;
detail.nFlaggedEpochs = sum(epochMask);
detail.globalEpochMask = epochMask;
detail.localChannelDetection = settings.localChannelDetection;
detail.localChannelIndices = settings.localChannelIndices;
detail.localProperties = localProperties;
detail.localPropertyNames = {'median-gradient', 'variance', ...
    'amplitude-range', 'deviation-from-channel-mean'};
detail.localMinZOptions = localMinZOptions;
detail.localChannelMask = channelMask;
detail.nFlaggedLocalChannelEpochs = sum(channelMask(:));
detail.exactVendorInterpolationRequested = ...
    settings.exactVendorInterpolation;
detail.fullPipelineEquivalent = false;
detail.contractScope = ...
    ['Released FASTER global epoch detector and optional v1.2.4 local ' ...
    'per-epoch channel detector called directly. These exact components ' ...
    'are embedded in NeuroFreq and do not claim full FASTER pipeline ' ...
    'ordering or equivalence.'];
end

function functionNames = nf_faster_function_names(settings)
functionNames = {'epoch_properties', 'min_z'};
if settings.localChannelDetection
    functionNames{end + 1} = 'single_epoch_channel_properties';
end
if settings.exactVendorInterpolation
    functionNames{end + 1} = 'h_epoch_interp_spl';
end
end

function [EEG, channelMask, epochMask, detail] = ...
        nf_detect_joint_probability(EEG, settings)
nf_require_functions({'pop_jointprob'}, 'jointprobability', 'EEGLAB');
EEG = nf_clear_rejection_fields(EEG, {'rejjp', 'rejjpE'});
EEG = nf_clear_stat_fields(EEG, {'jp', 'jpE'});
EEG = pop_jointprob(EEG, 1, settings.channelIndices, ...
    settings.localThresholdSD, settings.globalThresholdSD, ...
    double(settings.superpose), 0, settings.visualizationType);
EEG = nf_eeg_checkset_if_available(EEG);
channelMask = nf_get_channel_rejection_matrix(EEG, 'rejjpE', ...
    settings.channelIndices);
globalMask = nf_get_epoch_rejection_vector(EEG, 'rejjp', false);
epochMask = globalMask | any(channelMask, 1);
detail = struct();
detail.localThresholdSD = settings.localThresholdSD;
detail.globalThresholdSD = settings.globalThresholdSD;
detail.channelIndices = settings.channelIndices;
detail.superpose = logical(settings.superpose);
detail.visualizationType = settings.visualizationType;
detail.globalMask = globalMask;
detail.nFlaggedEpochs = sum(epochMask);
end

function [EEG, channelMask, epochMask, detail, functionNames] = ...
        nf_detect_eeglab_stats(EEG, settings)
channelMask = false(EEG.nbchan, EEG.trials);
epochMask = false(1, EEG.trials);
detail = struct();
functionNames = {};
for methodIndex = 1:numel(settings.methods)
    method = settings.methods{methodIndex};
    switch method
        case 'jointprobability'
            [EEG, methodChannelMask, methodEpochMask, methodDetail] = ...
                nf_detect_joint_probability(EEG, settings);
            functionNames{end + 1} = 'pop_jointprob';
        case 'kurtosis'
            nf_require_functions({'pop_rejkurt'}, 'eeglabstats', ...
                'EEGLAB');
            EEG = nf_clear_rejection_fields(EEG, ...
                {'rejkurt', 'rejkurtE'});
            EEG = nf_clear_stat_fields(EEG, {'kurt', 'kurtE'});
            EEG = pop_rejkurt(EEG, 1, settings.channelIndices, ...
                settings.localThresholdSD, settings.globalThresholdSD, ...
                double(settings.superpose), 0, ...
                settings.visualizationType);
            EEG = nf_eeg_checkset_if_available(EEG);
            methodChannelMask = nf_get_channel_rejection_matrix(EEG, ...
                'rejkurtE', settings.channelIndices);
            globalMask = nf_get_epoch_rejection_vector(EEG, ...
                'rejkurt', false);
            methodEpochMask = globalMask | any(methodChannelMask, 1);
            methodDetail = struct();
            methodDetail.localThresholdSD = ...
                settings.localThresholdSD;
            methodDetail.globalThresholdSD = ...
                settings.globalThresholdSD;
            methodDetail.channelIndices = settings.channelIndices;
            methodDetail.superpose = logical(settings.superpose);
            methodDetail.visualizationType = ...
                settings.visualizationType;
            methodDetail.globalMask = globalMask;
            methodDetail.nFlaggedEpochs = sum(methodEpochMask);
            functionNames{end + 1} = 'pop_rejkurt';
        otherwise
            error('nf_thresh:UnknownEEGLABStat', ...
                'Unknown EEGLAB epoch statistic: %s.', method);
    end
    detail.(nf_detector_field(method)) = methodDetail;
    channelMask = channelMask | methodChannelMask;
    epochMask = epochMask | methodEpochMask;
end
detail.methods = settings.methods;
detail.nFlaggedEpochs = sum(epochMask);
end

function [candidateRejected, actions] = nf_select_epoch_actions( ...
    channelMask, directEpochMask, frontalMask, excessMask, repair)
trials = size(channelMask, 2);
candidateRejected = false(1, trials);
actions = repmat({'clean'}, 1, trials);
for trialIndex = 1:trials
    applicable = {};
    if directEpochMask(trialIndex)
        applicable{end + 1} = repair.detectedEpochAction;
    end
    if excessMask(trialIndex)
        applicable{end + 1} = repair.excessChannelAction;
    end
    if frontalMask(trialIndex)
        applicable{end + 1} = repair.frontalAction;
    end
    if any(channelMask(:, trialIndex))
        applicable{end + 1} = repair.sparseChannelAction;
    end
    if isempty(applicable)
        action = 'clean';
    else
        action = nf_most_conservative_action(applicable);
    end
    actions{trialIndex} = action;
    candidateRejected(trialIndex) = strcmp(action, 'rejectepoch');
end
end

function action = nf_most_conservative_action(applicable)
if any(strcmp(applicable, 'rejectepoch'))
    action = 'rejectepoch';
elseif any(strcmp(applicable, 'interpolate'))
    action = 'interpolate';
else
    action = 'mark';
end
end

function merged = nf_merge_known_options(defaults, overrides, optionName)
if ~isstruct(overrides) || numel(overrides) ~= 1
    error('nf_thresh:InvalidOptionStruct', ...
        '%s must be one scalar structure.', optionName);
end
merged = defaults;
overrideFields = fieldnames(overrides);
knownFields = fieldnames(defaults);
for fieldIndex = 1:numel(overrideFields)
    fieldName = overrideFields{fieldIndex};
    if ~ismember(fieldName, knownFields)
        error('nf_thresh:UnknownOptionField', ...
            'Unknown field %s.%s.', optionName, fieldName);
    end
    merged.(fieldName) = overrides.(fieldName);
end
end

function merged = nf_merge_override_aliases(primary, alias, ...
    primaryName, aliasName)
if ~isstruct(primary) || numel(primary) ~= 1 || ...
        ~isstruct(alias) || numel(alias) ~= 1
    error('nf_thresh:InvalidOptionStruct', ...
        '%s and %s must be scalar structures.', primaryName, aliasName);
end
merged = primary;
aliasFields = fieldnames(alias);
for fieldIndex = 1:numel(aliasFields)
    fieldName = aliasFields{fieldIndex};
    if isfield(merged, fieldName)
        error('nf_thresh:ConflictingOptionAliases', ...
            ['%s.%s and its %s alias were both supplied. Use only one ' ...
            'spelling.'], primaryName, fieldName, aliasName);
    end
    merged.(fieldName) = alias.(fieldName);
end
end

function nf_validate_amplitude_options(EEG, settings)
if ~nf_is_positive_scalar(settings.thresholdUV)
    error('nf_thresh:InvalidAmplitudeOptions', ...
        'amplitudeOptions.thresholdUV must be positive.');
end
nf_validate_time_range(EEG, settings.timesSeconds, ...
    'amplitudeOptions.timesSeconds');
nf_validate_channel_indices(settings.channelIndices, EEG.nbchan, ...
    'amplitudeOptions.channelIndices');
end

function nf_validate_window_detector_options(EEG, settings, optionName)
if ~nf_is_positive_scalar(settings.windowSeconds)
    error('nf_thresh:InvalidWindowOptions', ...
        '%s.windowSeconds must be positive.', optionName);
end
if ~nf_is_positive_scalar(settings.stepSeconds)
    error('nf_thresh:InvalidWindowOptions', ...
        '%s.stepSeconds must be positive.', optionName);
end
nf_validate_time_range(EEG, settings.timesSeconds, ...
    [optionName '.timesSeconds']);
nf_validate_channel_indices(settings.channelIndices, EEG.nbchan, ...
    [optionName '.channelIndices']);
selectedDuration = diff(settings.timesSeconds);
if settings.windowSeconds > selectedDuration + 1 / EEG.srate
    error('nf_thresh:InvalidWindowOptions', ...
        '%s.windowSeconds exceeds the selected time range.', optionName);
end
end

function nf_validate_point_detector_options(EEG, settings, optionName)
nf_validate_time_range(EEG, settings.timesSeconds, ...
    [optionName '.timesSeconds']);
nf_validate_channel_indices(settings.channelIndices, EEG.nbchan, ...
    [optionName '.channelIndices']);
end

function nf_validate_eeglab_stat_options(EEG, settings, optionName)
nf_validate_channel_indices(settings.channelIndices, EEG.nbchan, ...
    [optionName '.channelIndices']);
if ~nf_is_positive_scalar(settings.localThresholdSD)
    error('nf_thresh:InvalidEEGLABStatOptions', ...
        '%s.localThresholdSD must be positive.', optionName);
end
if ~nf_is_positive_scalar(settings.globalThresholdSD)
    error('nf_thresh:InvalidEEGLABStatOptions', ...
        '%s.globalThresholdSD must be positive.', optionName);
end
if ~nf_is_logical_scalar(settings.superpose)
    error('nf_thresh:InvalidEEGLABStatOptions', ...
        '%s.superpose must be binary.', optionName);
end
if ~nf_is_logical_scalar(settings.visualizationType)
    error('nf_thresh:InvalidEEGLABStatOptions', ...
        '%s.visualizationType must be binary.', optionName);
end
end

function nf_validate_time_range(EEG, value, fieldName)
if ~nf_is_increasing_pair(value) || value(1) < EEG.xmin || ...
        value(2) > EEG.xmax
    error('nf_thresh:InvalidDetectorTimeRange', ...
        '%s must be an increasing pair inside the epoch.', fieldName);
end
end

function nf_validate_channel_indices(value, numberChannels, fieldName)
if ~(isnumeric(value) && isreal(value) && isvector(value) && ...
        ~isempty(value) && all(isfinite(value)) && ...
        all(value == round(value)) && all(value >= 1) && ...
        all(value <= numberChannels) && ...
        numel(unique(value)) == numel(value))
    error('nf_thresh:InvalidChannelIndices', ...
        '%s must contain unique valid channel indices.', fieldName);
end
end

function bands = nf_normalize_fft_bands(value, legacyRange, ...
    legacyPower, samplingRate)
if isempty(value)
    bands = struct();
    bands.frequencyRangeHz = reshape(legacyRange, 1, 2);
    bands.powerThresholdDb = reshape(legacyPower, 1, 2);
elseif isnumeric(value)
    if isvector(value) && numel(value) == 2
        value = reshape(value, 1, 2);
    end
    if size(value, 2) ~= 2 && size(value, 2) ~= 4
        error('nf_thresh:InvalidFFTBands', ...
            'Numeric fftBands must have two or four columns.');
    end
    bands = repmat(struct('frequencyRangeHz', [], ...
        'powerThresholdDb', []), 1, size(value, 1));
    for bandIndex = 1:size(value, 1)
        bands(bandIndex).frequencyRangeHz = value(bandIndex, 1:2);
        if size(value, 2) == 4
            bands(bandIndex).powerThresholdDb = ...
                value(bandIndex, 3:4);
        else
            bands(bandIndex).powerThresholdDb = legacyPower;
        end
    end
elseif isstruct(value)
    bands = repmat(struct('frequencyRangeHz', [], ...
        'powerThresholdDb', []), 1, numel(value));
    for bandIndex = 1:numel(value)
        current = value(bandIndex);
        if isfield(current, 'frequencyRangeHz')
            bands(bandIndex).frequencyRangeHz = ...
                current.frequencyRangeHz;
        elseif isfield(current, 'frequencyRange')
            bands(bandIndex).frequencyRangeHz = ...
                current.frequencyRange;
        else
            error('nf_thresh:InvalidFFTBands', ...
                ['Every fftBands structure needs frequencyRangeHz or ' ...
                'frequencyRange.']);
        end
        if isfield(current, 'powerThresholdDb')
            bands(bandIndex).powerThresholdDb = ...
                current.powerThresholdDb;
        elseif isfield(current, 'powerThreshold')
            bands(bandIndex).powerThresholdDb = ...
                current.powerThreshold;
        else
            bands(bandIndex).powerThresholdDb = legacyPower;
        end
        allowedFields = {'frequencyRangeHz', 'frequencyRange', ...
            'powerThresholdDb', 'powerThreshold'};
        if any(~ismember(fieldnames(current), allowedFields))
            error('nf_thresh:InvalidFFTBands', ...
                'An fftBands structure contains an unknown field.');
        end
    end
else
    error('nf_thresh:InvalidFFTBands', ...
        'fftBands must be numeric or a structure array.');
end
if isempty(bands)
    error('nf_thresh:InvalidFFTBands', ...
        'At least one FFT band is required when fft is selected.');
end
for bandIndex = 1:numel(bands)
    frequency = bands(bandIndex).frequencyRangeHz;
    power = bands(bandIndex).powerThresholdDb;
    if ~nf_is_increasing_pair(frequency) || frequency(1) < 0 || ...
            frequency(2) >= samplingRate / 2
        error('nf_thresh:InvalidFFTBands', ...
            'Every FFT frequency range must be below Nyquist.');
    end
    if ~nf_is_increasing_pair(power)
        error('nf_thresh:InvalidFFTBands', ...
            'Every FFT power threshold must be an increasing dB pair.');
    end
    bands(bandIndex).frequencyRangeHz = reshape(frequency, 1, 2);
    bands(bandIndex).powerThresholdDb = reshape(power, 1, 2);
end
end

function methods = nf_normalize_eeglab_stat_methods(value)
if ischar(value)
    value = {value};
elseif isstring(value)
    value = cellstr(value(:));
end
if ~iscell(value) || isempty(value)
    error('nf_thresh:InvalidEEGLABStatsMethods', ...
        'eeglabStatsOptions.methods must be nonempty text or a cell.');
end
methods = cell(1, numel(value));
for index = 1:numel(value)
    if ~nf_is_text(value{index})
        error('nf_thresh:InvalidEEGLABStatsMethods', ...
            'Every eeglabStatsOptions method must be scalar text.');
    end
    normalized = lower(regexprep(strtrim(char(value{index})), ...
        '[ _-]', ''));
    switch normalized
        case {'jointprobability', 'jointprob', 'probability'}
            methods{index} = 'jointprobability';
        case {'kurtosis', 'kurt'}
            methods{index} = 'kurtosis';
        otherwise
            error('nf_thresh:UnknownEEGLABStat', ...
                'Unknown EEGLAB epoch statistic: %s.', ...
                char(value{index}));
    end
end
methods = unique(methods, 'stable');
end

function functionNames = nf_eeglab_stat_function_names(methods)
functionNames = cell(1, numel(methods));
for index = 1:numel(methods)
    switch methods{index}
        case 'jointprobability'
            functionNames{index} = 'pop_jointprob';
        case 'kurtosis'
            functionNames{index} = 'pop_rejkurt';
        otherwise
            error('nf_thresh:UnknownEEGLABStat', ...
                'Unknown EEGLAB epoch statistic: %s.', methods{index});
    end
end
end

function detectors = nf_normalize_detectors(value)
if ischar(value)
    value = {value};
elseif isstring(value)
    value = cellstr(value(:));
end
if ~iscell(value) || isempty(value)
    error('nf_thresh:InvalidDetectors', ...
        'detectors must be nonempty text or a cell array of text.');
end
detectors = cell(1, numel(value));
for index = 1:numel(value)
    if ~nf_is_text(value{index})
        error('nf_thresh:InvalidDetectors', ...
            'Every detector name must be scalar text.');
    end
    normalized = lower(regexprep(strtrim(char(value{index})), ...
        '[^a-zA-Z0-9]', ''));
    switch normalized
        case {'amplitude', 'threshold', 'voltage'}
            detectors{index} = 'amplitude';
        case {'fft', 'spectral', 'muscle'}
            detectors{index} = 'fft';
        case {'peak2peak', 'peaktopeak', 'mwpp'}
            detectors{index} = 'peak2peak';
        case {'step', 'electrodepop', 'pop'}
            detectors{index} = 'step';
        case {'gradient', 'fasttransition', 'derivative'}
            detectors{index} = 'gradient';
        case {'flatline', 'flat'}
            detectors{index} = 'flatline';
        case {'clipping', 'clip', 'saturation'}
            detectors{index} = 'clipping';
        case 'faster'
            detectors{index} = 'faster';
        case {'jointprobability', 'jointprob'}
            detectors{index} = 'jointprobability';
        case {'eeglabstats', 'eeglabstatistics'}
            detectors{index} = 'eeglabstats';
        case 'none'
            detectors{index} = 'none';
        otherwise
            error('nf_thresh:UnknownDetector', ...
                'Unknown artifact detector: %s.', char(value{index}));
    end
end
detectors = unique(detectors, 'stable');
if ismember('none', detectors) && numel(detectors) > 1
    error('nf_thresh:InvalidDetectors', ...
        'Detector none cannot be combined with another detector.');
end
if isequal(detectors, {'none'})
    detectors = {};
end
end

function action = nf_normalize_channel_action(value, fieldName)
if ~nf_is_text(value)
    error('nf_thresh:InvalidRepairAction', ...
        '%s must be scalar text.', fieldName);
end
normalized = lower(regexprep(strtrim(char(value)), '[ _-]', ''));
switch normalized
    case {'interpolate', 'repair'}
        action = 'interpolate';
    case {'reject', 'rejectepoch', 'delete'}
        action = 'rejectepoch';
    case {'mark', 'retain', 'none'}
        action = 'mark';
    otherwise
        error('nf_thresh:InvalidRepairAction', ...
            ['%s must be interpolate, rejectepoch, or mark.'], ...
            fieldName);
end
end

function action = nf_normalize_detected_epoch_action(value)
if ~nf_is_text(value)
    error('nf_thresh:InvalidRepairAction', ...
        'repairOptions.detectedEpochAction must be scalar text.');
end
normalized = lower(regexprep(strtrim(char(value)), '[ _-]', ''));
switch normalized
    case {'reject', 'rejectepoch', 'delete'}
        action = 'rejectepoch';
    case {'mark', 'retain', 'none'}
        action = 'mark';
    otherwise
        error('nf_thresh:InvalidRepairAction', ...
            ['repairOptions.detectedEpochAction must be rejectepoch ' ...
            'or mark.']);
end
end

function action = nf_normalize_epoch_action(value)
if ~nf_is_text(value)
    error('nf_thresh:InvalidRepairAction', ...
        'repairOptions.epochAction must be scalar text.');
end
normalized = lower(regexprep(strtrim(char(value)), '[ _-]', ''));
switch normalized
    case {'reject', 'delete'}
        action = 'reject';
    case {'mark', 'retain', 'none'}
        action = 'mark';
    otherwise
        error('nf_thresh:InvalidRepairAction', ...
            'repairOptions.epochAction must be reject or mark.');
end
end

function [data, actualTimes] = nf_selected_data(EEG, channelIndices, ...
    times)
sampleTimes = EEG.xmin + (0:EEG.pnts - 1) / EEG.srate;
indices = find(sampleTimes >= times(1) - eps && ...
    sampleTimes <= times(2) + eps);
if numel(indices) < 2
    error('nf_thresh:EmptyDetectorRange', ...
        'The selected detector time range contains fewer than two samples.');
end
data = double(EEG.data(channelIndices, indices, :));
actualTimes = [sampleTimes(indices(1)) sampleTimes(indices(end))];
end

function samples = nf_seconds_to_samples(seconds, samplingRate, minimum)
samples = max(minimum, round(seconds * samplingRate));
end

function starts = nf_window_starts(numberSamples, windowSamples, ...
    stepSamples)
if windowSamples > numberSamples
    error('nf_thresh:DetectorWindowTooLong', ...
        'A detector window exceeds the selected sample range.');
end
lastStart = numberSamples - windowSamples + 1;
starts = 1:stepSamples:lastStart;
if starts(end) ~= lastStart
    starts(end + 1) = lastStart;
end
end

function matrix = nf_force_channel_trial_matrix(value, channels, trials)
if isequal(size(value), [channels trials])
    matrix = value;
elseif channels == 1 && numel(value) == trials
    matrix = reshape(value, 1, trials);
elseif trials == 1 && numel(value) == channels
    matrix = reshape(value, channels, 1);
else
    error('nf_thresh:InternalDetectorShape', ...
        'A native detector produced an unexpected array shape.');
end
end

function mask = nf_expand_local_mask(localMask, channelIndices, ...
    numberChannels, trials)
if ~isequal(size(localMask), [numel(channelIndices) trials])
    error('nf_thresh:InternalDetectorShape', ...
        'A detector produced an invalid local channel mask.');
end
mask = false(numberChannels, trials);
mask(channelIndices, :) = logical(localMask);
end

function found = nf_has_run(values, minimumLength)
values = reshape(logical(values), 1, []);
if minimumLength <= 1
    found = any(values);
    return
end
if numel(values) < minimumLength
    found = false;
    return
end
counts = conv(double(values), ones(1, minimumLength), 'valid');
found = any(counts >= minimumLength);
end

function vector = nf_expand_setting_vector(value, lengthNeeded, fieldName)
if isscalar(value)
    vector = repmat(value, 1, lengthNeeded);
elseif isvector(value) && numel(value) == lengthNeeded
    vector = reshape(value, 1, lengthNeeded);
else
    error('nf_thresh:InvalidFasterOptions', ...
        '%s must be scalar or have one value per FASTER measure.', ...
        fieldName);
end
end

function EEG = nf_clear_rejection_fields(EEG, fields)
if ~isfield(EEG, 'reject') || isempty(EEG.reject)
    EEG.reject = struct();
end
for index = 1:numel(fields)
    EEG.reject.(fields{index}) = [];
end
end

function EEG = nf_clear_stat_fields(EEG, fields)
if ~isfield(EEG, 'stats') || isempty(EEG.stats) || ...
        ~isstruct(EEG.stats)
    EEG.stats = struct();
end
for index = 1:numel(fields)
    EEG.stats.(fields{index}) = [];
end
end

function mask = nf_get_channel_rejection_matrix(EEG, fieldName, ...
    channelIndices)
if ~isfield(EEG, 'reject') || ~isfield(EEG.reject, fieldName)
    error('nf_thresh:MissingRejectionFlags', ...
        'EEGLAB did not create EEG.reject.%s.', fieldName);
end
value = logical(EEG.reject.(fieldName));
if isequal(size(value), [EEG.nbchan EEG.trials])
    mask = value;
    return
end
if EEG.nbchan ~= EEG.trials && ...
        isequal(size(value), [EEG.trials EEG.nbchan])
    mask = value';
    return
end
if isequal(size(value), [numel(channelIndices) EEG.trials])
    mask = false(EEG.nbchan, EEG.trials);
    mask(channelIndices, :) = value;
    return
end
if numel(channelIndices) ~= EEG.trials && ...
        isequal(size(value), [EEG.trials numel(channelIndices)])
    mask = false(EEG.nbchan, EEG.trials);
    mask(channelIndices, :) = value';
    return
end
error('nf_thresh:InvalidRejectionFlags', ...
    'EEG.reject.%s has unexpected dimensions.', fieldName);
end

function vector = nf_get_epoch_rejection_vector(EEG, fieldName, required)
if ~isfield(EEG, 'reject') || ~isfield(EEG.reject, fieldName) || ...
        isempty(EEG.reject.(fieldName))
    if required
        error('nf_thresh:MissingRejectionFlags', ...
            'EEGLAB did not create EEG.reject.%s.', fieldName);
    end
    vector = false(1, EEG.trials);
    return
end
value = logical(EEG.reject.(fieldName));
if numel(value) ~= EEG.trials
    error('nf_thresh:InvalidRejectionFlags', ...
        'EEG.reject.%s has an unexpected number of epochs.', fieldName);
end
vector = reshape(value, 1, EEG.trials);
end

function EEG = nf_eeg_checkset_if_available(EEG)
if exist('eeg_checkset', 'file') == 2
    EEG = eeg_checkset(EEG);
end
end

function nf_require_functions(functionNames, detector, provider)
missing = {};
for index = 1:numel(functionNames)
    fileType = exist(functionNames{index}, 'file');
    if ~ismember(fileType, [2 3 6])
        missing{end + 1} = functionNames{index};
    end
end
if ~isempty(missing)
    error('nf_thresh:MissingVendorDependency', ...
        ['Detector %s explicitly requires the installed %s functions: ' ...
        '%s. No fallback was used.'], detector, provider, ...
        strjoin(missing, ', '));
end
end

function contract = nf_vendor_contract(detector, provider, ...
    functionNames, level)
nf_require_functions(functionNames, detector, provider);
contract = nf_empty_provenance();
contract.detector = detector;
contract.provider = provider;
contract.contractLevel = level;
contract.strict = true;
contract.functions = functionNames;
contract.paths = cell(1, numel(functionNames));
for index = 1:numel(functionNames)
    contract.paths{index} = nf_unique_function_path( ...
        functionNames{index}, detector, provider);
end
nf_validate_provider_paths(provider, contract.paths, detector);
contract.sha256 = cell(1, numel(contract.paths));
for index = 1:numel(contract.paths)
    contract.sha256{index} = nf_file_sha256(contract.paths{index});
    if isempty(contract.sha256{index})
        error('nf_thresh:VendorHashFailed', ...
            ['Detector %s could not calculate a SHA-256 identity for ' ...
            '%s at %s.'], detector, functionNames{index}, ...
            contract.paths{index});
    end
end
contract.codeIdentityVerified = true;
contract.version = nf_provider_version(provider, contract.paths);
contract.definition = ...
    ['Installed vendor function called directly with unique path and ' ...
    'SHA-256 source-file identity; no fallback.'];
end

function path = nf_unique_function_path(functionName, detector, provider)
allPaths = which(functionName, '-all');
if ischar(allPaths) && size(allPaths, 1) > 1
    allPaths = cellstr(allPaths);
elseif ischar(allPaths)
    allPaths = {strtrim(allPaths)};
elseif isstring(allPaths)
    allPaths = cellstr(allPaths(:));
end
allPaths = cellfun(@strtrim, allPaths, 'UniformOutput', false);
allPaths = unique(allPaths, 'stable');
if isempty(allPaths)
    error('nf_thresh:MissingVendorDependency', ...
        'Detector %s cannot resolve %s function %s.', detector, ...
        provider, functionName);
end
if numel(allPaths) > 1
    error('nf_thresh:ShadowedVendorDependency', ...
        ['Detector %s found multiple %s implementations of %s. Remove ' ...
        'path shadowing before running a strict vendor contract:\n%s'], ...
        detector, provider, functionName, strjoin(allPaths, '\n'));
end
path = allPaths{1};
end

function nf_validate_provider_paths(provider, paths, detector)
lowerPaths = cellfun(@lower, paths, 'UniformOutput', false);
if strcmp(provider, 'FASTER')
    providerPresent = cellfun(@(value) contains(value, 'faster'), ...
        lowerPaths);
    providerDirectories = cellfun(@fileparts, lowerPaths, ...
        'UniformOutput', false);
    if numel(unique(providerDirectories)) ~= 1
        error('nf_thresh:MixedVendorDependencies', ...
            ['Detector %s resolved FASTER functions from more than one ' ...
            'installation directory:\n%s'], detector, ...
            strjoin(paths, '\n'));
    end
elseif strcmp(provider, 'EEGLAB')
    providerPresent = cellfun(@(value) contains(value, 'eeglab'), ...
        lowerPaths);
else
    providerPresent = true(size(lowerPaths));
end
if any(~providerPresent)
    error('nf_thresh:UnverifiedVendorDependency', ...
        ['Detector %s resolved a function outside an identifiable %s ' ...
        'installation. Resolved paths:\n%s'], detector, provider, ...
        strjoin(paths, '\n'));
end
end

function contract = nf_native_contract(detector, definition)
contract = nf_empty_provenance();
contract.detector = detector;
contract.provider = 'NeuroFreq';
contract.contractLevel = 'native';
contract.strict = true;
contract.functions = {['nf_thresh>nf_detect_' ...
    nf_detector_field(detector)]};
contract.paths = {mfilename('fullpath')};
contract.sha256 = {nf_file_sha256(mfilename('fullpath'))};
contract.codeIdentityVerified = ~isempty(contract.sha256{1});
contract.version = 'nf_thresh schema 3.0.0';
contract.definition = definition;
end

function contract = nf_empty_provenance()
contract = struct();
contract.detector = '';
contract.provider = '';
contract.contractLevel = '';
contract.strict = true;
contract.functions = {};
contract.paths = {};
contract.sha256 = {};
contract.codeIdentityVerified = false;
contract.version = '';
contract.definition = '';
end

function digest = nf_file_sha256(pathValue)
digest = '';
fileIdentifier = fopen(pathValue, 'rb');
if fileIdentifier < 0
    return
end
fileCleanup = onCleanup(@() fclose(fileIdentifier));
try
    engine = java.security.MessageDigest.getInstance('SHA-256');
    while true
        block = fread(fileIdentifier, 1048576, '*uint8');
        if isempty(block)
            break
        end
        engine.update(block);
    end
    rawDigest = typecast(engine.digest(), 'uint8');
    digest = lower(reshape(dec2hex(rawDigest, 2)', 1, []));
catch
    digest = '';
end
clear fileCleanup
end

function version = nf_provider_version(provider, paths)
version = 'not-reported';
if strcmp(provider, 'EEGLAB') && exist('eeg_getversion', 'file') == 2
    try
        version = eeg_getversion();
    catch
        version = 'installed-version-unavailable';
    end
elseif strcmp(provider, 'FASTER')
    joined = strjoin(paths, filesep);
    token = regexp(joined, ...
        '(?i)faster[^/\\]*', 'match', 'once');
    if ~isempty(token)
        version = token;
    end
end
end

function fieldName = nf_detector_field(detector)
fieldName = lower(regexprep(char(detector), '[^a-zA-Z0-9]', ''));
if isempty(fieldName)
    error('nf_thresh:InvalidDetectorName', ...
        'A detector name could not be represented as a structure field.');
end
end

function mask = nf_legacy_detector_mask(channelMasks, detector, ...
    channels, trials)
fieldName = nf_detector_field(detector);
if isfield(channelMasks, fieldName)
    mask = channelMasks.(fieldName);
else
    mask = false(channels, trials);
end
end

function masks = nf_detector_any_epoch_masks(channelMasks, epochMasks, ...
    trials)
masks = struct();
fields = fieldnames(channelMasks);
for index = 1:numel(fields)
    fieldName = fields{index};
    masks.(fieldName) = any(channelMasks.(fieldName), 1);
    if isfield(epochMasks, fieldName)
        masks.(fieldName) = masks.(fieldName) | ...
            reshape(epochMasks.(fieldName), 1, trials);
    end
end
epochFields = fieldnames(epochMasks);
for index = 1:numel(epochFields)
    fieldName = epochFields{index};
    if ~isfield(masks, fieldName)
        masks.(fieldName) = reshape(epochMasks.(fieldName), 1, trials);
    end
end
end

function EEG = nf_clear_threshold_marks(EEG)
if ~isfield(EEG, 'reject') || isempty(EEG.reject)
    EEG.reject = struct();
end
fields = {'rejthresh', 'rejthreshE', 'rejfreq', 'rejfreqE', ...
    'rejjp', 'rejjpE', 'rejkurt', 'rejkurtE', ...
    'rejglobal', 'rejglobalE'};
for index = 1:numel(fields)
    EEG.reject.(fields{index}) = [];
end
EEG.specdata = [];
EEG = nf_clear_stat_fields(EEG, {'jp', 'jpE', 'kurt', 'kurtE'});
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
    if ~iscell(requested)
        error('nf_thresh:InvalidFrontalChannels', ...
            'frontalChannels must contain scalar-text channel labels.');
    end
    for index = 1:numel(requested)
        if ~nf_is_text(requested{index})
            error('nf_thresh:InvalidFrontalChannels', ...
                ['frontalChannels must contain scalar-text channel ' ...
                'labels.']);
        end
        requested{index} = char(requested{index});
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
    frequencyRange, times, detectors)
if any(strcmp(detectors, 'amplitude')) && ...
        ~nf_is_positive_scalar(voltageThreshold)
    error('nf_thresh:InvalidVoltageThreshold', ...
        'voltageThreshold must be a positive scalar in microvolts.');
end
if any(strcmp(detectors, 'fft')) && ...
        ~nf_is_increasing_pair(powerThreshold)
    error('nf_thresh:InvalidPowerThreshold', ...
        'powerThreshold must be a finite increasing pair in dB.');
end
if any(strcmp(detectors, 'fft')) && ...
        (~nf_is_increasing_pair(frequencyRange) || ...
        frequencyRange(1) < 0 || frequencyRange(2) >= EEG.srate / 2)
    error('nf_thresh:InvalidFrequencyRange', ...
        'frequencyRange must be nonnegative and below Nyquist.');
end
timedDetectors = {'amplitude', 'peak2peak', 'step', 'gradient', ...
    'flatline', 'clipping'};
if any(ismember(detectors, timedDetectors)) && ...
        (~nf_is_increasing_pair(times) || times(1) < EEG.xmin || ...
        times(2) > EEG.xmax)
    error('nf_thresh:InvalidTimes', ...
        'times must be a finite increasing pair inside the epoch.');
end
end

function nf_validate_repair_limits(EEG, maxBadChannels, capable, ...
        requiresInterpolation)
if ~capable
    return
end
if ~nf_is_nonnegative_integer(maxBadChannels) || ...
        maxBadChannels >= EEG.nbchan
    error('nf_thresh:InvalidMaximum', ...
        'maxBadChannels must be a nonnegative integer below EEG.nbchan.');
end
if requiresInterpolation && maxBadChannels > EEG.nbchan - 3
    error('nf_thresh:InvalidMaximum', ...
        ['maxBadChannels must leave at least three donor channels for ' ...
        'the resolved interpolation actions.']);
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

function nf_validate_epoched_eeg(EEG, minimumTrials)
if nargin < 2
    minimumTrials = 2;
end
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
        ~isfinite(EEG.trials) || EEG.trials < minimumTrials || ...
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
        ['nf_thresh requires a consistent EEGLAB dataset with at least ' ...
        '%d trial(s).'], minimumTrials);
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
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value > 0;
end

function valid = nf_is_nonnegative_scalar(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value >= 0;
end

function valid = nf_is_positive_integer(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value >= 1 && value == round(value);
end

function valid = nf_is_increasing_pair(value)
valid = isnumeric(value) && isreal(value) && numel(value) == 2 && ...
    all(isfinite(value)) && value(1) < value(2);
end

function valid = nf_is_nonnegative_integer(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value >= 0 && value == round(value);
end

function valid = nf_is_logical_scalar(value)
valid = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
    isreal(value) && ismember(value, [0 1]);
end

function valid = nf_is_text(value)
valid = ischar(value) || (isstring(value) && isscalar(value));
end

function valid = nf_is_detector_list(value)
if nf_is_text(value)
    valid = true;
    return
end
if isstring(value)
    valid = ~isempty(value);
    return
end
if ~iscell(value) || isempty(value)
    valid = false;
    return
end
valid = true;
for index = 1:numel(value)
    valid = valid && nf_is_text(value{index});
end
end

function valid = nf_is_fft_band_input(value)
valid = isempty(value) || isnumeric(value) || isstruct(value);
end
