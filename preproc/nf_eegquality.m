function [quality, figureHandle] = nf_eegquality(EEG_preclean, EEG_postclean, varargin)
% NF_EEGQUALITY  Stage-aware EEG preprocessing quality control.
%
% [QUALITY, FIGUREHANDLE] = NF_EEGQUALITY(EEG_PRECLEAN, EEG_POSTCLEAN)
% compares two continuous or epoched processing stages. In the
% nf_preprocess integration these are the selected raw input before any
% filtering/artifact correction and the data after IC subtraction.
% Paired changes are calculated only when sample rate, points per epoch, and
% epoch count match. All windowed calculations respect epoch boundaries.
%
% Name/value inputs:
%   final        Final epoched EEG dataset. Default: [].
%   report       Preprocessing report. Default: recovered from FINAL or
%                EEG_POSTCLEAN when possible.
%   icaModel     Optional pre-subtraction EEG containing the fitted ICA
%                weights. Rejected-component sensor backprojections are
%                summarized as per-channel RMS microvolts.
%   plot         Create the dashboard. Default: true.
%   visible      Figure visibility, 'on' or 'off'. Default: 'on'.
%   thresholds   Structure overriding report-derived defaults stored in
%                QUALITY.thresholds. Voltage, muscle band, and line
%                frequency defaults are recovered from REPORT when present.
%   windows      Structure with artifactSeconds, psdSeconds, and
%                psdOverlap. Default: 1, 2, and 0.5.
%   maxWindows   Deterministic maximum windows per calculation. Default 500.
%   maxPcaSamples Maximum matched samples for covariance PCA. Default 100000.
%   frequencies  PSD display range in Hz. Default [1 80].
%   title        Optional dashboard title.
%   savePath     Optional exportgraphics-compatible output path.
%
% QUALITY preserves criterion-level measurements and flags. The dashboard
% does not collapse recording quality into a single score.

nf_validate_eeg(EEG_preclean, 'EEG_preclean');
nf_validate_eeg(EEG_postclean, 'EEG_postclean');

parser = inputParser;
parser.FunctionName = 'nf_eegquality';
addRequired(parser, 'EEG_preclean', @isstruct);
addRequired(parser, 'EEG_postclean', @isstruct);
addParameter(parser, 'final', [], @nf_optional_struct);
addParameter(parser, 'report', struct(), @isstruct);
addParameter(parser, 'icaModel', [], @nf_optional_struct);
addParameter(parser, 'plot', true, @nf_logical_scalar);
addParameter(parser, 'visible', 'on', @nf_visibility);
addParameter(parser, 'thresholds', struct(), @isstruct);
addParameter(parser, 'windows', struct(), @isstruct);
addParameter(parser, 'maxWindows', 500, @nf_positive_integer);
addParameter(parser, 'maxPcaSamples', 100000, @nf_positive_integer);
addParameter(parser, 'frequencies', [1 80], @nf_frequency_pair);
addParameter(parser, 'title', '', @nf_text_scalar);
addParameter(parser, 'savePath', '', @nf_text_scalar);
parse(parser, EEG_preclean, EEG_postclean, varargin{:});
options = parser.Results;

if ~isempty(options.final)
    nf_validate_eeg(options.final, 'final');
end
if ~isempty(options.icaModel)
    nf_validate_eeg(options.icaModel, 'icaModel');
end

report = nf_resolve_report(options.report, options.final, EEG_postclean);
defaultThresholds = nf_report_threshold_defaults(nf_default_thresholds(), report);
thresholds = nf_merge_options(defaultThresholds, options.thresholds, 'thresholds');
windows = nf_merge_options(nf_default_windows(), options.windows, 'windows');
nf_validate_thresholds(thresholds);
nf_validate_windows(windows);
filterActions = nf_parse_filter_actions(report, EEG_preclean, EEG_postclean);
lineEvaluationRequested = logical(thresholds.lineEvaluationEnabled);
precleanFilter = filterActions;
precleanFilter.available = false;
precleanFilter.effectiveLowpassHz = Inf;
precleanFilter.outputSampleRateHz = EEG_preclean.srate;
precleanSupport = nf_filter_support(precleanFilter, thresholds);
postcleanSupport = nf_filter_support(filterActions, thresholds);
finalFilter = filterActions;
if ~isempty(options.final)
    finalFilter.outputSampleRateHz = options.final.srate;
end
finalSupport = nf_filter_support(finalFilter, thresholds);
precleanThresholds = thresholds;
postcleanThresholds = thresholds;
finalThresholds = thresholds;
precleanThresholds.lineEvaluationEnabled = lineEvaluationRequested && ...
    precleanSupport.lineNoiseSupported;
postcleanThresholds.lineEvaluationEnabled = lineEvaluationRequested && ...
    postcleanSupport.lineNoiseSupported;
finalThresholds.lineEvaluationEnabled = lineEvaluationRequested && ...
    finalSupport.lineNoiseSupported;
filterSupport = struct();
filterSupport.requestedByThresholds = lineEvaluationRequested;
filterSupport.preclean = precleanSupport;
filterSupport.postclean = postcleanSupport;
filterSupport.final = finalSupport;
stageContext = nf_stage_context(report);

alignment = nf_align_channels(EEG_preclean, EEG_postclean);
locations = EEG_preclean.chanlocs(alignment.preIndices);
comparability = nf_comparability(EEG_preclean, EEG_postclean, alignment);

preWindows = nf_build_windows(EEG_preclean, windows.artifactSeconds, 0, options.maxWindows);
if comparability.paired
    postWindows = preWindows;
else
    postWindows = nf_build_windows(EEG_postclean, windows.artifactSeconds, 0, options.maxWindows);
end

preMetrics = nf_stage_metrics(EEG_preclean, alignment.preIndices, alignment.labels, ...
    preWindows, precleanThresholds, windows, options.maxWindows, options.frequencies, ...
    stageContext.preStage);
postMetrics = nf_stage_metrics(EEG_postclean, alignment.postIndices, alignment.labels, ...
    postWindows, postcleanThresholds, windows, options.maxWindows, options.frequencies, ...
    stageContext.postStage);

actions = nf_parse_actions(report, EEG_preclean, EEG_postclean, options.final, ...
    filterActions);
icaSpatialImpact = nf_ica_spatial_impact( ...
    options.icaModel, report, options.maxPcaSamples);
pcaMetrics = nf_fixed_pca(EEG_preclean, EEG_postclean, alignment, locations, ...
    comparability, options.maxPcaSamples);
pcaMetrics.earlierStage = stageContext.preStage;
pcaMetrics.laterStage = stageContext.postStage;
finalMetrics = nf_final_metrics(options.final, alignment.labelsLower, finalThresholds, ...
    windows, options.maxWindows, options.frequencies);

quality = struct();
quality.schemaVersion = '3.2.0';
quality.created = datestr(now, 30); %#ok<TNOW1,DATST>
quality.provenance = struct();
quality.provenance.function = 'nf_eegquality';
quality.provenance.diagnosticReference = 'Common-average reference over analyzed channels';
quality.provenance.unitAssumption = 'EEG.data amplitudes are expressed in microvolts';
quality.provenance.reportAvailable = ~isempty(fieldnames(report));
quality.provenance.reportAdjustedThresholdDefaults = defaultThresholds;
quality.provenance.explicitThresholdOverrides = options.thresholds;
quality.provenance.comparisonScope = stageContext.scope;
quality.provenance.filter = filterActions;
quality.provenance.filterSupport = filterSupport;
quality.provenance.lineNoiseEvaluation = struct();
quality.provenance.lineNoiseEvaluation.requested = lineEvaluationRequested;
quality.provenance.lineNoiseEvaluation.precleanEnabled = ...
    precleanThresholds.lineEvaluationEnabled;
quality.provenance.lineNoiseEvaluation.postcleanEnabled = ...
    postcleanThresholds.lineEvaluationEnabled;
quality.provenance.lineNoiseEvaluation.finalEnabled = ...
    finalThresholds.lineEvaluationEnabled;
quality.provenance.lineNoiseEvaluation.precleanReason = precleanSupport.reason;
quality.provenance.lineNoiseEvaluation.postcleanReason = postcleanSupport.reason;
quality.provenance.lineNoiseEvaluation.finalReason = finalSupport.reason;
quality.channelAlignment = alignment;
quality.comparability = comparability;
quality.thresholds = thresholds;
quality.stageThresholds = struct();
quality.stageThresholds.preclean = precleanThresholds;
quality.stageThresholds.postclean = postcleanThresholds;
quality.stageThresholds.final = finalThresholds;
quality.windows = windows;
quality.metrics = struct();
quality.metrics.preclean = preMetrics;
quality.metrics.postclean = postMetrics;
quality.metrics.final = finalMetrics;
quality.metrics.pca = pcaMetrics;
quality.metrics.icaSpatialImpact = icaSpatialImpact;
quality.metrics.epochDetectors = actions.epochs.detectorMetrics;
quality.actions = actions;
quality.change = nf_change_metrics(preMetrics, postMetrics, pcaMetrics, comparability);
quality.alerts = nf_quality_alerts(quality);

figureHandle = [];
if options.plot
    actionLocations = EEG_preclean.chanlocs;
    actionLabels = {EEG_preclean.chanlocs.labels};
    figureHandle = nf_quality_figure(quality, locations, actionLocations, ...
        actionLabels, char(options.visible), char(options.title));
    if ~isempty(options.savePath)
        exportgraphics(figureHandle, char(options.savePath), 'Resolution', 200);
    end
end

end

function defaults = nf_default_thresholds()
defaults = struct();
defaults.amplitudeMicrovolts = 125;
defaults.flatStandardDeviationMicrovolts = 0.5;
defaults.muscleBandHz = [20 40];
defaults.muscleReferenceBandHz = [1 15];
defaults.muscleRatioDb = -3;
defaults.lineFrequencyHz = 60;
defaults.lineHalfWidthHz = 1;
defaults.lineSideOffsetsHz = [2 5];
defaults.lineRatioDb = 4;
defaults.lineEvaluationEnabled = true;
defaults.minimumSpatialCorrelation = 0.15;
defaults.maximumNonfiniteFraction = 0;
defaults.severityCap = 4;
defaults.alertResidualWindowFraction = 0.20;
defaults.alertBadChannelFraction = 0.10;
defaults.alertRejectedComponentFraction = 0.50;
defaults.alertRejectedEpochFraction = 0.30;
end

function defaults = nf_report_threshold_defaults(defaults, report)
if ~isfield(report, 'options') || ~isstruct(report.options)
    return
end

options = report.options;
if isfield(options, 'voltageThreshold') && ...
        isnumeric(options.voltageThreshold) && ...
        isscalar(options.voltageThreshold) && ...
        isfinite(options.voltageThreshold) && ...
        options.voltageThreshold > 0
    defaults.amplitudeMicrovolts = options.voltageThreshold;
end

if isfield(options, 'muscleRange') && ...
        isnumeric(options.muscleRange) && ...
        numel(options.muscleRange) == 2 && ...
        all(isfinite(options.muscleRange)) && ...
        options.muscleRange(1) >= 0 && ...
        options.muscleRange(1) < options.muscleRange(2)
    defaults.muscleBandHz = double(options.muscleRange(:)');
end

if isfield(options, 'notch') && ...
        isnumeric(options.notch) && ...
        isscalar(options.notch) && ...
        isfinite(options.notch) && ...
        options.notch > 0
    defaults.lineFrequencyHz = options.notch;
end
end

function defaults = nf_default_windows()
defaults = struct();
defaults.artifactSeconds = 1;
defaults.psdSeconds = 2;
defaults.psdOverlap = 0.5;
end

function merged = nf_merge_options(defaults, supplied, optionName)
merged = defaults;
suppliedNames = fieldnames(supplied);
defaultNames = fieldnames(defaults);

for fieldIndex = 1:numel(suppliedNames)
    fieldName = suppliedNames{fieldIndex};
    if ~ismember(fieldName, defaultNames)
        error('nf_eegquality:UnknownOptionField', ...
            'Unknown %s field: %s.', optionName, fieldName);
    end
    merged.(fieldName) = supplied.(fieldName);
end
end

function nf_validate_thresholds(thresholds)
nf_require_positive(thresholds.amplitudeMicrovolts, 'amplitudeMicrovolts');
nf_require_positive(thresholds.flatStandardDeviationMicrovolts, ...
    'flatStandardDeviationMicrovolts');
nf_require_band(thresholds.muscleBandHz, 'muscleBandHz');
nf_require_band(thresholds.muscleReferenceBandHz, 'muscleReferenceBandHz');
nf_require_finite(thresholds.muscleRatioDb, 'muscleRatioDb');
nf_require_positive(thresholds.lineFrequencyHz, 'lineFrequencyHz');
nf_require_positive(thresholds.lineHalfWidthHz, 'lineHalfWidthHz');
nf_require_band(thresholds.lineSideOffsetsHz, 'lineSideOffsetsHz');
nf_require_finite(thresholds.lineRatioDb, 'lineRatioDb');
if ~nf_logical_scalar(thresholds.lineEvaluationEnabled)
    error('nf_eegquality:InvalidThreshold', ...
        'lineEvaluationEnabled must be a logical scalar.');
end

if thresholds.lineSideOffsetsHz(1) <= thresholds.lineHalfWidthHz
    error('nf_eegquality:InvalidThreshold', ...
        'lineSideOffsetsHz must begin above lineHalfWidthHz.');
end

if ~isnumeric(thresholds.minimumSpatialCorrelation) || ...
        ~isscalar(thresholds.minimumSpatialCorrelation) || ...
        ~isfinite(thresholds.minimumSpatialCorrelation) || ...
        thresholds.minimumSpatialCorrelation < 0 || ...
        thresholds.minimumSpatialCorrelation >= 1
    error('nf_eegquality:InvalidThreshold', ...
        'minimumSpatialCorrelation must be in [0, 1).');
end

if ~isnumeric(thresholds.maximumNonfiniteFraction) || ...
        ~isscalar(thresholds.maximumNonfiniteFraction) || ...
        ~isfinite(thresholds.maximumNonfiniteFraction) || ...
        thresholds.maximumNonfiniteFraction < 0 || ...
        thresholds.maximumNonfiniteFraction > 1
    error('nf_eegquality:InvalidThreshold', ...
        'maximumNonfiniteFraction must be in [0, 1].');
end

nf_require_positive(thresholds.severityCap, 'severityCap');
nf_require_fraction(thresholds.alertResidualWindowFraction, 'alertResidualWindowFraction');
nf_require_fraction(thresholds.alertBadChannelFraction, 'alertBadChannelFraction');
nf_require_fraction(thresholds.alertRejectedComponentFraction, ...
    'alertRejectedComponentFraction');
nf_require_fraction(thresholds.alertRejectedEpochFraction, 'alertRejectedEpochFraction');
end

function nf_validate_windows(windows)
nf_require_positive(windows.artifactSeconds, 'windows.artifactSeconds');
nf_require_positive(windows.psdSeconds, 'windows.psdSeconds');

if ~isnumeric(windows.psdOverlap) || ...
        ~isscalar(windows.psdOverlap) || ...
        ~isfinite(windows.psdOverlap) || ...
        windows.psdOverlap < 0 || ...
        windows.psdOverlap >= 1
    error('nf_eegquality:InvalidWindows', ...
        'windows.psdOverlap must be in [0, 1).');
end
end

function alignment = nf_align_channels(preclean, postclean)
preLabels = nf_channel_labels(preclean, 'EEG_preclean');
postLabels = nf_channel_labels(postclean, 'EEG_postclean');
preLower = lower(preLabels);
postLower = lower(postLabels);

if numel(unique(preLower)) ~= numel(preLower)
    error('nf_eegquality:DuplicateChannelLabels', ...
        'EEG_preclean contains duplicate channel labels when case is ignored.');
end

if numel(unique(postLower)) ~= numel(postLower)
    error('nf_eegquality:DuplicateChannelLabels', ...
        'EEG_postclean contains duplicate channel labels when case is ignored.');
end

[commonLower, preIndices, postIndices] = intersect(preLower, postLower, 'stable');
if numel(commonLower) < 3
    error('nf_eegquality:InsufficientCommonChannels', ...
        'At least three uniquely labeled common channels are required.');
end

alignment = struct();
alignment.labels = cellstr(preLabels(preIndices));
alignment.labelsLower = cellstr(commonLower);
alignment.preIndices = preIndices(:)';
alignment.postIndices = postIndices(:)';
alignment.preOnlyLabels = cellstr(preLabels(~ismember(preLower, commonLower)));
alignment.postOnlyLabels = cellstr(postLabels(~ismember(postLower, commonLower)));
alignment.nCommon = numel(commonLower);
alignment.nPreOnly = numel(alignment.preOnlyLabels);
alignment.nPostOnly = numel(alignment.postOnlyLabels);
end

function labels = nf_channel_labels(EEG, variableName)
labels = strings(1, EEG.nbchan);
for channelIndex = 1:EEG.nbchan
    if ~isfield(EEG.chanlocs(channelIndex), 'labels')
        error('nf_eegquality:MissingChannelLabels', ...
            '%s channel %d has no label field.', variableName, channelIndex);
    end
    label = string(EEG.chanlocs(channelIndex).labels);
    label = strtrim(label);
    if strlength(label) == 0
        error('nf_eegquality:EmptyChannelLabel', ...
            '%s channel %d has an empty label.', variableName, channelIndex);
    end
    labels(channelIndex) = label;
end
end

function comparability = nf_comparability(preclean, postclean, alignment)
comparability = struct();
comparability.sameSampleRate = preclean.srate == postclean.srate;
comparability.samePointsPerEpoch = preclean.pnts == postclean.pnts;
comparability.sameEpochCount = preclean.trials == postclean.trials;
comparability.sameEpochLimits = preclean.xmin == postclean.xmin && ...
    preclean.xmax == postclean.xmax;
comparability.commonChannels = alignment.nCommon;
comparability.preOnlyChannels = alignment.nPreOnly;
comparability.postOnlyChannels = alignment.nPostOnly;
comparability.paired = comparability.sameSampleRate && ...
    comparability.samePointsPerEpoch && comparability.sameEpochCount && ...
    comparability.sameEpochLimits;
comparability.assumption = ['Matching dimensions are treated as temporal correspondence. ' ...
    'The function cannot independently verify trial lineage.'];
reasons = {};

if ~comparability.sameSampleRate
    reasons{end + 1} = 'sample rates differ';
end
if ~comparability.samePointsPerEpoch
    reasons{end + 1} = 'points per epoch differ';
end
if ~comparability.sameEpochCount
    reasons{end + 1} = 'epoch counts differ';
end
if ~comparability.sameEpochLimits
    reasons{end + 1} = 'epoch time limits differ';
end

if comparability.paired
    comparability.reason = 'Paired changes are available.';
else
    comparability.reason = ['Paired changes withheld because ' strjoin(reasons, ', ') '.'];
end
end

function windows = nf_build_windows(EEG, windowSeconds, overlapFraction, maxWindows)
targetSamples = max(8, round(windowSeconds * EEG.srate));
targetSamples = min(targetSamples, EEG.pnts);
stepSamples = max(1, round(targetSamples * (1 - overlapFraction)));
epochValues = [];
startValues = [];
stopValues = [];
globalValues = [];
globalIndex = 0;

for epochIndex = 1:EEG.trials
    if EEG.pnts <= targetSamples
        starts = 1;
    else
        starts = 1:stepSamples:(EEG.pnts - targetSamples + 1);
    end
    for startIndex = 1:numel(starts)
        globalIndex = globalIndex + 1;
        epochValues(globalIndex) = epochIndex;
        startValues(globalIndex) = starts(startIndex);
        stopValues(globalIndex) = starts(startIndex) + targetSamples - 1;
        globalValues(globalIndex) = globalIndex;
    end
end

totalWindows = numel(globalValues);
if totalWindows > maxWindows
    selected = unique(round(linspace(1, totalWindows, maxWindows)), 'stable');
else
    selected = 1:totalWindows;
end

windows = struct();
windows.epoch = epochValues(selected);
windows.startSample = startValues(selected);
windows.stopSample = stopValues(selected);
windows.originalOrdinal = globalValues(selected);
windows.nAvailable = totalWindows;
windows.nAnalyzed = numel(selected);
windows.wasSubsampled = totalWindows > numel(selected);
windows.nominalSeconds = windowSeconds;
windows.samples = targetSamples;

if isfield(EEG, 'xmin') && isnumeric(EEG.xmin) && isscalar(EEG.xmin)
    epochStart = EEG.xmin;
else
    epochStart = 0;
end

windows.centerTimeWithinEpoch = epochStart + ...
    (((windows.startSample + windows.stopSample) ./ 2) - 1) ./ EEG.srate;
end

function metrics = nf_stage_metrics(EEG, channelIndices, channelLabels, windows, ...
        thresholds, windowOptions, maxWindows, frequencyRange, stageName)
metrics = struct();
metrics.available = true;
metrics.stage = stageName;
metrics.dataset = nf_dataset_description(EEG);
metrics.channelLabels = channelLabels;
metrics.channelIndices = channelIndices;
metrics.diagnosticReference = 'Common-average reference over analyzed channels';
metrics.artifact = nf_artifact_metrics(EEG, channelIndices, windows, thresholds);
metrics.psd = nf_psd_metrics(EEG, channelIndices, windowOptions.psdSeconds, ...
    windowOptions.psdOverlap, maxWindows, frequencyRange);
end

function description = nf_dataset_description(EEG)
description = struct();
description.setname = '';
if isfield(EEG, 'setname')
    description.setname = EEG.setname;
end
description.nbchan = EEG.nbchan;
description.pnts = EEG.pnts;
description.trials = EEG.trials;
description.srate = EEG.srate;
description.durationSeconds = (EEG.pnts * EEG.trials) / EEG.srate;
description.isContinuous = EEG.trials == 1;
end

function artifact = nf_artifact_metrics(EEG, channelIndices, windows, thresholds)
channelCount = numel(channelIndices);
windowCount = windows.nAnalyzed;
amplitude = nan(channelCount, windowCount);
standardDeviation = nan(channelCount, windowCount);
rmsValue = nan(channelCount, windowCount);
muscleRatioDb = nan(channelCount, windowCount);
lineRatioDb = nan(channelCount, windowCount);
spatialCorrelation = nan(channelCount, windowCount);
nonfiniteFraction = nan(channelCount, windowCount);

for windowIndex = 1:windowCount
    epochIndex = windows.epoch(windowIndex);
    startSample = windows.startSample(windowIndex);
    stopSample = windows.stopSample(windowIndex);
    data = double(reshape(EEG.data(channelIndices, startSample:stopSample, epochIndex), ...
        channelCount, []));
    [features, spectralAvailability] = nf_window_features(data, EEG.srate, thresholds);
    amplitude(:, windowIndex) = features.amplitude;
    standardDeviation(:, windowIndex) = features.standardDeviation;
    rmsValue(:, windowIndex) = features.rms;
    muscleRatioDb(:, windowIndex) = features.muscleRatioDb;
    lineRatioDb(:, windowIndex) = features.lineRatioDb;
    spatialCorrelation(:, windowIndex) = features.spatialCorrelation;
    nonfiniteFraction(:, windowIndex) = features.nonfiniteFraction;
end

amplitudeFlag = amplitude >= thresholds.amplitudeMicrovolts;
flatFlag = standardDeviation <= thresholds.flatStandardDeviationMicrovolts;
muscleFlag = muscleRatioDb >= thresholds.muscleRatioDb;
lineFlag = lineRatioDb >= thresholds.lineRatioDb;
correlationFlag = spatialCorrelation < thresholds.minimumSpatialCorrelation;
nonfiniteFlag = nonfiniteFraction > thresholds.maximumNonfiniteFraction;

muscleFlag(~isfinite(muscleRatioDb)) = false;
lineFlag(~isfinite(lineRatioDb)) = false;
correlationFlag(~isfinite(spatialCorrelation)) = true;

flagBits = zeros(channelCount, windowCount, 'uint8');
flagBits(amplitudeFlag) = bitor(flagBits(amplitudeFlag), uint8(1));
flagBits(flatFlag) = bitor(flagBits(flatFlag), uint8(2));
flagBits(muscleFlag) = bitor(flagBits(muscleFlag), uint8(4));
flagBits(lineFlag) = bitor(flagBits(lineFlag), uint8(8));
flagBits(correlationFlag) = bitor(flagBits(correlationFlag), uint8(16));
flagBits(nonfiniteFlag) = bitor(flagBits(nonfiniteFlag), uint8(32));
anyFlag = flagBits > 0;

severity = zeros(channelCount, windowCount);
severity = nf_add_severity(severity, amplitudeFlag, ...
    amplitude ./ thresholds.amplitudeMicrovolts);
severity = nf_add_severity(severity, flatFlag, ...
    thresholds.flatStandardDeviationMicrovolts ./ max(standardDeviation, eps));
severity = nf_add_severity(severity, muscleFlag, ...
    10 .^ ((muscleRatioDb - thresholds.muscleRatioDb) ./ 10));
severity = nf_add_severity(severity, lineFlag, ...
    10 .^ ((lineRatioDb - thresholds.lineRatioDb) ./ 10));

correlationSeverity = 1 + ...
    ((thresholds.minimumSpatialCorrelation - spatialCorrelation) ./ ...
    max(thresholds.minimumSpatialCorrelation, eps));
severity = nf_add_severity(severity, correlationFlag, correlationSeverity);
nonfiniteSeverity = 1 + nonfiniteFraction;
severity = nf_add_severity(severity, nonfiniteFlag, nonfiniteSeverity);
severity(~isfinite(severity)) = thresholds.severityCap;
severity = min(severity, thresholds.severityCap);

artifact = struct();
artifact.criteria = {'amplitude', 'flatness', 'muscle', 'lineNoise', ...
    'spatialCorrelation', 'nonfinite'};
artifact.bitValues = uint8([1 2 4 8 16 32]);
artifact.windows = windows;
artifact.measurements = struct();
artifact.measurements.amplitudeMicrovolts = amplitude;
artifact.measurements.standardDeviationMicrovolts = standardDeviation;
artifact.measurements.rmsMicrovolts = rmsValue;
artifact.measurements.muscleRatioDb = muscleRatioDb;
artifact.measurements.lineRatioDb = lineRatioDb;
artifact.measurements.spatialCorrelation = spatialCorrelation;
artifact.measurements.nonfiniteFraction = nonfiniteFraction;
artifact.flags = struct();
artifact.flags.amplitude = amplitudeFlag;
artifact.flags.flatness = flatFlag;
artifact.flags.muscle = muscleFlag;
artifact.flags.lineNoise = lineFlag;
artifact.flags.spatialCorrelation = correlationFlag;
artifact.flags.nonfinite = nonfiniteFlag;
artifact.flagBits = flagBits;
artifact.severity = severity;
artifact.severityDefinition = ['Maximum capped criterion threshold ratio; zero ' ...
    'when no criterion is violated'];
artifact.anyFlag = anyFlag;
artifact.channelFraction = mean(anyFlag, 2);
artifact.windowChannelFraction = mean(anyFlag, 1);
artifact.anyWindowFraction = mean(any(anyFlag, 1));
artifact.channelWindowFraction = mean(anyFlag(:));
artifact.medianRmsMicrovolts = median(rmsValue(:), 'omitnan');
artifact.criterionFraction = [mean(amplitudeFlag(:)), mean(flatFlag(:)), ...
    mean(muscleFlag(:)), mean(lineFlag(:)), mean(correlationFlag(:)), ...
    mean(nonfiniteFlag(:))];
artifact.spectralAvailability = spectralAvailability;
end

function severity = nf_add_severity(severity, flag, score)
score(~isfinite(score)) = 0;
candidate = zeros(size(score));
candidate(flag) = score(flag);
severity = max(severity, candidate);
end

function [features, availability] = nf_window_features(data, sampleRate, thresholds)
nonfinite = ~isfinite(data);
nonfiniteFraction = mean(nonfinite, 2);
unreferencedStandardDeviation = std(data, 0, 2, 'omitnan');
data = nf_common_average(data);
channelMedian = median(data, 2, 'omitnan');
centered = data - channelMedian;
amplitude = max(abs(centered), [], 2, 'omitnan');
rmsValue = sqrt(mean(centered .^ 2, 2, 'omitnan'));
filled = nf_fill_nonfinite(data);
spatialCorrelation = nf_spatial_correlation(filled);
[spectralPower, frequencies] = nf_unscaled_spectrum(filled, sampleRate);

muscleMask = frequencies >= thresholds.muscleBandHz(1) & ...
    frequencies <= thresholds.muscleBandHz(2);
muscleReferenceMask = frequencies >= thresholds.muscleReferenceBandHz(1) & ...
    frequencies <= thresholds.muscleReferenceBandHz(2);
muscleAvailable = any(muscleMask) && any(muscleReferenceMask);

if muscleAvailable
    musclePower = mean(spectralPower(:, muscleMask), 2);
    muscleReferencePower = mean(spectralPower(:, muscleReferenceMask), 2);
    muscleRatioDb = 10 .* log10(max(musclePower, realmin) ./ ...
        max(muscleReferencePower, realmin));
else
    muscleRatioDb = nan(size(data, 1), 1);
end

lineDistance = abs(frequencies - thresholds.lineFrequencyHz);
lineMask = lineDistance <= thresholds.lineHalfWidthHz;
sideMask = lineDistance >= thresholds.lineSideOffsetsHz(1) & ...
    lineDistance <= thresholds.lineSideOffsetsHz(2);
lineAvailable = logical(thresholds.lineEvaluationEnabled) && ...
    any(lineMask) && any(sideMask);

if lineAvailable
    linePower = mean(spectralPower(:, lineMask), 2);
    sidePower = mean(spectralPower(:, sideMask), 2);
    lineRatioDb = 10 .* log10(max(linePower, realmin) ./ max(sidePower, realmin));
else
    lineRatioDb = nan(size(data, 1), 1);
end

features = struct();
features.amplitude = amplitude;
features.standardDeviation = unreferencedStandardDeviation;
features.rms = rmsValue;
features.muscleRatioDb = muscleRatioDb;
features.lineRatioDb = lineRatioDb;
features.spatialCorrelation = spatialCorrelation;
features.nonfiniteFraction = nonfiniteFraction;
availability = struct();
availability.muscle = muscleAvailable;
availability.lineNoise = lineAvailable;
end

function data = nf_common_average(data)
reference = mean(data, 1, 'omitnan');
data = data - reference;
end

function data = nf_fill_nonfinite(data)
for channelIndex = 1:size(data, 1)
    row = data(channelIndex, :);
    valid = isfinite(row);
    if any(valid)
        replacement = mean(row(valid));
    else
        replacement = 0;
    end
    row(~valid) = replacement;
    data(channelIndex, :) = row;
end
end

function correlation = nf_spatial_correlation(data)
data = data - mean(data, 2);
energy = sqrt(sum(data .^ 2, 2));
valid = energy > eps;
normalized = zeros(size(data));
normalized(valid, :) = data(valid, :) ./ energy(valid);
correlationMatrix = normalized * normalized';
correlationMatrix(1:(size(correlationMatrix, 1) + 1):end) = NaN;
correlation = median(abs(correlationMatrix), 2, 'omitnan');
correlation(~valid) = 0;
correlation(~isfinite(correlation)) = 0;
end

function [power, frequencies] = nf_unscaled_spectrum(data, sampleRate)
sampleCount = size(data, 2);
if sampleCount < 2
    power = nan(size(data, 1), 1);
    frequencies = 0;
    return
end

data = nf_linear_detrend(data);
taper = nf_hann(sampleCount);
data = data .* taper;
nfft = max(8, 2 ^ nextpow2(sampleCount));
fourier = fft(data, nfft, 2);
oneSidedCount = floor(nfft / 2) + 1;
power = abs(fourier(:, 1:oneSidedCount)) .^ 2;
frequencies = (0:(oneSidedCount - 1)) .* (sampleRate / nfft);
end

function data = nf_linear_detrend(data)
sampleCount = size(data, 2);
data = data - mean(data, 2);
time = (0:(sampleCount - 1));
time = time - mean(time);
denominator = sum(time .^ 2);
if denominator > 0
    slopes = (data * time') ./ denominator;
    data = data - slopes * time;
end
end

function taper = nf_hann(sampleCount)
if sampleCount == 1
    taper = 1;
else
    taper = 0.5 - (0.5 .* cos((2 .* pi .* (0:(sampleCount - 1))) ./ ...
        (sampleCount - 1)));
end
end

function psd = nf_psd_metrics(EEG, channelIndices, segmentSeconds, overlapFraction, ...
        maxWindows, frequencyRange)
segments = nf_build_windows(EEG, segmentSeconds, overlapFraction, maxWindows);
segmentCount = segments.nAnalyzed;
channelCount = numel(channelIndices);
sampleCount = segments.samples;
nfft = max(8, 2 ^ nextpow2(sampleCount));
oneSidedCount = floor(nfft / 2) + 1;
frequencies = (0:(oneSidedCount - 1)) .* (EEG.srate / nfft);
frequencyMask = frequencies >= frequencyRange(1) & frequencies <= frequencyRange(2);

psd = struct();
psd.available = any(frequencyMask) && segmentCount > 0;
psd.segments = segments;
psd.units = 'dB microvolts^2/Hz';
psd.method = 'Segment-safe Hann-windowed Welch mean';
psd.intervalDefinition = ['10th-90th percentile across channel spectra after ' ...
    'averaging analyzed segments within channel'];

if ~psd.available
    psd.frequencies = [];
    psd.channelDb = [];
    psd.medianDb = [];
    psd.lowerDb = [];
    psd.upperDb = [];
    return
end

powerSum = zeros(channelCount, oneSidedCount);
validSegmentCount = zeros(channelCount, 1);

for segmentIndex = 1:segmentCount
    epochIndex = segments.epoch(segmentIndex);
    startSample = segments.startSample(segmentIndex);
    stopSample = segments.stopSample(segmentIndex);
    data = double(reshape(EEG.data(channelIndices, startSample:stopSample, epochIndex), ...
        channelCount, []));
    originalFinite = mean(isfinite(data), 2);
    data = nf_common_average(data);
    data = nf_fill_nonfinite(data);
    segmentPower = nf_periodogram(data, EEG.srate, nfft);
    validChannels = originalFinite > 0;
    powerSum(validChannels, :) = powerSum(validChannels, :) + segmentPower(validChannels, :);
    validSegmentCount(validChannels) = validSegmentCount(validChannels) + 1;
end

meanPower = nan(size(powerSum));
for channelIndex = 1:channelCount
    if validSegmentCount(channelIndex) > 0
        meanPower(channelIndex, :) = powerSum(channelIndex, :) ./ ...
            validSegmentCount(channelIndex);
    end
end

channelDb = 10 .* log10(max(meanPower(:, frequencyMask), realmin));
psd.frequencies = frequencies(frequencyMask);
psd.channelDb = channelDb;
psd.medianDb = nf_percentile(channelDb, 50);
psd.lowerDb = nf_percentile(channelDb, 10);
psd.upperDb = nf_percentile(channelDb, 90);
end

function power = nf_periodogram(data, sampleRate, nfft)
data = nf_linear_detrend(data);
taper = nf_hann(size(data, 2));
taperPower = sum(taper .^ 2);
data = data .* taper;
fourier = fft(data, nfft, 2);
oneSidedCount = floor(nfft / 2) + 1;
power = abs(fourier(:, 1:oneSidedCount)) .^ 2;
power = power ./ (sampleRate * max(taperPower, eps));

if mod(nfft, 2) == 0
    if oneSidedCount > 2
        power(:, 2:(end - 1)) = 2 .* power(:, 2:(end - 1));
    end
else
    if oneSidedCount > 1
        power(:, 2:end) = 2 .* power(:, 2:end);
    end
end
end

function values = nf_percentile(data, percentile)
values = nan(1, size(data, 2));
for columnIndex = 1:size(data, 2)
    column = data(:, columnIndex);
    column = sort(column(isfinite(column)));
    count = numel(column);
    if count == 0
        continue
    end
    if count == 1
        values(columnIndex) = column;
        continue
    end
    position = 1 + ((count - 1) .* percentile ./ 100);
    lowerIndex = floor(position);
    upperIndex = ceil(position);
    fraction = position - lowerIndex;
    values(columnIndex) = column(lowerIndex) + ...
        (fraction .* (column(upperIndex) - column(lowerIndex)));
end
end

function pca = nf_fixed_pca(preclean, postclean, alignment, locations, ...
        comparability, maxSamples)
pca = struct();
pca.available = false;
pca.ocularSelectionAvailable = false;
pca.reason = '';

if ~comparability.paired
    pca.reason = comparability.reason;
    return
end

totalSamples = preclean.pnts * preclean.trials;
sampleCount = min(totalSamples, maxSamples);
linearIndices = unique(round(linspace(1, totalSamples, sampleCount)), 'stable');
preData = nf_sample_data(preclean, alignment.preIndices, linearIndices);
postData = nf_sample_data(postclean, alignment.postIndices, linearIndices);
validSamples = all(isfinite(preData), 1) & all(isfinite(postData), 1);
preData = preData(:, validSamples);
postData = postData(:, validSamples);

if size(preData, 2) < 2
    pca.reason = 'Fewer than two matched finite samples were available.';
    return
end

preData = nf_common_average(preData);
postData = nf_common_average(postData);
preData = preData - mean(preData, 2);
postData = postData - mean(postData, 2);
preCovariance = (preData * preData') ./ (size(preData, 2) - 1);
postCovariance = (postData * postData') ./ (size(postData, 2) - 1);
preCovariance = (preCovariance + preCovariance') ./ 2;
postCovariance = (postCovariance + postCovariance') ./ 2;
[vectors, valueMatrix] = eig(preCovariance);
eigenvalues = real(diag(valueMatrix));
[eigenvalues, order] = sort(eigenvalues, 'descend');
vectors = real(vectors(:, order));
eigenvalues(eigenvalues < 0) = 0;
totalVariance = sum(eigenvalues);

if totalVariance <= eps
    pca.reason = 'Earlier-stage covariance contained no estimable variance.';
    return
end

explained = 100 .* eigenvalues ./ totalVariance;
tolerance = max(size(preCovariance)) .* eps(max(eigenvalues));
rankValue = sum(eigenvalues > tolerance);
leadingCount = min([8, rankValue, size(vectors, 2)]);

if leadingCount < 1
    pca.reason = 'No nonzero principal modes were available.';
    return
end

xCoordinates = nf_x_coordinates(locations);
xRange = max(xCoordinates) - min(xCoordinates);
coordinatesValid = all(isfinite(xCoordinates)) && isfinite(xRange) && xRange > 0;

if coordinatesValid
    normalizedX = (xCoordinates - min(xCoordinates)) ./ xRange;
    frontalWeights = 0.1 + (normalizedX .^ 2);
    frontality = zeros(1, leadingCount);
    selectionScore = zeros(1, leadingCount);
    for componentIndex = 1:leadingCount
        map = vectors(:, componentIndex);
        frontality(componentIndex) = sum((map .^ 2) .* frontalWeights) ./ ...
            max(sum(map .^ 2), eps);
        selectionScore(componentIndex) = frontality(componentIndex) .* ...
            sqrt(explained(componentIndex) ./ 100);
    end
    [~, candidateIndex] = max(selectionScore);
    pca.ocularSelectionAvailable = true;
    pca.ocularSelectionReason = ['Candidate selected among leading modes by ' ...
        'anterior X-weighted loading energy and explained variance.'];
else
    frontalWeights = nan(size(xCoordinates));
    frontality = nan(1, leadingCount);
    selectionScore = nan(1, leadingCount);
    candidateIndex = 1;
    pca.ocularSelectionReason = ['Valid, varying X coordinates were unavailable. ' ...
        'The first principal mode is shown and is not labeled ocular.'];
end

candidate = vectors(:, candidateIndex);
if coordinatesValid
    signedFrontalLoading = sum(candidate .* frontalWeights);
    if signedFrontalLoading < 0
        vectors(:, candidateIndex) = -vectors(:, candidateIndex);
        candidate = -candidate;
    end
end

preVariance = max(candidate' * preCovariance * candidate, 0);
postVariance = max(candidate' * postCovariance * candidate, 0);

pca.available = true;
pca.reason = 'Fixed earlier-stage PCA basis projected onto both matched stages.';
pca.nSamples = size(preData, 2);
pca.rank = rankValue;
pca.basis = vectors(:, 1:leadingCount);
pca.eigenvalues = eigenvalues(1:leadingCount);
pca.explainedPercent = explained(1:leadingCount);
pca.candidateIndex = candidateIndex;
pca.candidateBasisMap = candidate;
pca.candidateExplainedPercent = explained(candidateIndex);
pca.candidatePreVariance = preVariance;
pca.candidatePostVariance = postVariance;
pca.candidatePreRms = sqrt(preVariance);
pca.candidatePostRms = sqrt(postVariance);
pca.candidatePreScaledMap = candidate .* sqrt(preVariance);
pca.candidatePostScaledMap = candidate .* sqrt(postVariance);
pca.frontality = frontality;
pca.selectionScore = selectionScore;
pca.frontalWeights = frontalWeights;
end

function impact = nf_ica_spatial_impact(icaModel, report, maxSamples)
impact = struct();
impact.available = false;
impact.reason = '';
impact.ocularAvailable = false;
impact.ocularReason = '';
impact.channelLabels = {};
impact.locations = struct([]);
impact.rejectedComponents = [];
impact.ocularRejectedComponents = [];
impact.totalRejectedRmsMicrovolts = [];
impact.ocularRejectedRmsMicrovolts = [];
impact.originalRmsMicrovolts = [];
impact.totalRejectedToOriginalRmsPercent = [];
impact.totalFrontalPosteriorRmsRatio = NaN;
impact.ocularFrontalPosteriorRmsRatio = NaN;
impact.nSamples = 0;
impact.reconstructionEquation = ...
    'X_rejected = icawinv(:,R) * (icaweights * icasphere * X)(R,:)';
impact.scalingNote = ...
    ['Maps are RMS microvolts of the complete rejected-component sensor ' ...
    'backprojection. They retain component activation scale and covariance.'];

if isempty(icaModel)
    impact.reason = 'No pre-subtraction ICA model was supplied.';
    return
end
if isfield(report, 'steps') && isstruct(report.steps) && ...
        isfield(report.steps, 'ica') && isstruct(report.steps.ica) && ...
        isfield(report.steps.ica, 'applied') && ...
        nf_logical_scalar(report.steps.ica.applied) && ...
        ~logical(report.steps.ica.applied)
    impact.reason = 'ICA cleaning was explicitly skipped.';
    return
end
if isfield(report, 'ica') && isstruct(report.ica) && ...
        strcmpi(strtrim(nf_get_text_field(report.ica, 'method', '')), 'none')
    impact.reason = 'ICA cleaning method was none.';
    return
end
required = {'data', 'icaweights', 'icasphere', 'icawinv'};
for index = 1:numel(required)
    if ~isfield(icaModel, required{index}) || isempty(icaModel.(required{index}))
        impact.reason = sprintf('icaModel.%s was unavailable.', required{index});
        return
    end
end
if isfield(icaModel, 'icachansind') && ~isempty(icaModel.icachansind)
    channelIndices = double(icaModel.icachansind(:)');
else
    channelIndices = 1:size(icaModel.icaweights, 2);
end
if any(~isfinite(channelIndices)) || ...
        any(channelIndices ~= round(channelIndices)) || ...
        any(channelIndices < 1) || any(channelIndices > icaModel.nbchan) || ...
        numel(unique(channelIndices)) ~= numel(channelIndices)
    impact.reason = 'The ICA channel index vector was invalid.';
    return
end

weights = double(icaModel.icaweights);
sphere = double(icaModel.icasphere);
if ndims(weights) ~= 2 || ndims(sphere) ~= 2 || ...
        size(weights, 2) ~= size(sphere, 1)
    impact.reason = 'ICA weight and sphere dimensions are incompatible.';
    return
end
unmixing = weights * sphere;
mixing = double(icaModel.icawinv);
componentCount = size(unmixing, 1);
if any(~isfinite(unmixing(:))) || any(~isfinite(mixing(:)))
    impact.reason = 'ICA mixing or unmixing weights contain NaN or Inf values.';
    return
end
if size(unmixing, 2) ~= numel(channelIndices) || ...
        size(mixing, 1) ~= numel(channelIndices) || ...
        size(mixing, 2) ~= componentCount
    impact.reason = 'ICA mixing/unmixing dimensions did not match icachansind.';
    return
end

[rejected, rejectionAvailable, rejectionReason] = ...
    nf_ica_rejected_components(report);
if ~rejectionAvailable
    impact.reason = rejectionReason;
    return
end
if any(rejected < 1) || any(rejected > componentCount)
    impact.reason = 'The reported rejected-component indices exceed the ICA rank.';
    return
end
sampleTotal = icaModel.pnts * icaModel.trials;
sampleCount = min(sampleTotal, maxSamples);
sampleIndices = unique(round(linspace(1, sampleTotal, sampleCount)), 'stable');
data = nf_sample_data(icaModel, channelIndices, sampleIndices);
if any(~isfinite(data(:)))
    impact.reason = 'The sampled ICA-model data contain NaN or Inf values.';
    return
end
activations = unmixing * data;
if isempty(rejected)
    totalRemoved = zeros(numel(channelIndices), size(data, 2));
else
    totalRemoved = mixing(:, rejected) * activations(rejected, :);
end
[ocularRejected, ocularAvailable, ocularReason] = ...
    nf_ica_ocular_components(report, rejected, componentCount);
if ocularAvailable && ~isempty(ocularRejected)
    ocularRemoved = mixing(:, ocularRejected) * ...
        activations(ocularRejected, :);
else
    ocularRemoved = zeros(numel(channelIndices), size(data, 2));
end
data = data - mean(data, 2);
totalRemoved = totalRemoved - mean(totalRemoved, 2);
ocularRemoved = ocularRemoved - mean(ocularRemoved, 2);
originalRms = sqrt(mean(data .^ 2, 2));
totalRms = sqrt(mean(totalRemoved .^ 2, 2));
ocularRms = sqrt(mean(ocularRemoved .^ 2, 2));

impact.available = true;
impact.reason = ...
    'Rejected ICA activations were backprojected into sensor space.';
impact.ocularAvailable = ocularAvailable;
impact.ocularReason = ocularReason;
impact.rejectedComponents = rejected;
impact.ocularRejectedComponents = ocularRejected;
impact.nComponents = componentCount;
impact.nRejected = numel(rejected);
impact.nOcularRejected = numel(ocularRejected);
impact.nSamples = size(data, 2);
impact.originalRmsMicrovolts = originalRms;
impact.totalRejectedRmsMicrovolts = totalRms;
impact.ocularRejectedRmsMicrovolts = ocularRms;
impact.totalRejectedToOriginalRmsPercent = ...
    100 .* totalRms ./ max(originalRms, eps);
impact.channelIndices = channelIndices;
impact.channelLabels = {icaModel.chanlocs(channelIndices).labels};
impact.locations = icaModel.chanlocs(channelIndices);
impact.totalFrontalPosteriorRmsRatio = ...
    nf_frontal_posterior_rms_ratio(totalRms, impact.locations);
impact.ocularFrontalPosteriorRmsRatio = ...
    nf_frontal_posterior_rms_ratio(ocularRms, impact.locations);
end

function [rejected, available, reason] = nf_ica_rejected_components(report)
rejected = [];
available = false;
reason = 'Rejected-component provenance was unavailable.';
if ~isfield(report, 'ica') || ~isstruct(report.ica)
    return
end
if isfield(report.ica, 'components') && ...
        isstruct(report.ica.components) && ...
        isfield(report.ica.components, 'rejected')
    rejected = report.ica.components.rejected;
    available = true;
elseif isfield(report.ica, 'classification') && ...
        isstruct(report.ica.classification) && ...
        isfield(report.ica.classification, 'rejected')
    rejected = report.ica.classification.rejected;
    available = true;
end
if ~available
    return
end
if ~isnumeric(rejected) || ~isreal(rejected) || ...
        any(~isfinite(rejected(:))) || ...
        any(rejected(:) ~= round(rejected(:))) || ...
        any(rejected(:) < 1)
    rejected = [];
    available = false;
    reason = 'Rejected-component provenance contained invalid indices.';
    return
end
rejected = unique(reshape(rejected, 1, []), 'stable');
reason = 'Rejected-component provenance was recorded.';
end

function [ocular, available, reason] = ...
    nf_ica_ocular_components(report, rejected, componentCount)
ocular = [];
available = false;
reason = 'The selected classifier does not identify ocular subclasses.';
if ~isfield(report, 'ica') || ~isstruct(report.ica) || ...
        ~isfield(report.ica, 'classification') || ...
        ~isstruct(report.ica.classification)
    return
end
classification = report.ica.classification;
if isfield(classification, 'classNames')
    classNames = nf_to_cellstr(classification.classNames);
    eyeClass = find(strcmpi(strtrim(string(classNames)), 'eye'), 1);
else
    eyeClass = [];
end
explicitThresholds = isfield(classification, 'thresholds') && ...
    isnumeric(classification.thresholds) && ...
    ~isempty(classification.thresholds);
if ~isempty(eyeClass) && explicitThresholds && ...
        size(classification.thresholds, 1) >= eyeClass && ...
        size(classification.thresholds, 2) == 2
    eyeRange = double(classification.thresholds(eyeClass, :));
    if all(isfinite(eyeRange)) && eyeRange(1) < eyeRange(2) && ...
            isfield(classification, 'probabilities') && ...
            isnumeric(classification.probabilities) && ...
            size(classification.probabilities, 1) == componentCount && ...
            size(classification.probabilities, 2) >= eyeClass
        eyeProbability = classification.probabilities(:, eyeClass);
        ocular = find(eyeProbability > eyeRange(1) & ...
            eyeProbability < eyeRange(2));
        ocular = intersect(ocular, rejected, 'stable');
        available = true;
        reason = ...
            'Rejected components meeting the configured ICLabel Eye range.';
        return
    end
    reason = 'The explicit ICLabel Eye rejection range was disabled.';
    return
end
if ~isempty(eyeClass) && ...
        ~explicitThresholds && ...
        isfield(classification, 'winningClass') && ...
        isnumeric(classification.winningClass) && ...
        numel(classification.winningClass) == componentCount
        winningClass = reshape(classification.winningClass, 1, []);
        ocular = find(winningClass == eyeClass);
        ocular = intersect(ocular, rejected, 'stable');
        available = true;
        reason = 'Rejected components whose dominant ICLabel class was Eye.';
        return
end
fields = {'horizontalEyeMovement', 'verticalEyeMovement', 'blink'};
fieldAvailable = false;
candidate = [];
for index = 1:numel(fields)
    fieldName = fields{index};
    if isfield(classification, fieldName) && ...
            isnumeric(classification.(fieldName))
        fieldAvailable = true;
        candidate = [candidate reshape(classification.(fieldName), 1, [])];
    end
end
if fieldAvailable
    candidate = unique(round(candidate(isfinite(candidate))), 'stable');
    candidate = candidate(candidate >= 1 & candidate <= componentCount);
    ocular = intersect(candidate, rejected, 'stable');
    available = true;
    reason = ...
        'Rejected ADJUST/adjusted_ADJUST horizontal-eye, vertical-eye, or blink components.';
end
end

function sampled = nf_sample_data(EEG, channelIndices, linearIndices)
sampled = zeros(numel(channelIndices), numel(linearIndices));
for epochIndex = 1:EEG.trials
    firstLinear = ((epochIndex - 1) * EEG.pnts) + 1;
    lastLinear = epochIndex * EEG.pnts;
    selectedMask = linearIndices >= firstLinear & linearIndices <= lastLinear;
    if ~any(selectedMask)
        continue
    end
    localIndices = linearIndices(selectedMask) - firstLinear + 1;
    sampled(:, selectedMask) = double(reshape(EEG.data(channelIndices, ...
        localIndices, epochIndex), numel(channelIndices), []));
end
end

function ratio = nf_frontal_posterior_rms_ratio(values, locations)
ratio = NaN;
if isempty(values) || numel(values) ~= numel(locations)
    return
end
theta = nan(1, numel(locations));
radius = nan(1, numel(locations));
for index = 1:numel(locations)
    if isfield(locations(index), 'theta') && ...
            isnumeric(locations(index).theta) && ...
            isscalar(locations(index).theta) && ...
            isfinite(locations(index).theta)
        theta(index) = double(locations(index).theta);
    end
    if isfield(locations(index), 'radius') && ...
            isnumeric(locations(index).radius) && ...
            isscalar(locations(index).radius) && ...
            isfinite(locations(index).radius)
        radius(index) = double(locations(index).radius);
    end
end
frontal = abs(theta) <= 60 & radius >= 0.25;
posterior = abs(theta) >= 120 & radius >= 0.25;
if ~any(frontal) || ~any(posterior)
    return
end
frontalRms = sqrt(mean(double(values(frontal)) .^ 2));
posteriorRms = sqrt(mean(double(values(posterior)) .^ 2));
ratio = frontalRms ./ max(posteriorRms, eps);
end

function coordinates = nf_x_coordinates(locations)
coordinates = nan(numel(locations), 1);
for channelIndex = 1:numel(locations)
    if isfield(locations(channelIndex), 'X') && ...
            isnumeric(locations(channelIndex).X) && ...
            isscalar(locations(channelIndex).X) && ...
            isfinite(locations(channelIndex).X)
        coordinates(channelIndex) = locations(channelIndex).X;
    end
end
end

function report = nf_resolve_report(explicitReport, finalEEG, postclean)
if ~isempty(fieldnames(explicitReport))
    report = explicitReport;
    return
end

if ~isempty(finalEEG) && isfield(finalEEG, 'etc') && ...
        isfield(finalEEG.etc, 'preprocess') && isstruct(finalEEG.etc.preprocess)
    report = finalEEG.etc.preprocess;
    return
end

if isfield(postclean, 'etc') && isfield(postclean.etc, 'preprocess') && ...
        isstruct(postclean.etc.preprocess)
    report = postclean.etc.preprocess;
    return
end

report = struct();
end

function filter = nf_parse_filter_actions(report, preclean, postclean)
filter = struct();
filter.available = false;
filter.source = 'unavailable';
filter.requested = struct();
filter.requested.highpassHz = NaN;
filter.requested.lowpassHz = NaN;
filter.requested.notchHz = NaN;
filter.requested.notchHalfWidthHz = NaN;
filter.requested.targetRateHz = NaN;
filter.applied = struct();
filter.applied.highpass = false;
filter.applied.lowpass = false;
filter.applied.notch = false;
filter.applied.resample = false;
filter.skipped = {};
filter.inputSampleRateHz = preclean.srate;
filter.outputSampleRateHz = postclean.srate;
filter.effectiveHighpassHz = 0;
filter.effectiveLowpassHz = Inf;

hasRecordedFilter = isfield(report, 'steps') && isstruct(report.steps) && ...
    isfield(report.steps, 'filter') && isstruct(report.steps.filter) && ...
    ~isempty(fieldnames(report.steps.filter));
if hasRecordedFilter
    recorded = report.steps.filter;
    filter.available = true;
    filter.source = 'report.steps.filter';
    if isfield(recorded, 'requested') && isstruct(recorded.requested)
        filter.requested.highpassHz = nf_get_numeric_scalar( ...
            recorded.requested, 'highpassHz', NaN);
        filter.requested.lowpassHz = nf_get_numeric_scalar( ...
            recorded.requested, 'lowpassHz', NaN);
        filter.requested.notchHz = nf_get_numeric_scalar( ...
            recorded.requested, 'notchHz', NaN);
        filter.requested.notchHalfWidthHz = nf_get_numeric_scalar( ...
            recorded.requested, 'notchHalfWidthHz', NaN);
        filter.requested.targetRateHz = nf_get_numeric_scalar( ...
            recorded.requested, 'targetRateHz', NaN);
    end
    if isfield(recorded, 'applied') && isstruct(recorded.applied)
        filter.applied.highpass = nf_get_logical_field( ...
            recorded.applied, 'highpass', false);
        filter.applied.lowpass = nf_get_logical_field( ...
            recorded.applied, 'lowpass', false);
        filter.applied.notch = nf_get_logical_field( ...
            recorded.applied, 'notch', false);
        filter.applied.resample = nf_get_logical_field( ...
            recorded.applied, 'resample', false);
    end
    if isfield(recorded, 'skipped') && iscell(recorded.skipped)
        filter.skipped = recorded.skipped;
    end
    filter.inputSampleRateHz = nf_get_numeric_scalar( ...
        recorded, 'inputSrate', preclean.srate);
    filter.outputSampleRateHz = nf_get_numeric_scalar( ...
        recorded, 'outputSrate', postclean.srate);
elseif isfield(report, 'options') && isstruct(report.options)
    filter.available = true;
    filter.source = 'report.options (inferred)';
    filter.requested.highpassHz = nf_get_numeric_scalar( ...
        report.options, 'highpass', NaN);
    filter.requested.lowpassHz = nf_get_numeric_scalar( ...
        report.options, 'lowpass', NaN);
    filter.requested.notchHz = nf_get_numeric_scalar( ...
        report.options, 'notch', NaN);
    if nf_positive_setting(filter.requested.notchHz)
        filter.requested.notchHalfWidthHz = 2;
    end
    filter.requested.targetRateHz = nf_get_numeric_scalar( ...
        report.options, 'resample', NaN);
    filter.applied.highpass = nf_positive_setting( ...
        filter.requested.highpassHz);
    filter.applied.lowpass = nf_positive_setting( ...
        filter.requested.lowpassHz);
    filter.applied.notch = nf_positive_setting(filter.requested.notchHz);
    if filter.applied.notch && filter.applied.lowpass && ...
            filter.requested.lowpassHz <= filter.requested.notchHz - 2
        filter.applied.notch = false;
        filter.skipped = { ...
            'Notch inferred redundant because the low-pass is below its stop band.'};
    end
    filter.applied.resample = nf_positive_setting( ...
        filter.requested.targetRateHz) && ...
        filter.requested.targetRateHz ~= preclean.srate;
end

if isfield(report, 'options') && isstruct(report.options)
    if ~isfinite(filter.requested.highpassHz)
        filter.requested.highpassHz = nf_get_numeric_scalar( ...
            report.options, 'highpass', NaN);
    end
    if ~isfinite(filter.requested.lowpassHz)
        filter.requested.lowpassHz = nf_get_numeric_scalar( ...
            report.options, 'lowpass', NaN);
    end
    if ~isfinite(filter.requested.notchHz)
        filter.requested.notchHz = nf_get_numeric_scalar( ...
            report.options, 'notch', NaN);
    end
    if ~isfinite(filter.requested.targetRateHz)
        filter.requested.targetRateHz = nf_get_numeric_scalar( ...
            report.options, 'resample', NaN);
    end
end
if nf_positive_setting(filter.requested.notchHz) && ...
        ~isfinite(filter.requested.notchHalfWidthHz)
    filter.requested.notchHalfWidthHz = 2;
end

if filter.applied.highpass && ...
        nf_positive_setting(filter.requested.highpassHz)
    filter.effectiveHighpassHz = filter.requested.highpassHz;
end
if filter.applied.lowpass && ...
        nf_positive_setting(filter.requested.lowpassHz)
    filter.effectiveLowpassHz = filter.requested.lowpassHz;
end

filter.highpass = filter.effectiveHighpassHz;
filter.lowpass = filter.effectiveLowpassHz;
filter.resample = filter.outputSampleRateHz;
filter.notch = filter.requested.notchHz;
end

function support = nf_filter_support(filter, thresholds)
support = struct();
support.requiredLineUpperHz = thresholds.lineFrequencyHz + ...
    thresholds.lineSideOffsetsHz(2);
support.effectiveLowpassHz = filter.effectiveLowpassHz;
support.outputNyquistHz = filter.outputSampleRateHz ./ 2;
support.lineNoiseSupported = true;
support.reason = 'the recorded passband supports the requested line-noise bands';

if isfinite(filter.effectiveLowpassHz) && ...
        filter.effectiveLowpassHz < support.requiredLineUpperHz
    support.lineNoiseSupported = false;
    support.reason = sprintf( ...
        'the applied %.3g-Hz low-pass is below the %.3g-Hz line-analysis upper edge', ...
        filter.effectiveLowpassHz, support.requiredLineUpperHz);
elseif isfinite(support.outputNyquistHz) && ...
        support.outputNyquistHz <= support.requiredLineUpperHz
    support.lineNoiseSupported = false;
    support.reason = sprintf( ...
        'the %.3g-Hz output Nyquist limit does not exceed the %.3g-Hz line-analysis upper edge', ...
        support.outputNyquistHz, support.requiredLineUpperHz);
elseif ~filter.available
    support.reason = ['filter provenance was unavailable; line-noise support ' ...
        'is evaluated from each analysis window and Nyquist limit'];
end
end

function context = nf_stage_context(report)
context = struct();
context.preStage = 'selected raw input';
context.postStage = 'post-preprocessing';
context.scope = ['The two supplied stages are described independently. ' ...
    'Paired sample changes are shown only when their axes remain matched.'];
if isfield(report, 'quality') && isstruct(report.quality)
    context.preStage = nf_get_text_field( ...
        report.quality, 'earlierStage', context.preStage);
    context.postStage = nf_get_text_field( ...
        report.quality, 'laterStage', context.postStage);
    context.scope = nf_get_text_field( ...
        report.quality, 'comparisonScope', context.scope);
end
end

function actions = nf_parse_actions(report, preclean, postclean, finalEEG, filterActions)
actions = struct();
actions.reportAvailable = ~isempty(fieldnames(report));
actions.reportVersion = nf_get_text_field(report, 'schemaVersion', ...
    nf_get_text_field(report, 'version', 'unknown'));
actions.preset = nf_get_text_field(report, 'preset', 'unknown');
if isfield(report, 'steps') && isstruct(report.steps)
    actions.steps = report.steps;
else
    actions.steps = struct();
end
actions.filter = filterActions;
actions.channels = nf_empty_channel_actions();
actions.gedai = struct();
actions.gedai.available = false;
actions.ica = nf_empty_ica_actions();
actions.epochs = nf_empty_epoch_actions();

channelStepEnabled = true;
if isfield(report, 'steps') && isfield(report.steps, 'badChannels')
    if isstruct(report.steps.badChannels)
        channelStepEnabled = nf_get_logical_field( ...
            report.steps.badChannels, 'applied', true);
    else
        channelStepEnabled = logical(report.steps.badChannels);
    end
end
if isfield(report, 'channels') && ...
        isstruct(report.channels) && ~isempty(fieldnames(report.channels))
    channelReport = report.channels;
    actions.channels.available = true;
    actions.channels.detectionApplied = channelStepEnabled;
    if isfield(channelReport, 'detection') && isstruct(channelReport.detection)
        detection = channelReport.detection;
    else
        detection = channelReport;
    end
    if isfield(channelReport, 'finalization') && ...
            isstruct(channelReport.finalization)
        finalization = channelReport.finalization;
    else
        finalization = struct();
    end
    actions.channels.method = nf_get_text_field(detection, 'method', 'unknown');
    actions.channels.nOriginal = nf_get_numeric_scalar(detection, ...
        'nOriginal', NaN);
    if ~isfinite(actions.channels.nOriginal) && ...
            isfield(channelReport, 'selection') && ...
            isstruct(channelReport.selection)
        actions.channels.nOriginal = nf_get_numeric_scalar( ...
            channelReport.selection, 'nSelected', NaN);
    end

    if isfield(detection, 'artifact') && isstruct(detection.artifact)
        artifactLabels = nf_get_label_field(detection.artifact, 'labels');
        artifactIndices = nf_get_numeric_field( ...
            detection.artifact, 'indices', []);
        artifactCount = nf_get_numeric_scalar( ...
            detection.artifact, 'nDetected', numel(artifactLabels));
    else
        artifactLabels = nf_get_label_field(detection, 'labels');
        artifactIndices = nf_get_numeric_field(detection, 'indices', []);
        artifactCount = nf_get_numeric_scalar( ...
            detection, 'nRemoved', numel(artifactLabels));
    end
    actions.channels.categories.artifact = nf_make_channel_category( ...
        artifactLabels, artifactIndices, artifactCount);

    referenceLabels = {};
    referenceIndices = [];
    referenceCount = 0;
    if isfield(detection, 'reference') && isstruct(detection.reference)
        referenceRemoved = nf_get_logical_field( ...
            detection.reference, 'removedZeroReference', false);
        if referenceRemoved
            referenceLabel = nf_get_text_field( ...
                detection.reference, 'label', '');
            if ~isempty(referenceLabel)
                referenceLabels = {referenceLabel};
            end
            referenceIndex = nf_get_numeric_scalar( ...
                detection.reference, 'index', NaN);
            if isfinite(referenceIndex)
                referenceIndices = referenceIndex;
            end
            referenceCount = 1;
        end
    end
    actions.channels.categories.reference = nf_make_channel_category( ...
        referenceLabels, referenceIndices, referenceCount);

    if isfield(detection, 'removed') && isstruct(detection.removed)
        removedLabels = nf_get_label_field(detection.removed, 'labels');
        removedIndices = nf_get_numeric_field(detection.removed, 'indices', []);
        removedCount = nf_get_numeric_scalar( ...
            detection.removed, 'nRemoved', numel(removedLabels));
    else
        removedLabels = nf_union_channel_labels(artifactLabels, referenceLabels);
        removedIndices = unique([artifactIndices(:); referenceIndices(:)])';
        removedCount = numel(removedLabels);
    end
    actions.channels.removed = nf_make_channel_category( ...
        removedLabels, removedIndices, removedCount);

    missingLabels = nf_get_label_field( ...
        finalization, 'missingBeforeInterpolation');
    missingCount = nf_get_numeric_scalar(finalization, ...
        'nMissingBeforeInterpolation', numel(missingLabels));
    directInterpolation = false;
    if isfield(detection, 'interpolation') && ...
            isstruct(detection.interpolation)
        directInterpolation = nf_get_logical_field( ...
            detection.interpolation, 'performed', false);
    end
    actions.channels.missingAtFinalization = nf_make_channel_category( ...
        missingLabels, [], missingCount);
    actions.channels.interpolated = nf_get_logical_field(finalization, ...
        'interpolated', nf_get_logical_field( ...
        detection, 'interpolated', directInterpolation));
    actions.channels.globalInterpolationRequested = nf_get_logical_field( ...
        finalization, 'requestedGlobalInterpolation', ...
        actions.channels.interpolated);
    if actions.channels.interpolated
        if isempty(missingLabels)
            missingLabels = removedLabels;
            missingCount = removedCount;
            actions.channels.missingAtFinalization = ...
                nf_make_channel_category(missingLabels, [], missingCount);
        end
        actions.channels.categories.interpolated = ...
            nf_make_channel_category(missingLabels, [], missingCount);
    end
    actions.channels.cutoffUsed = nf_get_numeric_field(detection, ...
        'cutoffUsed', NaN);
    if isfield(detection, 'detector') && isstruct(detection.detector)
        actions.channels.metrics = detection.detector;
    elseif isfield(detection, 'metrics') && isstruct(detection.metrics)
        actions.channels.metrics = detection.metrics;
    end
end

if isfield(report, 'ica') && isstruct(report.ica)
    if isfield(report.ica, 'training') && isstruct(report.ica.training)
        preparationLabels = nf_get_label_field( ...
            report.ica.training, 'preparationBadChannelLabels');
        preparationIndices = nf_get_numeric_field( ...
            report.ica.training, 'preparationBadChannelIndices', []);
        actions.channels.categories.icaPreparation = ...
            nf_make_channel_category(preparationLabels, preparationIndices, ...
            numel(preparationLabels));
    end
    if isfield(report.ica, 'rank') && isstruct(report.ica.rank)
        dependentLabels = nf_get_label_field( ...
            report.ica.rank, 'dependentChannelLabels');
        dependentIndices = nf_get_numeric_field( ...
            report.ica.rank, 'dependentChannelIndices', []);
        actions.channels.categories.rankDependent = ...
            nf_make_channel_category(dependentLabels, dependentIndices, ...
            numel(dependentLabels));
    end
end

actions.channels = nf_add_prep_channel_actions( ...
    actions.channels, report, preclean, finalEEG);
actions.channels = nf_finalize_channel_actions(actions.channels, preclean);

precleanRecorded = isfield(report, 'preclean') && ...
    isstruct(report.preclean) && ...
    ~isempty(fieldnames(report.preclean));
if precleanRecorded
    actions.gedai = report.preclean;
    actions.gedai.available = nf_get_logical_field( ...
        report.preclean, 'applied', true);
elseif isfield(report, 'gedai') && isstruct(report.gedai) && ...
        ~isempty(fieldnames(report.gedai))
    actions.gedai = report.gedai;
    actions.gedai.available = nf_get_logical_field( ...
        report.gedai, 'applied', true);
elseif isfield(report, 'steps') && isfield(report.steps, 'preclean')
    actions.gedai.available = logical(report.steps.preclean);
end

if isfield(report, 'ica') && isstruct(report.ica) && ...
        ~isempty(fieldnames(report.ica))
    actions.ica = nf_parse_ica_actions(report.ica);
end

if isfield(report, 'epochs') && isstruct(report.epochs) && ...
        (isfield(report.epochs, 'threshold') || ...
        isfield(report.epochs, 'nBefore') || ...
        isfield(report.epochs, 'rejectionReason') || ...
        isfield(report.epochs, 'rejectedMask'))
    actions.epochs = nf_parse_epoch_actions(report.epochs);
end
actions.localInterpolation = struct();
actions.localInterpolation.nEpochs = actions.epochs.nLocallyRepaired;
actions.localInterpolation.nChannelEpochs = ...
    actions.epochs.nChannelEpochRepairs;
actions.localInterpolation.mask = actions.epochs.localInterpolationMask;

actions.input = nf_dataset_description(preclean);
actions.postclean = nf_dataset_description(postclean);
if ~isempty(finalEEG)
    actions.final = nf_dataset_description(finalEEG);
else
    actions.final = struct();
end
end

function channels = nf_empty_channel_actions()
channels = struct();
channels.available = false;
channels.detectionApplied = false;
channels.method = '';
channels.labels = {};
channels.indices = [];
channels.nDetected = 0;
channels.nOriginal = NaN;
channels.interpolated = false;
channels.globalInterpolationRequested = false;
channels.nInterpolated = 0;
channels.nRemoved = 0;
channels.cutoffUsed = NaN;
channels.metrics = struct();
channels.prep = nf_empty_prep_channel_actions();
channels.categoryCountsAreExclusive = false;
channels.categoryNote = ['Interpolated channels can also belong to artifact, ' ...
    'reference, PREP, ICA-preparation, or rank-dependent categories.'];
channels.categories = struct();
channels.categories.artifact = nf_make_channel_category({}, [], 0);
channels.categories.reference = nf_make_channel_category({}, [], 0);
channels.categories.icaPreparation = nf_make_channel_category({}, [], 0);
channels.categories.rankDependent = nf_make_channel_category({}, [], 0);
channels.categories.interpolated = nf_make_channel_category({}, [], 0);
channels.categories.prepDetected = nf_make_channel_category({}, [], 0);
channels.removed = nf_make_channel_category({}, [], 0);
channels.missingAtFinalization = nf_make_channel_category({}, [], 0);
end

function category = nf_make_channel_category(labels, indices, reportedCount)
category = struct();
category.labels = nf_unique_channel_labels(nf_to_cellstr(labels));
if isnumeric(indices)
    indices = indices(isfinite(indices));
    category.indices = unique(round(indices(:)'));
else
    category.indices = [];
end

if isnumeric(reportedCount) && isscalar(reportedCount) && ...
        isfinite(reportedCount) && reportedCount >= 0
    category.n = max([round(reportedCount), numel(category.labels), ...
        numel(category.indices)]);
else
    category.n = max(numel(category.labels), numel(category.indices));
end
end

function prep = nf_empty_prep_channel_actions()
prep = struct();
prep.available = false;
prep.source = '';
prep.mappingNote = '';
prep.detected = nf_make_channel_category({}, [], 0);
prep.interpolated = nf_make_channel_category({}, [], 0);
prep.removed = nf_make_channel_category({}, [], 0);
prep.stillNoisy = nf_make_channel_category({}, [], 0);
end

function channels = nf_add_prep_channel_actions( ...
        channels, report, preclean, finalEEG)
if ~isfield(report, 'preclean') || ~isstruct(report.preclean) || ...
        ~isfield(report.preclean, 'nativeNoiseDetectionSummary') || ...
        ~isstruct(report.preclean.nativeNoiseDetectionSummary)
    return
end

summary = report.preclean.nativeNoiseDetectionSummary;
rawLabels = {preclean.chanlocs.labels};
if isempty(finalEEG)
    finalLabels = {};
else
    finalLabels = {finalEEG.chanlocs.labels};
end
interpolatedIndices = nf_reported_channel_indices( ...
    nf_get_numeric_field(summary, 'interpolatedChannelNumbers', []));
removedIndices = nf_reported_channel_indices( ...
    nf_get_numeric_field(summary, 'removedChannelNumbers', []));
stillNoisyIndices = nf_reported_channel_indices( ...
    nf_get_numeric_field(summary, 'stillNoisyChannelNumbers', []));
detectedIndices = unique([interpolatedIndices, removedIndices, ...
    stillNoisyIndices], 'stable');

prep = nf_empty_prep_channel_actions();
prep.available = true;
prep.source = ...
    'report.preclean.nativeNoiseDetectionSummary (official PREP prepPipeline)';
prep.mappingNote = ['Reported PREP acquisition-channel numbers are mapped to ' ...
    'selected-raw labels. Final-montage labels are used only when raw labels ' ...
    'are unavailable; unresolved numbers remain explicitly recorded.'];
prep.detected = nf_prep_channel_category( ...
    detectedIndices, rawLabels, finalLabels);
prep.interpolated = nf_prep_channel_category( ...
    interpolatedIndices, rawLabels, finalLabels);
prep.removed = nf_prep_channel_category( ...
    removedIndices, rawLabels, finalLabels);
prep.stillNoisy = nf_prep_channel_category( ...
    stillNoisyIndices, rawLabels, finalLabels);
channels.prep = prep;

channels.categories.artifact = nf_merge_channel_categories( ...
    channels.categories.artifact, prep.detected);
channels.categories.interpolated = nf_merge_channel_categories( ...
    channels.categories.interpolated, prep.interpolated);
channels.categories.prepDetected = prep.detected;
channels.removed = nf_merge_channel_categories( ...
    channels.removed, prep.removed);
if ~isempty(finalLabels)
    missingLabels = prep.removed.labels(~ismember( ...
        lower(string(prep.removed.labels)), lower(string(finalLabels))));
    prepMissing = nf_make_channel_category( ...
        missingLabels, [], numel(missingLabels));
    channels.missingAtFinalization = nf_merge_channel_categories( ...
        channels.missingAtFinalization, prepMissing);
end
channels.available = true;
channels.detectionApplied = true;
channels.interpolated = channels.interpolated || ...
    prep.interpolated.n > 0;
channels.globalInterpolationRequested = ...
    channels.globalInterpolationRequested || prep.interpolated.n > 0;
if isempty(channels.method) || strcmpi(channels.method, 'none') || ...
        strcmpi(channels.method, 'unknown')
    channels.method = 'PREP prepPipeline';
elseif isempty(strfind(lower(channels.method), 'prep'))
    channels.method = [channels.method ' + PREP prepPipeline'];
end
end

function indices = nf_reported_channel_indices(value)
if ~isnumeric(value) || ~isreal(value)
    indices = [];
    return
end
indices = double(value(:)');
indices = indices(isfinite(indices) & indices == round(indices) & indices >= 1);
indices = unique(indices, 'stable');
end

function category = nf_prep_channel_category(indices, rawLabels, finalLabels)
indices = nf_reported_channel_indices(indices);
labels = {};
rawMappedIndices = [];
rawMappedCount = 0;
finalMappedCount = 0;
unresolvedIndices = [];
for index = 1:numel(indices)
    channelIndex = indices(index);
    if channelIndex <= numel(rawLabels)
        labels{end + 1} = rawLabels{channelIndex};
        rawMappedIndices(end + 1) = channelIndex;
        rawMappedCount = rawMappedCount + 1;
    elseif isempty(rawLabels) && channelIndex <= numel(finalLabels)
        labels{end + 1} = finalLabels{channelIndex};
        finalMappedCount = finalMappedCount + 1;
    else
        unresolvedIndices(end + 1) = channelIndex;
    end
end
category = nf_make_channel_category( ...
    labels, rawMappedIndices, numel(indices));
category.reportedIndices = indices;
category.unresolvedIndices = unresolvedIndices;
category.nMappedFromRaw = rawMappedCount;
category.nMappedFromFinalFallback = finalMappedCount;
end

function merged = nf_merge_channel_categories(first, second)
labels = nf_union_channel_labels(first.labels, second.labels);
indices = unique([first.indices(:); second.indices(:)])';
reportedCount = max([first.n, second.n, numel(labels), numel(indices)]);
merged = nf_make_channel_category(labels, indices, reportedCount);
end

function channels = nf_finalize_channel_actions(channels, preclean)
categories = channels.categories;
channels.nDetected = categories.artifact.n;
channels.nInterpolated = categories.interpolated.n;
channels.nRemoved = channels.removed.n;
channels.labels = nf_union_channel_labels(categories.artifact.labels, ...
    categories.reference.labels, categories.icaPreparation.labels, ...
    categories.rankDependent.labels, categories.interpolated.labels);
channels.indices = unique([categories.artifact.indices(:); ...
    categories.reference.indices(:)])';

if ~isfinite(channels.nOriginal) || channels.nOriginal < 1
    channels.nOriginal = preclean.nbchan + channels.nRemoved;
end
if ~channels.available
    categoryCounts = [categories.artifact.n, categories.reference.n, ...
        categories.icaPreparation.n, categories.rankDependent.n, ...
        categories.interpolated.n, categories.prepDetected.n];
    channels.available = any(categoryCounts > 0);
end
end

function labels = nf_union_channel_labels(varargin)
labels = {};
for inputIndex = 1:nargin
    labels = [labels nf_to_cellstr(varargin{inputIndex})];
end
labels = nf_unique_channel_labels(labels);
end

function labels = nf_unique_channel_labels(labels)
if isempty(labels)
    labels = {};
    return
end

labels = cellfun(@strtrim, labels, 'UniformOutput', false);
labels = labels(~cellfun(@isempty, labels));
if isempty(labels)
    return
end
[~, firstIndices] = unique(lower(string(labels)), 'stable');
labels = labels(sort(firstIndices));
end

function ica = nf_empty_ica_actions()
ica = struct();
ica.available = false;
ica.method = '';
ica.algorithm = '';
ica.rank = NaN;
ica.nComponents = 0;
ica.rejected = [];
ica.nRejected = 0;
ica.nRetained = 0;
ica.rejectedFraction = NaN;
ica.classNames = {'Brain', 'Muscle', 'Eye', 'Heart', 'Line noise', ...
    'Channel noise', 'Other'};
ica.retainedByClass = [];
ica.rejectedByClass = [];
ica.probabilities = [];
ica.aggressive = false;
end

function ica = nf_parse_ica_actions(icaReport)
ica = nf_empty_ica_actions();
ica.available = true;
ica.method = nf_get_text_field(icaReport, 'method', 'unknown');
ica.algorithm = nf_get_text_field(icaReport, 'algorithm', 'unknown');
if strcmpi(strtrim(ica.method), 'none')
    ica.available = false;
    ica.algorithm = 'none';
    return
end
if isfield(icaReport, 'rank') && isstruct(icaReport.rank)
    ica.rank = nf_get_numeric_field(icaReport.rank, 'afterRepair', NaN);
else
    ica.rank = nf_get_numeric_field(icaReport, 'rank', NaN);
end
if isfield(icaReport, 'components') && isstruct(icaReport.components)
    componentReport = icaReport.components;
else
    componentReport = icaReport;
end
ica.nComponents = nf_get_numeric_field(componentReport, 'nComponents', 0);
ica.rejected = nf_get_numeric_field(componentReport, 'rejected', []);
ica.nRejected = nf_get_numeric_field(componentReport, 'nRejected', ...
    numel(ica.rejected));
ica.nRetained = max(0, ica.nComponents - ica.nRejected);
ica.aggressive = nf_get_logical_field(componentReport, 'aggressive', false);
if ica.nComponents > 0
    ica.rejectedFraction = ica.nRejected ./ ica.nComponents;
end

if isfield(icaReport, 'classification') && ...
        isstruct(icaReport.classification) && ...
        isfield(icaReport.classification, 'probabilities') && ...
        isnumeric(icaReport.classification.probabilities)
    probabilities = icaReport.classification.probabilities;
    if isfield(icaReport.classification, 'classNames')
        ica.classNames = nf_to_cellstr(icaReport.classification.classNames);
    end
elseif isfield(icaReport, 'iclabelProbabilities') && ...
        isnumeric(icaReport.iclabelProbabilities)
    probabilities = icaReport.iclabelProbabilities;
else
    probabilities = [];
end
ica.probabilities = probabilities;

if size(probabilities, 1) == ica.nComponents && size(probabilities, 2) >= 7
    [~, topClass] = max(probabilities(:, 1:7), [], 2);
    rejectedMask = false(ica.nComponents, 1);
    rejected = ica.rejected;
    rejected = rejected(rejected >= 1 & rejected <= ica.nComponents);
    rejectedMask(rejected) = true;
    retainedByClass = zeros(1, 7);
    rejectedByClass = zeros(1, 7);
    for classIndex = 1:7
        retainedByClass(classIndex) = sum(topClass == classIndex & ~rejectedMask);
        rejectedByClass(classIndex) = sum(topClass == classIndex & rejectedMask);
    end
    ica.retainedByClass = retainedByClass;
    ica.rejectedByClass = rejectedByClass;
end
end

function epochs = nf_empty_epoch_actions()
epochs = struct();
epochs.available = false;
epochs.nCandidate = NaN;
epochs.nRetained = NaN;
epochs.nRejected = NaN;
epochs.rejectedFraction = NaN;
epochs.nCleanRetained = NaN;
epochs.nLocallyRepaired = 0;
epochs.nChannelEpochRepairs = 0;
epochs.nVoltageFlagged = 0;
epochs.nSpectralFlagged = 0;
epochs.reasonNames = {};
epochs.reasonCounts = [];
epochs.dispositionNames = {'Clean retained', 'Locally repaired', 'Rejected'};
epochs.dispositionCounts = [];
epochs.rejectedMask = [];
epochs.retainedIndices = [];
epochs.localBadCount = [];
epochs.rejectionReason = {};
epochs.voltageFlags = [];
epochs.spectralFlags = [];
epochs.badChannelEpoch = [];
epochs.localInterpolationMask = [];
epochs.activeDetectors = {};
epochs.byDetectorChannel = struct();
epochs.byDetectorEpoch = struct();
epochs.detectorAnyEpoch = struct();
epochs.detectors = repmat(nf_empty_detector_metric(), 1, 0);
epochs.detectorMetrics = repmat(nf_empty_detector_metric(), 1, 0);
epochs.detectorEpochCounts = [];
epochs.detectorChannelEpochCounts = [];
end

function epochs = nf_parse_epoch_actions(epochReport)
epochs = nf_empty_epoch_actions();
epochs.available = true;
wrapper = epochReport;
if isfield(epochReport, 'threshold') && isstruct(epochReport.threshold)
    epochReport = epochReport.threshold;
end
epochs.nCandidate = nf_get_numeric_field(epochReport, 'nOriginal', ...
    nf_get_numeric_field(epochReport, 'nBefore', NaN));
epochs.nRetained = nf_get_numeric_field(epochReport, 'nRetained', NaN);
epochs.nRejected = nf_get_numeric_field(epochReport, 'nRejected', NaN);

if isfield(epochReport, 'rejectionReason')
    reasons = nf_to_cellstr(epochReport.rejectionReason);
else
    reasons = {};
end
epochs.rejectionReason = reasons;
if isfield(epochReport, 'masks') && isstruct(epochReport.masks)
    masks = epochReport.masks;
else
    masks = struct();
end
epochs.rejectedMask = nf_get_array_field(masks, 'rejected', ...
    nf_get_array_field(wrapper, 'rejectedMask', ...
    nf_get_array_field(epochReport, 'rejectedMask', [])));
epochs.retainedIndices = nf_get_array_field(epochReport, 'retainedIndices', []);
epochs.localBadCount = nf_get_array_field(epochReport, 'localBadCount', []);
epochs.voltageFlags = nf_get_array_field(masks, 'voltage', ...
    nf_get_array_field(epochReport, 'voltageFlags', []));
epochs.spectralFlags = nf_get_array_field(masks, 'spectral', ...
    nf_get_array_field(epochReport, 'spectralFlags', []));
epochs.badChannelEpoch = nf_get_array_field(masks, 'anyArtifact', ...
    nf_get_array_field(epochReport, 'badChannelEpoch', []));
epochs.localInterpolationMask = nf_get_array_field(masks, ...
    'localInterpolation', []);
epochs = nf_parse_epoch_detector_attribution(epochs, epochReport, masks);

if isnan(epochs.nCandidate) && ~isempty(reasons)
    epochs.nCandidate = numel(reasons);
end
if isnan(epochs.nRejected) && isfield(epochReport, 'rejectedMask')
    epochs.nRejected = sum(logical(epochReport.rejectedMask));
end
if isnan(epochs.nRetained) && isfinite(epochs.nCandidate) && isfinite(epochs.nRejected)
    epochs.nRetained = epochs.nCandidate - epochs.nRejected;
end

repairedMask = strcmpi(reasons, 'locally interpolated') | ...
    strcmpi(reasons, 'retained-local-interpolation');
epochs.nLocallyRepaired = sum(repairedMask);
if isfield(epochReport, 'nLocallyRepairedEpochs')
    epochs.nLocallyRepaired = epochReport.nLocallyRepairedEpochs;
end
if isfield(epochReport, 'nLocallyRepairedChannelEpochs')
    epochs.nChannelEpochRepairs = epochReport.nLocallyRepairedChannelEpochs;
elseif isfield(epochReport, 'localBadCount') && isnumeric(epochReport.localBadCount) && ...
        numel(epochReport.localBadCount) == numel(repairedMask)
    epochs.nChannelEpochRepairs = sum(epochReport.localBadCount(repairedMask));
end

if ~isempty(epochs.voltageFlags)
    epochs.nVoltageFlagged = sum(any(logical(epochs.voltageFlags), 1));
end
if ~isempty(epochs.spectralFlags)
    epochs.nSpectralFlagged = sum(any(logical(epochs.spectralFlags), 1));
end
amplitudeIndex = find(strcmp(epochs.activeDetectors, 'amplitude'), 1);
if isempty(epochs.voltageFlags) && ~isempty(amplitudeIndex)
    epochs.nVoltageFlagged = ...
        epochs.detectorMetrics(amplitudeIndex).nFlaggedEpochs;
end
fftIndex = find(strcmp(epochs.activeDetectors, 'fft'), 1);
if isempty(epochs.spectralFlags) && ~isempty(fftIndex)
    epochs.nSpectralFlagged = ...
        epochs.detectorMetrics(fftIndex).nFlaggedEpochs;
end

if ~isempty(reasons)
    [reasonNames, ~, reasonIndex] = unique(reasons, 'stable');
    reasonCounts = zeros(1, numel(reasonNames));
    for reasonIndexValue = 1:numel(reasonNames)
        reasonCounts(reasonIndexValue) = sum(reasonIndex == reasonIndexValue);
    end
    epochs.reasonNames = reasonNames;
    epochs.reasonCounts = reasonCounts;
end

if isfinite(epochs.nRetained)
    epochs.nCleanRetained = max(0, epochs.nRetained - epochs.nLocallyRepaired);
end
if isfinite(epochs.nCandidate) && epochs.nCandidate > 0 && isfinite(epochs.nRejected)
    epochs.rejectedFraction = epochs.nRejected ./ epochs.nCandidate;
end
if isfinite(epochs.nCleanRetained) && isfinite(epochs.nRejected)
    epochs.dispositionCounts = [epochs.nCleanRetained, epochs.nLocallyRepaired, ...
        epochs.nRejected];
end
end

function metric = nf_empty_detector_metric()
metric = struct();
metric.name = '';
metric.displayName = '';
metric.fieldName = '';
metric.masksAvailable = false;
metric.channelMaskAvailable = false;
metric.epochMaskAvailable = false;
metric.masksConsistent = true;
metric.nFlaggedEpochs = NaN;
metric.nFlaggedChannelEpochs = NaN;
metric.nDirectEpochFlags = NaN;
metric.nChannelsFlagged = NaN;
end

function epochs = nf_parse_epoch_detector_attribution( ...
        epochs, epochReport, masks)
activeDetectors = {};
if isfield(epochReport, 'activeDetectors')
    activeDetectors = nf_to_cellstr(epochReport.activeDetectors);
elseif isfield(epochReport, 'detectors')
    activeDetectors = nf_to_cellstr(epochReport.detectors);
end
if isfield(masks, 'byDetectorChannel') && ...
        isstruct(masks.byDetectorChannel)
    epochs.byDetectorChannel = masks.byDetectorChannel;
else
    epochs.byDetectorChannel = struct();
end
if isfield(masks, 'byDetectorEpoch') && ...
        isstruct(masks.byDetectorEpoch)
    epochs.byDetectorEpoch = masks.byDetectorEpoch;
else
    epochs.byDetectorEpoch = struct();
end
if isempty(activeDetectors)
    activeDetectors = [fieldnames(epochs.byDetectorChannel); ...
        fieldnames(epochs.byDetectorEpoch)]';
end

normalizedDetectors = {};
for index = 1:numel(activeDetectors)
    detector = nf_normalize_qc_detector_name(activeDetectors{index});
    if isempty(detector) || strcmp(detector, 'none')
        continue
    end
    normalizedDetectors{end + 1} = detector;
end
epochs.activeDetectors = unique(normalizedDetectors, 'stable');
epochs.detectors = repmat(nf_empty_detector_metric(), ...
    1, numel(epochs.activeDetectors));
epochs.detectorAnyEpoch = struct();

for index = 1:numel(epochs.activeDetectors)
    detector = epochs.activeDetectors{index};
    fieldName = nf_qc_detector_field(detector);
    channelMask = nf_detector_mask(epochs.byDetectorChannel, fieldName);
    epochMask = nf_detector_mask(epochs.byDetectorEpoch, fieldName);
    metric = nf_empty_detector_metric();
    metric.name = detector;
    metric.displayName = nf_detector_display_name(detector);
    metric.fieldName = fieldName;
    metric.channelMaskAvailable = ~isempty(channelMask);
    metric.epochMaskAvailable = ~isempty(epochMask);
    metric.masksAvailable = ...
        metric.channelMaskAvailable || metric.epochMaskAvailable;
    if metric.channelMaskAvailable
        channelAnyEpoch = any(channelMask, 1);
        metric.nFlaggedChannelEpochs = sum(channelMask(:));
        metric.nChannelsFlagged = sum(any(channelMask, 2));
    else
        channelAnyEpoch = [];
    end
    if metric.epochMaskAvailable
        epochMask = reshape(epochMask, 1, []);
        metric.nDirectEpochFlags = sum(epochMask);
    end
    if metric.channelMaskAvailable && metric.epochMaskAvailable && ...
            numel(channelAnyEpoch) ~= numel(epochMask)
        metric.masksConsistent = false;
        anyEpoch = [];
    elseif metric.channelMaskAvailable && metric.epochMaskAvailable
        anyEpoch = channelAnyEpoch | epochMask;
    elseif metric.channelMaskAvailable
        anyEpoch = channelAnyEpoch;
    elseif metric.epochMaskAvailable
        anyEpoch = epochMask;
    else
        anyEpoch = [];
    end
    if ~isempty(anyEpoch)
        metric.nFlaggedEpochs = sum(anyEpoch);
    end
    epochs.detectorAnyEpoch.(fieldName) = anyEpoch;
    epochs.detectors(index) = metric;
end
epochs.detectorMetrics = epochs.detectors;
if isempty(epochs.detectors)
    epochs.detectorEpochCounts = [];
    epochs.detectorChannelEpochCounts = [];
else
    epochs.detectorEpochCounts = [epochs.detectors.nFlaggedEpochs];
    epochs.detectorChannelEpochCounts = ...
        [epochs.detectors.nFlaggedChannelEpochs];
end
end

function mask = nf_detector_mask(source, fieldName)
mask = [];
if ~isstruct(source) || ~isfield(source, fieldName)
    return
end
candidate = source.(fieldName);
if ~(isnumeric(candidate) || islogical(candidate)) || ...
        ~isreal(candidate) || ndims(candidate) > 2
    return
end
mask = logical(candidate);
end

function detector = nf_normalize_qc_detector_name(value)
detector = lower(regexprep(strtrim(char(string(value))), ...
    '[^a-zA-Z0-9]', ''));
switch detector
    case {'threshold', 'voltage'}
        detector = 'amplitude';
    case {'spectral', 'muscle'}
        detector = 'fft';
    case {'peaktopeak', 'mwpp'}
        detector = 'peak2peak';
    case {'electrodepop', 'pop'}
        detector = 'step';
    case {'fasttransition', 'derivative'}
        detector = 'gradient';
    case 'flat'
        detector = 'flatline';
    case {'clip', 'saturation'}
        detector = 'clipping';
    case 'jointprob'
        detector = 'jointprobability';
    case 'eeglabstatistics'
        detector = 'eeglabstats';
end
end

function fieldName = nf_qc_detector_field(detector)
fieldName = lower(regexprep(char(detector), '[^a-zA-Z0-9]', ''));
end

function displayName = nf_detector_display_name(detector)
switch detector
    case 'amplitude'
        displayName = 'amplitude';
    case 'fft'
        displayName = 'FFT';
    case 'peak2peak'
        displayName = 'peak-to-peak';
    case 'step'
        displayName = 'step/pop';
    case 'gradient'
        displayName = 'gradient';
    case 'flatline'
        displayName = 'flatline';
    case 'clipping'
        displayName = 'clipping';
    case 'faster'
        displayName = 'FASTER';
    case 'eeglabstats'
        displayName = 'EEGLAB stats';
    case 'jointprobability'
        displayName = 'joint probability';
    otherwise
        displayName = detector;
end
end

function final = nf_final_metrics(finalEEG, baseLabelsLower, thresholds, windows, ...
        maxWindows, frequencyRange)
final = struct();
final.available = false;
final.reason = 'No final dataset was supplied.';

if isempty(finalEEG)
    return
end

nf_validate_eeg(finalEEG, 'final');
finalLabels = nf_channel_labels(finalEEG, 'final');
finalLower = lower(finalLabels);
baseLabelsLower = string(baseLabelsLower);
if numel(unique(finalLower)) ~= numel(finalLower)
    final.reason = 'Final dataset contains duplicate channel labels.';
    return
end

[sharedLabels, baseIndices, finalIndices] = intersect(baseLabelsLower, finalLower, 'stable');
finalWindows = nf_build_windows(finalEEG, windows.artifactSeconds, 0, maxWindows);
stage = nf_stage_metrics(finalEEG, 1:finalEEG.nbchan, cellstr(finalLabels), ...
    finalWindows, thresholds, windows, maxWindows, frequencyRange, ...
    'final retained dataset');
final = stage;
final.available = true;
final.baseChannelIndices = baseIndices;
final.sharedFinalChannelIndices = finalIndices;
final.sharedChannelLabels = cellstr(sharedLabels);
final.nSharedWithComparedStages = numel(sharedLabels);
final.reason = 'Final-stage descriptive metrics are available.';
end

function change = nf_change_metrics(pre, post, pca, comparability)
change = struct();
change.available = comparability.paired;
change.reason = comparability.reason;
change.spatialDescriptiveAvailable = ...
    numel(pre.artifact.channelFraction) == ...
    numel(post.artifact.channelFraction);
if change.spatialDescriptiveAvailable
    change.channelArtifactFraction = post.artifact.channelFraction - ...
        pre.artifact.channelFraction;
else
    change.channelArtifactFraction = [];
end

if ~comparability.paired
    change.windowChannelFraction = [];
    change.anyWindowFraction = NaN;
    change.channelWindowFraction = NaN;
    change.medianRmsPercent = NaN;
    change.psdMedianDb = [];
    change.pcaCandidateVariancePercent = NaN;
    return
end

change.windowChannelFraction = post.artifact.windowChannelFraction - ...
    pre.artifact.windowChannelFraction;
change.anyWindowFraction = post.artifact.anyWindowFraction - ...
    pre.artifact.anyWindowFraction;
change.channelWindowFraction = post.artifact.channelWindowFraction - ...
    pre.artifact.channelWindowFraction;
change.medianRmsPercent = 100 .* ...
    (post.artifact.medianRmsMicrovolts - pre.artifact.medianRmsMicrovolts) ./ ...
    max(pre.artifact.medianRmsMicrovolts, eps);

if pre.psd.available && post.psd.available && ...
        isequal(pre.psd.frequencies, post.psd.frequencies)
    change.psdFrequencies = pre.psd.frequencies;
    change.psdMedianDb = post.psd.medianDb - pre.psd.medianDb;
else
    change.psdFrequencies = [];
    change.psdMedianDb = [];
end

if pca.available && pca.candidatePreVariance > 0
    change.pcaCandidateVariancePercent = 100 .* ...
        (pca.candidatePostVariance - pca.candidatePreVariance) ./ ...
        pca.candidatePreVariance;
else
    change.pcaCandidateVariancePercent = NaN;
end
end

function alerts = nf_quality_alerts(quality)
alerts = struct();
alerts.available = true;
alerts.comparativeAvailable = quality.comparability.paired;
alerts.items = struct('code', {}, 'message', {}, 'value', {}, 'threshold', {});
if alerts.comparativeAvailable
    alerts.reason = ['Absolute alerts use transparent thresholds. Matched stage ' ...
        'changes are also available and are not collapsed into a score.'];
else
    alerts.reason = ['Absolute alerts use transparent thresholds. ' ...
        quality.comparability.reason];
end

if quality.metrics.final.available
    residualMetrics = quality.metrics.final;
else
    residualMetrics = quality.metrics.postclean;
end
alerts.residualStage = residualMetrics.stage;
residualWindowFraction = residualMetrics.artifact.anyWindowFraction;
thresholdValue = quality.thresholds.alertResidualWindowFraction;
if residualWindowFraction > thresholdValue
    alerts.items(end + 1) = nf_alert_item('residual_windows', ...
        sprintf('Artifact criteria occur in %.1f%% of %s windows.', ...
        100 .* residualWindowFraction, residualMetrics.stage), ...
        residualWindowFraction, thresholdValue);
end

channelActions = quality.actions.channels;
if channelActions.available && channelActions.nDetected > 0
    if isfinite(channelActions.nOriginal) && channelActions.nOriginal > 0
        channelDenominator = channelActions.nOriginal;
    else
        channelDenominator = quality.actions.input.nbchan + ...
            channelActions.nRemoved;
    end
    badFraction = channelActions.nDetected ./ max(channelDenominator, 1);
    thresholdValue = quality.thresholds.alertBadChannelFraction;
    if badFraction > thresholdValue
        alerts.items(end + 1) = nf_alert_item('bad_channels', ...
            sprintf('%.1f%% of input channels were artifact channels.', ...
            100 .* badFraction), ...
            badFraction, thresholdValue);
    end
end

icaActions = quality.actions.ica;
if icaActions.available && isfinite(icaActions.rejectedFraction)
    thresholdValue = quality.thresholds.alertRejectedComponentFraction;
    if icaActions.rejectedFraction > thresholdValue
        alerts.items(end + 1) = nf_alert_item('rejected_components', ...
            sprintf('%.1f%% of independent components were rejected.', ...
            100 .* icaActions.rejectedFraction), icaActions.rejectedFraction, ...
            thresholdValue);
    end
end

epochActions = quality.actions.epochs;
if epochActions.available && isfinite(epochActions.rejectedFraction)
    thresholdValue = quality.thresholds.alertRejectedEpochFraction;
    if epochActions.rejectedFraction > thresholdValue
        alerts.items(end + 1) = nf_alert_item('rejected_epochs', ...
            sprintf('%.1f%% of candidate epochs were rejected.', ...
            100 .* epochActions.rejectedFraction), epochActions.rejectedFraction, ...
            thresholdValue);
    end
end
end

function item = nf_alert_item(code, message, value, threshold)
item = struct();
item.code = code;
item.message = message;
item.value = value;
item.threshold = threshold;
end

function figureHandle = nf_quality_figure(quality, locations, ...
    actionLocations, actionLabels, visibility, suppliedTitle)
figureHandle = figure('Color', 'w', 'Visible', visibility, ...
    'Position', [40 40 1700 1050]);
layout = tiledlayout(figureHandle, 4, 4, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

nf_plot_ledger(nexttile(layout, 1), quality);
nf_plot_psd(nexttile(layout, 2, [1 2]), quality.metrics.preclean.psd, ...
    quality.metrics.postclean.psd, quality.metrics.preclean.stage, ...
    quality.metrics.postclean.stage);
nf_plot_ica(nexttile(layout, 4), quality.actions.ica);
nf_plot_heatmap(nexttile(layout, 5, [1 2]), quality.metrics.preclean.artifact, ...
    quality.channelAlignment.labels, quality.thresholds.severityCap, ...
    [quality.metrics.preclean.stage ' artifact severity']);
nf_plot_heatmap(nexttile(layout, 7, [1 2]), quality.metrics.postclean.artifact, ...
    quality.channelAlignment.labels, quality.thresholds.severityCap, ...
    [quality.metrics.postclean.stage ' artifact severity']);

preBurden = quality.metrics.preclean.artifact.channelFraction;
postBurden = quality.metrics.postclean.artifact.channelFraction;
burdenLimit = max([preBurden(:); postBurden(:)]);
burdenLimit = max(burdenLimit, 0.01);
nf_plot_topography(nexttile(layout, 9), preBurden, locations, ...
    [quality.metrics.preclean.stage ' spatial burden'], [0 burdenLimit]);
nf_plot_topography(nexttile(layout, 10), postBurden, locations, ...
    [quality.metrics.postclean.stage ' spatial burden'], [0 burdenLimit]);

if quality.change.spatialDescriptiveAvailable
    changeLimit = max(abs(quality.change.channelArtifactFraction));
    changeLimit = max(changeLimit, 0.01);
    nf_plot_topography(nexttile(layout, 11), ...
        quality.change.channelArtifactFraction, locations, ...
        'Later - earlier stage burden', [-changeLimit changeLimit]);
else
    nf_plot_message(nexttile(layout, 11), ...
        'Spatial burden change unavailable', quality.change.reason);
end

nf_plot_channel_actions(nexttile(layout, 12), actionLocations, ...
    actionLabels, quality.actions.channels);

icaImpact = quality.metrics.icaSpatialImpact;
if icaImpact.available
    totalLimit = max(icaImpact.totalRejectedRmsMicrovolts);
    totalLimit = max(totalLimit, eps);
    totalTitle = sprintf('All %d rejected ICs: sensor RMS (uV)', ...
        icaImpact.nRejected);
    nf_plot_topography(nexttile(layout, 13), ...
        icaImpact.totalRejectedRmsMicrovolts, icaImpact.locations, ...
        totalTitle, [0 totalLimit]);
    if icaImpact.ocularAvailable
        ocularLimit = max(icaImpact.ocularRejectedRmsMicrovolts);
        ocularLimit = max(ocularLimit, eps);
        ocularTitle = sprintf('%d ocular rejected ICs: sensor RMS (uV)', ...
            icaImpact.nOcularRejected);
        nf_plot_topography(nexttile(layout, 14), ...
            icaImpact.ocularRejectedRmsMicrovolts, ...
            icaImpact.locations, ocularTitle, [0 ocularLimit]);
    else
        ratioLimit = max(icaImpact.totalRejectedToOriginalRmsPercent);
        ratioLimit = max(ratioLimit, eps);
        nf_plot_topography(nexttile(layout, 14), ...
            icaImpact.totalRejectedToOriginalRmsPercent, ...
            icaImpact.locations, ...
            'Rejected IC RMS / pre-subtraction RMS (%)', [0 ratioLimit]);
    end
else
    nf_plot_message(nexttile(layout, 13), ...
        'Rejected-IC scalp contribution unavailable', icaImpact.reason);
    nf_plot_message(nexttile(layout, 14), ...
        'Ocular-removal scalp contribution unavailable', icaImpact.reason);
end

nf_plot_epochs(nexttile(layout, 15), quality.actions.epochs);
nf_plot_alerts(nexttile(layout, 16), quality.alerts);

if isempty(suppliedTitle)
    reportTitle = quality.metrics.postclean.dataset.setname;
    if isempty(reportTitle)
        reportTitle = 'EEG preprocessing quality';
    end
else
    reportTitle = suppliedTitle;
end
title(layout, reportTitle, 'Interpreter', 'none', 'FontWeight', 'bold');
end

function nf_plot_ledger(axisHandle, quality)
axis(axisHandle, 'off');
pre = quality.metrics.preclean.dataset;
post = quality.metrics.postclean.dataset;
channels = quality.actions.channels;
ica = quality.actions.ica;
epochs = quality.actions.epochs;
filter = quality.actions.filter;

if filter.available
    filterText = sprintf('filter: %.3g-%s Hz; rate %.3g->%.3g Hz', ...
        filter.effectiveHighpassHz, ...
        nf_upper_frequency_text(filter.effectiveLowpassHz), ...
        filter.inputSampleRateHz, ...
        filter.outputSampleRateHz);
else
    filterText = 'filter provenance: unavailable';
end

if filter.applied.notch
    notchText = sprintf('notch: %.3g Hz applied', filter.requested.notchHz);
elseif nf_positive_setting(filter.requested.notchHz)
    notchText = sprintf('notch: %.3g Hz not applied (%d skip notes)', ...
        filter.requested.notchHz, numel(filter.skipped));
else
    notchText = 'notch: disabled or unavailable';
end

if isfield(quality.actions.gedai, 'available') && quality.actions.gedai.available
    precleanMethod = nf_get_text_field(quality.actions.gedai, ...
        'method', 'GEDAI');
    sensaiScore = nf_get_numeric_scalar(quality.actions.gedai, ...
        'sensaiScore', NaN);
    artifactFraction = NaN;
    if isfield(quality.actions.gedai, 'artifacts') && ...
            isstruct(quality.actions.gedai.artifacts)
        artifactFraction = nf_get_numeric_scalar( ...
            quality.actions.gedai.artifacts, 'nonzeroFraction', NaN);
    end
    if isfinite(sensaiScore) && isfinite(artifactFraction)
        gedaiText = sprintf( ...
            '%s: SENS-AI %.3g; modeled samples %.1f%%', ...
            precleanMethod, ...
            sensaiScore, 100 .* artifactFraction);
    elseif isfinite(sensaiScore)
        gedaiText = sprintf( ...
            '%s: SENS-AI %.3g', precleanMethod, sensaiScore);
    else
        removedRms = nf_get_numeric_scalar( ...
            quality.actions.gedai, 'removedRMSMicrovolts', NaN);
        if isfinite(removedRms)
            gedaiText = sprintf( ...
                '%s preclean: removed RMS %.3g uV', ...
                precleanMethod, removedRms);
        else
            gedaiText = sprintf('%s preclean: recorded', precleanMethod);
        end
    end
else
    gedaiText = 'GEDAI/preclean: unavailable';
end

if quality.comparability.paired
    comparisonText = 'paired changes: yes';
else
    comparisonText = 'paired changes: WITHHELD';
end
stageText = sprintf('dashboard stages: %s -> %s', ...
    quality.metrics.preclean.stage, quality.metrics.postclean.stage);
precleanLineEnabled = ...
    quality.metrics.preclean.artifact.spectralAvailability.lineNoise;
postcleanLineEnabled = ...
    quality.metrics.postclean.artifact.spectralAvailability.lineNoise;
if precleanLineEnabled && postcleanLineEnabled
    lineText = 'line criterion: enabled at both displayed stages';
elseif precleanLineEnabled
    lineText = 'line criterion: raw stage only (later passband excludes it)';
elseif postcleanLineEnabled
    lineText = 'line criterion: later stage only';
else
    lineText = 'line criterion: disabled or unsupported';
end
if quality.metrics.final.available
    finalText = sprintf('final flagged channel-windows: %.1f%%', ...
        100 .* quality.metrics.final.artifact.channelWindowFraction);
else
    finalText = 'final flagged channel-windows: unavailable';
end
icaImpact = quality.metrics.icaSpatialImpact;
if icaImpact.available
    icaImpactText = sprintf( ...
        'rejected-IC sensor RMS: median %.3g uV; max %.3g uV', ...
        median(icaImpact.totalRejectedRmsMicrovolts), ...
        max(icaImpact.totalRejectedRmsMicrovolts));
    if icaImpact.ocularAvailable && ...
            isfinite(icaImpact.ocularFrontalPosteriorRmsRatio)
        ocularSpatialText = sprintf( ...
            'ocular removal frontal/posterior RMS ratio: %.2f', ...
            icaImpact.ocularFrontalPosteriorRmsRatio);
    elseif isfinite(icaImpact.totalFrontalPosteriorRmsRatio)
        ocularSpatialText = sprintf( ...
            'all rejected ICs frontal/posterior RMS ratio: %.2f', ...
            icaImpact.totalFrontalPosteriorRmsRatio);
    else
        ocularSpatialText = ...
            'rejected-IC frontal/posterior ratio: unavailable';
    end
else
    icaImpactText = 'rejected-IC sensor RMS: unavailable';
    ocularSpatialText = ...
        'rejected-IC frontal/posterior ratio: unavailable';
end
if isempty(channels.method)
    channelMethodText = 'channel actions: unavailable';
else
    channelMethodText = sprintf('channel actions: %s', channels.method);
end
prepText = '';
if channels.prep.available
    prepText = sprintf( ...
        'PREP channels: detected %d; interpolated %d; removed %d; still noisy %d', ...
        channels.prep.detected.n, channels.prep.interpolated.n, ...
        channels.prep.removed.n, channels.prep.stillNoisy.n);
end
detectorText = nf_epoch_detector_ledger_text(epochs);

lines = {
    comparisonText
    stageText
    sprintf('common channels: %d', quality.channelAlignment.nCommon)
    sprintf('sample rate: %.3g -> %.3g Hz', pre.srate, post.srate)
    sprintf('points x epochs: %d x %d -> %d x %d', ...
        pre.pnts, pre.trials, post.pnts, post.trials)
    filterText
    notchText
    lineText
    gedaiText
    channelMethodText
    prepText
    sprintf('channels: artifact %d; reference %d', ...
        channels.categories.artifact.n, channels.categories.reference.n)
    sprintf('final montage: missing %d; interpolated %d', ...
        channels.missingAtFinalization.n, ...
        channels.categories.interpolated.n)
    sprintf('ICA channel exclusions: preparation %d; rank %d', ...
        channels.categories.icaPreparation.n, ...
        channels.categories.rankDependent.n)
    sprintf('ICs: %d retained; %d rejected', ica.nRetained, ica.nRejected)
    icaImpactText
    ocularSpatialText
    sprintf('epochs: %s candidate; %s retained; %s rejected', ...
        nf_number_text(epochs.nCandidate), nf_number_text(epochs.nRetained), ...
        nf_number_text(epochs.nRejected))
    detectorText
    sprintf('local repairs: %d epochs / %d channel-epochs', ...
        epochs.nLocallyRepaired, epochs.nChannelEpochRepairs)
    sprintf('matched flagged channel-windows: %.1f%% -> %.1f%%', ...
        100 .* quality.metrics.preclean.artifact.channelWindowFraction, ...
        100 .* quality.metrics.postclean.artifact.channelWindowFraction)
    finalText
    };
lines = lines(~cellfun(@isempty, lines));

text(axisHandle, 0, 1, strjoin(lines, newline), 'VerticalAlignment', 'top', ...
    'FontName', 'FixedWidth', 'FontSize', 9, 'Interpreter', 'none');
title(axisHandle, 'Preprocessing ledger');
end

function textValue = nf_upper_frequency_text(value)
if isnumeric(value) && isscalar(value) && isfinite(value)
    textValue = sprintf('%.3g', value);
else
    textValue = 'open';
end
end

function textValue = nf_number_text(value)
if isnumeric(value) && isscalar(value) && isfinite(value)
    textValue = sprintf('%d', round(value));
else
    textValue = 'NA';
end
end

function textValue = nf_epoch_detector_ledger_text(epochs)
textValue = '';
if isempty(epochs.detectorMetrics)
    return
end
parts = cell(1, numel(epochs.detectorMetrics));
for index = 1:numel(epochs.detectorMetrics)
    metric = epochs.detectorMetrics(index);
    parts{index} = sprintf('%s %s', metric.displayName, ...
        nf_number_text(metric.nFlaggedEpochs));
end
textValue = ['detector flagged epochs: ' strjoin(parts, '; ')];
end

function nf_plot_psd(axisHandle, pre, post, preStage, postStage)
if ~pre.available || ~post.available
    nf_plot_message(axisHandle, 'Power spectral density', ...
        'One or both stage spectra were unavailable.');
    return
end

hold(axisHandle, 'on');
nf_plot_psd_band(axisHandle, pre, [0.80 0.25 0.20]);
nf_plot_psd_band(axisHandle, post, [0.10 0.35 0.70]);
xlabel(axisHandle, 'Frequency (Hz)');
ylabel(axisHandle, 'PSD (dB microV^2/Hz)');
title(axisHandle, 'Segment-safe spectrum: median and 10th-90th percentile');
legend(axisHandle, {[preStage ' range'], [preStage ' median'], ...
    [postStage ' range'], [postStage ' median']}, ...
    'Location', 'best', 'Box', 'off');
grid(axisHandle, 'on');
end

function nf_plot_psd_band(axisHandle, psd, color)
xValues = psd.frequencies;
patchX = [xValues, fliplr(xValues)];
patchY = [psd.lowerDb, fliplr(psd.upperDb)];
patch(axisHandle, patchX, patchY, color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none');
plot(axisHandle, xValues, psd.medianDb, 'Color', color, 'LineWidth', 1.5);
end

function nf_plot_ica(axisHandle, ica)
if ~ica.available
    nf_plot_message(axisHandle, 'Independent components', ...
        'No ICA action report was available.');
    return
end

if ~isempty(ica.retainedByClass)
    values = [ica.retainedByClass; ica.rejectedByClass];
    bar(axisHandle, values, 'stacked');
    set(axisHandle, 'XTick', [1 2], 'XTickLabel', {'Retained', 'Rejected'});
    legend(axisHandle, ica.classNames, 'Location', 'eastoutside', ...
        'FontSize', 6, 'Box', 'off');
else
    values = [ica.nRetained, ica.nRejected];
    bar(axisHandle, values, 'FaceColor', [0.35 0.55 0.75]);
    set(axisHandle, 'XTick', [1 2], 'XTickLabel', {'Retained', 'Rejected'});
end
ylabel(axisHandle, 'Components');
title(axisHandle, sprintf('ICA: %s / %s', ica.algorithm, ica.method), ...
    'Interpreter', 'none');
grid(axisHandle, 'on');
end

function nf_plot_heatmap(axisHandle, artifact, labels, severityCap, plotTitle)
imagesc(axisHandle, 1:artifact.windows.nAnalyzed, 1:numel(labels), artifact.severity);
set(axisHandle, 'YDir', 'normal');
caxis(axisHandle, [0 severityCap]);
colorbar(axisHandle);
xlabel(axisHandle, 'Analyzed window ordinal');
ylabel(axisHandle, 'Channel');
title(axisHandle, sprintf('%s (%d/%d windows)', plotTitle, ...
    artifact.windows.nAnalyzed, artifact.windows.nAvailable));

maximumTicks = 10;
if numel(labels) <= maximumTicks
    tickIndices = 1:numel(labels);
else
    tickIndices = unique(round(linspace(1, numel(labels), maximumTicks)));
end
set(axisHandle, 'YTick', tickIndices, 'YTickLabel', labels(tickIndices), ...
    'TickLabelInterpreter', 'none');
end

function nf_plot_topography(axisHandle, values, locations, plotTitle, limits)
if nf_topography_available(locations)
    try
        axes(axisHandle);
        topoplot(values, locations, 'electrodes', 'on', 'style', 'map', ...
            'conv', 'on');
        caxis(axisHandle, limits);
        colorbar(axisHandle);
        title(axisHandle, plotTitle);
        return
    catch
        cla(axisHandle);
    end
end

bar(axisHandle, values, 'FaceColor', [0.25 0.50 0.70], 'EdgeColor', 'none');
ylim(axisHandle, limits);
xlabel(axisHandle, 'Channel index');
title(axisHandle, [plotTitle ' (topography unavailable)']);
grid(axisHandle, 'on');
end

function available = nf_topography_available(locations)
available = exist('topoplot', 'file') == 2;
if ~available
    return
end

hasCartesian = true;
hasPolar = true;
for channelIndex = 1:numel(locations)
    cartesianValid = isfield(locations(channelIndex), 'X') && ...
        isfield(locations(channelIndex), 'Y') && ...
        isnumeric(locations(channelIndex).X) && ...
        isnumeric(locations(channelIndex).Y) && ...
        isscalar(locations(channelIndex).X) && ...
        isscalar(locations(channelIndex).Y) && ...
        isfinite(locations(channelIndex).X) && ...
        isfinite(locations(channelIndex).Y);
    polarValid = isfield(locations(channelIndex), 'theta') && ...
        isfield(locations(channelIndex), 'radius') && ...
        isnumeric(locations(channelIndex).theta) && ...
        isnumeric(locations(channelIndex).radius) && ...
        isscalar(locations(channelIndex).theta) && ...
        isscalar(locations(channelIndex).radius) && ...
        isfinite(locations(channelIndex).theta) && ...
        isfinite(locations(channelIndex).radius);
    hasCartesian = hasCartesian && cartesianValid;
    hasPolar = hasPolar && polarValid;
end
available = hasCartesian || hasPolar;
end

function nf_plot_channel_actions(axisHandle, locations, labels, actions)
categoryNames = {'interpolated', 'artifact', 'reference', ...
    'icaPreparation', 'rankDependent', 'prepDetected'};
displayNames = {'Interpolated', 'Artifact', 'Reference', ...
    'ICA preparation', 'Rank dependent', 'PREP detected'};
markers = {'+', 'o', 'd', 's', 'x', 'p'};
colors = [0.15 0.60 0.30; 0.85 0.15 0.10; 0.55 0.20 0.70; ...
    0.95 0.50 0.10; 0.10 0.40 0.85; 0.05 0.05 0.05];
categoryIndices = cell(1, numel(categoryNames));
reportedCategoryEntries = 0;
mappedCategoryEntries = 0;
for categoryIndex = 1:numel(categoryNames)
    category = actions.categories.(categoryNames{categoryIndex});
    categoryIndices{categoryIndex} = nf_channel_category_indices( ...
        category, labels);
    if ~strcmp(categoryNames{categoryIndex}, 'prepDetected')
        reportedCategoryEntries = reportedCategoryEntries + category.n;
        mappedCategoryEntries = mappedCategoryEntries + ...
            numel(categoryIndices{categoryIndex});
    end
end
unmappedCategoryEntries = max(0, ...
    reportedCategoryEntries - mappedCategoryEntries);

[xCoordinates, yCoordinates, coordinatesAvailable] = ...
    nf_topoplot_coordinates(locations);
if nf_topography_available(locations) && coordinatesAvailable
    try
        axes(axisHandle);
        topoplot(zeros(numel(locations), 1), locations, 'style', 'blank', ...
            'electrodes', 'on');
        hold(axisHandle, 'on');
        legendHandles = gobjects(0, 1);
        legendNames = {};
        for categoryIndex = 1:numel(categoryNames)
            indices = categoryIndices{categoryIndex};
            if isempty(indices)
                continue
            end
            handle = plot(axisHandle, xCoordinates(indices), ...
                yCoordinates(indices), 'LineStyle', 'none', ...
                'Marker', markers{categoryIndex}, ...
                'MarkerEdgeColor', colors(categoryIndex, :), ...
                'MarkerSize', 8, 'LineWidth', 1.8);
            legendHandles(end + 1, 1) = handle;
            legendNames{end + 1} = displayNames{categoryIndex};
        end
        if ~isempty(legendHandles)
            legend(axisHandle, legendHandles, legendNames, ...
                'Location', 'southoutside', 'FontSize', 6, 'Box', 'off');
        end
        title(axisHandle, sprintf( ...
            'Channel-action provenance (%d entries off montage)', ...
            unmappedCategoryEntries));
        return
    catch
        cla(axisHandle);
    end
end

hold(axisHandle, 'on');
legendHandles = gobjects(0, 1);
legendNames = {};
for categoryIndex = 1:numel(categoryNames)
    indices = categoryIndices{categoryIndex};
    if isempty(indices)
        continue
    end
    yValues = categoryIndex .* ones(size(indices));
    handle = plot(axisHandle, indices, yValues, 'LineStyle', 'none', ...
        'Marker', markers{categoryIndex}, ...
        'MarkerEdgeColor', colors(categoryIndex, :), ...
        'MarkerSize', 8, 'LineWidth', 1.8);
    legendHandles(end + 1, 1) = handle;
    legendNames{end + 1} = displayNames{categoryIndex};
end
xlim(axisHandle, [0.5 max(numel(labels) + 0.5, 1.5)]);
ylim(axisHandle, [0.5 numel(categoryNames) + 0.5]);
set(axisHandle, 'YTick', 1:numel(categoryNames), ...
    'YTickLabel', displayNames, 'TickLabelInterpreter', 'none');
xlabel(axisHandle, 'Channel index');
title(axisHandle, sprintf( ...
    'Channel actions (%d entries off montage; no topography)', ...
    unmappedCategoryEntries));
grid(axisHandle, 'on');
if ~isempty(legendHandles)
    legend(axisHandle, legendHandles, legendNames, ...
        'Location', 'southoutside', 'FontSize', 6, 'Box', 'off');
end
end

function indices = nf_channel_category_indices(category, labels)
indices = [];
if ~isempty(category.labels)
    indices = find(ismember(lower(string(labels)), ...
        lower(string(category.labels))));
elseif ~isempty(category.indices)
    indices = category.indices;
    indices = indices(indices >= 1 & indices <= numel(labels));
end
indices = unique(round(indices(:)'));
end

function [xCoordinates, yCoordinates, available] = nf_topoplot_coordinates(locations)
xCoordinates = nan(1, numel(locations));
yCoordinates = nan(1, numel(locations));
available = false;

polarValid = true;
for channelIndex = 1:numel(locations)
    polarValid = polarValid && ...
        isfield(locations(channelIndex), 'theta') && ...
        isfield(locations(channelIndex), 'radius') && ...
        isnumeric(locations(channelIndex).theta) && ...
        isnumeric(locations(channelIndex).radius) && ...
        isscalar(locations(channelIndex).theta) && ...
        isscalar(locations(channelIndex).radius) && ...
        isfinite(locations(channelIndex).theta) && ...
        isfinite(locations(channelIndex).radius);
end
if polarValid
    for channelIndex = 1:numel(locations)
        thetaRadians = locations(channelIndex).theta .* pi ./ 180;
        xCoordinates(channelIndex) = ...
            locations(channelIndex).radius .* sin(thetaRadians);
        yCoordinates(channelIndex) = ...
            locations(channelIndex).radius .* cos(thetaRadians);
    end
    available = true;
    return
end

cartesianValid = true;
for channelIndex = 1:numel(locations)
    cartesianValid = cartesianValid && ...
        isfield(locations(channelIndex), 'X') && ...
        isfield(locations(channelIndex), 'Y') && ...
        isnumeric(locations(channelIndex).X) && ...
        isnumeric(locations(channelIndex).Y) && ...
        isscalar(locations(channelIndex).X) && ...
        isscalar(locations(channelIndex).Y) && ...
        isfinite(locations(channelIndex).X) && ...
        isfinite(locations(channelIndex).Y);
end
if ~cartesianValid
    return
end

for channelIndex = 1:numel(locations)
    xCoordinates(channelIndex) = -locations(channelIndex).Y;
    yCoordinates(channelIndex) = locations(channelIndex).X;
end
maximumRadius = max(hypot(xCoordinates, yCoordinates));
if maximumRadius > 0
    xCoordinates = 0.5 .* xCoordinates ./ maximumRadius;
    yCoordinates = 0.5 .* yCoordinates ./ maximumRadius;
end
available = true;
end

function nf_plot_epochs(axisHandle, epochs)
if ~epochs.available
    nf_plot_message(axisHandle, 'Epoch disposition', ...
        'No complete epoch action report was available.');
    return
end

if ~isempty(epochs.badChannelEpoch) && ~isempty(epochs.rejectedMask)
    actionImage = double(logical(epochs.badChannelEpoch));
    if isequal(size(epochs.localInterpolationMask), size(actionImage))
        actionImage(logical(epochs.localInterpolationMask)) = 2;
    end
    rejectionStrip = 3 .* double(logical(epochs.rejectedMask(:)'));
    actionImage = [actionImage; rejectionStrip];
    imagesc(axisHandle, actionImage);
    set(axisHandle, 'YDir', 'normal');
    caxis(axisHandle, [-0.5 3.5]);
    colormap(axisHandle, [1 1 1; 0.95 0.65 0.20; 0.25 0.65 0.40; 0.85 0.15 0.10]);
    colorbarHandle = colorbar(axisHandle);
    colorbarHandle.Ticks = [0 1 2 3];
    colorbarHandle.TickLabels = {'Clean', 'Flagged', ...
        'Locally interpolated', 'Rejected'};
    xlabel(axisHandle, 'Original epoch');
    ylabel(axisHandle, 'Channels + reject strip');
    actionTitle = sprintf('Threshold actions: %d repaired; %d rejected', ...
        epochs.nLocallyRepaired, epochs.nRejected);
    detectorTitle = nf_epoch_detector_title_lines(epochs);
    title(axisHandle, [{actionTitle}, detectorTitle], ...
        'Interpreter', 'none');
    return
end

if isempty(epochs.dispositionCounts)
    nf_plot_message(axisHandle, 'Epoch disposition', ...
        'No complete epoch action report was available.');
    return
end

bar(axisHandle, epochs.dispositionCounts, 'FaceColor', [0.35 0.60 0.45]);
set(axisHandle, 'XTick', 1:numel(epochs.dispositionNames), ...
    'XTickLabel', epochs.dispositionNames, 'XTickLabelRotation', 25);
ylabel(axisHandle, 'Epochs');
actionTitle = sprintf('Epochs: voltage flags %d; spectral flags %d', ...
    epochs.nVoltageFlagged, epochs.nSpectralFlagged);
detectorTitle = nf_epoch_detector_title_lines(epochs);
title(axisHandle, [{actionTitle}, detectorTitle], 'Interpreter', 'none');
grid(axisHandle, 'on');
end

function lines = nf_epoch_detector_title_lines(epochs)
lines = {};
if isempty(epochs.detectorMetrics)
    return
end
parts = cell(1, numel(epochs.detectorMetrics));
for index = 1:numel(epochs.detectorMetrics)
    metric = epochs.detectorMetrics(index);
    parts{index} = sprintf('%s %s', ...
        nf_detector_short_name(metric.name), ...
        nf_number_text(metric.nFlaggedEpochs));
end
detectorsPerLine = 5;
for firstIndex = 1:detectorsPerLine:numel(parts)
    lastIndex = min(firstIndex + detectorsPerLine - 1, numel(parts));
    if firstIndex == 1
        prefix = 'flagged epochs: ';
    else
        prefix = '';
    end
    lines{end + 1} = [prefix strjoin(parts(firstIndex:lastIndex), ' | ')];
end
end

function name = nf_detector_short_name(detector)
switch detector
    case 'amplitude'
        name = 'amp';
    case 'peak2peak'
        name = 'p2p';
    case 'gradient'
        name = 'grad';
    case 'flatline'
        name = 'flat';
    case 'clipping'
        name = 'clip';
    case 'jointprobability'
        name = 'jointP';
    case 'eeglabstats'
        name = 'EEGLAB';
    otherwise
        name = detector;
end
end

function nf_plot_alerts(axisHandle, alerts)
axis(axisHandle, 'off');
if ~alerts.available
    textValue = alerts.reason;
elseif isempty(alerts.items)
    textValue = sprintf('No threshold-based alerts.\nResidual source: %s.', ...
        alerts.residualStage);
else
    messages = cell(1, numel(alerts.items));
    for itemIndex = 1:numel(alerts.items)
        messages{itemIndex} = ['- ' alerts.items(itemIndex).message];
    end
    textValue = strjoin(messages, newline);
end
text(axisHandle, 0, 1, textValue, 'VerticalAlignment', 'top', ...
    'Interpreter', 'none', 'FontSize', 9);
title(axisHandle, 'Transparent QC alerts');
end

function nf_plot_message(axisHandle, plotTitle, message)
axis(axisHandle, 'off');
text(axisHandle, 0.5, 0.5, message, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'Interpreter', 'none');
title(axisHandle, plotTitle, 'Interpreter', 'none');
end

function value = nf_get_numeric_field(source, fieldName, defaultValue)
if isfield(source, fieldName) && isnumeric(source.(fieldName))
    value = source.(fieldName);
else
    value = defaultValue;
end
end

function value = nf_get_numeric_scalar(source, fieldName, defaultValue)
value = defaultValue;
if ~isfield(source, fieldName)
    return
end
candidate = source.(fieldName);
if isnumeric(candidate) && isreal(candidate) && isscalar(candidate) && ...
        isfinite(candidate)
    value = double(candidate);
end
end

function valid = nf_positive_setting(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value > 0;
end

function value = nf_get_array_field(source, fieldName, defaultValue)
if isfield(source, fieldName) && ...
        (isnumeric(source.(fieldName)) || islogical(source.(fieldName)))
    value = source.(fieldName);
else
    value = defaultValue;
end
end

function value = nf_get_text_field(source, fieldName, defaultValue)
if isfield(source, fieldName) && ...
        (ischar(source.(fieldName)) || ...
        (isstring(source.(fieldName)) && isscalar(source.(fieldName))))
    value = char(source.(fieldName));
else
    value = defaultValue;
end
end

function value = nf_get_logical_field(source, fieldName, defaultValue)
value = defaultValue;
if ~isfield(source, fieldName)
    return
end
candidate = source.(fieldName);
if islogical(candidate) && isscalar(candidate)
    value = candidate;
elseif isnumeric(candidate) && isreal(candidate) && isscalar(candidate) && ...
        isfinite(candidate) && ismember(candidate, [0 1])
    value = logical(candidate);
else
    value = defaultValue;
end
end

function value = nf_get_label_field(source, fieldName)
if isfield(source, fieldName)
    value = nf_to_cellstr(source.(fieldName));
else
    value = {};
end
end

function values = nf_to_cellstr(value)
if isempty(value)
    values = {};
elseif ischar(value)
    values = {value};
elseif isstring(value)
    values = cellstr(value);
elseif iscell(value)
    values = cell(size(value));
    for valueIndex = 1:numel(value)
        values{valueIndex} = char(string(value{valueIndex}));
    end
else
    values = {};
end
values = values(:)';
end

function nf_validate_eeg(EEG, variableName)
requiredFields = {'data', 'srate', 'nbchan', 'pnts', 'trials', ...
    'chanlocs', 'xmin', 'xmax'};
for fieldIndex = 1:numel(requiredFields)
    fieldName = requiredFields{fieldIndex};
    if ~isfield(EEG, fieldName)
        error('nf_eegquality:InvalidEEG', '%s.%s is required.', variableName, fieldName);
    end
end

if ~isnumeric(EEG.data) || ~isreal(EEG.data)
    error('nf_eegquality:InvalidEEG', '%s.data must be real numeric EEG data.', variableName);
end
if EEG.nbchan < 3 || EEG.nbchan ~= size(EEG.data, 1)
    error('nf_eegquality:InvalidEEG', ...
        '%s must contain at least three channels and a consistent nbchan.', variableName);
end
if EEG.pnts < 2 || EEG.pnts ~= size(EEG.data, 2)
    error('nf_eegquality:InvalidEEG', '%s.pnts does not match EEG.data.', variableName);
end
if EEG.trials < 1 || EEG.trials ~= size(EEG.data, 3)
    error('nf_eegquality:InvalidEEG', '%s.trials does not match EEG.data.', variableName);
end
if ~isnumeric(EEG.srate) || ~isscalar(EEG.srate) || ...
        ~isfinite(EEG.srate) || EEG.srate <= 0
    error('nf_eegquality:InvalidEEG', '%s.srate must be positive.', variableName);
end
if numel(EEG.chanlocs) ~= EEG.nbchan
    error('nf_eegquality:InvalidEEG', ...
        '%s.chanlocs must contain one entry per channel.', variableName);
end
if ~isnumeric(EEG.xmin) || ~isscalar(EEG.xmin) || ~isfinite(EEG.xmin) || ...
        ~isnumeric(EEG.xmax) || ~isscalar(EEG.xmax) || ...
        ~isfinite(EEG.xmax) || EEG.xmin >= EEG.xmax
    error('nf_eegquality:InvalidEEG', ...
        '%s must contain finite increasing xmin/xmax values.', variableName);
end
end

function valid = nf_optional_struct(value)
valid = isempty(value) || isstruct(value);
end

function valid = nf_logical_scalar(value)
valid = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
    ismember(value, [0 1]);
end

function valid = nf_visibility(value)
valid = nf_text_scalar(value);
if valid
    valid = ismember(lower(char(value)), {'on', 'off'});
end
end

function valid = nf_positive_integer(value)
valid = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value >= 1 && value == round(value);
end

function valid = nf_frequency_pair(value)
valid = isnumeric(value) && numel(value) == 2 && all(isfinite(value)) && ...
    value(1) >= 0 && value(1) < value(2);
end

function valid = nf_text_scalar(value)
valid = ischar(value) || (isstring(value) && isscalar(value));
end

function nf_require_positive(value, name)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
    error('nf_eegquality:InvalidOption', '%s must be a positive scalar.', name);
end
end

function nf_require_finite(value, name)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    error('nf_eegquality:InvalidOption', '%s must be a finite scalar.', name);
end
end

function nf_require_fraction(value, name)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value < 0 || value > 1
    error('nf_eegquality:InvalidOption', '%s must be in [0, 1].', name);
end
end

function nf_require_band(value, name)
if ~isnumeric(value) || numel(value) ~= 2 || any(~isfinite(value)) || ...
        value(1) < 0 || value(1) >= value(2)
    error('nf_eegquality:InvalidOption', ...
        '%s must be a finite increasing two-element frequency band.', name);
end
end
