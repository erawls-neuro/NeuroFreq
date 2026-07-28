function [EEG, info, EEG_preRemoval] = nf_cleanic(EEG, method, aggressive, varargin)
% NF_CLEANIC  Fit one ICA decomposition and classify artifactual components.
%
% [EEG, INFO, EEG_PREREMOVAL] = NF_CLEANIC(EEG, METHOD, AGGRESSIVE, ...)
%
% METHOD is 'iclabel' for the adult preset or 'made'/'adjustedadjust' for
% the child preset. Both methods use the same training preparation and ICA;
% only the component classifier differs.
% The child path verifies MADE's required retained theta/radius scalp zones
% before adjusted_ADJUST so sparse or damaged montages cannot fail silently.
%
% ICA training follows the MADE intent: a 1-Hz filtered continuous copy is
% divided into one-second chunks, channels artifactual in more than 20% of
% chunks are removed, +/-1000-uV and 20-40-Hz FFT masks are unioned, and
% the union is rejected once. Rank-dependent channels are removed before
% ICA so adjusted_ADJUST receives the square decomposition it requires.
%
% Name/value inputs:
%   algorithm             'runica' (default) or optional 'runamica15'
%   randomSeed            1 (controls runica; AMICA uses internal seeding)
%   trainingHighpass      1 Hz
%   trainingEpochLength   1 second
%   trainingVoltage       1000 microvolts
%   trainingPower         [-100 30] dB
%   trainingFrequencies   [20 40] Hz
%   badChannelFraction    0.20
%   minimumTrainingEpochs 10
%   minimumSamplesPerRankSquared 20
%   iclabelThresholds     [] uses dominant-class NeuroFreq behavior, or 7x2
%   adjustReportFile      writable output path; temporary when empty
%   amicaMaxIterations    2000
%   amicaThreads          4
%   amicaProcesses        1

nf_validate_continuous_eeg(EEG);
EEG_preRemoval = [];

if nargin < 2 || isempty(method)
    method = 'iclabel';
end
if nargin < 3 || isempty(aggressive)
    aggressive = false;
end

parser = inputParser;
parser.FunctionName = 'nf_cleanic';
addParameter(parser, 'algorithm', 'runica', @nf_is_text);
addParameter(parser, 'randomSeed', 1, @nf_is_nonnegative_integer);
addParameter(parser, 'trainingHighpass', 1, @nf_is_positive_scalar);
addParameter(parser, 'trainingEpochLength', 1, @nf_is_positive_scalar);
addParameter(parser, 'trainingVoltage', 1000, @nf_is_positive_scalar);
addParameter(parser, 'trainingPower', [-100 30], @nf_is_increasing_pair);
addParameter(parser, 'trainingFrequencies', [20 40], @nf_is_increasing_pair);
addParameter(parser, 'badChannelFraction', 0.20, @nf_is_fraction);
addParameter(parser, 'minimumTrainingEpochs', 10, @nf_is_positive_integer);
addParameter(parser, 'minimumSamplesPerRankSquared', 20, ...
    @nf_is_positive_scalar);
addParameter(parser, 'iclabelThresholds', [], @nf_is_iclabel_thresholds);
addParameter(parser, 'adjustReportFile', '', @nf_is_text_or_empty);
addParameter(parser, 'amicaMaxIterations', 2000, @nf_is_positive_integer);
addParameter(parser, 'amicaThreads', 4, @nf_is_positive_integer);
addParameter(parser, 'amicaProcesses', 1, @nf_is_positive_integer);
addParameter(parser, 'runicaStop', 1e-7, @nf_is_positive_scalar);
parse(parser, varargin{:});
options = parser.Results;
options.trainingPower = reshape(options.trainingPower, 1, 2);
options.trainingFrequencies = reshape(options.trainingFrequencies, 1, 2);

method = nf_normalize_method(method);
algorithm = lower(char(options.algorithm));
nf_validate_options(EEG, method, algorithm, aggressive, options);
nf_preflight(method, algorithm, options);
EEG = nf_normalize_locations(EEG);

if strcmp(algorithm, 'runica')
    previousRandomState = rng;
    randomCleanup = onCleanup(@() rng(previousRandomState));
    rng(options.randomSeed, 'twister');
end

EEG = nf_clear_ica(EEG);
training = EEG;
training = pop_eegfiltnew(training, 'locutoff', options.trainingHighpass);
training = eeg_checkset(training);
training.datfile = '';
training.filepath = '';
training = nf_ensure_event_fields(training);
training = eeg_regepochs(training, ...
    'recurrence', options.trainingEpochLength, ...
    'limits', [0 options.trainingEpochLength], ...
    'rmbase', NaN, 'eventtype', 'nf_ica_training');
training = eeg_checkset(training);

[training, initialMasks] = nf_training_artifact_masks(training, options);
initialChannelFraction = mean(initialMasks.any, 2);
initialChannelLabels = {training.chanlocs.labels};
preparationBadIndices = find(initialChannelFraction > options.badChannelFraction);
preparationBadLabels = {training.chanlocs(preparationBadIndices).labels};

if numel(preparationBadIndices) >= training.nbchan - 2
    error('nf_cleanic:TooManyTrainingChannels', ...
        'ICA preparation would leave fewer than three channels.');
end
if ~isempty(preparationBadIndices)
    training = pop_select(training, 'nochannel', preparationBadIndices);
    EEG = pop_select(EEG, 'nochannel', preparationBadIndices);
    training = eeg_checkset(training);
    EEG = eeg_checkset(EEG);
    EEG = nf_append_bad_channels(EEG, preparationBadLabels);
end

[training, finalMasks] = nf_training_artifact_masks(training, options);
preRankChannelLabels = {training.chanlocs.labels};
preRankChannelArtifactFraction = mean(finalMasks.any, 2);
trainingRejectedMask = any(finalMasks.any, 1);
if all(trainingRejectedMask)
    error('nf_cleanic:NoCleanTrainingData', ...
        'Every ICA training epoch failed voltage or spectral criteria.');
end
training = pop_rejepoch(training, trainingRejectedMask, 0);
training = eeg_checkset(training);
if training.trials < options.minimumTrainingEpochs
    error('nf_cleanic:InsufficientTrainingData', ...
        ['Only %d clean one-second ICA training epochs remain; the configured ' ...
        'minimum is %d.'], training.trials, options.minimumTrainingEpochs);
end

[rankBeforeRepair, independentIndices, dependentIndices, singularValues, ...
    rankTolerance] = ...
    nf_rank_repair_plan(training);
dependentLabels = {training.chanlocs(dependentIndices).labels};
if ~isempty(dependentIndices)
    training = pop_select(training, 'channel', independentIndices);
    EEG = pop_select(EEG, 'channel', independentIndices);
    training = eeg_checkset(training);
    EEG = eeg_checkset(EEG);
end

[rankAfterRepair, ~, remainingDependent] = nf_rank_repair_plan(training);
if rankAfterRepair ~= training.nbchan || ~isempty(remainingDependent)
    error('nf_cleanic:RankRepairFailed', ...
        'ICA training data remain rank deficient after channel repair.');
end

trainingSamples = training.pnts * training.trials;
samplesPerRankSquared = trainingSamples / (rankAfterRepair ^ 2);
if samplesPerRankSquared < options.minimumSamplesPerRankSquared
    error('nf_cleanic:InsufficientICASamples', ...
        ['ICA training retained %d samples for rank %d (%.2f x rank^2); ' ...
        'the configured minimum is %.2f x rank^2.'], ...
        trainingSamples, rankAfterRepair, samplesPerRankSquared, ...
        options.minimumSamplesPerRankSquared);
end

training = nf_clear_ica(training);
EEG = nf_clear_ica(EEG);
[training, algorithmInfo] = nf_run_ica(training, algorithm, options);
nf_validate_decomposition(training);

EEG.icachansind = training.icachansind;
EEG.icasphere = training.icasphere;
EEG.icaweights = training.icaweights;
EEG.icawinv = training.icawinv;
EEG.icaact = [];
EEG = eeg_checkset(EEG, 'ica');

componentCount = size(EEG.icaweights, 1);
if strcmp(method, 'iclabel')
    [EEG, rejected, classification] = nf_classify_iclabel( ...
        EEG, logical(aggressive), options.iclabelThresholds);
else
    [rejected, classification] = nf_classify_adjusted_adjust( ...
        EEG, options.trainingEpochLength, options.adjustReportFile);
end

rejected = rejected(:)';
if isempty(rejected)
    rejected = zeros(1, 0);
end
if any(~isfinite(rejected)) || any(rejected ~= round(rejected)) || ...
        any(rejected < 1) || any(rejected > componentCount) || ...
        numel(unique(rejected)) ~= numel(rejected)
    error('nf_cleanic:InvalidClassifierOutput', ...
        'The component classifier returned invalid component indices.');
end
if numel(rejected) == componentCount
    error('nf_cleanic:AllComponentsRejected', ...
        'The classifier marked every independent component for removal.');
end

EEG.reject.gcompreject = false(1, componentCount);
EEG.reject.gcompreject(rejected) = true;
preRemovalSummary = nf_dataset_summary(EEG);
if nargout >= 3
    EEG_preRemoval = EEG;
end
if ~isempty(rejected)
    EEG = pop_subcomp(EEG, rejected, 0);
    EEG = eeg_checkset(EEG);
end

info = struct();
info.schemaVersion = '2.0.0';
info.method = method;
info.algorithm = algorithm;
info.algorithmDetails = algorithmInfo;
info.randomness.requestedSeed = options.randomSeed;
if strcmp(algorithm, 'runica')
    info.randomSeed = options.randomSeed;
    info.randomness.algorithmSeedControlled = true;
    info.randomness.control = 'MATLAB rng twister before runica';
else
    info.randomSeed = [];
    info.randomness.algorithmSeedControlled = false;
    info.randomness.control = ...
        'AMICA binary initialization is plugin/version specific';
end
info.preRemovalDataset = preRemovalSummary;
info.training.highpassHz = options.trainingHighpass;
info.training.epochLengthSeconds = options.trainingEpochLength;
info.training.voltageThresholdMicrovolts = options.trainingVoltage;
info.training.powerThresholdDb = options.trainingPower;
info.training.frequencyRangeHz = options.trainingFrequencies;
info.training.badChannelFractionLimit = options.badChannelFraction;
info.training.nInitialEpochs = numel(trainingRejectedMask);
info.training.rejectedEpochMask = trainingRejectedMask;
info.training.nRejectedEpochs = sum(trainingRejectedMask);
info.training.nRetainedEpochs = training.trials;
info.training.nRetainedSamples = trainingSamples;
info.training.minimumSamplesPerRankSquared = ...
    options.minimumSamplesPerRankSquared;
info.training.samplesPerRankSquared = samplesPerRankSquared;
info.training.initialChannelArtifactFraction = initialChannelFraction;
info.training.initialChannelLabels = initialChannelLabels;
info.training.preparationBadChannelIndices = preparationBadIndices;
info.training.preparationBadChannelLabels = preparationBadLabels;
info.training.preRankChannelArtifactFraction = ...
    preRankChannelArtifactFraction;
info.training.preRankChannelLabels = preRankChannelLabels;
info.training.finalChannelArtifactFraction = ...
    preRankChannelArtifactFraction(independentIndices);
info.training.finalChannelLabels = {training.chanlocs.labels};
info.rank.beforeRepair = rankBeforeRepair;
info.rank.afterRepair = rankAfterRepair;
info.rank.singularValues = singularValues;
info.rank.tolerance = rankTolerance;
info.rank.dependentChannelIndices = dependentIndices;
info.rank.dependentChannelLabels = dependentLabels;
info.rank.finalChannels = training.nbchan;
info.components.nComponents = componentCount;
info.components.rejected = rejected;
info.components.nRejected = numel(rejected);
info.components.nRetained = componentCount - numel(rejected);
info.components.aggressive = logical(aggressive);
info.classification = classification;
info.outputChannelLabels = {EEG.chanlocs.labels};

if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
EEG.etc.nf_cleanic = info;
if ~isfield(EEG.etc, 'nf_cleanic_history')
    EEG.etc.nf_cleanic_history = {};
elseif ~iscell(EEG.etc.nf_cleanic_history)
    EEG.etc.nf_cleanic_history = {EEG.etc.nf_cleanic_history};
end
historyEntry = struct();
historyEntry.schemaVersion = info.schemaVersion;
historyEntry.method = info.method;
historyEntry.algorithm = info.algorithm;
historyEntry.randomSeed = info.randomSeed;
historyEntry.trainingEpochsRetained = info.training.nRetainedEpochs;
historyEntry.rank = info.rank.afterRepair;
historyEntry.components = info.components;
EEG.etc.nf_cleanic_history{end + 1} = historyEntry;

clear randomCleanup

end

function [training, masks] = nf_training_artifact_masks(training, options)
training = nf_clear_rejection_marks(training);
training = pop_eegthresh(training, 1, 1:training.nbchan, ...
    -options.trainingVoltage, options.trainingVoltage, ...
    training.xmin, training.xmax, 0, 0);
training = eeg_checkset(training);
voltage = nf_get_rejection_matrix(training, 'rejthreshE');

training.specdata = [];
training = pop_rejspec(training, 1, 'elecrange', 1:training.nbchan, ...
    'method', 'fft', 'threshold', options.trainingPower, ...
    'freqlimits', options.trainingFrequencies, 'specdata', [], ...
    'eegplotplotallrej', 0, 'eegplotreject', 0);
training = eeg_checkset(training);
spectral = nf_get_rejection_matrix(training, 'rejfreqE');
training.specdata = [];

masks = struct();
masks.voltage = voltage;
masks.spectral = spectral;
masks.any = voltage | spectral;
end

function matrix = nf_get_rejection_matrix(EEG, fieldName)
if ~isfield(EEG, 'reject') || ~isfield(EEG.reject, fieldName)
    error('nf_cleanic:MissingRejectionFlags', ...
        'EEGLAB did not create EEG.reject.%s.', fieldName);
end
matrix = logical(EEG.reject.(fieldName));
if EEG.trials ~= EEG.nbchan && ...
        isequal(size(matrix), [EEG.trials EEG.nbchan])
    matrix = matrix';
end
if ~isequal(size(matrix), [EEG.nbchan EEG.trials])
    error('nf_cleanic:InvalidRejectionFlags', ...
        'EEG.reject.%s has unexpected dimensions.', fieldName);
end
end

function EEG = nf_clear_rejection_marks(EEG)
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

function [rankValue, independentIndices, dependentIndices, singularValues, ...
    tolerance] = ...
    nf_rank_repair_plan(EEG)
data = reshape(EEG.data, EEG.nbchan, []);
dataClass = class(data);
maximumSamples = 100000;
if size(data, 2) > maximumSamples
    sampleIndices = unique(round(linspace(1, size(data, 2), maximumSamples)));
    data = data(:, sampleIndices);
end
data = double(data);
data = data - mean(data, 2);
if any(~isfinite(data(:)))
    error('nf_cleanic:NonfiniteTrainingData', ...
        'ICA training data contain NaN or Inf values.');
end

singularValues = svd(data, 'econ');
if isempty(singularValues) || singularValues(1) <= 0
    rankValue = 0;
    tolerance = 0;
else
    if isa(dataClass, 'single')
        precisionEpsilon = double(eps('single'));
    else
        precisionEpsilon = eps('double');
    end
    doubleTolerance = max(size(data)) * eps('double') * singularValues(1);
    sourceTolerance = 10 * precisionEpsilon * singularValues(1);
    tolerance = max(doubleTolerance, sourceTolerance);
    rankValue = sum(singularValues > tolerance);
end
if rankValue < 2
    error('nf_cleanic:RankDeficient', ...
        'ICA training data have rank below two.');
end

[~, ~, pivot] = qr(data', 'vector');
independentIndices = sort(pivot(1:rankValue));
dependentIndices = setdiff(1:EEG.nbchan, independentIndices, 'stable');
end

function [training, details] = nf_run_ica(training, algorithm, options)
details = struct();
details.started = datestr(now, 30); %#ok<TNOW1,DATST>
details.nChannels = training.nbchan;
details.nSamples = training.pnts * training.trials;
details.samplesPerChannelSquared = details.nSamples / (training.nbchan ^ 2);

if strcmp(algorithm, 'runamica15')
    amicaOutputDirectory = tempname;
    amicaCleanup = onCleanup(@() nf_remove_temporary_directory(amicaOutputDirectory));
    training = pop_runica(training, 'icatype', 'runamica15', ...
        'outdir', amicaOutputDirectory, ...
        'pcakeep', training.nbchan, ...
        'maxiter', options.amicaMaxIterations, ...
        'max_threads', options.amicaThreads, ...
        'numprocs', options.amicaProcesses);
    details.maxIterations = options.amicaMaxIterations;
    details.maxThreads = options.amicaThreads;
    details.processes = options.amicaProcesses;
    details.temporaryOutputRemoved = true;
    details.algorithmSeedControlled = false;
    details.seedControl = ...
        'AMICA binary initialization is plugin/version specific';
    clear amicaCleanup
elseif strcmp(algorithm, 'runica')
    training = pop_runica(training, 'icatype', 'runica', ...
        'extended', 1, 'stop', options.runicaStop, 'interrupt', 'off');
    details.extended = 1;
    details.stop = options.runicaStop;
    details.algorithmSeedControlled = true;
    details.randomSeed = options.randomSeed;
    details.seedControl = 'MATLAB rng twister';
else
    error('nf_cleanic:UnknownAlgorithm', ...
        'algorithm must be ''runamica15'' or ''runica''.');
end
training = eeg_checkset(training, 'ica');
details.finished = datestr(now, 30); %#ok<TNOW1,DATST>
end

function [EEG, rejected, classification] = nf_classify_iclabel( ...
    EEG, aggressive, thresholds)
EEG = pop_iclabel(EEG, 'default');
probabilities = EEG.etc.ic_classification.ICLabel.classifications;
classNames = EEG.etc.ic_classification.ICLabel.classes;
if size(probabilities, 2) ~= 7
    error('nf_cleanic:InvalidICLabelOutput', ...
        'ICLabel did not return seven class probabilities.');
end

[winningProbability, winningClass] = max(probabilities, [], 2);
if isempty(thresholds)
    if aggressive
        rejected = find(winningClass ~= 1);
        rule = 'dominant-class; retain Brain only';
    else
        rejected = find(~ismember(winningClass, [1 7]));
        rule = 'dominant-class; retain Brain and Other';
    end
else
    EEG = pop_icflag(EEG, thresholds);
    rejected = find(logical(EEG.reject.gcompreject));
    rule = 'explicit ICLabel probability thresholds';
end

classification = struct();
classification.name = 'ICLabel';
classification.classNames = classNames;
classification.probabilities = probabilities;
classification.winningClass = winningClass;
classification.winningProbability = winningProbability;
classification.rule = rule;
classification.thresholds = thresholds;
classification.rejectedByWinningClass = nf_class_counts( ...
    winningClass(rejected), numel(classNames));
end

function [rejected, classification] = nf_classify_adjusted_adjust( ...
    EEG, epochLength, requestedReportFile)
if size(EEG.icaweights, 1) ~= size(EEG.icaweights, 2) || ...
        size(EEG.icaweights, 1) ~= EEG.nbchan
    error('nf_cleanic:AdjustedAdjustRequiresSquareICA', ...
        'adjusted_ADJUST requires a square full-channel ICA decomposition.');
end
zoneCoverage = nf_adjusted_adjust_zone_coverage(EEG.chanlocs);

classifierEEG = nf_ensure_event_fields(EEG);
classifierEEG = eeg_regepochs(classifierEEG, 'recurrence', epochLength, ...
    'limits', [0 epochLength], 'rmbase', NaN, ...
    'eventtype', 'nf_adjust_training');
classifierEEG = eeg_checkset(classifierEEG, 'ica');
classifierEEG.setname = 'nf_adjusted_adjust';
classifierEEG.filename = 'nf_adjusted_adjust.set';
classifierEEG.filepath = '';
classifierEEG.datfile = '';

temporaryReport = isempty(requestedReportFile);
if temporaryReport
    reportFile = [tempname '_adjust_report.txt'];
else
    reportFile = nf_absolute_path(char(requestedReportFile));
    nf_validate_report_path(reportFile);
end
reportCleanup = onCleanup(@() nf_remove_temporary_report(reportFile, temporaryReport));

workingDirectory = tempname;
[created, message] = mkdir(workingDirectory);
if ~created
    error('nf_cleanic:TemporaryDirectoryFailed', ...
        'Could not create the adjusted_ADJUST working directory: %s', message);
end
previousDirectory = pwd;
workingCleanup = onCleanup(@() nf_restore_and_remove_directory( ...
    previousDirectory, workingDirectory));
cd(workingDirectory);

figuresBefore = findall(0, 'Type', 'figure');
figureCleanup = onCleanup(@() nf_close_new_figures(figuresBefore));
classifierRandomState = rng;
classifierWarningState = warning;
classifierPathState = path;
environmentCleanup = onCleanup(@() nf_restore_environment( ...
    classifierRandomState, classifierWarningState, classifierPathState));
[rejected, horizontal, vertical, blink, discontinuity] = ...
    adjusted_ADJUST(classifierEEG, reportFile);
clear environmentCleanup
figuresAfter = findall(0, 'Type', 'figure');
newFigures = setdiff(figuresAfter, figuresBefore);
for index = 1:numel(newFigures)
    if ishghandle(newFigures(index))
        close(newFigures(index));
    end
end
clear figureCleanup
clear workingCleanup

classification = struct();
classification.name = 'MADE adjusted_ADJUST';
classification.horizontalEyeMovement = horizontal(:)';
classification.verticalEyeMovement = vertical(:)';
classification.blink = blink(:)';
classification.discontinuity = discontinuity(:)';
classification.montageZoneCoverage = zoneCoverage;
classification.reportWasTemporary = temporaryReport;
if temporaryReport
    classification.reportFile = '';
else
    classification.reportFile = reportFile;
end

clear reportCleanup
end

function coverage = nf_adjusted_adjust_zone_coverage(chanlocs)
theta = double([chanlocs.theta]);
radius = double([chanlocs.radius]);

horizontalLeft = theta > -62 & theta < -35 & radius > 0.5;
horizontalRight = theta > 35 & theta < 62 & radius > 0.5;
center = abs(theta) > 35 & abs(theta) < 109 & radius < 0.45;
backLeft = theta <= -109 & radius < 0.55;
backRight = theta >= 109 & radius < 0.55;
blinkLeft = theta > -60 & theta < 0 & ...
    radius > 0.45 & radius < 0.60;
blinkRight = theta > 0 & theta < 60 & ...
    radius > 0.45 & radius < 0.60;
blinkCenterPrimary = abs(theta) < 20 & radius > 0.45;
blinkCenterFallback = theta == 0 & radius > 0.39;
if any(blinkCenterPrimary)
    blinkCenter = blinkCenterPrimary;
    blinkCenterSource = 'primary';
else
    blinkCenter = blinkCenterFallback;
    blinkCenterSource = 'fallback';
end

requiredNames = {'horizontal-left-eye', 'horizontal-right-eye', ...
    'blink-left-eye', 'blink-right-eye', 'blink-center-eye', ...
    'shared-center', 'shared-back-left', 'shared-back-right'};
requiredCounts = [sum(horizontalLeft), sum(horizontalRight), ...
    sum(blinkLeft), sum(blinkRight), sum(blinkCenter), sum(center), ...
    sum(backLeft), sum(backRight)];
emptyRequired = requiredNames(requiredCounts == 0);
if ~isempty(emptyRequired)
    error('nf_cleanic:AdjustedAdjustMontageCoverage', ...
        ['The retained ICA montage leaves MADE adjusted_ADJUST scalp ' ...
        'zone(s) empty: %s. Restore the affected channels or use a ' ...
        'supported denser montage; MADE otherwise divides by zero and can ' ...
        'silently miss eye components.'], strjoin(emptyRequired, ', '));
end

horizontalSpatialEye = radius > 0.51 & radius < 0.60 & ...
    abs(theta) > 27 & abs(theta) < 58;
blinkSpatialEye = radius > 0.45 & radius < 0.54 & abs(theta) < 17;
spatialOuterRing = radius > 0.51;
spatialHigherThreshold = radius > 0.35 & radius < 0.5 & ...
    abs(theta) < 20;

coverage = struct();
coverage.source = 'MADE adjusted_ADJUST exact theta/radius inequalities';
coverage.horizontal.leftEye = sum(horizontalLeft);
coverage.horizontal.rightEye = sum(horizontalRight);
coverage.horizontal.center = sum(center);
coverage.horizontal.backLeft = sum(backLeft);
coverage.horizontal.backRight = sum(backRight);
coverage.blink.leftEye = sum(blinkLeft);
coverage.blink.rightEye = sum(blinkRight);
coverage.blink.centerEye = sum(blinkCenter);
coverage.blink.centerEyeSource = blinkCenterSource;
coverage.blink.center = sum(center);
coverage.blink.backLeft = sum(backLeft);
coverage.blink.backRight = sum(backRight);
coverage.secondary.horizontalSpatialEye = sum(horizontalSpatialEye);
coverage.secondary.blinkSpatialEye = sum(blinkSpatialEye);
coverage.secondary.spatialOuterRing = sum(spatialOuterRing);
coverage.secondary.spatialHigherThreshold = sum(spatialHigherThreshold);
coverage.emptyRequiredZones = emptyRequired;

secondaryNames = {'horizontal-spatial-eye', 'blink-spatial-eye', ...
    'spatial-outer-ring', 'spatial-higher-threshold'};
secondaryCounts = [coverage.secondary.horizontalSpatialEye, ...
    coverage.secondary.blinkSpatialEye, coverage.secondary.spatialOuterRing, ...
    coverage.secondary.spatialHigherThreshold];
coverage.emptySecondaryZones = secondaryNames(secondaryCounts == 0);
if ~isempty(coverage.emptySecondaryZones)
    warning('nf_cleanic:SparseAdjustedAdjustSpatialZones', ...
        ['The retained montage leaves secondary adjusted_ADJUST spatial ' ...
        'zone(s) empty: %s. Core classification can run, but spatial ' ...
        'confirmation may be less informative.'], ...
        strjoin(coverage.emptySecondaryZones, ', '));
end
end

function absolutePath = nf_absolute_path(pathValue)
directory = fileparts(pathValue);
isAbsolute = ~isempty(regexp(pathValue, ...
    '^(?:[A-Za-z]:[\\/]|[\\/]{1,2})', 'once'));
if isempty(directory) || ~isAbsolute
    absolutePath = fullfile(pwd, pathValue);
else
    absolutePath = pathValue;
end
end

function nf_restore_and_remove_directory(previousDirectory, workingDirectory)
try
    cd(previousDirectory);
catch
end
nf_remove_temporary_directory(workingDirectory);
end

function nf_restore_environment(randomState, warningState, pathState)
try
    rng(randomState);
catch
end
try
    path(pathState);
catch
end
try
    warning(warningState);
catch
end
end

function counts = nf_class_counts(classes, classCount)
counts = zeros(1, classCount);
for index = 1:classCount
    counts(index) = sum(classes == index);
end
end

function summary = nf_dataset_summary(EEG)
summary = struct();
summary.nbchan = EEG.nbchan;
summary.channelLabels = {EEG.chanlocs.labels};
summary.pntsPerEpoch = EEG.pnts;
summary.trials = EEG.trials;
summary.srate = EEG.srate;
summary.aggregateDataSeconds = (EEG.pnts * EEG.trials) / EEG.srate;
summary.epochLimitsSeconds = [EEG.xmin EEG.xmax];
summary.nComponents = size(EEG.icaweights, 1);
end

function nf_validate_decomposition(EEG)
if isempty(EEG.icaweights) || isempty(EEG.icasphere) || ...
        isempty(EEG.icawinv) || isempty(EEG.icachansind)
    error('nf_cleanic:ICAFailed', ...
        'ICA did not return a complete EEGLAB decomposition.');
end
values = [EEG.icaweights(:); EEG.icasphere(:); EEG.icawinv(:)];
if any(~isfinite(values))
    error('nf_cleanic:ICAFailed', ...
        'ICA returned NaN or Inf values.');
end
if size(EEG.icaweights, 1) ~= EEG.nbchan || ...
        size(EEG.icaweights, 2) ~= EEG.nbchan || ...
        numel(EEG.icachansind) ~= EEG.nbchan
    error('nf_cleanic:NonSquareICA', ...
        'The final ICA is not square across all retained channels.');
end
end

function EEG = nf_append_bad_channels(EEG, labels)
if isempty(labels)
    return
end
if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
if ~isfield(EEG.etc, 'badchanlabels') || isempty(EEG.etc.badchanlabels)
    EEG.etc.badchanlabels = {};
end
existingLabels = EEG.etc.badchanlabels;
if isstring(existingLabels)
    existingLabels = cellstr(existingLabels);
elseif ischar(existingLabels)
    existingLabels = {existingLabels};
end
combined = [existingLabels(:)' labels(:)'];
[~, uniqueIndices] = unique(lower(string(combined)), 'stable');
EEG.etc.badchanlabels = combined(sort(uniqueIndices));
EEG.etc.badchans = numel(EEG.etc.badchanlabels);
if isfield(EEG.etc, 'ogchan') && isstruct(EEG.etc.ogchan) && ...
        isfield(EEG.etc.ogchan, 'labels')
    originalLabels = lower(string({EEG.etc.ogchan.labels}));
    [present, originalIndices] = ismember( ...
        lower(string(EEG.etc.badchanlabels)), originalLabels);
    if all(present)
        EEG.etc.badchanindices = originalIndices;
    elseif isfield(EEG.etc, 'badchanindices')
        EEG.etc = rmfield(EEG.etc, 'badchanindices');
    end
elseif isfield(EEG.etc, 'badchanindices')
    EEG.etc = rmfield(EEG.etc, 'badchanindices');
end
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

function method = nf_normalize_method(method)
method = lower(char(method));
if ismember(method, {'made', 'adjustedadjust', 'adjustedadust', ...
        'adjusted_adjust', 'adjusted_adust', 'adjusted-adjust', ...
        'adjusted-adust'})
    method = 'adjustedadjust';
elseif ~strcmp(method, 'iclabel')
    error('nf_cleanic:UnknownMethod', ...
        'method must be ''iclabel'' or MADE adjusted_ADJUST.');
end
end

function nf_preflight(method, algorithm, options)
required = {'pop_eegfiltnew', 'eeg_regepochs', 'pop_eegthresh', ...
    'pop_rejspec', 'pop_rejepoch', 'pop_runica', 'pop_subcomp', ...
    'convertlocs'};
for index = 1:numel(required)
    if exist(required{index}, 'file') ~= 2
        error('nf_cleanic:MissingDependency', ...
            'Required EEGLAB function %s was not found.', required{index});
    end
end
if strcmp(algorithm, 'runamica15') && ...
        (exist('pop_runamica', 'file') ~= 2 || exist('runamica15', 'file') ~= 2)
    error('nf_cleanic:MissingAMICA', ...
        'The AMICA plugin with pop_runamica.m and runamica15.m is required.');
end
if strcmp(method, 'iclabel')
    if exist('pop_iclabel', 'file') ~= 2
        error('nf_cleanic:MissingICLabel', ...
            'The ICLabel EEGLAB plugin is required for the adult preset.');
    end
    if ~isempty(options.iclabelThresholds) && exist('pop_icflag', 'file') ~= 2
        error('nf_cleanic:MissingICLabel', ...
            'pop_icflag.m is required when iclabelThresholds are supplied.');
    end
elseif exist('adjusted_ADJUST', 'file') ~= 2
    error('nf_cleanic:MissingAdjustedAdjust', ...
        'adjusted_ADJUST.m from the MADE pipeline is required for the child preset.');
end
if strcmp(method, 'adjustedadjust')
    companions = {'compute_GD_feat', 'computeSED_NOnorm', 'computeSAD', ...
        'kurt', 'trim_and_mean', 'trim_and_max', 'EM', ...
        'MARA_extract_time_freq_features', 'beall_horizontal', ...
        'beall_blink_detection', 'Spatial_Info_eyes', 'spectopo', ...
        'fitlm', 'findpeaks'};
    for index = 1:numel(companions)
        if exist(companions{index}, 'file') ~= 2
            error('nf_cleanic:MissingAdjustedAdjustCompanion', ...
                ['The MADE adjusted_ADJUST companion %s.m was not found. ' ...
                'Add the complete adjusted_adjust_scripts folder to the path.'], ...
                companions{index});
        end
    end
end
end

function EEG = nf_ensure_event_fields(EEG)
if ~isfield(EEG, 'event')
    EEG.event = [];
end
if ~isfield(EEG, 'urevent')
    EEG.urevent = [];
end
if ~isfield(EEG, 'setname') || ...
        ~(ischar(EEG.setname) || ...
        (isstring(EEG.setname) && isscalar(EEG.setname))) || ...
        isempty(strtrim(char(EEG.setname)))
    EEG.setname = 'nf_dataset';
else
    EEG.setname = char(EEG.setname);
end
end

function nf_validate_options(EEG, method, algorithm, aggressive, options)
if ~ismember(algorithm, {'runamica15', 'runica'})
    error('nf_cleanic:UnknownAlgorithm', ...
        'algorithm must be ''runamica15'' or ''runica''.');
end
if ~nf_is_logical_scalar(aggressive)
    error('nf_cleanic:InvalidAggressiveFlag', ...
        'aggressive must be a logical scalar.');
end
if options.trainingHighpass >= EEG.srate / 2
    error('nf_cleanic:InvalidTrainingFilter', ...
        'trainingHighpass must be below Nyquist.');
end
if options.trainingFrequencies(1) < 0 || ...
        options.trainingFrequencies(2) >= EEG.srate / 2
    error('nf_cleanic:InvalidTrainingFrequencies', ...
        'trainingFrequencies must be nonnegative and below Nyquist.');
end
if strcmp(method, 'adjustedadjust') && ~isempty(options.iclabelThresholds)
    error('nf_cleanic:ClassifierOptionMismatch', ...
        'iclabelThresholds apply only to the ICLabel method.');
end
if strcmp(method, 'adjustedadjust') && logical(aggressive)
    error('nf_cleanic:ClassifierOptionMismatch', ...
        'aggressive applies only to the ICLabel method.');
end
if strcmp(method, 'adjustedadjust') && EEG.srate < 100
    error('nf_cleanic:ChildSamplingRate', ...
        ['The MADE adjusted_ADJUST feature extractor requires a sampling ' ...
        'rate of at least 100 Hz.']);
end
if strcmp(method, 'adjustedadjust') && ...
        options.trainingEpochLength * EEG.srate < 100
    error('nf_cleanic:ChildEpochLength', ...
        ['trainingEpochLength must yield at least 100 samples per MADE ' ...
        'adjusted_ADJUST feature epoch.']);
end
end

function nf_validate_report_path(reportFile)
[directory, ~, ~] = fileparts(reportFile);
if isempty(directory)
    directory = pwd;
end
if exist(directory, 'dir') ~= 7
    error('nf_cleanic:InvalidAdjustReportPath', ...
        'The adjusted_ADJUST report directory does not exist.');
end
testFile = [tempname(directory) '.tmp'];
fileId = fopen(testFile, 'w');
if fileId < 0
    error('nf_cleanic:InvalidAdjustReportPath', ...
        'The adjusted_ADJUST report directory is not writable.');
end
fclose(fileId);
delete(testFile);
end

function nf_remove_temporary_report(reportFile, temporaryReport)
if temporaryReport && exist(reportFile, 'file') == 2
    try
        delete(reportFile);
    catch
    end
end
end

function nf_close_new_figures(figuresBefore)
figuresAfter = findall(0, 'Type', 'figure');
newFigures = setdiff(figuresAfter, figuresBefore);
for index = 1:numel(newFigures)
    if ishghandle(newFigures(index))
        close(newFigures(index));
    end
end
end

function nf_remove_temporary_directory(directory)
if exist(directory, 'dir') == 7
    try
        rmdir(directory, 's');
    catch
    end
end
end

function EEG = nf_normalize_locations(EEG)
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
    error('nf_cleanic:MissingChannelCoordinates', ...
        ['Every ICA channel needs finite nonzero X/Y/Z or finite ' ...
        'theta/radius coordinates for component classification.']);
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
labels = cell(1, numel(chanlocs));
for index = 1:numel(chanlocs)
    if ~isfield(chanlocs, 'labels') || ...
            ~(ischar(chanlocs(index).labels) || ...
            (isstring(chanlocs(index).labels) && ...
            isscalar(chanlocs(index).labels))) || ...
            isempty(strtrim(char(chanlocs(index).labels)))
        error('nf_cleanic:InvalidChannelLabels', ...
            'Every ICA channel must have a nonempty scalar-text label.');
    end
    labels{index} = char(chanlocs(index).labels);
    if ~nf_has_xyz(chanlocs(index)) || ...
            ~nf_has_topography(chanlocs(index))
        error('nf_cleanic:MissingChannelCoordinates', ...
            'Channel %s lacks complete classifier geometry.', labels{index});
    end
    coordinates(index, :) = double([chanlocs(index).X ...
        chanlocs(index).Y chanlocs(index).Z]);
    coordinates(index, :) = coordinates(index, :) ./ ...
        norm(coordinates(index, :));
end
if numel(unique(lower(string(labels)))) ~= numel(labels)
    error('nf_cleanic:InvalidChannelLabels', ...
        'ICA channel labels must be unique, ignoring case.');
end
for first = 1:size(coordinates, 1) - 1
    distances = sqrt(sum((coordinates(first + 1:end, :) - ...
        coordinates(first, :)) .^ 2, 2));
    if any(distances < 1e-10)
        second = first + find(distances < 1e-10, 1);
        error('nf_cleanic:DuplicateChannelCoordinates', ...
            'Channels %s and %s occupy the same scalp coordinate.', ...
            labels{first}, labels{second});
    end
end
end

function nf_validate_continuous_eeg(EEG)
if ~isstruct(EEG) || numel(EEG) ~= 1 || ~isfield(EEG, 'data') || ...
        ~isfield(EEG, 'srate') || ~isfield(EEG, 'nbchan') || ...
        ~isfield(EEG, 'pnts') || ~isfield(EEG, 'trials') || ...
        ~isfield(EEG, 'chanlocs') || ~isfield(EEG, 'xmin') || ...
        ~isfield(EEG, 'xmax')
    error('nf_cleanic:InvalidEEG', ...
        'EEG must be one valid EEGLAB dataset structure.');
end
if ~isnumeric(EEG.trials) || ~isscalar(EEG.trials) || ...
        ~isfinite(EEG.trials) || EEG.trials ~= 1 || ...
        (~ismatrix(EEG.data) && size(EEG.data, 3) ~= 1)
    error('nf_cleanic:EpochedInput', ...
        'nf_cleanic requires continuous data so filtering does not create epoch-edge artifacts.');
end
if ~isnumeric(EEG.nbchan) || ~isscalar(EEG.nbchan) || ...
        ~isfinite(EEG.nbchan) || EEG.nbchan < 3 || ...
        EEG.nbchan ~= round(EEG.nbchan) || ...
        ~isnumeric(EEG.pnts) || ~isscalar(EEG.pnts) || ...
        ~isfinite(EEG.pnts) || EEG.pnts < 2 || EEG.pnts ~= round(EEG.pnts) || ...
        ~isnumeric(EEG.srate) || ~isscalar(EEG.srate) || ...
        ~isfinite(EEG.srate) || EEG.srate <= 0 || ...
        ~isnumeric(EEG.data) || ~isreal(EEG.data) || isempty(EEG.data) || ...
        size(EEG.data, 1) ~= EEG.nbchan || size(EEG.data, 2) ~= EEG.pnts || ...
        numel(EEG.chanlocs) ~= EEG.nbchan || any(~isfinite(EEG.data(:)))
    error('nf_cleanic:InvalidEEG', ...
        'EEG data, channel count, or channel locations are inconsistent.');
end
if nf_contains_boundary_event(EEG)
    error('nf_cleanic:BoundaryUnsupported', ...
        ['nf_cleanic filters and trains ICA on a continuous stream and cannot ' ...
        'cross EEGLAB boundary events. Segment the recording first.']);
end
end

function present = nf_contains_boundary_event(EEG)
present = false;
if ~isfield(EEG, 'event') || isempty(EEG.event) || ...
        ~isfield(EEG.event, 'type')
    return
end
for index = 1:numel(EEG.event)
    value = EEG.event(index).type;
    if (ischar(value) || (isstring(value) && isscalar(value))) && ...
            strcmpi(strtrim(char(value)), 'boundary')
        present = true;
        return
    end
end
end

function valid = nf_is_text(value)
valid = ischar(value) || (isstring(value) && isscalar(value));
end

function valid = nf_is_text_or_empty(value)
valid = isempty(value) || nf_is_text(value);
end

function valid = nf_is_positive_scalar(value)
valid = isnumeric(value) && isscalar(value) && isfinite(value) && value > 0;
end

function valid = nf_is_positive_integer(value)
valid = nf_is_positive_scalar(value) && value == round(value);
end

function valid = nf_is_nonnegative_integer(value)
valid = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value >= 0 && value == round(value);
end

function valid = nf_is_increasing_pair(value)
valid = isnumeric(value) && numel(value) == 2 && all(isfinite(value)) && ...
    value(1) < value(2);
end

function valid = nf_is_fraction(value)
valid = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value >= 0 && value < 1;
end

function valid = nf_is_logical_scalar(value)
valid = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
    ismember(value, [0 1]);
end

function valid = nf_is_iclabel_thresholds(value)
if isempty(value)
    valid = true;
    return
end
valid = isnumeric(value) && isequal(size(value), [7 2]);
if ~valid
    return
end
for row = 1:7
    pair = value(row, :);
    if all(isnan(pair))
        continue
    end
    if any(~isfinite(pair)) || pair(1) < 0 || pair(2) > 1 || ...
            pair(1) > pair(2)
        valid = false;
        return
    end
end
end
