function [EEG, info, EEG_preRemoval] = nf_cleanic(EEG, method, aggressive, varargin)
% NF_CLEANIC  Fit ICA and classify artifactual components.
%
% [EEG, INFO, EEG_PREREMOVAL] = NF_CLEANIC(EEG, METHOD, AGGRESSIVE, ...)
%
% METHOD selects 'iclabel', 'adjustedadjust'/'made', 'adjust', 'mara',
% 'faster', or 'none'. Except for 'none', every method uses the same
% training preparation and ICA; only the component classifier differs.
% The MADE path verifies required retained theta/radius scalp zones
% before adjusted_ADJUST so sparse or damaged montages cannot fail silently.
% The 'none' method returns without preparing or fitting ICA.
%
% ICA training uses a separately filtered copy. By default, that copy is
% divided into fixed one-second chunks following the MADE intent. When
% trainingEvents and trainingEpochLimits are supplied, the filtered copy is
% instead epoched around those events and the same rejection rules operate
% directly on the task epochs. The main working EEG remains continuous until
% the fitted components have been classified and subtracted. Event time is
% locked to zero without per-epoch baseline or mean removal. Overlapping
% training epochs are rejected to prevent duplicate sample weighting.
% ICA rank is checked from covariance eigenvalues using an absolute 1e-7
% threshold. Pivoted-QR channel repair runs only when this check finds a
% deficient training matrix.
%
% Name/value inputs:
%   algorithm             'runica' (default) or optional 'runamica15'
%   randomSeed            1 (controls runica; AMICA uses internal seeding)
%   trainingHighpass      1 Hz
%   trainingEpochLength   1 second
%   trainingVoltage       1000 microvolts
%   trainingPower         [-100 30] dB
%   trainingFrequencies   [20 40] Hz
%   trainingEvents        [] uses fixed chunks; otherwise event type(s)
%                         used to select the ICA fitting samples
%   trainingEpochLimits   [] with fixed chunks; [start end] seconds with
%                         trainingEvents
%   badChannelFraction    0.20
%   minimumTrainingEpochs 10
%   minimumSamplesPerRankSquared 20
%   iclabelThresholds     [] uses dominant-class NeuroFreq behavior, or 7x2
%   adjustReportFile      writable output path; temporary when empty
%   adjustedReportFile    compatibility alias for adjustReportFile
%   adjustOptions         scalar struct with fields:
%                           epochLength - 5 seconds for ADJUST; the configured
%                               trainingEpochLength for adjusted_ADJUST
%                           reportFile - temporary when empty
%   maraOptions           scalar struct with field:
%                           artifactProbabilityThreshold - [] uses MARA's
%                               released classifier decision; otherwise [0,1]
%                         Official MARA additionally requires at least 100 Hz
%                         input and enough data for its first 15-second feature
%                         window (approximately 16 seconds).
%   fasterOptions         scalar struct with fields:
%                           eogChannels - [] by default; numeric indices in
%                               the nf_cleanic input montage, remapped by
%                               unique label after ICA channel removals
%                           spectralBand - [] disables spectral-slope property
%                           measure - [1 1 1 1 1]
%                           z - [3 3 3 3 3]
%   amicaMaxIterations    2000
%   amicaThreads          4
%   amicaProcesses        1

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
addParameter(parser, 'trainingEvents', [], @nf_is_events);
addParameter(parser, 'trainingEpochLimits', [], @nf_is_limits_or_empty);
addParameter(parser, 'badChannelFraction', 0.20, @nf_is_fraction);
addParameter(parser, 'minimumTrainingEpochs', 10, @nf_is_positive_integer);
addParameter(parser, 'minimumSamplesPerRankSquared', 20, ...
    @nf_is_positive_scalar);
addParameter(parser, 'iclabelThresholds', [], @nf_is_iclabel_thresholds);
addParameter(parser, 'adjustReportFile', '', @nf_is_text_or_empty);
addParameter(parser, 'adjustedReportFile', '', @nf_is_text_or_empty);
addParameter(parser, 'adjustOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'maraOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'fasterOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'amicaMaxIterations', 2000, @nf_is_positive_integer);
addParameter(parser, 'amicaThreads', 4, @nf_is_positive_integer);
addParameter(parser, 'amicaProcesses', 1, @nf_is_positive_integer);
addParameter(parser, 'runicaStop', 1e-7, @nf_is_positive_scalar);
parse(parser, varargin{:});
options = parser.Results;
options.trainingPower = reshape(options.trainingPower, 1, 2);
options.trainingFrequencies = reshape(options.trainingFrequencies, 1, 2);
if ~isempty(options.trainingEpochLimits)
    options.trainingEpochLimits = ...
        reshape(options.trainingEpochLimits, 1, 2);
end

method = nf_normalize_method(method);
if strcmp(method, 'none')
    nf_validate_eeg_for_skip(EEG);
else
    nf_validate_continuous_eeg(EEG);
end
algorithm = lower(char(options.algorithm));
options = nf_resolve_classifier_options(options, method);
nf_validate_options(EEG, method, algorithm, aggressive, options);
if strcmp(method, 'none')
    [EEG, info, EEG_preRemoval] = nf_skip_ica(EEG, method);
    return
end
if ~isempty(options.trainingEvents)
    EEG = nf_prepare_training_events(EEG);
end
classifierProvenance = nf_preflight(method, algorithm, options);
if strcmp(method, 'iclabel')
    classifierProvenance.parameters.aggressive = logical(aggressive);
end
removalProvenance = nf_removal_provenance();
EEG = nf_normalize_locations(EEG);
fasterEogMapping = struct();
if strcmp(method, 'faster')
    fasterEogMapping = nf_capture_faster_eog_request( ...
        EEG, options.fasterOptions.eogChannels);
end

if strcmp(algorithm, 'runica')
    previousRandomState = rng;
    randomCleanup = onCleanup(@() rng(previousRandomState));
    rng(options.randomSeed, 'twister');
end

EEG = nf_clear_ica(EEG);
training = EEG;
continuousTrainingSampleCount = training.pnts;
training = pop_eegfiltnew(training, 'locutoff', options.trainingHighpass);
training = eeg_checkset(training);
training.datfile = '';
training.filepath = '';
training = nf_ensure_event_fields(training);
trainingCandidateEventPositions = [];
trainingCandidateEventIndices = [];
trainingCandidateEventLatencies = [];
trainingRejectedEventPositions = [];
trainingRejectedEventIndices = [];
trainingRejectedEventLatencies = [];
trainingRetainedEventPositions = [];
trainingRetainedEventIndices = [];
trainingRetainedEventLatencies = [];
trainingOverlap = nf_empty_overlap_info();
if isempty(options.trainingEvents)
    training = eeg_regepochs(training, ...
        'recurrence', options.trainingEpochLength, ...
        'limits', [0 options.trainingEpochLength], ...
        'rmbase', NaN, 'eventtype', 'nf_ica_training');
    trainingEpochSource = 'fixed-length';
    trainingEventTypes = {};
    trainingEpochLimits = [0 options.trainingEpochLength];
    trainingEpochLengthSeconds = options.trainingEpochLength;
else
    sourceEventLatencies = double([training.event.latency]);
    candidateEventIndices = nf_matching_event_indices( ...
        training, options.trainingEvents);
    [training, trainingCandidateEventPositions] = pop_epoch(training, ...
        nf_pop_epoch_events(options.trainingEvents), ...
        options.trainingEpochLimits, ...
        'eventindices', candidateEventIndices, ...
        'epochinfo', 'yes');
    trainingCandidateEventPositions = ...
        reshape(trainingCandidateEventPositions, 1, []);
    trainingCandidateEventIndices = ...
        candidateEventIndices(trainingCandidateEventPositions);
    trainingCandidateEventLatencies = ...
        sourceEventLatencies(trainingCandidateEventIndices);
    trainingOverlap = nf_training_epoch_overlap( ...
        sourceEventLatencies, trainingCandidateEventIndices, ...
        options.trainingEpochLimits, training.srate);
    if trainingOverlap.nOverlappingEpochs > 0
        error('nf_cleanic:OverlappingTrainingEpochs', ...
            ['The requested ICA-training epochs overlap in source time. ' ...
            'Overlapping epochs duplicate samples and inflate ICA sample ' ...
            'sufficiency. Shorten epochLimits/trainingEpochLimits or use ' ...
            'epochStage=''afterica'' with fixed ICA-training epochs.']);
    end
    trainingEpochSource = 'requested-event-locked';
    trainingEventTypes = nf_pop_epoch_events(options.trainingEvents);
    trainingEpochLimits = options.trainingEpochLimits;
    trainingEpochLengthSeconds = training.pnts / training.srate;
end
training = eeg_checkset(training);
initialTrainingEpochCount = training.trials;
initialTrainingSampleCount = training.pnts * training.trials;

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
if strcmp(trainingEpochSource, 'requested-event-locked')
    trainingRejectedEventPositions = ...
        trainingCandidateEventPositions(trainingRejectedMask);
    trainingRejectedEventIndices = ...
        trainingCandidateEventIndices(trainingRejectedMask);
    trainingRejectedEventLatencies = ...
        trainingCandidateEventLatencies(trainingRejectedMask);
    trainingRetainedEventPositions = ...
        trainingCandidateEventPositions(~trainingRejectedMask);
    trainingRetainedEventIndices = ...
        trainingCandidateEventIndices(~trainingRejectedMask);
    trainingRetainedEventLatencies = ...
        trainingCandidateEventLatencies(~trainingRejectedMask);
end
training = pop_rejepoch(training, trainingRejectedMask, 0);
training = eeg_checkset(training);
if training.trials < options.minimumTrainingEpochs
    error('nf_cleanic:InsufficientTrainingData', ...
        ['Only %d clean ICA training epochs remain from the %s preparation; ' ...
        'the configured minimum is %d.'], ...
        training.trials, trainingEpochSource, ...
        options.minimumTrainingEpochs);
end

[rankBeforeRepair, eigenvaluesBeforeRepair, rankTolerance] = ...
    nf_getrank(training);
if rankBeforeRepair < 2
    error('nf_cleanic:RankDeficient', ...
        'ICA training data have rank below two.');
end

fprintf( ...
    ['[NeuroFreq] ICA rank check: %d/%d channels using covariance ' ...
    'eigenvalues > %.1e.\n'], ...
    rankBeforeRepair, ...
    training.nbchan, ...
    rankTolerance);

independentIndices = 1:training.nbchan;
dependentIndices = [];
dependentLabels = {};
rankRepairApplied = false;
if rankBeforeRepair < training.nbchan
    [independentIndices, dependentIndices] = ...
        nf_rank_repair_plan(training, rankBeforeRepair);
    dependentLabels = {training.chanlocs(dependentIndices).labels};
    rankRepairApplied = true;
    fprintf( ...
        ['[NeuroFreq] Rank deficiency detected; removing %d dependent ' ...
        'ICA channel(s) before decomposition.\n'], ...
        numel(dependentIndices));
    training = pop_select(training, 'channel', independentIndices);
    EEG = pop_select(EEG, 'channel', independentIndices);
    training = eeg_checkset(training);
    EEG = eeg_checkset(EEG);

    [rankAfterRepair, eigenvaluesAfterRepair] = nf_getrank(training);
    if rankAfterRepair ~= training.nbchan
        error('nf_cleanic:RankRepairFailed', ...
            'ICA training data remain rank deficient after channel repair.');
    end
else
    rankAfterRepair = rankBeforeRepair;
    eigenvaluesAfterRepair = eigenvaluesBeforeRepair;
    fprintf( ...
        '[NeuroFreq] ICA data are full rank; channel rank repair skipped.\n');
end
if strcmp(method, 'faster')
    [options.fasterOptions.eogChannels, fasterEogMapping] = ...
        nf_remap_faster_eog_channels(EEG, fasterEogMapping);
    classifierProvenance.parameters.requestedEogChannelIndices = ...
        fasterEogMapping.requestedInputIndices;
    classifierProvenance.parameters.requestedEogChannelLabels = ...
        fasterEogMapping.requestedLabels;
    classifierProvenance.parameters.eogChannels = ...
        fasterEogMapping.retainedIndices;
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
elseif strcmp(method, 'adjustedadjust')
    [rejected, classification] = nf_classify_adjusted_adjust( ...
        EEG, options.adjustOptions.epochLength, ...
        options.adjustOptions.reportFile);
elseif strcmp(method, 'adjust')
    [rejected, classification] = nf_classify_adjust( ...
        EEG, options.adjustOptions);
elseif strcmp(method, 'mara')
    [rejected, classification] = nf_classify_mara( ...
        EEG, options.maraOptions);
elseif strcmp(method, 'faster')
    [rejected, classification] = nf_classify_faster( ...
        EEG, options.fasterOptions);
    classification.eogChannelMapping = fasterEogMapping;
else
    error('nf_cleanic:UnknownMethod', ...
        'No classifier implementation is registered for %s.', method);
end
classification.provenance = classifierProvenance;

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
classification = nf_normalize_classification( ...
    classification, method, rejected, componentCount);

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
info.schemaVersion = '2.1.0';
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
info.training.source = trainingEpochSource;
info.training.highpassHz = options.trainingHighpass;
info.training.epochLengthSeconds = trainingEpochLengthSeconds;
info.training.samplesPerEpoch = training.pnts;
info.training.dataSecondsPerEpoch = training.pnts / training.srate;
info.training.requestedFixedEpochLengthSeconds = ...
    options.trainingEpochLength;
info.training.events = trainingEventTypes;
info.training.requestedLimitsSeconds = trainingEpochLimits;
info.training.candidateEventIndices = ...
    reshape(trainingCandidateEventIndices, 1, []);
info.training.candidateEventPositions = ...
    reshape(trainingCandidateEventPositions, 1, []);
info.training.candidateEventLatenciesSamples = ...
    reshape(trainingCandidateEventLatencies, 1, []);
info.training.candidateEventLatenciesSeconds = ...
    (info.training.candidateEventLatenciesSamples - 1) / EEG.srate;
info.training.rejectedEventIndices = ...
    reshape(trainingRejectedEventIndices, 1, []);
info.training.rejectedEventPositions = ...
    reshape(trainingRejectedEventPositions, 1, []);
info.training.rejectedEventLatenciesSamples = ...
    reshape(trainingRejectedEventLatencies, 1, []);
info.training.rejectedEventLatenciesSeconds = ...
    (info.training.rejectedEventLatenciesSamples - 1) / EEG.srate;
info.training.retainedEventIndices = ...
    reshape(trainingRetainedEventIndices, 1, []);
info.training.retainedEventPositions = ...
    reshape(trainingRetainedEventPositions, 1, []);
info.training.retainedEventLatenciesSamples = ...
    reshape(trainingRetainedEventLatencies, 1, []);
info.training.retainedEventLatenciesSeconds = ...
    (info.training.retainedEventLatenciesSamples - 1) / EEG.srate;
info.training.acceptedEventIndices = info.training.candidateEventIndices;
info.training.acceptedEventPositions = info.training.candidateEventPositions;
info.training.timeLockedToEvent = ...
    strcmp(trainingEpochSource, 'requested-event-locked');
if info.training.timeLockedToEvent
    info.training.acceptedEventSemantics = ...
        ['Complete candidate epochs accepted by pop_epoch before ' ...
        'ICA-training artifact rejection; retainedEventIndices identifies ' ...
        'fitted epochs.'];
else
    info.training.acceptedEventSemantics = ...
        ['Not applicable to fixed-length eeg_regepochs preparation; event ' ...
        'identity ledgers are empty.'];
end
if info.training.timeLockedToEvent
    info.training.timeLockLatencySeconds = 0;
else
    info.training.timeLockLatencySeconds = NaN;
end
info.training.perEpochMeanRemoved = false;
info.training.baselineRemoval = 'none';
info.training.overlap = trainingOverlap;
info.training.voltageThresholdMicrovolts = options.trainingVoltage;
info.training.powerThresholdDb = options.trainingPower;
info.training.frequencyRangeHz = options.trainingFrequencies;
info.training.badChannelFractionLimit = options.badChannelFraction;
info.training.nInitialEpochs = initialTrainingEpochCount;
info.training.nInitialSamples = initialTrainingSampleCount;
info.training.continuousSourceSamples = continuousTrainingSampleCount;
info.training.selectedSourceSampleFraction = ...
    initialTrainingSampleCount / continuousTrainingSampleCount;
info.training.rejectedEpochMask = trainingRejectedMask;
info.training.nRejectedEpochs = sum(trainingRejectedMask);
info.training.nRetainedEpochs = training.trials;
info.training.nRetainedSamples = trainingSamples;
info.training.retainedSourceSampleFraction = ...
    trainingSamples / continuousTrainingSampleCount;
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
info.rank.method = ...
    'covariance eigenvalue count using an absolute 1e-7 threshold';
info.rank.eigenvaluesBeforeRepair = eigenvaluesBeforeRepair;
info.rank.eigenvaluesAfterRepair = eigenvaluesAfterRepair;
info.rank.repairApplied = rankRepairApplied;
info.rank.repairMethod = ...
    'pivoted QR channel selection, invoked only after deficient rank';
info.rank.singularValues = [];
info.rank.tolerance = rankTolerance;
info.rank.dependentChannelIndices = dependentIndices;
info.rank.dependentChannelLabels = dependentLabels;
info.rank.finalChannels = training.nbchan;
info.components.nComponents = componentCount;
info.components.rejected = rejected;
info.components.nRejected = numel(rejected);
info.components.nRetained = componentCount - numel(rejected);
if strcmp(method, 'iclabel')
    info.components.aggressive = logical(aggressive);
else
    info.components.aggressive = [];
end
info.classification = classification;
info.classification.dataScope = ...
    ['Full continuous post-preclean data carrying the decomposition fitted ' ...
    'to the reported ICA-training scope.'];
info.componentSubtractionDataScope = ...
    ['Full continuous post-preclean data; requested final epochs are ' ...
    'materialized after component subtraction.'];
info.provenance.classifier = classifierProvenance;
removalProvenance.executed = ~isempty(rejected);
info.provenance.removal = removalProvenance;
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
historyEntry.trainingSource = info.training.source;
historyEntry.trainingEpochsRetained = info.training.nRetainedEpochs;
historyEntry.rank = info.rank.afterRepair;
historyEntry.components = info.components;
historyEntry.classifierProvenance = info.provenance.classifier;
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

function [rankValue, eigenvalues, tolerance] = nf_getrank(EEG)
tmpdata = double(EEG.data(:, :));
if any(~isfinite(tmpdata(:)))
    error('nf_cleanic:NonfiniteTrainingData', ...
        'ICA training data contain NaN or Inf values.');
end
covarianceMatrix = cov(tmpdata', 1);
[~, diagonalMatrix] = eig(covarianceMatrix);
tolerance = 1e-7;
eigenvalues = real(diag(diagonalMatrix));
rankValue = sum(eigenvalues > tolerance);
end

function [independentIndices, dependentIndices] = ...
    nf_rank_repair_plan(EEG, rankValue)
data = reshape(EEG.data, EEG.nbchan, []);
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

if ~isnumeric(rankValue) || ~isscalar(rankValue) || ...
        ~isfinite(rankValue) || rankValue ~= round(rankValue) || ...
        rankValue < 2 || rankValue >= EEG.nbchan
    error('nf_cleanic:InvalidRepairRank', ...
        'Rank repair requires an integer rank from 2 through nbchan - 1.');
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
    training = pop_runamica(training, ...
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

function [rejected, classification] = nf_classify_adjust(EEG, options)
if size(EEG.icaweights, 1) ~= size(EEG.icaweights, 2) || ...
        size(EEG.icaweights, 1) ~= EEG.nbchan
    error('nf_cleanic:AdjustRequiresSquareICA', ...
        'ADJUST requires a square full-channel ICA decomposition.');
end

classifierEEG = nf_prepare_adjust_eeg( ...
    EEG, options.epochLength, 'nf_adjust_training', 'nf_adjust');
[reportFile, temporaryReport, reportCleanup] = ...
    nf_prepare_adjust_report(options.reportFile);
[~, environmentCleanup] = ...
    nf_prepare_vendor_environment();

[rejected, horizontal, vertical, blink, discontinuity, ...
    spatialVarianceDifferenceThreshold, spatialVarianceDifference, ...
    temporalKurtosisThreshold, temporalKurtosisMedian, temporalKurtosis, ...
    spatialEyeDifferenceThreshold, spatialEyeDifferenceMedian, ...
    spatialEyeDifference, spatialAverageDifferenceThreshold, ...
    spatialAverageDifferenceMedian, spatialAverageDifference, ...
    discontinuitySpatialThreshold, discontinuitySpatialMedian, ...
    discontinuitySpatial, maximumEpochVarianceThreshold, ...
    maximumEpochVarianceMedian, maximumEpochVariance, ...
    maximumDerivativeThreshold, maximumDerivative] = ...
    ADJUST(classifierEEG, reportFile);

clear environmentCleanup

classification = struct();
classification.name = 'ADJUST';
classification.classNames = {'horizontalEyeMovement', ...
    'verticalEyeMovement', 'blink', 'discontinuity'};
classification.horizontalEyeMovement = horizontal(:)';
classification.verticalEyeMovement = vertical(:)';
classification.blink = blink(:)';
classification.discontinuity = discontinuity(:)';
classification.rule = 'official ADJUST artifact-component union';
classification.thresholds.spatialVarianceDifference = ...
    spatialVarianceDifferenceThreshold;
classification.thresholds.temporalKurtosis = temporalKurtosisThreshold;
classification.thresholds.spatialEyeDifference = ...
    spatialEyeDifferenceThreshold;
classification.thresholds.spatialAverageDifference = ...
    spatialAverageDifferenceThreshold;
classification.thresholds.discontinuitySpatial = ...
    discontinuitySpatialThreshold;
classification.thresholds.maximumEpochVariance = ...
    maximumEpochVarianceThreshold;
classification.thresholds.maximumDerivative = maximumDerivativeThreshold;
classification.features.spatialVarianceDifference = ...
    spatialVarianceDifference;
classification.features.temporalKurtosisMedian = temporalKurtosisMedian;
classification.features.temporalKurtosis = temporalKurtosis;
classification.features.spatialEyeDifferenceMedian = ...
    spatialEyeDifferenceMedian;
classification.features.spatialEyeDifference = spatialEyeDifference;
classification.features.spatialAverageDifferenceMedian = ...
    spatialAverageDifferenceMedian;
classification.features.spatialAverageDifference = spatialAverageDifference;
classification.features.discontinuitySpatialMedian = ...
    discontinuitySpatialMedian;
classification.features.discontinuitySpatial = discontinuitySpatial;
classification.features.maximumEpochVarianceMedian = ...
    maximumEpochVarianceMedian;
classification.features.maximumEpochVariance = maximumEpochVariance;
classification.features.maximumDerivative = maximumDerivative;
classification.epochLengthSeconds = options.epochLength;
classification.reportWasTemporary = temporaryReport;
if temporaryReport
    classification.reportFile = '';
else
    classification.reportFile = reportFile;
end

clear reportCleanup
end

function [rejected, classification] = nf_classify_mara(EEG, options)
[officialRejected, maraInformation] = MARA(EEG);
if ~isstruct(maraInformation) || ...
        ~isfield(maraInformation, 'posterior_artefactprob')
    error('nf_cleanic:InvalidMARAOutput', ...
        ['MARA did not return posterior_artefactprob. Its released MARA.m ' ...
        'suppresses some internal errors; inspect the preceding output.']);
end

probabilities = double(maraInformation.posterior_artefactprob(:)');
componentCount = size(EEG.icaweights, 1);
if numel(probabilities) ~= componentCount || ...
        any(~isfinite(probabilities)) || ...
        any(probabilities < 0) || any(probabilities > 1)
    error('nf_cleanic:InvalidMARAOutput', ...
        'MARA returned invalid artifact probabilities.');
end

if isempty(options.artifactProbabilityThreshold)
    rejected = officialRejected;
    rule = 'official MARA binary classifier decision';
else
    rejected = find( ...
        probabilities >= options.artifactProbabilityThreshold);
    rule = 'NeuroFreq threshold applied to MARA artifact probabilities';
end

classification = struct();
classification.name = 'MARA';
classification.classNames = {'artifact', 'accept'};
classification.artifactProbability = probabilities;
classification.probabilities = [probabilities(:) 1 - probabilities(:)];
classification.normalizedFeatures = [];
if isfield(maraInformation, 'normfeats')
    classification.normalizedFeatures = maraInformation.normfeats;
end
classification.officialRejected = officialRejected(:)';
classification.thresholds.artifactProbability = ...
    options.artifactProbabilityThreshold;
classification.rule = rule;
end

function mapping = nf_capture_faster_eog_request(EEG, requestedIndices)
mapping = struct();
mapping.requestedInputIndices = reshape(requestedIndices, 1, []);
mapping.requestedLabels = {};
mapping.retainedIndices = zeros(1, 0);
mapping.retainedLabels = {};
mapping.indicesChanged = false;

if isempty(mapping.requestedInputIndices)
    return
end
mapping.requestedLabels = ...
    {EEG.chanlocs(mapping.requestedInputIndices).labels};
if numel(unique(lower(string(mapping.requestedLabels)))) ~= ...
        numel(mapping.requestedLabels)
    error('nf_cleanic:DuplicateFASTERRequestedEOGLabels', ...
        ['fasterOptions.eogChannels resolves to duplicate channel labels ' ...
        'in the nf_cleanic input montage.']);
end
end

function [retainedIndices, mapping] = ...
    nf_remap_faster_eog_channels(EEG, mapping)
retainedIndices = zeros(1, 0);
if isempty(mapping.requestedInputIndices)
    return
end

retainedLabels = {EEG.chanlocs.labels};
[present, retainedIndices] = ismember( ...
    lower(string(mapping.requestedLabels)), lower(string(retainedLabels)));
if any(~present)
    removedLabels = mapping.requestedLabels(~present);
    error('nf_cleanic:FASTERRequestedEOGRemoved', ...
        ['FASTER EOG channel(s) %s were removed during ICA preparation or ' ...
        'rank repair. Supply surviving EOG channels or revise the channel/' ...
        'rank-cleaning settings; numeric input indices cannot be reused ' ...
        'after their channels are removed.'], strjoin(removedLabels, ', '));
end

retainedIndices = reshape(double(retainedIndices), 1, []);
if numel(unique(retainedIndices)) ~= numel(retainedIndices)
    error('nf_cleanic:DuplicateFASTERRemappedEOGChannels', ...
        'FASTER EOG labels did not remap to unique retained channels.');
end
mapping.retainedIndices = retainedIndices;
mapping.retainedLabels = retainedLabels(retainedIndices);
mapping.indicesChanged = ...
    ~isequal(mapping.requestedInputIndices, retainedIndices);
end

function [rejected, classification] = nf_classify_faster(EEG, options)
effectiveMeasure = logical(options.measure);
if isempty(options.spectralBand)
    spectralBand = [];
    effectiveMeasure(2) = false;
else
    spectralBand = options.spectralBand;
end
if isempty(options.eogChannels)
    effectiveMeasure(5) = false;
end
if ~any(effectiveMeasure)
    error('nf_cleanic:NoFASTERComponentMeasures', ...
        'At least one usable FASTER component property must be enabled.');
end

properties = component_properties( ...
    EEG, options.eogChannels, spectralBand);
if ~isequal(size(properties), [size(EEG.icaweights, 1) 5]) || ...
        any(~isfinite(properties(:)))
    error('nf_cleanic:InvalidFASTEROutput', ...
        'FASTER component_properties returned invalid output.');
end

rejectionOptions = struct();
rejectionOptions.measure = effectiveMeasure;
rejectionOptions.z = options.z;
rejectedMask = logical(min_z(properties, rejectionOptions));
if numel(rejectedMask) ~= size(EEG.icaweights, 1)
    error('nf_cleanic:InvalidFASTEROutput', ...
        'FASTER min_z returned an invalid component mask.');
end
rejected = find(rejectedMask);

classification = struct();
classification.name = 'FASTER';
classification.propertyNames = {'medianGradient', 'spectralSlope', ...
    'spatialKurtosis', 'hurstExponent', 'eogCorrelation'};
classification.properties = properties;
classification.requestedMeasure = logical(options.measure);
classification.effectiveMeasure = effectiveMeasure;
classification.thresholds.z = options.z;
classification.eogChannels = options.eogChannels;
classification.spectralBandHz = spectralBand;
classification.rule = ...
    'official FASTER component_properties followed by official min_z';
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
[rejected, horizontal, vertical, blink, discontinuity, ...
    spatialVarianceDifferenceThreshold, spatialVarianceDifference, ...
    temporalKurtosisThreshold, temporalKurtosisMedian, temporalKurtosis, ...
    spatialEyeDifferenceThreshold, spatialEyeDifferenceMedian, ...
    spatialEyeDifference, spatialAverageDifferenceThreshold, ...
    spatialAverageDifferenceMedian, spatialAverageDifference, ...
    discontinuitySpatialThreshold, discontinuitySpatialMedian, ...
    discontinuitySpatial, maximumEpochVarianceThreshold, ...
    maximumEpochVarianceMedian, maximumEpochVariance, ...
    maximumDerivativeThreshold, maximumDerivative] = ...
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
classification.classNames = {'horizontalEyeMovement', ...
    'verticalEyeMovement', 'blink', 'discontinuity'};
classification.horizontalEyeMovement = horizontal(:)';
classification.verticalEyeMovement = vertical(:)';
classification.blink = blink(:)';
classification.discontinuity = discontinuity(:)';
classification.montageZoneCoverage = zoneCoverage;
classification.rule = ...
    'official MADE adjusted_ADJUST artifact-component union';
classification.thresholds.spatialVarianceDifference = ...
    spatialVarianceDifferenceThreshold;
classification.thresholds.temporalKurtosis = temporalKurtosisThreshold;
classification.thresholds.spatialEyeDifference = ...
    spatialEyeDifferenceThreshold;
classification.thresholds.spatialAverageDifference = ...
    spatialAverageDifferenceThreshold;
classification.thresholds.discontinuitySpatial = ...
    discontinuitySpatialThreshold;
classification.thresholds.maximumEpochVariance = ...
    maximumEpochVarianceThreshold;
classification.thresholds.maximumDerivative = maximumDerivativeThreshold;
classification.features.spatialVarianceDifference = ...
    spatialVarianceDifference;
classification.features.temporalKurtosisMedian = temporalKurtosisMedian;
classification.features.temporalKurtosis = temporalKurtosis;
classification.features.spatialEyeDifferenceMedian = ...
    spatialEyeDifferenceMedian;
classification.features.spatialEyeDifference = spatialEyeDifference;
classification.features.spatialAverageDifferenceMedian = ...
    spatialAverageDifferenceMedian;
classification.features.spatialAverageDifference = spatialAverageDifference;
classification.features.discontinuitySpatialMedian = ...
    discontinuitySpatialMedian;
classification.features.discontinuitySpatial = discontinuitySpatial;
classification.features.maximumEpochVarianceMedian = ...
    maximumEpochVarianceMedian;
classification.features.maximumEpochVariance = maximumEpochVariance;
classification.features.maximumDerivative = maximumDerivative;
classification.epochLengthSeconds = epochLength;
classification.reportWasTemporary = temporaryReport;
if temporaryReport
    classification.reportFile = '';
else
    classification.reportFile = reportFile;
end

clear reportCleanup
end

function classifierEEG = nf_prepare_adjust_eeg( ...
    EEG, epochLength, eventType, setName)
classifierEEG = nf_ensure_event_fields(EEG);
classifierEEG = eeg_regepochs(classifierEEG, ...
    'recurrence', epochLength, ...
    'limits', [0 epochLength], ...
    'rmbase', NaN, ...
    'eventtype', eventType);
classifierEEG = eeg_checkset(classifierEEG, 'ica');
classifierEEG.setname = setName;
classifierEEG.filename = [setName '.set'];
classifierEEG.filepath = '';
classifierEEG.datfile = '';
end

function [reportFile, temporaryReport, reportCleanup] = ...
    nf_prepare_adjust_report(requestedReportFile)
temporaryReport = isempty(requestedReportFile);
if temporaryReport
    reportFile = [tempname '_adjust_report.txt'];
else
    reportFile = nf_absolute_path(char(requestedReportFile));
    nf_validate_report_path(reportFile);
end
reportCleanup = onCleanup(@() nf_remove_temporary_report( ...
    reportFile, temporaryReport));
end

function [workingDirectory, environmentCleanup] = ...
    nf_prepare_vendor_environment()
workingDirectory = tempname;
[created, message] = mkdir(workingDirectory);
if ~created
    error('nf_cleanic:TemporaryDirectoryFailed', ...
        'Could not create the classifier working directory: %s', message);
end
previousDirectory = pwd;
figuresBefore = findall(0, 'Type', 'figure');
randomState = rng;
warningState = warning;
pathState = path;
environmentCleanup = onCleanup(@() nf_restore_vendor_environment( ...
    previousDirectory, workingDirectory, figuresBefore, ...
    randomState, warningState, pathState));
cd(workingDirectory);
end

function nf_restore_vendor_environment( ...
    previousDirectory, workingDirectory, figuresBefore, ...
    randomState, warningState, pathState)
nf_close_new_figures(figuresBefore);
nf_restore_environment(randomState, warningState, pathState);
nf_restore_and_remove_directory(previousDirectory, workingDirectory);
end

function classification = nf_normalize_classification( ...
    classification, method, rejected, componentCount)
classification.schemaVersion = '1.0.0';
classification.method = method;
classification.componentCount = componentCount;
classification.rejected = rejected(:)';
classification.rejectedMask = false(1, componentCount);
classification.rejectedMask(rejected) = true;
classification.nRejected = numel(rejected);
classification.nRetained = componentCount - numel(rejected);
classification.removalOwner = 'NeuroFreq';
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
if isfield(EEG, 'icaweights') && ~isempty(EEG.icaweights)
    summary.nComponents = size(EEG.icaweights, 1);
else
    summary.nComponents = 0;
end
end

function [EEG, info, EEG_preRemoval] = nf_skip_ica(EEG, method)
EEG_preRemoval = EEG;
info = struct();
info.schemaVersion = '2.1.0';
info.method = method;
info.algorithm = 'none';
info.algorithmDetails = struct();
info.randomSeed = [];
info.randomness.requestedSeed = [];
info.randomness.algorithmSeedControlled = false;
info.randomness.control = 'not applicable';
info.preRemovalDataset = nf_dataset_summary(EEG);
info.training = struct();
info.training.source = 'none';
info.training.highpassHz = NaN;
info.training.epochLengthSeconds = NaN;
info.training.samplesPerEpoch = 0;
info.training.dataSecondsPerEpoch = 0;
info.training.events = {};
info.training.requestedLimitsSeconds = [];
info.training.candidateEventIndices = [];
info.training.candidateEventPositions = [];
info.training.candidateEventLatenciesSamples = [];
info.training.candidateEventLatenciesSeconds = [];
info.training.rejectedEventIndices = [];
info.training.rejectedEventPositions = [];
info.training.rejectedEventLatenciesSamples = [];
info.training.rejectedEventLatenciesSeconds = [];
info.training.retainedEventIndices = [];
info.training.retainedEventPositions = [];
info.training.retainedEventLatenciesSamples = [];
info.training.retainedEventLatenciesSeconds = [];
info.training.acceptedEventIndices = [];
info.training.acceptedEventPositions = [];
info.training.acceptedEventSemantics = 'ICA training was skipped.';
info.training.timeLockedToEvent = false;
info.training.timeLockLatencySeconds = NaN;
info.training.perEpochMeanRemoved = false;
info.training.baselineRemoval = 'none';
info.training.overlap = nf_empty_overlap_info();
info.training.nInitialEpochs = 0;
info.training.nRejectedEpochs = 0;
info.training.nRetainedEpochs = 0;
info.training.nInitialSamples = 0;
info.training.nRetainedSamples = 0;
info.rank = struct();
componentCount = info.preRemovalDataset.nComponents;
info.components.nComponents = componentCount;
info.components.rejected = zeros(1, 0);
info.components.nRejected = 0;
info.components.nRetained = componentCount;
info.components.aggressive = [];
info.classification = struct();
info.classification.schemaVersion = '1.0.0';
info.classification.method = method;
info.classification.name = 'None';
info.classification.componentCount = componentCount;
info.classification.rejected = zeros(1, 0);
info.classification.rejectedMask = false(1, componentCount);
info.classification.nRejected = 0;
info.classification.nRetained = componentCount;
info.classification.rule = 'ICA fitting and component removal skipped';
info.classification.removalOwner = 'none';
info.classification.dataScope = 'not applicable';
info.classification.provenance.contractLevel = 'native-control';
info.classification.provenance.provider = 'NeuroFreq';
info.classification.provenance.release = '';
info.classification.provenance.functions = struct([]);
info.provenance.classifier = info.classification.provenance;
info.provenance.removal.contractLevel = 'not-run';
info.provenance.removal.provider = '';
info.componentSubtractionDataScope = 'not applicable';
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
historyEntry.trainingSource = info.training.source;
historyEntry.trainingEpochsRetained = [];
historyEntry.rank = [];
historyEntry.components = info.components;
historyEntry.classifierProvenance = info.provenance.classifier;
EEG.etc.nf_cleanic_history{end + 1} = historyEntry;
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

function options = nf_resolve_classifier_options(options, method)
if ismember(method, {'adjust', 'adjustedadjust'})
    options.adjustOptions = nf_resolve_adjust_options(options, method);
elseif strcmp(method, 'mara')
    maraDefaults = struct();
    maraDefaults.artifactProbabilityThreshold = [];
    options.maraOptions = nf_merge_option_struct( ...
        maraDefaults, options.maraOptions, 'maraOptions');
    nf_validate_mara_options(options.maraOptions);
elseif strcmp(method, 'faster')
    fasterDefaults = struct();
    fasterDefaults.eogChannels = [];
    fasterDefaults.spectralBand = [];
    fasterDefaults.measure = true(1, 5);
    fasterDefaults.z = 3 * ones(1, 5);
    options.fasterOptions = nf_merge_option_struct( ...
        fasterDefaults, options.fasterOptions, 'fasterOptions');
    nf_validate_faster_classifier_options(options.fasterOptions);
    options.fasterOptions.eogChannels = ...
        reshape(options.fasterOptions.eogChannels, 1, []);
    options.fasterOptions.spectralBand = ...
        reshape(options.fasterOptions.spectralBand, 1, []);
    options.fasterOptions.measure = ...
        reshape(logical(options.fasterOptions.measure), 1, 5);
    options.fasterOptions.z = reshape(options.fasterOptions.z, 1, 5);
end
end

function adjustOptions = nf_resolve_adjust_options(options, method)
if ~isempty(options.adjustReportFile) && ...
        ~isempty(options.adjustedReportFile) && ...
        ~nf_same_report_file(options.adjustReportFile, ...
        options.adjustedReportFile)
    error('nf_cleanic:ConflictingAdjustReportFiles', ...
        ['adjustReportFile and adjustedReportFile specify different ' ...
        'paths. Supply only one compatibility argument.']);
end
legacyReportFile = options.adjustReportFile;
if isempty(legacyReportFile)
    legacyReportFile = options.adjustedReportFile;
end

suppliedOptions = options.adjustOptions;
adjustDefaults = struct();
if strcmp(method, 'adjust')
    adjustDefaults.epochLength = 5;
else
    adjustDefaults.epochLength = options.trainingEpochLength;
end
adjustDefaults.reportFile = legacyReportFile;
adjustOptions = nf_merge_option_struct( ...
    adjustDefaults, suppliedOptions, 'adjustOptions');
if isempty(adjustOptions.reportFile)
    adjustOptions.reportFile = legacyReportFile;
end
nf_validate_adjust_options(adjustOptions);

if ~isempty(legacyReportFile) && ...
        isfield(suppliedOptions, 'reportFile') && ...
        ~isempty(suppliedOptions.reportFile) && ...
        ~nf_same_report_file(legacyReportFile, suppliedOptions.reportFile)
    error('nf_cleanic:ConflictingAdjustReportFiles', ...
        ['The compatibility report-file argument and ' ...
        'adjustOptions.reportFile specify different paths.']);
end
end

function equal = nf_same_report_file(firstValue, secondValue)
firstPath = nf_absolute_path(char(firstValue));
secondPath = nf_absolute_path(char(secondValue));
equal = nf_paths_equal(firstPath, secondPath);
end

function resolved = nf_merge_option_struct(defaults, supplied, optionName)
resolved = defaults;
names = fieldnames(supplied);
allowed = fieldnames(defaults);
for index = 1:numel(names)
    name = names{index};
    if ~ismember(name, allowed)
        error('nf_cleanic:UnknownClassifierOption', ...
            'Unknown %s field: %s.', optionName, name);
    end
    resolved.(name) = supplied.(name);
end
end

function nf_validate_adjust_options(options)
if ~nf_is_positive_scalar(options.epochLength)
    error('nf_cleanic:InvalidAdjustOptions', ...
        'adjustOptions.epochLength must be a positive finite scalar.');
end
if ~nf_is_text_or_empty(options.reportFile)
    error('nf_cleanic:InvalidAdjustOptions', ...
        'adjustOptions.reportFile must be empty or scalar text.');
end
end

function nf_validate_mara_options(options)
threshold = options.artifactProbabilityThreshold;
if ~isempty(threshold) && ...
        (~isnumeric(threshold) || ~isscalar(threshold) || ...
        ~isfinite(threshold) || threshold < 0 || threshold > 1)
    error('nf_cleanic:InvalidMARAOptions', ...
        ['maraOptions.artifactProbabilityThreshold must be empty or a ' ...
        'finite scalar in [0, 1].']);
end
end

function nf_validate_faster_classifier_options(options)
if ~isempty(options.eogChannels)
    channels = options.eogChannels;
    if ~isnumeric(channels) || ~isvector(channels) || ...
            any(~isfinite(channels)) || any(channels < 1) || ...
            any(channels ~= round(channels)) || ...
            numel(unique(channels)) ~= numel(channels)
        error('nf_cleanic:InvalidFASTEROptions', ...
            ['fasterOptions.eogChannels must contain unique positive ' ...
            'integer channel indices.']);
    end
end
if ~isempty(options.spectralBand) && ...
        ~nf_is_increasing_pair(options.spectralBand)
    error('nf_cleanic:InvalidFASTEROptions', ...
        ['fasterOptions.spectralBand must be empty or a two-value ' ...
        'increasing frequency range.']);
end
if ~isnumeric(options.measure) && ~islogical(options.measure)
    error('nf_cleanic:InvalidFASTEROptions', ...
        'fasterOptions.measure must contain five binary values.');
end
if numel(options.measure) ~= 5 || ...
        any(~ismember(options.measure, [0 1]))
    error('nf_cleanic:InvalidFASTEROptions', ...
        'fasterOptions.measure must contain five binary values.');
end
if ~isnumeric(options.z) || numel(options.z) ~= 5 || ...
        any(~isfinite(options.z)) || any(options.z <= 0)
    error('nf_cleanic:InvalidFASTEROptions', ...
        'fasterOptions.z must contain five positive finite values.');
end
end

function method = nf_normalize_method(method)
method = lower(char(method));
if ismember(method, {'made', 'adjustedadjust', 'adjustedadust', ...
        'adjusted_adjust', 'adjusted_adust', 'adjusted-adjust', ...
        'adjusted-adust'})
    method = 'adjustedadjust';
elseif ismember(method, {'adjust', 'mara', 'faster', 'none'})
    return
elseif ~strcmp(method, 'iclabel')
    error('nf_cleanic:UnknownMethod', ...
        ['method must be ''iclabel'', ''adjustedadjust'', ''adjust'', ' ...
        '''mara'', ''faster'', or ''none''.']);
end
end

function provenance = nf_preflight(method, algorithm, options)
required = {'pop_eegfiltnew', 'eeg_regepochs', 'pop_eegthresh', ...
    'pop_rejspec', 'pop_rejepoch', 'pop_subcomp', ...
    'convertlocs'};
if ~isempty(options.trainingEvents)
    required{end + 1} = 'pop_epoch';
end
if strcmp(algorithm, 'runica')
    required{end + 1} = 'pop_runica';
end
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
            'The ICLabel EEGLAB plugin is required.');
    end
    if ~isempty(options.iclabelThresholds) && exist('pop_icflag', 'file') ~= 2
        error('nf_cleanic:MissingICLabel', ...
            'pop_icflag.m is required when iclabelThresholds are supplied.');
    end
elseif strcmp(method, 'adjustedadjust')
    if exist('adjusted_ADJUST', 'file') ~= 2
        error('nf_cleanic:MissingAdjustedAdjust', ...
            ['adjusted_ADJUST.m from the official MADE distribution is ' ...
            'required.']);
    end
elseif strcmp(method, 'adjust')
    if exist('ADJUST', 'file') ~= 2
        error('nf_cleanic:MissingAdjust', ...
            'ADJUST.m from the official ADJUST plugin is required.');
    end
elseif strcmp(method, 'mara')
    if exist('MARA', 'file') ~= 2
        error('nf_cleanic:MissingMARA', ...
            'MARA.m from the official MARA plugin is required.');
    end
    maraRequired = {'classify', 'pwelch'};
    for index = 1:numel(maraRequired)
        if exist(maraRequired{index}, 'file') ~= 2
            error('nf_cleanic:MissingMARADependency', ...
                ['MARA requires %s.m and its associated MATLAB ' ...
                'toolbox.'], maraRequired{index});
        end
    end
    hasOptimizationFitter = exist('lsqcurvefit', 'file') == 2;
    hasCurveFittingFitter = exist('fit', 'file') == 2 && ...
        exist('fitoptions', 'file') == 2 && ...
        exist('fittype', 'file') == 2;
    hasStatisticsFitter = exist('NonLinearModel', 'class') == 8;
    if ~hasOptimizationFitter && ~hasCurveFittingFitter && ...
            ~hasStatisticsFitter
        error('nf_cleanic:MissingMARAFitter', ...
            ['Official MARA requires one released nonlinear-fit path: ' ...
            'lsqcurvefit, the fit/fitoptions/fittype set, or NonLinearModel.']);
    end
elseif strcmp(method, 'faster')
    fasterVendorRequired = {'component_properties', 'min_z', ...
        'hurst_exponent'};
    for index = 1:numel(fasterVendorRequired)
        if exist(fasterVendorRequired{index}, 'file') ~= 2
            error('nf_cleanic:MissingFASTER', ...
                ['%s.m from the official FASTER distribution was not ' ...
                'found.'], fasterVendorRequired{index});
        end
    end
    fasterDependencies = {'eeg_getica', 'pwelch', 'kurt', 'nanmean'};
    for index = 1:numel(fasterDependencies)
        if exist(fasterDependencies{index}, 'file') ~= 2
            error('nf_cleanic:MissingFASTERDependency', ...
                ['Official FASTER component_properties requires %s.m, ' ...
                'which was not found on the MATLAB path.'], ...
                fasterDependencies{index});
        end
    end
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
    madeCompanions = {'MARA_extract_time_freq_features', 'beall_horizontal', ...
        'beall_blink_detection', 'Spatial_Info_eyes'};
    nf_assert_vendor_directory('adjusted_ADJUST', madeCompanions);
    if exist('ADJUST', 'file') ~= 2
        error('nf_cleanic:MissingAdjustedAdjustCompanion', ...
            ['The official ADJUST plugin is required because MADE ships ' ...
            'adjusted_ADJUST changes but reuses ADJUST feature functions.']);
    end
    adjustCompanions = {'compute_GD_feat', 'computeSED_NOnorm', ...
        'computeSAD', 'trim_and_mean', 'trim_and_max', 'EM'};
    nf_assert_vendor_directory('ADJUST', adjustCompanions);
end
if strcmp(method, 'adjust')
    companions = {'compute_GD_feat', 'computeSED_NOnorm', 'computeSAD', ...
        'kurt', 'trim_and_mean', 'trim_and_max', 'EM'};
    for index = 1:numel(companions)
        if exist(companions{index}, 'file') ~= 2
            error('nf_cleanic:MissingAdjustCompanion', ...
                ['The official ADJUST companion %s.m was not found. Add ' ...
                'the complete ADJUST plugin folder to the path.'], ...
                companions{index});
        end
    end
    vendorCompanions = {'compute_GD_feat', 'computeSED_NOnorm', ...
        'computeSAD', 'trim_and_mean', 'trim_and_max', 'EM'};
    nf_assert_vendor_directory('ADJUST', vendorCompanions);
end
if strcmp(method, 'mara')
    maraData = {'fv_training_MARA.mat', 'inv_matrix_icbm152.mat'};
    maraPath = nf_resolved_path('MARA');
    maraDirectory = fileparts(maraPath);
    for index = 1:numel(maraData)
        dataPath = fullfile(maraDirectory, maraData{index});
        if exist(dataPath, 'file') ~= 2
            error('nf_cleanic:MissingMARAData', ...
                ['The MARA data file %s was not found beside MARA.m. Add ' ...
                'the complete official MARA plugin folder to the path.'], ...
                maraData{index});
        end
    end
end
if strcmp(method, 'faster')
    nf_assert_vendor_directory('component_properties', ...
        {'min_z', 'hurst_exponent'});
end
provenance = nf_classifier_provenance(method, options);
end

function nf_assert_vendor_directory(entryName, companionNames)
entryPath = which(entryName);
entryDirectory = fileparts(entryPath);
for index = 1:numel(companionNames)
    companionPath = which(companionNames{index});
    companionDirectory = fileparts(companionPath);
    if ~nf_paths_equal(entryDirectory, companionDirectory)
        error('nf_cleanic:MixedVendorCode', ...
            ['%s resolves to %s, but companion %s resolves outside that ' ...
            'vendor folder at %s. Reorder the MATLAB path so one complete ' ...
            'official distribution supplies this classifier.'], ...
            entryName, entryPath, companionNames{index}, companionPath);
    end
end
end

function equal = nf_paths_equal(firstPath, secondPath)
if ispc
    equal = strcmpi(firstPath, secondPath);
else
    equal = strcmp(firstPath, secondPath);
end
end

function provenance = nf_classifier_provenance(method, options)
provenance = struct();
provenance.contractLevel = 'vendor-exact-classifier';
provenance.decisionProvider = 'NeuroFreq';
provenance.removalProvider = 'NeuroFreq through EEGLAB pop_subcomp';

if strcmp(method, 'iclabel')
    provider = 'ICLabel';
    names = {'pop_iclabel'};
    if ~isempty(options.iclabelThresholds)
        names{end + 1} = 'pop_icflag';
        provenance.decisionProvider = 'ICLabel pop_icflag';
    end
    source = 'https://github.com/sccn/ICLabel';
    parameters = struct();
    parameters.thresholds = options.iclabelThresholds;
elseif strcmp(method, 'adjustedadjust')
    provider = 'MADE adjusted_ADJUST';
    names = {'adjusted_ADJUST', 'MARA_extract_time_freq_features', ...
        'beall_horizontal', 'beall_blink_detection', 'Spatial_Info_eyes', ...
        'compute_GD_feat', 'computeSED_NOnorm', 'computeSAD', ...
        'trim_and_mean', 'trim_and_max', 'EM'};
    source = 'https://github.com/ChildDevLab/MADE-EEG-preprocessing-pipeline';
    parameters = options.adjustOptions;
    provenance.decisionProvider = 'MADE adjusted_ADJUST';
elseif strcmp(method, 'adjust')
    provider = 'ADJUST';
    names = {'ADJUST', 'compute_GD_feat', 'computeSED_NOnorm', ...
        'computeSAD', 'trim_and_mean', 'trim_and_max', 'EM'};
    source = 'https://www.nitrc.org/projects/adjust/';
    parameters = options.adjustOptions;
    provenance.decisionProvider = 'ADJUST';
elseif strcmp(method, 'mara')
    provider = 'MARA';
    names = {'MARA'};
    source = 'https://github.com/irenne/MARA';
    parameters = options.maraOptions;
    if isempty(options.maraOptions.artifactProbabilityThreshold)
        provenance.decisionProvider = 'MARA';
    end
elseif strcmp(method, 'faster')
    provider = 'FASTER';
    names = {'component_properties', 'min_z', 'hurst_exponent'};
    source = 'https://sourceforge.net/projects/faster/';
    parameters = options.fasterOptions;
    provenance.decisionProvider = 'FASTER min_z';
else
    error('nf_cleanic:UnknownMethod', ...
        'Cannot construct provenance for method %s.', method);
end

records = repmat(struct( ...
    'name', '', ...
    'path', '', ...
    'sha256', ''), 1, numel(names));
for index = 1:numel(names)
    records(index).name = names{index};
    records(index).path = nf_resolved_path(names{index});
    records(index).sha256 = nf_file_sha256(records(index).path);
    if isempty(records(index).sha256)
        error('nf_cleanic:VendorHashFailed', ...
            ['Could not calculate a SHA-256 identity for %s at %s. ' ...
            'The strict classifier contract was not satisfied.'], ...
            names{index}, records(index).path);
    end
end

provenance.provider = provider;
provenance.source = source;
provenance.functions = records;
provenance.release = nf_detect_vendor_release(provider, records(1).path);
provenance.parameters = parameters;
provenance.checkedAt = datestr(now, 30); %#ok<TNOW1,DATST>
end

function provenance = nf_removal_provenance()
provenance = struct();
provenance.contractLevel = 'vendor-primitive';
provenance.provider = 'EEGLAB';
provenance.entryPoint = 'pop_subcomp';
provenance.decisionProvider = 'NeuroFreq';
provenance.path = nf_resolved_path('pop_subcomp');
provenance.sha256 = nf_file_sha256(provenance.path);
if isempty(provenance.sha256)
    error('nf_cleanic:VendorHashFailed', ...
        ['Could not calculate a SHA-256 identity for pop_subcomp at %s. ' ...
        'The strict removal contract was not satisfied.'], ...
        provenance.path);
end
provenance.executed = false;
end

function pathValue = nf_resolved_path(functionName)
resolved = which(functionName, '-all');
if isempty(resolved)
    error('nf_cleanic:MissingDependency', ...
        'Required function %s was not found.', functionName);
end
if ischar(resolved)
    if size(resolved, 1) > 1
        resolved = cellstr(resolved);
    else
        resolved = {resolved};
    end
elseif isstring(resolved)
    resolved = cellstr(resolved(:));
elseif ~iscell(resolved)
    resolved = {char(resolved)};
end
resolved = resolved(~cellfun(@isempty, resolved));
normalized = string(resolved);
if ispc
    normalized = lower(normalized);
end
[~, uniqueIndices] = unique(normalized, 'stable');
resolved = resolved(sort(uniqueIndices));
if numel(resolved) ~= 1
    error('nf_cleanic:AmbiguousDependency', ...
        ['%s resolves to multiple files. The vendor contract requires ' ...
        'one unambiguous implementation. Resolved files: %s'], ...
        functionName, strjoin(resolved, ' | '));
end
pathValue = resolved{1};
if exist(pathValue, 'file') ~= 2
    error('nf_cleanic:InvalidDependencyPath', ...
        '%s did not resolve to a readable source file.', functionName);
end
end

function release = nf_detect_vendor_release(provider, entryPath)
release = '';
directory = fileparts(entryPath);
if strcmp(provider, 'ADJUST')
    pluginFile = fullfile(directory, 'eegplugin_adjust.m');
    release = nf_extract_release(pluginFile, ...
        'ADJUST([0-9]+(?:\.[0-9]+)+)');
elseif strcmp(provider, 'FASTER')
    guiFile = fullfile(directory, 'FASTER_GUI.m');
    release = nf_extract_release(guiFile, ...
        'FASTER_version\s*=\s*''([^'']+)''');
elseif strcmp(provider, 'ICLabel')
    pluginFiles = dir(fullfile(directory, 'eegplugin_iclabel.m'));
    if ~isempty(pluginFiles)
        release = nf_extract_release( ...
            fullfile(directory, pluginFiles(1).name), ...
            '[Vv]ersion[^0-9]*([0-9]+(?:\.[0-9]+)+)');
    end
end
if isempty(release)
    [~, directoryName] = fileparts(directory);
    token = regexp(directoryName, ...
        '([0-9]+(?:\.[0-9]+)+(?:[A-Za-z][0-9A-Za-z.-]*)?)', ...
        'tokens', 'once');
    if ~isempty(token)
        release = token{1};
    end
end
end

function release = nf_extract_release(filePath, expression)
release = '';
if exist(filePath, 'file') ~= 2
    return
end
try
    contents = fileread(filePath);
    token = regexp(contents, expression, 'tokens', 'once');
    if ~isempty(token)
        release = token{1};
    end
catch
    release = '';
end
end

function value = nf_file_sha256(filePath)
value = '';
if isempty(filePath) || exist(filePath, 'file') ~= 2
    return
end
fileId = fopen(filePath, 'r');
if fileId < 0
    return
end
fileCleanup = onCleanup(@() fclose(fileId));
try
    bytes = fread(fileId, Inf, '*uint8');
    digest = java.security.MessageDigest.getInstance('SHA-256');
    digest.update(bytes);
    digestBytes = typecast(digest.digest(), 'uint8');
    value = lower(reshape(dec2hex(digestBytes, 2).', 1, []));
catch
    value = '';
end
clear fileCleanup
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

function events = nf_pop_epoch_events(events)
if isstring(events)
    events = cellstr(events);
elseif ischar(events)
    events = {events};
elseif isnumeric(events)
    events = num2cell(events);
end
end

function EEG = nf_prepare_training_events(EEG)
if ~isfield(EEG, 'event') || isempty(EEG.event) || ...
        ~isstruct(EEG.event) || ~isfield(EEG.event, 'type') || ...
        ~isfield(EEG.event, 'latency')
    error('nf_cleanic:InvalidTrainingEvents', ...
        ['Event-locked ICA training requires EEG.event structures with ' ...
        'type and latency fields.']);
end

types = {EEG.event.type};
containsText = any(cellfun(@(value) ischar(value) || ...
    (isstring(value) && isscalar(value)), types));
latencies = nan(1, numel(EEG.event));
for eventIndex = 1:numel(EEG.event)
    value = types{eventIndex};
    if isstring(value) && isscalar(value)
        value = char(value);
    end
    if containsText && isnumeric(value) && isscalar(value) && ...
            isfinite(value)
        value = num2str(value);
    end
    if containsText
        validType = ischar(value) && isrow(value);
    else
        validType = isnumeric(value) && isreal(value) && ...
            isscalar(value) && isfinite(value);
    end
    if ~validType
        error('nf_cleanic:InvalidTrainingEvents', ...
            ['EEG.event.type must use consistently representable scalar ' ...
            'text or finite numeric values.']);
    end
    EEG.event(eventIndex).type = value;

    latency = EEG.event(eventIndex).latency;
    if ~isnumeric(latency) || ~isreal(latency) || ~isscalar(latency) || ...
            ~isfinite(latency) || latency < 0.5 || ...
            latency > EEG.pnts + 0.5
        error('nf_cleanic:InvalidTrainingEvents', ...
            'EEG.event(%d).latency is outside the continuous sample range.', ...
            eventIndex);
    end
    latencies(eventIndex) = double(latency);
end

[~, eventOrder] = sort(latencies);
EEG.event = EEG.event(eventOrder);
end

function indices = nf_matching_event_indices(EEG, requested)
requestedValues = nf_pop_epoch_events(requested);
if isempty(requestedValues) || ~isfield(EEG, 'event') || ...
        isempty(EEG.event) || ~isfield(EEG.event, 'type')
    error('nf_cleanic:MissingTrainingEvents', ...
        'No valid events are available for event-locked ICA training.');
end

availableValues = {EEG.event.type};
availableAreText = ischar(availableValues{1}) || ...
    (isstring(availableValues{1}) && isscalar(availableValues{1}));
matched = cell(1, numel(requestedValues));
for requestIndex = 1:numel(requestedValues)
    requestedValue = requestedValues{requestIndex};
    if availableAreText
        if ischar(requestedValue)
            normalized = deblank(requestedValue);
        elseif isstring(requestedValue) && isscalar(requestedValue)
            normalized = deblank(char(requestedValue));
        elseif isnumeric(requestedValue) && isscalar(requestedValue) && ...
                isfinite(requestedValue)
            normalized = num2str(requestedValue);
        else
            error('nf_cleanic:InvalidTrainingEventType', ...
                'Training event types must be scalar text or numbers.');
        end
        availableText = cell(size(availableValues));
        for eventIndex = 1:numel(availableValues)
            value = availableValues{eventIndex};
            if isstring(value) && isscalar(value)
                value = char(value);
            end
            if ~ischar(value)
                error('nf_cleanic:MixedTrainingEventTypes', ...
                    'EEG.event.type must use one consistent storage type.');
            end
            availableText{eventIndex} = deblank(value);
        end
        matched{requestIndex} = find(strcmp(availableText, normalized));
    else
        if ischar(requestedValue) || ...
                (isstring(requestedValue) && isscalar(requestedValue))
            normalized = str2double(strtrim(char(requestedValue)));
        else
            normalized = requestedValue;
        end
        if ~isnumeric(normalized) || ~isscalar(normalized) || ...
                ~isfinite(normalized)
            error('nf_cleanic:InvalidTrainingEventType', ...
                'A requested event type cannot be matched to numeric EEG events.');
        end
        numericAvailable = nan(1, numel(availableValues));
        for eventIndex = 1:numel(availableValues)
            value = availableValues{eventIndex};
            if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
                error('nf_cleanic:MixedTrainingEventTypes', ...
                    'EEG.event.type must use one consistent storage type.');
            end
            numericAvailable(eventIndex) = value;
        end
        matched{requestIndex} = find(numericAvailable == normalized);
    end
    if isempty(matched{requestIndex})
        error('nf_cleanic:TrainingEventTypeNotFound', ...
            'A requested ICA-training event type does not occur in EEG.event.');
    end
end
indices = sort(unique([matched{:}]));
end

function info = nf_empty_overlap_info()
info = struct();
info.checked = false;
info.hasOverlap = false;
info.nAcceptedEpochs = 0;
info.nOverlappingEpochs = 0;
info.overlappingAcceptedPositions = [];
info.minimumEventSeparationSeconds = NaN;
info.requestedEpochDurationSeconds = NaN;
end

function info = nf_training_epoch_overlap(sourceLatencies, acceptedIndices, ...
        limits, samplingRate)
info = nf_empty_overlap_info();
info.checked = true;
acceptedIndices = reshape(double(acceptedIndices), 1, []);
info.nAcceptedEpochs = numel(acceptedIndices);
info.requestedEpochDurationSeconds = diff(limits);
if isempty(acceptedIndices)
    return
end
if any(~isfinite(acceptedIndices)) || ...
        any(acceptedIndices ~= round(acceptedIndices)) || ...
        any(acceptedIndices < 1) || ...
        any(acceptedIndices > numel(sourceLatencies))
    error('nf_cleanic:InvalidAcceptedEventIndices', ...
        'pop_epoch returned event indices outside the source event table.');
end

centers = sourceLatencies(acceptedIndices);
if any(~isfinite(centers))
    error('nf_cleanic:InvalidTrainingEventLatency', ...
        'An accepted ICA-training event has a nonfinite source latency.');
end
[centers, order] = sort(centers);
starts = centers + limits(1) * samplingRate;
ends = centers + limits(2) * samplingRate;
if numel(centers) > 1
    info.minimumEventSeparationSeconds = ...
        min(diff(centers)) / samplingRate;
end

overlapping = false(1, numel(centers));
runningEnd = ends(1);
tolerance = 10 * eps(max(abs([starts ends])) + 1);
for index = 2:numel(centers)
    if starts(index) < runningEnd - tolerance
        overlapping(index) = true;
    end
    runningEnd = max(runningEnd, ends(index));
end

overlappingOriginalOrder = false(size(overlapping));
overlappingOriginalOrder(order) = overlapping;
info.hasOverlap = any(overlapping);
info.nOverlappingEpochs = sum(overlapping);
info.overlappingAcceptedPositions = find(overlappingOriginalOrder);
end

function nf_validate_options(EEG, method, algorithm, aggressive, options)
if ~nf_is_logical_scalar(aggressive)
    error('nf_cleanic:InvalidAggressiveFlag', ...
        'aggressive must be a logical scalar.');
end
if ~isempty(options.trainingEvents) && ...
        isempty(options.trainingEpochLimits)
    error('nf_cleanic:MissingTrainingEpochLimits', ...
        'trainingEpochLimits are required when trainingEvents are supplied.');
end
if isempty(options.trainingEvents) && ...
        ~isempty(options.trainingEpochLimits)
    error('nf_cleanic:UnexpectedTrainingEpochLimits', ...
        'trainingEpochLimits apply only when trainingEvents are supplied.');
end
if ~isempty(options.trainingEvents) && ...
        (~isfield(EEG, 'event') || isempty(EEG.event))
    error('nf_cleanic:MissingTrainingEvents', ...
        'The input dataset has no events for event-locked ICA training.');
end
if strcmp(method, 'none')
    if ~isempty(options.trainingEvents)
        error('nf_cleanic:TrainingEventsWithoutICA', ...
            'trainingEvents cannot be used when method=''none''.');
    end
    return
end
if ~strcmp(method, 'iclabel') && logical(aggressive)
    error('nf_cleanic:ClassifierOptionMismatch', ...
        'aggressive applies only to the ICLabel method.');
end
if ~strcmp(method, 'iclabel') && ~isempty(options.iclabelThresholds)
    error('nf_cleanic:ClassifierOptionMismatch', ...
        'iclabelThresholds apply only to the ICLabel method.');
end
if ~ismember(algorithm, {'runamica15', 'runica'})
    error('nf_cleanic:UnknownAlgorithm', ...
        'algorithm must be ''runamica15'' or ''runica''.');
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
if strcmp(method, 'adjustedadjust') && EEG.srate < 100
    error('nf_cleanic:MADESamplingRate', ...
        ['The MADE adjusted_ADJUST feature extractor requires a sampling ' ...
        'rate of at least 100 Hz.']);
end
if strcmp(method, 'adjustedadjust') && ...
        options.trainingEpochLength * EEG.srate < 100
    error('nf_cleanic:MADEEpochLength', ...
        ['trainingEpochLength must yield at least 100 samples per MADE ' ...
        'adjusted_ADJUST feature epoch.']);
end
if ismember(method, {'adjust', 'adjustedadjust'}) && ...
        options.adjustOptions.epochLength * EEG.srate < 100
    error('nf_cleanic:ClassifierEpochLength', ...
        ['adjustOptions.epochLength must yield at least 100 samples for ' ...
        'ADJUST feature extraction.']);
end
if strcmp(method, 'mara')
    nf_validate_mara_input(EEG);
end
if strcmp(method, 'faster')
    if any(options.fasterOptions.eogChannels > EEG.nbchan)
        error('nf_cleanic:InvalidFASTEROptions', ...
            'fasterOptions.eogChannels contains an index above EEG.nbchan.');
    end
    if ~isempty(options.fasterOptions.spectralBand) && ...
            (options.fasterOptions.spectralBand(1) < 0 || ...
            options.fasterOptions.spectralBand(2) > EEG.srate / 2)
        error('nf_cleanic:InvalidFASTEROptions', ...
            ['fasterOptions.spectralBand must be nonnegative and no ' ...
            'higher than Nyquist.']);
    end
end
end

function nf_validate_mara_input(EEG)
if EEG.srate < 100
    error('nf_cleanic:MARASamplingRate', ...
        ['Official MARA expects input sampled at least at 100 Hz before its ' ...
        'released 100-200-Hz internal downsampling and spectral features.']);
end

downsampleFactor = max(floor(EEG.srate / 100), 1);
maraSamplingRate = round(EEG.srate / downsampleFactor);
if maraSamplingRate / 2 < 39
    error('nf_cleanic:MARANyquist', ...
        ['MARA''s internally downsampled data must retain its released ' ...
        '33-39-Hz spectral feature band.']);
end

maraSampleCount = floor((EEG.pnts - 1) / downsampleFactor) + 1;
minimumMaraSamples = 16 * maraSamplingRate;
if maraSampleCount < minimumMaraSamples
    requiredSeconds = ...
        (minimumMaraSamples * downsampleFactor) / EEG.srate;
    error('nf_cleanic:MARADuration', ...
        ['Official MARA''s 15-second local-skewness feature starts at one ' ...
        'second and requires at least %d internally downsampled samples ' ...
        '(approximately %.3g input seconds); only %d are available.'], ...
        minimumMaraSamples, requiredSeconds, maraSampleCount);
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

function nf_validate_eeg_for_skip(EEG)
required = {'data', 'srate', 'nbchan', 'pnts', 'trials', ...
    'chanlocs', 'xmin', 'xmax'};
if ~isstruct(EEG) || numel(EEG) ~= 1
    error('nf_cleanic:InvalidEEG', ...
        'EEG must be one valid EEGLAB dataset structure.');
end
for index = 1:numel(required)
    if ~isfield(EEG, required{index})
        error('nf_cleanic:InvalidEEG', ...
            'EEG is missing required field %s.', required{index});
    end
end
if ~isnumeric(EEG.data) || isempty(EEG.data) || ...
        ~isnumeric(EEG.srate) || ~isscalar(EEG.srate) || ...
        ~isfinite(EEG.srate) || EEG.srate <= 0 || ...
        ~isnumeric(EEG.nbchan) || ~isscalar(EEG.nbchan) || ...
        EEG.nbchan ~= size(EEG.data, 1) || ...
        numel(EEG.chanlocs) ~= EEG.nbchan || ...
        ~isfield(EEG.chanlocs, 'labels')
    error('nf_cleanic:InvalidEEG', ...
        'EEG data, sampling rate, and channel metadata are inconsistent.');
end
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

function valid = nf_is_events(value)
valid = isempty(value) || ischar(value) || isstring(value) || ...
    isnumeric(value) || iscell(value);
end

function valid = nf_is_limits_or_empty(value)
valid = isempty(value) || nf_is_increasing_pair(value);
end

function valid = nf_is_scalar_struct(value)
valid = isstruct(value) && isscalar(value);
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
