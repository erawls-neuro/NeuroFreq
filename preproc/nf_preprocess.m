function [EEG, report, quality, figureHandle] = nf_preprocess(EEG, varargin)
% NF_PREPROCESS  Complete NeuroFreq adult/child EEGLAB preprocessing.
%
% [EEG, REPORT, QUALITY, FIGUREHANDLE] = NF_PREPROCESS(EEG, ...)
%
% The two presets share filtering, FASTER, GEDAI, ICA training, epoch
% handling, 125-uV thresholding, FFT muscle screening, interpolation, and
% final reference. Their sole preset difference is the IC classifier:
%
%   adult : ICLabel
%   child : MADE adjusted_ADJUST
%
% Common name/value inputs:
%   preset                 'adult' (default) or 'child'
%   events                 event type(s) for pop_epoch; [] makes fixed epochs
%   epochLimits            [start end] seconds when events are supplied
%   continuousEpochLength  1 second when events are empty
%   baseline               false (default) or [start end] milliseconds
%   lowpass                45 Hz; retains the complete 20-40-Hz muscle band
%   highpass               0.3 Hz
%   notch                  60 Hz; skipped when redundant below low-pass
%   resample               250 Hz
%   eegChannels            [] auto-selects channels typed EEG/empty
%   channelMethod          'faster' (default) (may add other options later)
%   maxBadChannels         floor(10%% of selected channels), used as QC limit
%   badChannelReference    [] automatic, channel label, or index
%   cleanHighpass          1-Hz diagnostic copy for cleanrawdata channels
%   precleanMethod         'gedai' (default) or 'asr'
%   gedaiOptions           GEDAI v1.7 options; default uses the standard
%                          precomputed 10-5 database. Use referenceMatrixType
%                          'interpolated' for nonstandard coordinate montages.
%   icaMethod              '' uses the preset classifier
%   icaAlgorithm           'runica' (default) or optional 'runamica15'
%   randomSeed             1
%   minimumSamplesPerRankSquared 20 for production ICA sufficiency
%   iclabelThresholds      [] uses dominant ICLabel class, or a 7x2 matrix
%   voltageThreshold       125 microvolts
%   powerThreshold         [-100 30] dB
%   muscleRange            [20 40] Hz
%   thresholdTimes         [] uses the full final epoch
%   localInterp            true
%   maxLocalBad            floor(10%% of selected channels)
%   frontalChannels        {} uses validated frontopolar labels/coordinates
%   globalInterpolation    true
%   rereference            true for final common-average reference
%   qualityCompute         false
%   qualityPlot            false
%   qualityVisible         'on' or 'off'
%   save                   false (default), true, an output directory, or
%                          an exact .set path. Existing outputs are never
%                          overwritten.
%   log                    false (default); true records Command Window
%                          output and full failure diagnostics for this job.
%
% Save behavior:
%
%   save=false             Return the processed EEG without writing a .set.
%   save=true              Write <source>_preprocessed.set beside the source.
%   save=<directory>       Write <source>_preprocessed.set in that directory.
%   save=<file.set>        Write the exact requested .set path.
%
% When save and qualityPlot are both true, the QC figure is written as .fig
% and .pdf beside the .set. When save and log are both true, the job log is
% written there as well. With log=true and save=false, the log is written
% beside the source dataset. QC export failures are recorded without
% preventing the dataset and log from being preserved. Output bundles are
% reserved with a per-output bundle lock, and every artifact is checked
% immediately before writing. Existing artifacts are never overwritten. A
% pre-existing MATLAB diary is restored. A relative active DiaryFile is
% resolved against the starting working directory before redirection. If
% MATLAB or a plugin disables or redirects the job diary, nf_preprocess
% restarts it at the next checkpoint, records the possible capture gap, and
% continues.
%
% Future priority:
%   BIDS-aware input and derivative export will be integrated after the
%   dedicated vetted-raw nf_bidsify platform is complete. nf_preprocess
%   currently performs no BIDS inference, entity naming, metadata writing,
%   or directory construction.
%
% For compatibility, an event specification may be supplied as the first
% positional argument, followed by name/value inputs.

nf_validate_input_eeg(EEG);
EEG = nf_ensure_event_fields(EEG);
EEG = nf_normalize_event_types(EEG);
varargin = nf_normalize_legacy_events(varargin);

parser = inputParser;
parser.FunctionName = 'nf_preprocess';
parser.KeepUnmatched = false;
addParameter(parser, 'preset', 'adult', @nf_is_text);
addParameter(parser, 'events', [], @nf_is_events);
addParameter(parser, 'epochLimits', [], @nf_is_limits_or_empty);
addParameter(parser, 'continuousEpochLength', 1, @nf_is_positive_scalar);
addParameter(parser, 'baseline', false, @nf_is_baseline);
addParameter(parser, 'lowpass', 45, @nf_is_positive_scalar);
addParameter(parser, 'highpass', 0.3, @nf_is_nonnegative_scalar);
addParameter(parser, 'notch', 60, @nf_is_nonnegative_scalar);
addParameter(parser, 'resample', 250, @nf_is_positive_scalar);
addParameter(parser, 'eegChannels', [], @nf_is_channel_specification);
addParameter(parser, 'channelMethod', 'faster', @nf_is_text);
addParameter(parser, 'maxBadChannels', floor(EEG.nbchan / 5), ...
    @nf_is_nonnegative_integer);
addParameter(parser, 'badChannelReference', [], @nf_is_reference);
addParameter(parser, 'fasterOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'cleanCorrelation', 0.8, @nf_is_correlation);
addParameter(parser, 'cleanHighpass', 1, @nf_is_positive_scalar);
addParameter(parser, 'precleanMethod', 'gedai', @nf_is_text);
addParameter(parser, 'gedaiOptions', struct(), @nf_is_gedai_options);
addParameter(parser, 'icaMethod', '', @nf_is_text);
addParameter(parser, 'icaAlgorithm', 'runica', @nf_is_text);
addParameter(parser, 'aggressiveICA', false, @nf_is_logical_scalar);
addParameter(parser, 'randomSeed', 1, @nf_is_nonnegative_integer);
addParameter(parser, 'iclabelThresholds', [], @nf_is_iclabel_thresholds);
addParameter(parser, 'adjustReportFile', '', @nf_is_text_or_empty);
addParameter(parser, 'icaTrainingHighpass', 1, @nf_is_positive_scalar);
addParameter(parser, 'icaTrainingEpochLength', 1, @nf_is_positive_scalar);
addParameter(parser, 'icaTrainingVoltage', 1000, @nf_is_positive_scalar);
addParameter(parser, 'icaTrainingPower', [-100 30], @nf_is_increasing_pair);
addParameter(parser, 'icaTrainingFrequencies', [20 40], @nf_is_increasing_pair);
addParameter(parser, 'icaBadChannelFraction', 0.20, @nf_is_fraction);
addParameter(parser, 'minimumTrainingEpochs', 10, @nf_is_positive_integer);
addParameter(parser, 'minimumSamplesPerRankSquared', 20, ...
    @nf_is_positive_scalar);
addParameter(parser, 'amicaMaxIterations', 2000, @nf_is_positive_integer);
addParameter(parser, 'amicaThreads', 4, @nf_is_positive_integer);
addParameter(parser, 'amicaProcesses', 1, @nf_is_positive_integer);
addParameter(parser, 'runicaStop', 1e-7, @nf_is_positive_scalar);
addParameter(parser, 'voltageThreshold', 125, @nf_is_positive_scalar);
addParameter(parser, 'powerThreshold', [-100 30], @nf_is_increasing_pair);
addParameter(parser, 'muscleRange', [20 40], @nf_is_increasing_pair);
addParameter(parser, 'thresholdTimes', [], @nf_is_limits_or_empty);
addParameter(parser, 'localInterp', true, @nf_is_logical_scalar);
addParameter(parser, 'maxLocalBad', floor(EEG.nbchan / 10), ...
    @nf_is_nonnegative_integer);
addParameter(parser, 'frontalChannels', {}, @nf_is_labels);
addParameter(parser, 'globalInterpolation', true, @nf_is_logical_scalar);
addParameter(parser, 'interpolationMethod', 'sphericalKang', @nf_is_text);
addParameter(parser, 'rereference', true, @nf_is_logical_scalar);
addParameter(parser, 'qualityCompute', false, @nf_is_logical_scalar);
addParameter(parser, 'qualityPlot', false, @nf_is_logical_scalar);
addParameter(parser, 'qualityVisible', 'on', @nf_is_visibility);
addParameter(parser, 'qualityOptions', {}, @iscell);
addParameter(parser, 'save', false, @nf_is_save_request);
addParameter(parser, 'log', false, @nf_is_logical_scalar);
parse(parser, varargin{:});
options = parser.Results;
usingDefaultMaxBad = ismember('maxBadChannels', parser.UsingDefaults);
usingDefaultMaxLocal = ismember('maxLocalBad', parser.UsingDefaults);

options.preset = lower(char(options.preset));
options.channelMethod = lower(char(options.channelMethod));
options.precleanMethod = lower(char(options.precleanMethod));
options.icaAlgorithm = lower(char(options.icaAlgorithm));
options.icaMethod = nf_resolve_ica_method(options.preset, options.icaMethod);
options.interpolationMethod = nf_normalize_interpolation_method( ...
    options.interpolationMethod);
options.qualityCompute = logical(options.qualityCompute || options.qualityPlot);
options.log = logical(options.log);
options.epochLimits = nf_row_pair(options.epochLimits);
if isnumeric(options.baseline)
    options.baseline = nf_row_pair(options.baseline);
end
options.icaTrainingPower = nf_row_pair(options.icaTrainingPower);
options.icaTrainingFrequencies = nf_row_pair(options.icaTrainingFrequencies);
options.powerThreshold = nf_row_pair(options.powerThreshold);
options.muscleRange = nf_row_pair(options.muscleRange);
options.thresholdTimes = nf_row_pair(options.thresholdTimes);
options.gedaiOptions = nf_resolve_gedai_configuration(options.gedaiOptions);

artifactPlan = nf_resolve_artifact_plan(EEG, options);
% Retain both onCleanup objects until the job reaches its finalization path.
bundleCleanup = [];
[artifactPlan, bundleCleanup] = ...
    nf_prepare_artifact_destination(artifactPlan);
diaryCleanup = [];
diaryState = struct();
if artifactPlan.log
    diaryState = nf_start_job_log(artifactPlan, options);
    diaryCleanup = onCleanup(@() nf_restore_diary(diaryState));
end
datasetCommitted = false;
qualityFigureOwned = false;
qualityPdfOwned = false;

try
    nf_validate_options(EEG, options);
    nf_check_job_log(artifactPlan, 'option validation');
    if isempty(options.events)
        options.eventValidation = struct();
    else
        options.eventValidation = nf_validate_requested_events(EEG, options.events);
    end
    nf_preflight_dependencies(options);

    report = struct();
    report.schemaVersion = '3.0.0';
    report.started = datestr(now, 30); %#ok<TNOW1,DATST>
    report.preset = options.preset;
    if strcmp(options.channelMethod, 'faster') && ...
            strcmp(options.precleanMethod, 'gedai')
        report.pipelineMode = 'standard-preset';
    else
        report.pipelineMode = 'customized-common-stages';
    end
    report.presetDefinition.commonPipeline = ...
        'Filtering + FASTER + GEDAI + common ICA preparation + voltage/FFT cleaning';
    report.presetDefinition.standardVoltageThresholdMicrovolts = 125;
    report.presetDefinition.classifier = options.icaMethod;
    report.presetDefinition.onlyPresetDifference = ...
        'adult uses ICLabel; child uses MADE adjusted_ADJUST';
    report.options = nf_report_options(options);
    report.input = nf_dataset_summary(EEG);
    report.steps = struct();
    report.channels = struct();
    report.gedai = struct();
    report.ica = struct();
    report.epochs = struct();
    report.quality = struct('computed', false, 'plotted', false, 'error', '');
    report.persistence = nf_initial_persistence_report(artifactPlan);

    [EEG, channelSelection] = nf_select_eeg_channels(EEG, options.eegChannels);
    EEG = nf_normalize_channel_locations(EEG);
    report.channels.selection = channelSelection;
    nf_validate_interpolation_locations(EEG.chanlocs);
    if usingDefaultMaxBad
        options.maxBadChannels = floor(EEG.nbchan / 5);
    end
    if usingDefaultMaxLocal
        options.maxLocalBad = floor(EEG.nbchan / 5);
    end
    if strcmp(options.precleanMethod, 'gedai')
        options.gedaiOptions = nf_validate_gedai_montage( ...
            options.gedaiOptions, EEG.chanlocs);
    end
    report.options = nf_report_options(options);
    if options.maxBadChannels >= EEG.nbchan - 2
        error('nf_preprocess:InvalidBadChannelLimit', ...
            'maxBadChannels must leave at least three selected EEG channels.');
    end
    if options.localInterp && options.maxLocalBad > EEG.nbchan - 3
        error('nf_preprocess:InvalidLocalChannelLimit', ...
            'maxLocalBad must leave at least three interpolation donor channels.');
    end

    originalMontage = EEG.chanlocs;
    if ~isfield(EEG, 'etc') || isempty(EEG.etc)
        EEG.etc = struct();
    end
    EEG.etc.ogchan = originalMontage;

    [EEG, filterInfo] = nf_filter(EEG, options.lowpass, options.highpass, ...
        options.notch, options.resample);
    report.steps.filter = filterInfo;
    nf_check_job_log(artifactPlan, 'filtering');

    if strcmp(options.channelMethod, 'none')
        channelInfo = struct();
        channelInfo.method = 'none';
        channelInfo.removed.labels = {};
        channelInfo.removed.nRemoved = 0;
    else
        [EEG, channelInfo] = nf_badchans(EEG, options.maxBadChannels, ...
            false, options.channelMethod, ...
            'reference', options.badChannelReference, ...
            'fasterOptions', options.fasterOptions, ...
            'cleanCorrelation', options.cleanCorrelation, ...
            'cleanHighpass', options.cleanHighpass, ...
            'interpolationMethod', options.interpolationMethod);
    end
    report.channels.detection = channelInfo;
    report.steps.badChannels.applied = ~strcmp(options.channelMethod, 'none');
    report.steps.badChannels.method = options.channelMethod;
    nf_check_job_log(artifactPlan, 'bad-channel detection');

    if strcmp(options.precleanMethod, 'gedai')
        gedaiConfiguration = nf_subset_gedai_reference( ...
            options.gedaiOptions, originalMontage, EEG.chanlocs);
        if ischar(gedaiConfiguration.referenceMatrixType) && ...
                strcmp(gedaiConfiguration.referenceMatrixType, 'precomputed')
            nf_validate_gedai_standard_labels(EEG.chanlocs);
        end
        [EEG, gedaiInfo] = nf_run_gedai(EEG, gedaiConfiguration);
    else
        gedaiInfo = struct();
        gedaiInfo.applied = false;
    end
    EEG = eeg_checkset(EEG);
    report.gedai = gedaiInfo;
    report.steps.gedai.applied = gedaiInfo.applied;
    nf_check_job_log(artifactPlan, 'GEDAI');
    if options.qualityCompute
        EEG_preclean = EEG;
    else
        EEG_preclean = [];
    end

    [EEG, icaInfo] = nf_cleanic(EEG, ...
        options.icaMethod, options.aggressiveICA, ...
        'algorithm', options.icaAlgorithm, ...
        'randomSeed', options.randomSeed, ...
        'trainingHighpass', options.icaTrainingHighpass, ...
        'trainingEpochLength', options.icaTrainingEpochLength, ...
        'trainingVoltage', options.icaTrainingVoltage, ...
        'trainingPower', options.icaTrainingPower, ...
        'trainingFrequencies', options.icaTrainingFrequencies, ...
        'badChannelFraction', options.icaBadChannelFraction, ...
        'minimumTrainingEpochs', options.minimumTrainingEpochs, ...
        'minimumSamplesPerRankSquared', ...
        options.minimumSamplesPerRankSquared, ...
        'iclabelThresholds', options.iclabelThresholds, ...
        'adjustReportFile', options.adjustReportFile, ...
        'amicaMaxIterations', options.amicaMaxIterations, ...
        'amicaThreads', options.amicaThreads, ...
        'amicaProcesses', options.amicaProcesses, ...
        'runicaStop', options.runicaStop);
    EEG = eeg_checkset(EEG);
    if options.qualityCompute
        EEG_postclean = EEG;
    else
        EEG_postclean = [];
    end
    report.ica = icaInfo;
    report.steps.ica.applied = true;
    report.steps.ica.method = options.icaMethod;
    report.steps.ica.algorithm = options.icaAlgorithm;
    nf_check_job_log(artifactPlan, 'ICA cleaning');

    [EEG, epochInfo] = nf_make_final_epochs(EEG, options);
    report.epochs.creation = epochInfo;
    report.steps.epoch = epochInfo;
    nf_check_job_log(artifactPlan, 'epoch creation');

    if isnumeric(options.baseline) && numel(options.baseline) == 2
        [baselineWindow, baselineAdjustment] = nf_fit_epoch_window( ...
            options.baseline, [EEG.xmin EEG.xmax] * 1000, ...
            1000 / EEG.srate, 'baseline');
        EEG = pop_rmbase(EEG, baselineWindow);
        EEG = eeg_checkset(EEG);
        report.epochs.baselineApplied = true;
        report.epochs.baselineRequestedMilliseconds = options.baseline;
        report.epochs.baselineMilliseconds = baselineWindow;
        report.epochs.baselineEndpointAdjusted = baselineAdjustment;
    else
        report.epochs.baselineApplied = false;
        report.epochs.baselineMilliseconds = [];
    end

    thresholdTimes = options.thresholdTimes;
    if isempty(thresholdTimes)
        thresholdTimes = [EEG.xmin EEG.xmax];
        thresholdAdjustment = false;
    else
        [thresholdTimes, thresholdAdjustment] = nf_fit_epoch_window( ...
            thresholdTimes, [EEG.xmin EEG.xmax], 1 / EEG.srate, ...
            'thresholdTimes');
    end
    [EEG, ~, thresholdInfo] = nf_thresh(EEG, ...
        options.voltageThreshold, options.powerThreshold, ...
        options.muscleRange, thresholdTimes, options.localInterp, ...
        options.maxLocalBad, options.frontalChannels, ...
        'interpolationMethod', options.interpolationMethod);
    report.epochs.threshold = thresholdInfo;
    report.epochs.threshold.requestedTimesSeconds = options.thresholdTimes;
    report.epochs.threshold.appliedTimesSeconds = thresholdTimes;
    report.epochs.threshold.endpointAdjusted = thresholdAdjustment;
    report.steps.threshold.applied = true;
    report.steps.threshold.nRejected = thresholdInfo.nRejected;
    report.steps.threshold.nLocallyRepaired = ...
        thresholdInfo.nLocallyRepairedEpochs;
    nf_check_job_log(artifactPlan, 'epoch thresholding');

    [EEG, interpolationInfo] = nf_finalize_montage( ...
        EEG, originalMontage, options, channelInfo);
    report.channels.finalization = interpolationInfo;
    report.steps.finalization = interpolationInfo;
    EEG = eeg_checkset(EEG);
    nf_check_job_log(artifactPlan, 'montage finalization');

    quality = struct();
    figureHandle = [];
    if options.qualityCompute
        try
            [quality, figureHandle] = nf_eegquality(EEG_preclean, ...
                EEG_postclean, 'final', EEG, 'report', report, ...
                'plot', logical(options.qualityPlot), ...
                'visible', char(options.qualityVisible), ...
                options.qualityOptions{:});
            report.quality.computed = true;
            report.quality.plotted = logical(options.qualityPlot);
            report.quality.alerts = quality.alerts;
        catch qualityException
            report.quality.error = qualityException.message;
            report.quality.exceptionIdentifier = qualityException.identifier;
            try
                qualityDiagnostic = getReport(qualityException, ...
                    'extended', 'hyperlinks', 'off');
            catch
                qualityDiagnostic = qualityException.message;
            end
            report.quality.diagnostic = qualityDiagnostic;
            quality = struct();
            quality.error = qualityException.message;
            quality.exceptionIdentifier = qualityException.identifier;
            quality.diagnostic = qualityDiagnostic;
            nf_nonfatal_warning('nf_preprocess:QualityFailed', ...
                'Preprocessing completed, but EEG quality evaluation failed: %s', ...
                qualityException.message);
            nf_nonfatal_warning( ...
                'nf_preprocess:QualityDiagnostic', ...
                'QC diagnostic:\n%s', ...
                qualityDiagnostic);
        end
    end

    nf_check_job_log(artifactPlan, 'quality evaluation');
    [report, quality] = nf_export_quality_artifacts( ...
        report, quality, figureHandle, artifactPlan);
    if artifactPlan.exportQuality
        qualityFigureOwned = report.persistence.quality.fig.saved;
        qualityPdfOwned = report.persistence.quality.pdf.saved;
    end
    nf_check_job_log(artifactPlan, 'quality export');
    report.processingFinished = datestr(now, 30); %#ok<TNOW1,DATST>
    report.finished = report.processingFinished;
    report.output = nf_dataset_summary(EEG);
    if artifactPlan.saveDataset
        report.output.filename = artifactPlan.datasetFilename;
        report.output.filepath = artifactPlan.outputDirectory;
    end
    report.software = nf_software_summary();
    report.persistence = nf_prepare_persistence_report( ...
        report.persistence, artifactPlan);

    if ~isfield(EEG, 'etc') || isempty(EEG.etc)
        EEG.etc = struct();
    end
    EEG = nf_remove_duplicate_helper_ledgers(EEG);
    EEG.etc.preprocess = report;
    if ~isfield(EEG.etc, 'nf_preprocess_history')
        EEG.etc.nf_preprocess_history = {};
    elseif ~iscell(EEG.etc.nf_preprocess_history)
        EEG.etc.nf_preprocess_history = {EEG.etc.nf_preprocess_history};
    end
    historyIndex = numel(EEG.etc.nf_preprocess_history) + 1;
    EEG.etc.nf_preprocess_history{historyIndex} = ...
        nf_compact_history(report);

    if artifactPlan.saveDataset
        EEG = nf_save_dataset(EEG, artifactPlan);
        datasetCommitted = true;
        report.persistence = nf_mark_dataset_saved( ...
            report.persistence);
        report.output = nf_dataset_summary(EEG);
        EEG.etc.preprocess = report;
        EEG.etc.nf_preprocess_history{historyIndex} = ...
            nf_compact_history(report);
    end
    nf_check_job_log(artifactPlan, 'dataset save');
    completionTime = datestr(now, 30); %#ok<TNOW1,DATST>
    if artifactPlan.log
        report.persistence = nf_mark_log_completed( ...
            report.persistence, completionTime);
    end
    report.persistenceFinished = completionTime;
    report.finished = report.persistenceFinished;
    report.persistence.finished = report.persistenceFinished;
    report.persistence.snapshotPhase = 'persistence-complete';
    EEG.etc.preprocess = report;
    EEG.etc.nf_preprocess_history{historyIndex} = ...
        nf_compact_history(report);
    nf_write_log_success(artifactPlan, report);
    if artifactPlan.log
        nf_finish_job_log(artifactPlan, diaryState);
        clear diaryCleanup
    end
    if artifactPlan.saveDataset || artifactPlan.log
        clear bundleCleanup
    end

catch preprocessingException
    if artifactPlan.saveDataset && ~datasetCommitted
        nf_rollback_failed_quality_artifacts( ...
            artifactPlan, ...
            qualityFigureOwned, ...
            qualityPdfOwned);
    end
    nf_write_log_failure(artifactPlan, preprocessingException);
    if artifactPlan.log
        try
            nf_finish_job_log(artifactPlan, diaryState);
        catch logFinalizationException
            try
                preprocessingException = addCause( ...
                    preprocessingException, logFinalizationException);
            catch
                nf_nonfatal_warning( ...
                    'nf_preprocess:FailureLogFinalizationFailed', ...
                    ['The original preprocessing error is being rethrown, ' ...
                    'but the failure log could not be finalized: %s'], ...
                    logFinalizationException.message);
            end
        end
        clear diaryCleanup
    end
    rethrow(preprocessingException)
end

end

function plan = nf_resolve_artifact_plan(EEG, options)
plan = struct();
plan.schemaVersion = '2.0.0';
plan.saveDataset = nf_save_requested(options.save);
plan.log = logical(options.log);
plan.exportQuality = logical(plan.saveDataset && options.qualityPlot);
explicitSavePath = nf_is_text(options.save);
sourceLocationRequired = ...
    (plan.log && ~plan.saveDataset) || ...
    (plan.saveDataset && ~explicitSavePath);
plan.sourceDirectory = nf_eeg_source_directory( ...
    EEG, sourceLocationRequired);
plan.sourceStem = nf_eeg_source_stem(EEG);
plan.outputDirectory = '';
plan.datasetFilename = '';
plan.datasetPath = '';
plan.legacyFdtPath = '';
plan.legacyDatPath = '';
plan.qualityFigurePath = '';
plan.qualityPdfPath = '';
plan.logPath = '';
plan.bundleLockPath = '';

[plan.outputDirectory, plan.datasetFilename] = ...
    nf_resolve_save_destination( ...
        options.save, ...
        plan.sourceDirectory, ...
        plan.sourceStem);
if plan.saveDataset
    plan.datasetPath = fullfile( ...
        plan.outputDirectory, ...
        plan.datasetFilename);
    [~, outputStem] = fileparts(plan.datasetFilename);
    if plan.exportQuality
        plan.qualityFigurePath = fullfile( ...
            plan.outputDirectory, ...
            [outputStem '_qc.fig']);
        plan.qualityPdfPath = fullfile( ...
            plan.outputDirectory, ...
            [outputStem '_qc.pdf']);
    end
    if plan.log
        plan.logPath = fullfile( ...
            plan.outputDirectory, ...
            [outputStem '_log.txt']);
    end
elseif plan.log
    plan.outputDirectory = plan.sourceDirectory;
    plan.logPath = fullfile( ...
        plan.outputDirectory, ...
        [plan.sourceStem '_preprocess_log.txt']);
end
if plan.saveDataset
    [~, datasetStem] = fileparts(plan.datasetFilename);
    plan.legacyFdtPath = fullfile( ...
        plan.outputDirectory, [datasetStem '.fdt']);
    plan.legacyDatPath = fullfile( ...
        plan.outputDirectory, [datasetStem '.dat']);
end
if plan.saveDataset || plan.log
    if plan.saveDataset
        lockStem = datasetStem;
    else
        [~, lockStem] = fileparts(plan.logPath);
    end
    plan.bundleLockPath = fullfile(plan.outputDirectory, ...
        ['.' lockStem '.nf_preprocess.lock']);
end
end

function [plan, cleanup] = nf_prepare_artifact_destination(plan)
cleanup = [];
if ~plan.saveDataset && ~plan.log
    return
end

nf_assert_new_artifact(plan.datasetPath);
nf_assert_new_artifact(plan.legacyFdtPath);
nf_assert_new_artifact(plan.legacyDatPath);
nf_assert_new_artifact(plan.qualityFigurePath);
nf_assert_new_artifact(plan.qualityPdfPath);
nf_assert_new_artifact(plan.logPath);

if exist(plan.outputDirectory, 'file') == 2
    error('nf_preprocess:OutputDirectoryIsFile', ...
        'The output directory is an existing file: %s', ...
        plan.outputDirectory);
end
if exist(plan.outputDirectory, 'dir') ~= 7
    [created, message] = mkdir(plan.outputDirectory);
    if ~created && exist(plan.outputDirectory, 'dir') ~= 7
        error('nf_preprocess:CannotCreateOutputDirectory', ...
            'Could not create output directory %s: %s', ...
            plan.outputDirectory, message);
    end
end
nf_assert_directory_writable(plan.outputDirectory);
cleanup = nf_acquire_bundle_lock(plan.bundleLockPath);
nf_assert_new_artifact(plan.datasetPath);
nf_assert_new_artifact(plan.legacyFdtPath);
nf_assert_new_artifact(plan.legacyDatPath);
nf_assert_new_artifact(plan.qualityFigurePath);
nf_assert_new_artifact(plan.qualityPdfPath);
nf_assert_new_artifact(plan.logPath);
end

function cleanup = nf_acquire_bundle_lock(lockPath)
if exist(lockPath, 'file') == 2 || exist(lockPath, 'dir') == 7
    error('nf_preprocess:OutputBundleLocked', ...
        ['Another nf_preprocess job has reserved this output bundle: %s. ' ...
        'If no job is active, remove the stale lock directory manually.'], ...
        lockPath);
end
[created, message, messageIdentifier] = mkdir(lockPath);
if ~created || ~isempty(message) || ~isempty(messageIdentifier)
    error('nf_preprocess:OutputBundleLockFailed', ...
        'Could not reserve output bundle %s: %s', lockPath, message);
end
cleanup = onCleanup(@() nf_release_bundle_lock(lockPath));
end

function nf_release_bundle_lock(lockPath)
if exist(lockPath, 'dir') ~= 7
    return
end
try
    rmdir(lockPath);
catch releaseException
    nf_nonfatal_warning( ...
        'nf_preprocess:OutputBundleLockReleaseFailed', ...
        'Could not release output bundle lock %s: %s', ...
        lockPath, releaseException.message);
end
end

function [directory, filename] = nf_resolve_save_destination( ...
        saveRequest, sourceDirectory, sourceStem)
directory = '';
filename = '';
if ~nf_save_requested(saveRequest)
    return
end

defaultFilename = [sourceStem '_preprocessed.set'];
if nf_is_text(saveRequest)
    requestedPath = strtrim(char(saveRequest));
    requestedPath = nf_absolute_path(requestedPath);
    [requestedDirectory, requestedName, requestedExtension] = ...
        fileparts(requestedPath);
    if strcmpi(requestedExtension, '.set')
        if isempty(requestedName)
            error('nf_preprocess:InvalidSavePath', ...
                'An exact .set path must include a filename stem.');
        end
        if isempty(requestedDirectory)
            requestedDirectory = pwd;
        end
        directory = requestedDirectory;
        filename = [requestedName '.set'];
    else
        if exist(requestedPath, 'file') == 2
            error('nf_preprocess:InvalidSavePath', ...
                ['The supplied save path is an existing file without a ' ...
                '.set extension: %s'], requestedPath);
        end
        directory = requestedPath;
        filename = defaultFilename;
    end
else
    directory = sourceDirectory;
    filename = defaultFilename;
end
end

function directory = nf_eeg_source_directory(EEG, required)
if nargin < 2
    required = false;
end
directory = pwd;
if ~isfield(EEG, 'filepath') || ~nf_is_text(EEG.filepath)
    if required
        error('nf_preprocess:MissingSourceDirectory', ...
            ['EEG.filepath is required to resolve the requested output. ' ...
            'Supply an explicit save directory or .set path.']);
    end
    return
end
candidate = strtrim(char(EEG.filepath));
if isempty(candidate)
    if required
        error('nf_preprocess:MissingSourceDirectory', ...
            ['EEG.filepath is empty. Supply an explicit save directory ' ...
            'or .set path.']);
    end
    return
end
candidate = nf_absolute_path(candidate);
if exist(candidate, 'dir') == 7
    directory = candidate;
elseif required
    error('nf_preprocess:InvalidSourceDirectory', ...
        ['EEG.filepath does not name an existing directory: %s. Supply ' ...
        'an explicit save directory or .set path.'], candidate);
end
end

function stem = nf_eeg_source_stem(EEG)
stem = '';
if isfield(EEG, 'filename') && nf_is_text(EEG.filename)
    filename = strtrim(char(EEG.filename));
    if ~isempty(filename)
        [~, stem] = fileparts(filename);
    end
end
if isempty(stem) && isfield(EEG, 'setname') && nf_is_text(EEG.setname)
    stem = strtrim(char(EEG.setname));
end
if isempty(stem)
    stem = 'nf_dataset';
end
stem = regexprep(stem, '[^A-Za-z0-9_+\-]+', '-');
stem = regexprep(stem, '^[-_]+|[-_]+$', '');
if isempty(stem)
    stem = 'nf_dataset';
end
end

function pathValue = nf_absolute_path(pathValue)
pathValue = char(pathValue);
if isempty(pathValue)
    return
end
if pathValue(1) == '~' && ...
        (isscalar(pathValue) || ismember(pathValue(2), '/\'))
    homeDirectory = getenv('HOME');
    if isempty(homeDirectory)
        homeDirectory = getenv('USERPROFILE');
    end
    if isempty(homeDirectory)
        error('nf_preprocess:CannotExpandHomeDirectory', ...
            'Could not expand the home-directory path: %s', pathValue);
    end
    if isscalar(pathValue)
        pathValue = homeDirectory;
    else
        pathValue = fullfile(homeDirectory, pathValue(3:end));
    end
end
if ~nf_is_absolute_path(pathValue)
    pathValue = fullfile(pwd, pathValue);
end
try
    fileObject = javaObject('java.io.File', pathValue);
    pathValue = char(fileObject.getCanonicalPath());
catch
    pathValue = nf_strip_trailing_separators(pathValue);
end
end

function pathValue = nf_strip_trailing_separators(pathValue)
pathValue = char(pathValue);
while numel(pathValue) > 1 && ismember(pathValue(end), '/\')
    if ispc && ~isempty(regexp(pathValue, ...
            '^[A-Za-z]:[\\/]$', 'once'))
        break
    end
    pathValue(end) = [];
end
end

function absolute = nf_is_absolute_path(pathValue)
absolute = false;
if isempty(pathValue)
    return
end
if pathValue(1) == filesep
    absolute = true;
    return
end
if ispc
    absolute = ~isempty(regexp(pathValue, ...
        '^[A-Za-z]:[\\/]|^\\\\', 'once'));
end
end

function equal = nf_paths_equal(first, second)
first = nf_absolute_path(first);
second = nf_absolute_path(second);
first = nf_strip_trailing_separators(first);
second = nf_strip_trailing_separators(second);
if ispc
    equal = strcmpi(first, second);
else
    equal = strcmp(first, second);
end
end

function nf_assert_new_artifact(pathValue)
if isempty(pathValue)
    return
end
if exist(pathValue, 'file') == 2 || exist(pathValue, 'dir') == 7
    error('nf_preprocess:OutputExists', ...
        ['Output already exists and nf_preprocess never overwrites ' ...
        'artifacts: %s'], pathValue);
end
end

function nf_assert_directory_writable(directory)
probePath = [tempname(directory) '.nf-write-test'];
[fileIdentifier, message] = fopen(probePath, 'w');
if fileIdentifier < 0
    error('nf_preprocess:OutputDirectoryNotWritable', ...
        'Output directory %s is not writable: %s', directory, message);
end
closeStatus = fclose(fileIdentifier);
if exist(probePath, 'file') == 2
    delete(probePath);
end
if closeStatus ~= 0
    error('nf_preprocess:OutputDirectoryNotWritable', ...
        'A test file could not be closed cleanly in %s.', directory);
end
end

function persistence = nf_initial_persistence_report(plan)
persistence = struct();
persistence.schemaVersion = plan.schemaVersion;
persistence.overwrite = false;
persistence.outputDirectory = plan.outputDirectory;
persistence.bundleLockPath = plan.bundleLockPath;
persistence.concurrentReservation = ...
    logical(plan.saveDataset || plan.log);
persistence.dataset = struct();
persistence.dataset.requested = plan.saveDataset;
persistence.dataset.path = plan.datasetPath;
persistence.dataset.saveMode = 'onefile';
persistence.dataset.protectedLegacySidecars = ...
    {plan.legacyFdtPath, plan.legacyDatPath};
persistence.dataset.saved = false;
if plan.saveDataset
    persistence.dataset.status = 'pending';
else
    persistence.dataset.status = 'not-requested';
end
persistence.log = struct();
persistence.log.requested = plan.log;
persistence.log.path = plan.logPath;
persistence.log.capture = 'MATLAB Command Window output';
persistence.log.saved = false;
if plan.log
    persistence.log.status = 'active';
else
    persistence.log.status = 'not-requested';
end
persistence.quality = struct();
persistence.quality.requested = plan.exportQuality;
persistence.quality.fig = nf_artifact_status( ...
    plan.qualityFigurePath, plan.exportQuality);
persistence.quality.pdf = nf_artifact_status( ...
    plan.qualityPdfPath, plan.exportQuality);
end

function status = nf_artifact_status(pathValue, requested)
status = struct();
status.path = pathValue;
status.requested = requested;
status.saved = false;
status.error = '';
if requested
    status.status = 'pending';
else
    status.status = 'not-requested';
end
end

function [report, quality] = nf_export_quality_artifacts( ...
        report, quality, figureHandle, plan)
if ~plan.exportQuality
    return
end
if isempty(figureHandle) || numel(figureHandle) ~= 1 || ...
        ~ishghandle(figureHandle) || ...
        ~strcmpi(get(figureHandle, 'Type'), 'figure')
    message = ['The QC figure was not created. The dataset and job log ' ...
        'will still be preserved.'];
    report.persistence.quality.fig.status = 'failed';
    report.persistence.quality.fig.error = message;
    report.persistence.quality.pdf.status = 'failed';
    report.persistence.quality.pdf.error = message;
    quality = nf_record_quality_export(quality, plan, false, false, ...
        message, message);
    nf_nonfatal_warning( ...
        'nf_preprocess:QualityExportFailed', ...
        '%s', ...
        message);
    return
end

[figSaved, figError] = nf_save_quality_figure( ...
    figureHandle, plan.qualityFigurePath);
[pdfSaved, pdfError] = nf_save_quality_pdf( ...
    figureHandle, plan.qualityPdfPath);
report.persistence.quality.fig.saved = figSaved;
report.persistence.quality.pdf.saved = pdfSaved;
if figSaved
    report.persistence.quality.fig.status = 'saved';
else
    report.persistence.quality.fig.status = 'failed';
    report.persistence.quality.fig.error = figError;
    nf_nonfatal_warning( ...
        'nf_preprocess:QualityFigureSaveFailed', ...
        ['The QC .fig could not be saved, but the dataset and log will ' ...
        'still be preserved: %s'], figError);
end
if pdfSaved
    report.persistence.quality.pdf.status = 'saved';
else
    report.persistence.quality.pdf.status = 'failed';
    report.persistence.quality.pdf.error = pdfError;
    nf_nonfatal_warning( ...
        'nf_preprocess:QualityPdfSaveFailed', ...
        ['The QC PDF could not be saved, but the dataset and log will ' ...
        'still be preserved: %s'], pdfError);
end
quality = nf_record_quality_export(quality, plan, ...
    figSaved, pdfSaved, figError, pdfError);
end

function [saved, errorMessage] = nf_save_quality_figure( ...
        figureHandle, outputPath)
saved = false;
errorMessage = '';
attempted = false;
try
    nf_assert_new_artifact(outputPath);
    attempted = true;
    savefig(figureHandle, outputPath);
    if ~nf_is_nonempty_file(outputPath)
        error('nf_preprocess:QualityFigureSaveFailed', ...
            'savefig returned without creating a nonempty file at %s.', ...
            outputPath);
    end
    saved = true;
catch figureException
    errorMessage = figureException.message;
    if attempted
        nf_delete_partial_artifact(outputPath);
    end
end
end

function [saved, errorMessage] = nf_save_quality_pdf( ...
        figureHandle, outputPath)
saved = false;
errorMessage = '';
attempted = false;
try
    nf_assert_new_artifact(outputPath);
    attempted = true;
    if nf_function_available('exportgraphics')
        try
            exportgraphics(figureHandle, outputPath, ...
                'ContentType', 'vector');
        catch exportException
            nf_delete_partial_artifact(outputPath);
            try
                print(figureHandle, outputPath, ...
                    '-dpdf', '-vector', '-bestfit');
            catch printException
                error('nf_preprocess:QualityPdfSaveFailed', ...
                    ['exportgraphics failed (%s), and the print fallback ' ...
                    'also failed (%s).'], ...
                    exportException.message, printException.message);
            end
        end
    else
        print(figureHandle, outputPath, ...
            '-dpdf', '-vector', '-bestfit');
    end
    if ~nf_is_nonempty_file(outputPath)
        error('nf_preprocess:QualityPdfSaveFailed', ...
            'PDF export returned without creating a nonempty file at %s.', ...
            outputPath);
    end
    saved = true;
catch pdfException
    errorMessage = pdfException.message;
    if attempted
        nf_delete_partial_artifact(outputPath);
    end
end
end

function nf_delete_partial_artifact(pathValue)
if isempty(pathValue)
    return
end
if exist(pathValue, 'file') == 2
    try
        delete(pathValue);
    catch deletionException
        nf_nonfatal_warning( ...
            'nf_preprocess:PartialArtifactCleanupFailed', ...
            'Could not remove partial artifact %s: %s', ...
            pathValue, ...
            deletionException.message);
    end
end
end

function nf_nonfatal_warning(identifier, formatSpec, varargin)
try
    message = sprintf(formatSpec, varargin{:});
catch
    message = char(formatSpec);
end
try
    fprintf(2, '[NeuroFreq warning] %s: %s\n', ...
        identifier, ...
        message);
catch
end
end

function present = nf_is_nonempty_file(pathValue)
present = false;
if exist(pathValue, 'file') ~= 2
    return
end
details = dir(pathValue);
present = isscalar(details) && details.bytes > 0;
end

function quality = nf_record_quality_export(quality, plan, ...
        figSaved, pdfSaved, figError, pdfError)
if ~isstruct(quality) || numel(quality) ~= 1
    quality = struct();
end
quality.export = struct();
quality.export.figPath = plan.qualityFigurePath;
quality.export.pdfPath = plan.qualityPdfPath;
quality.export.figSaved = figSaved;
quality.export.pdfSaved = pdfSaved;
quality.export.figError = figError;
quality.export.pdfError = pdfError;
end

function persistence = nf_prepare_persistence_report( ...
        persistence, plan)
if plan.saveDataset
    persistence.snapshotPhase = 'set-serialization';
    persistence.finalStatusLocation = ...
        ['Returned report/EEG and, when requested, the completed job log. ' ...
        'The report embedded in the .set records the serialization boundary.'];
    persistence.dataset.saved = true;
    persistence.dataset.status = 'committed-by-containing-set';
    persistence.dataset.evidence = ...
        ['A readable .set containing this ledger is the successful output ' ...
        'of the guarded pop_saveset serialization.'];
elseif plan.log
    persistence.snapshotPhase = 'pre-log-finalization';
    persistence.finalStatusLocation = ...
        'Returned report/EEG and the completed job log.';
else
    persistence.snapshotPhase = 'processing-complete';
    persistence.finalStatusLocation = 'Returned report/EEG.';
end
if plan.log
    if plan.saveDataset
        persistence.log.status = 'active-at-dataset-save';
    else
        persistence.log.status = 'active-before-log-finalization';
    end
end
if plan.exportQuality
    persistence.quality.complete = ...
        persistence.quality.fig.saved && ...
        persistence.quality.pdf.saved;
else
    persistence.quality.complete = false;
end
end

function persistence = nf_mark_dataset_saved(persistence)
persistence.dataset.saved = true;
persistence.dataset.status = 'saved';
persistence.dataset.savedAt = datestr(now, 30); %#ok<TNOW1,DATST>
end

function persistence = nf_mark_log_completed( ...
        persistence, completedAt)
persistence.log.saved = true;
persistence.log.status = 'completed';
persistence.log.savedAt = completedAt;
end

function EEG = nf_save_dataset(EEG, plan)
datasetWriteAttempted = false;
try
    nf_assert_new_artifact(plan.datasetPath);
    nf_assert_new_artifact(plan.legacyFdtPath);
    nf_assert_new_artifact(plan.legacyDatPath);
    datasetWriteAttempted = true;
    EEG = pop_saveset(EEG, ...
        'filename', plan.datasetFilename, ...
        'filepath', plan.outputDirectory, ...
        'savemode', 'onefile', ...
        'check', 'on');
    if ~nf_is_nonempty_file(plan.datasetPath)
        error('nf_preprocess:DatasetSaveFailed', ...
            ['pop_saveset returned without creating a nonempty .set ' ...
            'file at %s.'], plan.datasetPath);
    end
catch saveException
    if datasetWriteAttempted
        fprintf(2, ['Dataset save failed; removing partial dataset files ' ...
            'created by this job.\n']);
        nf_delete_partial_artifact(plan.datasetPath);
        nf_delete_partial_artifact(plan.legacyFdtPath);
        nf_delete_partial_artifact(plan.legacyDatPath);
    end
    rethrow(saveException)
end
end

function nf_rollback_failed_quality_artifacts( ...
        plan, figureOwned, pdfOwned)
if ~figureOwned && ~pdfOwned
    return
end
fprintf(2, ['Preprocessing failed before dataset commit; removing QC files ' ...
    'created by this job so a no-overwrite rerun remains possible.\n']);
if figureOwned
    nf_delete_partial_artifact(plan.qualityFigurePath);
end
if pdfOwned
    nf_delete_partial_artifact(plan.qualityPdfPath);
end
end

function state = nf_start_job_log(plan, options)
state = struct();
state.previousStatus = 'off';
previousFile = '';
try
    previousStatus = get(0, 'Diary');
    if nf_is_text(previousStatus)
        state.previousStatus = char(previousStatus);
    end
catch
end
try
    previousFile = get(0, 'DiaryFile');
catch
end
if nf_is_text(previousFile)
    state.previousFile = char(previousFile);
else
    state.previousFile = '';
end
state.jobLogPath = plan.logPath;
if strcmpi(state.previousStatus, 'on') && ...
        ~isempty(state.previousFile) && ...
        ~nf_is_absolute_path(state.previousFile)
    try
        state.previousFile = nf_absolute_path(state.previousFile);
    catch pathException
        nf_nonfatal_warning( ...
            'nf_preprocess:PriorDiaryPathResolutionFailed', ...
            ['Could not resolve the prior relative DiaryFile %s: %s. ' ...
            'Restoration will use the reported relative path.'], ...
            state.previousFile, ...
            pathException.message);
    end
end
if strcmpi(state.previousStatus, 'on')
    try
        fprintf('\n[NeuroFreq] Redirecting active diary to job log: %s\n', ...
            plan.logPath);
    catch
    end
    try
        diary off
    catch diaryException
        nf_nonfatal_warning( ...
            'nf_preprocess:PriorDiaryStopFailed', ...
            'Could not stop the prior diary before job logging: %s', ...
            diaryException.message);
    end
end

try
    reportedOptions = nf_report_options(options);
    optionsText = evalc('disp(reportedOptions)');
catch optionsException
    optionsText = sprintf( ...
        '<Resolved options could not be formatted: %s>\n', ...
        optionsException.message);
end
try
    header = sprintf( ...
        ['\n============================================================\n' ...
        'NeuroFreq nf_preprocess job log\n' ...
        'Started: %s\n' ...
        'Source: %s\n' ...
        'Output directory: %s\n' ...
        'Capture: MATLAB Command Window output (source echo disabled)\n' ...
        '============================================================\n\n' ...
        'Resolved options:\n' ...
        '%s\n'], ...
        datestr(now, 31), ... %#ok<TNOW1,DATST>
        plan.sourceStem, ...
        plan.outputDirectory, ...
        optionsText);
catch headerException
    header = sprintf( ...
        ['\nNeuroFreq nf_preprocess job log\n' ...
        'Header formatting failed: %s\n'], ...
        headerException.message);
end
try
    diary(plan.logPath)
    [active, observedStatus, observedFile] = ...
        nf_job_diary_status(plan);
    if ~active
        error('nf_preprocess:LogInitializationFailed', ...
            ['MATLAB reported diary status %s and DiaryFile %s after ' ...
            'the requested initialization at %s.'], ...
            observedStatus, ...
            observedFile, ...
            plan.logPath);
    end
    count = fprintf(1, '%s', header);
    if count < numel(header)
        error('nf_preprocess:LogInitializationFailed', ...
            'MATLAB did not print the complete job-log header.');
    end
catch logException
    try
        currentStatus = get(0, 'Diary');
        if nf_is_text(currentStatus) && strcmpi(char(currentStatus), 'on')
            diary off
        end
    catch
    end
    initializationMarker = sprintf( ...
        ['\n============================================================\n' ...
        'LOG CAPTURE INITIALIZATION FAILED\n' ...
        'Expected DiaryFile: %s\n' ...
        'Initialization error: %s\n' ...
        'Logging will retry at the first checkpoint.\n' ...
        '============================================================\n'], ...
        plan.logPath, ...
        logException.message);
    fallbackText = [header initializationMarker];
    appended = nf_append_job_log_text(plan.logPath, fallbackText);
    try
        fprintf(1, '%s', fallbackText);
    catch
    end
    if appended
        nf_nonfatal_warning( ...
            'nf_preprocess:LogInitializationRecovered', ...
            ['MATLAB diary initialization failed, but the header was ' ...
            'written directly to %s. Logging will retry at the first ' ...
            'checkpoint. Cause: %s'], ...
            plan.logPath, ...
            logException.message);
    else
        nf_nonfatal_warning( ...
            'nf_preprocess:LogInitializationFailed', ...
            ['MATLAB diary initialization failed, and the job-log header ' ...
            'could not be written to %s. Logging will retry at the first ' ...
            'checkpoint. Cause: %s'], ...
            plan.logPath, ...
            logException.message);
    end
end
end

function nf_finish_job_log(plan, state)
try
    [active, ~, ~] = nf_job_diary_status(plan);
    if active
        diary off
    end
catch completionException
    nf_nonfatal_warning( ...
        'nf_preprocess:LogCompletionFailed', ...
        'Could not stop the completed job diary cleanly: %s', ...
        completionException.message);
end
try
    if ~nf_is_nonempty_file(plan.logPath)
        nf_nonfatal_warning( ...
            'nf_preprocess:LogCompletionFailed', ...
            'The completed job log is missing or empty: %s', ...
            plan.logPath);
    end
catch completionException
    nf_nonfatal_warning( ...
        'nf_preprocess:LogCompletionFailed', ...
        'Could not verify the completed job log %s: %s', ...
        plan.logPath, ...
        completionException.message);
end
nf_restore_diary(state);
end

function nf_restore_diary(state)
try
    if nf_diary_state_matches(state)
        return
    end
catch
end
try
    if strcmpi(get(0, 'Diary'), 'on')
        diary off
    end
    if strcmpi(state.previousStatus, 'on')
        if isempty(state.previousFile)
            diary on
        else
            diary(state.previousFile)
        end
        fprintf(['\n[NeuroFreq] Resumed diary after nf_preprocess job ' ...
            'log: %s\n'], state.jobLogPath);
    elseif ~isempty(state.previousFile)
        set(0, 'DiaryFile', state.previousFile);
    end
    if ~nf_diary_state_matches(state)
        error('nf_preprocess:DiaryRestoreFailed', ...
            ['MATLAB diary state did not match its pre-job status and ' ...
            'target after restoration.']);
    end
catch restoreException
    nf_nonfatal_warning( ...
        'nf_preprocess:DiaryRestoreFailed', ...
        'The previous MATLAB diary could not be restored: %s', ...
        restoreException.message);
end
end

function matches = nf_diary_state_matches(state)
currentStatus = get(0, 'Diary');
currentFile = get(0, 'DiaryFile');
statusMatches = strcmpi(currentStatus, state.previousStatus);
if isempty(state.previousFile)
    fileMatches = true;
elseif ~nf_is_text(currentFile) || isempty(strtrim(char(currentFile)))
    fileMatches = false;
else
    fileMatches = nf_paths_equal(currentFile, state.previousFile);
end
matches = statusMatches && fileMatches;
end

function [active, observedStatus, observedFile] = ...
        nf_job_diary_status(plan)
active = false;
observedStatus = 'unknown';
observedFile = '';
try
    currentStatus = get(0, 'Diary');
    if nf_is_text(currentStatus)
        observedStatus = char(currentStatus);
    end
catch
end
try
    currentFile = get(0, 'DiaryFile');
    if nf_is_text(currentFile)
        observedFile = char(currentFile);
    end
catch
end
if ~strcmpi(observedStatus, 'on')
    return
end
active = nf_diary_file_matches(observedFile, plan.logPath);
end

function matches = nf_diary_file_matches(observedFile, expectedFile)
matches = false;
if isempty(observedFile) || isempty(expectedFile)
    return
end
try
    matches = nf_paths_equal(observedFile, expectedFile);
catch
    matches = false;
end
end

function active = nf_ensure_job_diary(plan, stage)
[active, observedStatus, observedFile] = ...
    nf_job_diary_status(plan);
if active
    return
end
observedDirectory = pwd;
try
    currentStatus = get(0, 'Diary');
    if nf_is_text(currentStatus) && strcmpi(char(currentStatus), 'on')
        diary off
    end
    diary(plan.logPath)
    [verified, recoveredStatus, recoveredFile] = ...
        nf_job_diary_status(plan);
    if ~strcmpi(recoveredStatus, 'on')
        error('nf_preprocess:LogRecoveryFailed', ...
            'MATLAB did not report an active diary after restart.');
    end
    if ~verified
        error('nf_preprocess:LogRecoveryWrongTarget', ...
            ['MATLAB restarted the diary but reported DiaryFile %s ' ...
            'instead of %s.'], ...
            recoveredFile, ...
            plan.logPath);
    end
    active = true;
catch recoveryException
    active = false;
    try
        currentStatus = get(0, 'Diary');
        if nf_is_text(currentStatus) && strcmpi(char(currentStatus), 'on')
            diary off
        end
    catch
    end
    marker = sprintf( ...
        ['\n============================================================\n' ...
        'LOG CAPTURE RECOVERY FAILED\n' ...
        'Stage: %s\n' ...
        'Observed status: %s\n' ...
        'Observed DiaryFile: %s\n' ...
        'Expected DiaryFile: %s\n' ...
        'Working directory: %s\n' ...
        'Recovery error: %s\n' ...
        'Further Command Window output may be absent from this log.\n' ...
        '============================================================\n'], ...
        stage, ...
        observedStatus, ...
        observedFile, ...
        plan.logPath, ...
        observedDirectory, ...
        recoveryException.message);
    appended = nf_append_job_log_text(plan.logPath, marker);
    try
        fprintf(2, '%s', marker);
    catch
    end
    if ~appended
        nf_nonfatal_warning( ...
            'nf_preprocess:LogRecoveryFailed', ...
            ['The job diary could not be restarted, and the recovery ' ...
            'diagnostic could not be appended to %s.'], ...
            plan.logPath);
    end
    return
end
marker = sprintf( ...
    ['\n============================================================\n' ...
    'LOG CAPTURE RESTARTED\n' ...
    'Stage: %s\n' ...
    'Observed status before restart: %s\n' ...
    'Observed DiaryFile before restart: %s\n' ...
    'Expected DiaryFile: %s\n' ...
    'Working directory: %s\n' ...
    'Output produced during the interruption may be missing.\n' ...
    '============================================================\n'], ...
    stage, ...
    observedStatus, ...
    observedFile, ...
    plan.logPath, ...
    observedDirectory);
markerWritten = false;
try
    count = fprintf(2, '%s', marker);
    markerWritten = count >= numel(marker);
catch markerException
    nf_nonfatal_warning( ...
        'nf_preprocess:LogRecoveryMarkerFailed', ...
        'Could not print the log-recovery marker during %s: %s', ...
        stage, ...
        markerException.message);
end
if ~markerWritten
    try
        diary off
    catch
    end
    active = false;
    appended = nf_append_job_log_text(plan.logPath, marker);
    if ~appended
        nf_nonfatal_warning( ...
            'nf_preprocess:LogRecoveryMarkerFailed', ...
            ['The job diary restarted, but its recovery marker could not ' ...
            'be recorded in %s.'], ...
            plan.logPath);
    end
end
end

function appended = nf_append_job_log_text(pathValue, textValue)
appended = false;
fileIdentifier = -1;
try
    [fileIdentifier, message] = fopen(pathValue, 'a');
    if fileIdentifier < 0
        error('nf_preprocess:LogAppendFailed', ...
            'Could not open %s for append: %s', ...
            pathValue, ...
            message);
    end
    count = fprintf(fileIdentifier, '%s', textValue);
    if count < numel(textValue)
        error('nf_preprocess:LogAppendFailed', ...
            'MATLAB did not append the complete log text to %s.', ...
            pathValue);
    end
    closeStatus = fclose(fileIdentifier);
    fileIdentifier = -1;
    if closeStatus ~= 0
        error('nf_preprocess:LogAppendFailed', ...
            'MATLAB could not close the job log cleanly: %s', ...
            pathValue);
    end
    appended = true;
catch
    if fileIdentifier >= 0
        try
            fclose(fileIdentifier);
        catch
        end
    end
end
end

function nf_emit_job_log_text(plan, textValue, stage, errorStream)
active = false;
try
    active = nf_ensure_job_diary(plan, stage);
catch ensureException
    nf_nonfatal_warning( ...
        'nf_preprocess:LogRecoveryFailed', ...
        'Could not verify or restart the job diary during %s: %s', ...
        stage, ...
        ensureException.message);
end
if errorStream
    stream = 2;
else
    stream = 1;
end
if ~active
    appended = nf_append_job_log_text(plan.logPath, textValue);
    if ~appended
        nf_nonfatal_warning( ...
            'nf_preprocess:LogAppendFailed', ...
            'Could not append job-log text during %s to %s.', ...
            stage, ...
            plan.logPath);
    end
end
screenWritten = false;
try
    count = fprintf(stream, '%s', textValue);
    screenWritten = count >= numel(textValue);
catch outputException
    nf_nonfatal_warning( ...
        'nf_preprocess:LogScreenWriteFailed', ...
        'Could not print job-log text during %s: %s', ...
        stage, ...
        outputException.message);
end
if active && ~screenWritten
    try
        diary off
    catch
    end
    appended = nf_append_job_log_text(plan.logPath, textValue);
    if ~appended
        nf_nonfatal_warning( ...
            'nf_preprocess:LogAppendFailed', ...
            'Could not append job-log text during %s to %s.', ...
            stage, ...
            plan.logPath);
    end
end
end

function nf_check_job_log(plan, stage)
if ~plan.log
    return
end
try
    checkpoint = sprintf('[NeuroFreq] Checkpoint: %s | %s\n', ...
        stage, ...
        datestr(now, 31)); %#ok<TNOW1,DATST>
    nf_emit_job_log_text(plan, checkpoint, stage, false);
catch checkpointException
    nf_nonfatal_warning( ...
        'nf_preprocess:LogCheckpointFailed', ...
        'Could not record the %s checkpoint: %s', ...
        stage, ...
        checkpointException.message);
end
end

function nf_write_log_success(plan, report)
if ~plan.log
    return
end
try
    textValue = sprintf( ...
        ['\n============================================================\n' ...
        'nf_preprocess completed successfully\n' ...
        'Finished: %s\n'], ...
        datestr(now, 31)); %#ok<TNOW1,DATST>
    if plan.saveDataset
        textValue = sprintf( ...
            '%sDataset: %s\n', ...
            textValue, ...
            plan.datasetPath);
    end
    if plan.exportQuality
        textValue = sprintf( ...
            '%sQC FIG: %s [%s]\n', ...
            textValue, ...
            plan.qualityFigurePath, ...
            report.persistence.quality.fig.status);
        textValue = sprintf( ...
            '%sQC PDF: %s [%s]\n', ...
            textValue, ...
            plan.qualityPdfPath, ...
            report.persistence.quality.pdf.status);
    end
    textValue = sprintf( ...
        ['%sLog: %s\n' ...
        '============================================================\n'], ...
        textValue, ...
        plan.logPath);
    nf_emit_job_log_text( ...
        plan, ...
        textValue, ...
        'success reporting', ...
        false);
catch logException
    nf_nonfatal_warning( ...
        'nf_preprocess:SuccessLogWriteFailed', ...
        'Could not write the final success block to the job log: %s', ...
        logException.message);
end
end

function nf_write_log_failure(plan, preprocessingException)
if ~plan.log
    return
end
try
    try
        diagnostic = getReport(preprocessingException, ...
            'extended', 'hyperlinks', 'off');
    catch
        diagnostic = preprocessingException.message;
    end
    textValue = sprintf( ...
        ['\n============================================================\n' ...
        'nf_preprocess FAILED\n' ...
        'Finished: %s\n' ...
        '%s\n' ...
        'Log preserved at: %s\n' ...
        '============================================================\n'], ...
        datestr(now, 31), ... %#ok<TNOW1,DATST>
        diagnostic, ...
        plan.logPath);
    nf_emit_job_log_text( ...
        plan, ...
        textValue, ...
        'failure reporting', ...
        true);
catch logException
    nf_nonfatal_warning( ...
        'nf_preprocess:FailureLogWriteFailed', ...
        'Could not append the final failure diagnostic to the job log: %s', ...
        logException.message);
end
end

function [EEG, info] = nf_run_gedai(EEG, configuration)
gedaiArguments = {configuration.artifactThresholdType, ...
    configuration.epochSizeInCycles, configuration.lowcutFrequency, ...
    configuration.referenceMatrixType, configuration.parallel, ...
    configuration.visualizeArtifacts, ...
    configuration.enovaThresholdPerEpoch, ...
    configuration.enovaThresholdPerChannel, configuration.signalType, ...
    configuration.smoothingWindowSeconds, ...
    configuration.outputReferenceChannel};

warningState = warning;
pathState = path;
workingDirectory = pwd;
randomState = rng;
environmentCleanup = onCleanup(@() nf_restore_matlab_environment( ...
    warningState, pathState, workingDirectory, randomState));
[EEG, artifacts, sensai, sensaiByBand, thresholdByBand, meanEnova, ...
    enovaByEpoch, command, enovaByBand, enovaByChannel] = ...
    GEDAI(EEG, gedaiArguments{:});
clear environmentCleanup
EEG = eeg_checkset(EEG);

info = struct();
info.applied = true;
info.configuration = nf_compact_gedai_configuration(configuration);
info.sensaiScore = sensai;
info.sensaiScorePerBand = sensaiByBand;
info.artifactThresholdPerBand = thresholdByBand;
info.meanEnova = meanEnova;
info.enovaPerEpoch = enovaByEpoch;
info.enovaPerBand = enovaByBand;
info.enovaPerChannel = enovaByChannel;
info.command = command;
info.artifacts = nf_artifact_summary(artifacts);
info.nativeLedgerRemoved = false;
if isfield(EEG, 'etc') && isfield(EEG.etc, 'GEDAI')
    EEG.etc = rmfield(EEG.etc, 'GEDAI');
    info.nativeLedgerRemoved = true;
end
end

function nf_restore_matlab_environment(warningState, pathState, ...
    workingDirectory, randomState)
try
    cd(workingDirectory);
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
try
    rng(randomState);
catch
end
end

function configuration = nf_resolve_gedai_configuration(supplied)
configuration = nf_default_gedai_configuration();
if isstruct(supplied)
    suppliedFields = fieldnames(supplied);
    allowedFields = fieldnames(configuration);
    for index = 1:numel(suppliedFields)
        fieldName = suppliedFields{index};
        if ~ismember(fieldName, allowedFields)
            error('nf_preprocess:UnknownGEDAIOption', ...
                'Unknown gedaiOptions field: %s.', fieldName);
        end
        configuration.(fieldName) = supplied.(fieldName);
    end
elseif iscell(supplied)
    if ~ismember(numel(supplied), [10 11])
        error('nf_preprocess:InvalidGEDAIOptions', ...
            'A positional gedaiOptions cell must contain 10 or 11 values.');
    end
    fields = {'artifactThresholdType', 'epochSizeInCycles', ...
        'lowcutFrequency', 'referenceMatrixType', 'parallel', ...
        'visualizeArtifacts', 'enovaThresholdPerEpoch', ...
        'enovaThresholdPerChannel', 'signalType', ...
        'smoothingWindowSeconds', 'outputReferenceChannel'};
    for index = 1:numel(supplied)
        configuration.(fields{index}) = supplied{index};
    end
end
configuration = nf_normalize_gedai_configuration(configuration);
end

function configuration = nf_normalize_gedai_configuration(configuration)
configuration.artifactThresholdType = lower(strtrim( ...
    nf_require_text(configuration.artifactThresholdType, ...
    'gedaiOptions.artifactThresholdType')));
if ~ismember(configuration.artifactThresholdType, {'auto-', 'auto', 'auto+'})
    error('nf_preprocess:InvalidGEDAIOptions', ...
        'GEDAI artifactThresholdType must be auto-, auto, or auto+.');
end
nf_require_positive(configuration.epochSizeInCycles, ...
    'gedaiOptions.epochSizeInCycles');
nf_require_positive(configuration.lowcutFrequency, ...
    'gedaiOptions.lowcutFrequency');
reference = configuration.referenceMatrixType;
if ischar(reference) || (isstring(reference) && isscalar(reference))
    reference = lower(strtrim(char(reference)));
    if ~ismember(reference, {'precomputed', 'interpolated', 'warped'})
        error('nf_preprocess:InvalidGEDAIOptions', ...
            ['GEDAI referenceMatrixType must be precomputed, interpolated, ' ...
            'warped, or a custom covariance matrix.']);
    end
    configuration.referenceMatrixType = reference;
else
    nf_validate_covariance(reference, 'gedaiOptions.referenceMatrixType');
end
if ~nf_is_logical_scalar(configuration.parallel) || ...
        ~nf_is_logical_scalar(configuration.visualizeArtifacts)
    error('nf_preprocess:InvalidGEDAIOptions', ...
        'GEDAI parallel and visualizeArtifacts must be logical scalars.');
end
configuration.parallel = logical(configuration.parallel);
configuration.visualizeArtifacts = logical(configuration.visualizeArtifacts);
nf_require_nonnegative_or_inf(configuration.enovaThresholdPerEpoch, ...
    'gedaiOptions.enovaThresholdPerEpoch');
nf_require_nonnegative_or_inf(configuration.enovaThresholdPerChannel, ...
    'gedaiOptions.enovaThresholdPerChannel');
configuration.signalType = lower(strtrim(nf_require_text( ...
    configuration.signalType, 'gedaiOptions.signalType')));
if ~strcmp(configuration.signalType, 'eeg')
    error('nf_preprocess:InvalidGEDAIOptions', ...
        'nf_preprocess supports GEDAI signalType=''eeg'' only.');
end
nf_require_positive_or_inf(configuration.smoothingWindowSeconds, ...
    'gedaiOptions.smoothingWindowSeconds');
outputReference = strtrim(nf_require_text( ...
    configuration.outputReferenceChannel, ...
    'gedaiOptions.outputReferenceChannel'));
if isempty(outputReference) || strcmpi(outputReference, 'avgref')
    outputReference = 'AvgRef';
elseif strcmpi(outputReference, 'rest')
    outputReference = 'REST';
end
configuration.outputReferenceChannel = outputReference;
if isnumeric(configuration.referenceMatrixType) && strcmp(outputReference, 'REST')
    error('nf_preprocess:UnsupportedGEDAIReference', ...
        ['GEDAI v1.7 does not apply a valid REST transform when a custom ' ...
        'reference covariance is supplied. Use AvgRef or a named channel.']);
end
end

function configuration = nf_validate_gedai_montage(configuration, chanlocs)
reference = configuration.referenceMatrixType;
if isnumeric(reference) && ~isequal(size(reference), ...
        [numel(chanlocs) numel(chanlocs)])
    error('nf_preprocess:GEDAIReferenceSize', ...
        ['A custom GEDAI reference covariance must be sized to the selected ' ...
        'EEG montage (%d x %d).'], numel(chanlocs), numel(chanlocs));
end
outputReference = configuration.outputReferenceChannel;
if ~ismember(lower(outputReference), {'avgref', 'rest'})
    labels = {chanlocs.labels};
    if sum(strcmpi(labels, outputReference)) ~= 1
        error('nf_preprocess:GEDAIOutputReference', ...
            'GEDAI output reference channel %s is not in the selected montage.', ...
            outputReference);
    end
end
end

function nf_validate_gedai_standard_labels(chanlocs)
gedaiDirectory = fileparts(which('GEDAI'));
databasePath = fullfile(gedaiDirectory, 'auxiliaries', ...
    'fsavLEADFIELD_4_GEDAI.mat');
try
    database = load(databasePath, 'leadfield4GEDAI');
catch databaseException
    error('nf_preprocess:GEDAIStandardDatabaseUnreadable', ...
        'GEDAI''s standard database could not be read: %s', ...
        databaseException.message);
end
if ~isfield(database, 'leadfield4GEDAI') || ...
        ~isstruct(database.leadfield4GEDAI) || ...
        ~isfield(database.leadfield4GEDAI, 'electrodes') || ...
        ~isstruct(database.leadfield4GEDAI.electrodes) || ...
        ~isfield(database.leadfield4GEDAI.electrodes, 'Name')
    error('nf_preprocess:GEDAIStandardDatabaseInvalid', ...
        ['GEDAI''s standard database does not contain the expected ' ...
        'leadfield4GEDAI.electrodes.Name labels.']);
end
databaseLabels = strtrim(string( ...
    {database.leadfield4GEDAI.electrodes.Name}));
inputLabels = strtrim(string({chanlocs.labels}));
[present, databaseIndices] = ismember(lower(inputLabels), ...
    lower(databaseLabels));
if any(~present)
    missing = strjoin(cellstr(inputLabels(~present)), ', ');
    error('nf_preprocess:GEDAIStandardLabelMismatch', ...
        ['GEDAI''s precomputed standard database has no exact label match ' ...
        'for: %s. For a non-10-5 montage, set ' ...
        'gedaiOptions.referenceMatrixType=''interpolated'' so GEDAI maps ' ...
        'the standard leadfield to the supplied coordinates.'], missing);
end
if numel(unique(databaseIndices)) ~= numel(databaseIndices)
    error('nf_preprocess:GEDAIStandardLabelMismatch', ...
        'The selected montage does not map uniquely to GEDAI''s standard database.');
end
end

function configuration = nf_subset_gedai_reference(configuration, ...
    originalMontage, retainedMontage)
originalLabels = lower(string({originalMontage.labels}));
retainedLabels = lower(string({retainedMontage.labels}));
[present, retainedIndices] = ismember(retainedLabels, originalLabels);
if any(~present)
    error('nf_preprocess:GEDAIReferenceOrder', ...
        'The retained channel montage cannot be mapped to the GEDAI reference montage.');
end
if isnumeric(configuration.referenceMatrixType)
    originalSize = size(configuration.referenceMatrixType);
    configuration.referenceMatrixType = ...
        configuration.referenceMatrixType(retainedIndices, retainedIndices);
    configuration.referenceMatrixOriginalSize = originalSize;
    configuration.referenceMatrixSubsetIndices = retainedIndices;
end
outputReference = configuration.outputReferenceChannel;
if ~ismember(lower(outputReference), {'avgref', 'rest'}) && ...
        ~any(strcmpi({retainedMontage.labels}, outputReference))
    error('nf_preprocess:GEDAIOutputReferenceRemoved', ...
        ['GEDAI output reference channel %s was removed by bad-channel ' ...
        'screening. Choose AvgRef/REST or another retained reference.'], ...
        outputReference);
end
end

function compact = nf_compact_gedai_configuration(configuration)
compact = configuration;
if isnumeric(configuration.referenceMatrixType)
    matrix = double(configuration.referenceMatrixType);
    descriptor = struct();
    descriptor.type = 'custom covariance';
    descriptor.size = size(matrix);
    descriptor.trace = trace(matrix);
    descriptor.frobeniusNorm = norm(matrix, 'fro');
    compact.referenceMatrixType = descriptor;
end
end

function configuration = nf_default_gedai_configuration()
configuration = struct();
configuration.artifactThresholdType = 'auto';
configuration.epochSizeInCycles = 12;
configuration.lowcutFrequency = 0.5;
configuration.referenceMatrixType = 'precomputed';
configuration.parallel = false;
configuration.visualizeArtifacts = false;
configuration.enovaThresholdPerEpoch = Inf;
configuration.enovaThresholdPerChannel = Inf;
configuration.signalType = 'eeg';
configuration.smoothingWindowSeconds = Inf;
configuration.outputReferenceChannel = 'AvgRef';
end

function summary = nf_artifact_summary(artifacts)
summary = struct();
summary.available = false;
summary.nonzeroFraction = NaN;
summary.channelNonzeroFraction = [];
if isstruct(artifacts) && isfield(artifacts, 'data') && ~isempty(artifacts.data)
    data = reshape(artifacts.data, size(artifacts.data, 1), []);
elseif isnumeric(artifacts) && ~isempty(artifacts)
    data = reshape(artifacts, size(artifacts, 1), []);
else
    return
end
summary.available = true;
summary.nonzeroFraction = nnz(data) / numel(data);
summary.channelNonzeroFraction = sum(data ~= 0, 2) ./ size(data, 2);
end

function [EEG, info] = nf_make_final_epochs(EEG, options)
EEG = nf_ensure_event_fields(EEG);
if ~isempty(options.events)
    survivingValidation = nf_validate_requested_events(EEG, options.events);
    [EEG, acceptedIndices] = pop_epoch(EEG, ...
        survivingValidation.popEpochEvents, ...
        options.epochLimits, 'epochinfo', 'yes');
    mode = 'event-locked';
    eventTypes = options.events;
else
    survivingValidation = struct();
    acceptedIndices = [];
    EEG = eeg_regepochs(EEG, 'recurrence', options.continuousEpochLength, ...
        'limits', [0 options.continuousEpochLength], 'rmbase', NaN, ...
        'eventtype', 'nf_fixed_epoch');
    mode = 'fixed-length';
    eventTypes = {};
end
EEG = eeg_checkset(EEG);
if EEG.trials < 2
    error('nf_preprocess:InsufficientEpochs', ...
        'Final epoching produced fewer than two complete epochs.');
end
if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
EEG.etc.nf_epoch_ids = 1:EEG.trials;

info = struct();
info.mode = mode;
info.events = eventTypes;
info.inputEventValidation = options.eventValidation;
info.preEpochEventValidation = survivingValidation;
info.acceptedEventIndices = acceptedIndices;
info.nAcceptedEvents = numel(acceptedIndices);
info.requestedLimitsSeconds = options.epochLimits;
info.actualLimitsSeconds = [EEG.xmin EEG.xmax];
info.nCreated = EEG.trials;
info.pnts = EEG.pnts;
end

function [window, adjusted] = nf_fit_epoch_window(window, actualLimits, ...
    sampleInterval, optionName)
window = reshape(window, 1, 2);
actualLimits = reshape(actualLimits, 1, 2);
adjusted = false;
tolerance = sampleInterval + 10 * eps(max(abs(actualLimits)) + 1);
if window(1) < actualLimits(1)
    if actualLimits(1) - window(1) <= tolerance
        window(1) = actualLimits(1);
        adjusted = true;
    else
        error('nf_preprocess:EpochWindowOutsideData', ...
            '%s begins before the actual epoch.', optionName);
    end
end
if window(2) > actualLimits(2)
    if window(2) - actualLimits(2) <= tolerance
        window(2) = actualLimits(2);
        adjusted = true;
    else
        error('nf_preprocess:EpochWindowOutsideData', ...
            '%s ends after the actual epoch.', optionName);
    end
end
if window(1) >= window(2)
    error('nf_preprocess:EpochWindowOutsideData', ...
        '%s does not span at least two ordered time points.', optionName);
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

function events = nf_pop_epoch_events(events)
if isstring(events)
    events = cellstr(events);
elseif ischar(events)
    events = {events};
elseif isnumeric(events)
    events = num2cell(events);
end
end

function [EEG, info] = nf_finalize_montage(EEG, originalMontage, options, ...
    channelInfo)
hadIcaBeforeFinalization = isfield(EEG, 'icaweights') && ...
    ~isempty(EEG.icaweights);
currentLabels = lower(string({EEG.chanlocs.labels}));
originalLabels = lower(string({originalMontage.labels}));
missingMask = ~ismember(originalLabels, currentLabels);
referenceMask = false(1, numel(originalMontage));
referenceRemoved = false;
referenceLabel = '';
if isfield(channelInfo, 'reference') && isstruct(channelInfo.reference) && ...
        isfield(channelInfo.reference, 'removedZeroReference') && ...
        logical(channelInfo.reference.removedZeroReference)
    referenceLabel = channelInfo.reference.label;
    referenceMask = strcmpi({originalMontage.labels}, referenceLabel);
    referenceRemoved = sum(referenceMask) == 1;
end
artifactMissingMask = missingMask & ~referenceMask;
missingLabels = {originalMontage(artifactMissingMask).labels};

info = struct();
info.requestedGlobalInterpolation = logical(options.globalInterpolation);
info.missingBeforeInterpolation = missingLabels;
info.nMissingBeforeInterpolation = numel(missingLabels);
info.referenceRemoved = referenceRemoved;
info.referenceLabel = referenceLabel;
info.referenceRestored = false;
info.interpolationMethod = char(options.interpolationMethod);
info.interpolated = false;
info.rereferenced = false;

if options.globalInterpolation && ~isempty(missingLabels)
    interpolationMontage = originalMontage(~referenceMask);
    EEG = eeg_interp(EEG, interpolationMontage, ...
        char(options.interpolationMethod));
    EEG = eeg_checkset(EEG);
    info.interpolated = true;
end
if options.rereference
    if referenceRemoved
        referenceLocation = originalMontage(referenceMask);
        EEG = pop_reref(EEG, [], 'refloc', referenceLocation, ...
            'refica', 'remove');
        info.referenceRestored = true;
    else
        EEG = pop_reref(EEG, [], 'refica', 'remove');
    end
    EEG = eeg_checkset(EEG);
    info.rereferenced = true;
    info.reference = 'common average';
else
    info.reference = 'unchanged';
end
EEG = nf_clear_final_ica(EEG);
[EEG, orderRestored] = nf_restore_channel_order(EEG, originalMontage);
info.icaDecompositionRemoved = hadIcaBeforeFinalization;
info.originalChannelOrderRestored = orderRestored;
info.outputLabels = {EEG.chanlocs.labels};
end

function [EEG, restored] = nf_restore_channel_order(EEG, originalMontage)
expected = lower(string({originalMontage.labels}));
current = lower(string({EEG.chanlocs.labels}));
[present, originalIndices] = ismember(current, expected);
if any(~present) || numel(unique(originalIndices)) ~= numel(originalIndices)
    error('nf_preprocess:MontageRestorationFailed', ...
        'The final channel set cannot be mapped uniquely to the original montage.');
end
[~, order] = sort(originalIndices);
restored = ~isequal(order, 1:numel(order));
if restored
    EEG.data = EEG.data(order, :, :);
    EEG.chanlocs = EEG.chanlocs(order);
    if isfield(EEG, 'reject') && isstruct(EEG.reject)
        rejectionFields = {'rejmanualE', 'rejthreshE', 'rejconstE', ...
            'rejfreqE', 'rejjpE', 'rejkurtE', 'rejglobalE'};
        for index = 1:numel(rejectionFields)
            fieldName = rejectionFields{index};
            if isfield(EEG.reject, fieldName)
                value = EEG.reject.(fieldName);
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
end

function EEG = nf_clear_final_ica(EEG)
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

function [EEG, info] = nf_select_eeg_channels(EEG, requested)
labels = {EEG.chanlocs.labels};
nf_require_unique_labels(labels);

if isempty(requested)
    selected = true(1, EEG.nbchan);
    typesAvailable = isfield(EEG.chanlocs, 'type');
    if typesAvailable
        for index = 1:EEG.nbchan
            channelType = EEG.chanlocs(index).type;
            if isempty(channelType)
                selected(index) = true;
            else
                selected(index) = strcmpi(strtrim(char(channelType)), 'eeg');
            end
        end
        if ~any(selected)
            error('nf_preprocess:NoEEGChannels', ...
                'No channels typed EEG or empty were found. Supply eegChannels explicitly.');
        end
    end
    selectionSource = 'chanlocs.type (EEG or empty)';
else
    selected = false(1, EEG.nbchan);
    indices = nf_resolve_channels(EEG, requested);
    selected(indices) = true;
    selectionSource = 'eegChannels option';
end

if sum(selected) < 3
    error('nf_preprocess:TooFewEEGChannels', ...
        'At least three scalp EEG channels are required.');
end
excludedLabels = labels(~selected);
if any(~selected)
    EEG = pop_select(EEG, 'channel', find(selected));
    EEG = eeg_checkset(EEG);
end

info = struct();
info.source = selectionSource;
info.selectedOriginalIndices = find(selected);
info.selectedLabels = labels(selected);
info.excludedLabels = excludedLabels;
info.nSelected = sum(selected);
info.nExcluded = sum(~selected);
end

function indices = nf_resolve_channels(EEG, requested)
if islogical(requested)
    if numel(requested) ~= EEG.nbchan
        error('nf_preprocess:InvalidEEGChannels', ...
            'A logical eegChannels mask must contain EEG.nbchan elements.');
    end
    indices = find(requested);
elseif isnumeric(requested)
    indices = requested(:)';
else
    if ischar(requested)
        requested = {requested};
    elseif isstring(requested)
        requested = cellstr(requested);
    end
    labels = lower(string({EEG.chanlocs.labels}));
    [present, indices] = ismember(lower(string(requested)), labels);
    if any(~present)
        missing = strjoin(cellstr(string(requested(~present))), ', ');
        error('nf_preprocess:InvalidEEGChannels', ...
            'Requested EEG channels were not found: %s', missing);
    end
end
if isempty(indices) || any(~isfinite(indices)) || any(indices ~= round(indices)) || ...
        any(indices < 1) || any(indices > EEG.nbchan) || ...
        numel(unique(indices)) ~= numel(indices)
    error('nf_preprocess:InvalidEEGChannels', ...
        'eegChannels must resolve to unique valid channel indices.');
end
end

function method = nf_resolve_ica_method(preset, supplied)
if ~ismember(preset, {'adult', 'child'})
    error('nf_preprocess:UnknownPreset', ...
        'preset must be ''adult'' or ''child''.');
end
if strcmp(preset, 'adult')
    expected = 'iclabel';
else
    expected = 'adjustedadjust';
end
if isempty(supplied)
    method = expected;
    return
end
method = lower(char(supplied));
if ismember(method, {'made', 'adjusted_adjust', 'adjusted_adust', ...
        'adjusted-adjust', 'adjusted-adust', 'adjustedadjust', ...
        'adjustedadust'})
    method = 'adjustedadjust';
end
if ~strcmp(method, expected)
    error('nf_preprocess:PresetClassifierInvariant', ...
        ['The %s preset requires %s. adult and child share every other ' ...
        'stage and differ only in their component classifier.'], ...
        preset, expected);
end
end

function nf_validate_options(EEG, options)
if mod(numel(options.qualityOptions), 2) ~= 0
    error('nf_preprocess:InvalidQualityOptions', ...
        'qualityOptions must contain nf_eegquality name/value pairs.');
end
for optionIndex = 1:2:numel(options.qualityOptions)
    if ~nf_is_text(options.qualityOptions{optionIndex})
        error('nf_preprocess:InvalidQualityOptions', ...
            'Every qualityOptions parameter name must be scalar text.');
    end
end
reservedQualityOptions = {'final', 'report', 'plot', 'visible'};
for optionIndex = 1:numel(reservedQualityOptions)
    if nf_quality_options_contains(options.qualityOptions, ...
            reservedQualityOptions{optionIndex})
        error('nf_preprocess:ReservedQualityOption', ...
            ['qualityOptions cannot contain %s because nf_preprocess ' ...
            'controls that nf_eegquality input.'], ...
            reservedQualityOptions{optionIndex});
    end
end
if nf_quality_options_contains(options.qualityOptions, 'savePath')
    error('nf_preprocess:ReservedQualitySavePath', ...
        ['qualityOptions cannot contain savePath. nf_preprocess owns all ' ...
        'persistence so requested .fig, .pdf, .set, and log artifacts ' ...
        'remain in one guarded output bundle.']);
end
if EEG.trials ~= 1
    error('nf_preprocess:EpochedInput', ...
        ['nf_preprocess requires continuous input so filtering, GEDAI, and ' ...
        'ICA training do not cross or create epoch boundaries.']);
end
if options.highpass >= options.lowpass
    error('nf_preprocess:InvalidPassband', ...
        'highpass must be below lowpass.');
end
if options.lowpass >= EEG.srate / 2 || options.lowpass >= options.resample / 2
    error('nf_preprocess:InvalidLowpass', ...
        'lowpass must be below both input and output Nyquist frequencies.');
end
if options.highpass >= EEG.srate / 2 || ...
        options.highpass >= options.resample / 2
    error('nf_preprocess:InvalidHighpass', ...
        'highpass must be below both input and output Nyquist frequencies.');
end
if options.muscleRange(1) < 0 || ...
        options.muscleRange(2) > options.lowpass || ...
        options.muscleRange(2) >= options.resample / 2
    error('nf_preprocess:InvalidMuscleBand', ...
        'The complete muscleRange must be retained by lowpass and resampling.');
end
if options.icaTrainingFrequencies(1) < 0 || ...
        options.icaTrainingFrequencies(2) > options.lowpass || ...
        options.icaTrainingFrequencies(2) >= options.resample / 2
    error('nf_preprocess:InvalidICATrainingBand', ...
        'The complete ICA training frequency range must be retained.');
end
if ~isempty(options.events) && isempty(options.epochLimits)
    error('nf_preprocess:MissingEpochLimits', ...
        'epochLimits are required when events are supplied.');
end
if isempty(options.events) && ~isempty(options.epochLimits)
    error('nf_preprocess:UnexpectedEpochLimits', ...
        'epochLimits apply only when events are supplied.');
end
if ~isempty(options.events) && ...
        (~isfield(EEG, 'event') || isempty(EEG.event))
    error('nf_preprocess:MissingEvents', ...
        'The dataset contains no events for event-locked epoching.');
end
if ~isempty(options.epochLimits)
    finalLimits = options.epochLimits;
else
    finalLimits = [0 options.continuousEpochLength];
end
if isnumeric(options.baseline) && numel(options.baseline) == 2
    epochMilliseconds = finalLimits * 1000;
    if options.baseline(1) < epochMilliseconds(1) || ...
            options.baseline(2) > epochMilliseconds(2)
        error('nf_preprocess:InvalidBaseline', ...
            'baseline must fall inside epochLimits.');
    end
end
if ~isempty(options.thresholdTimes) && ...
        (options.thresholdTimes(1) < finalLimits(1) || ...
        options.thresholdTimes(2) > finalLimits(2))
    error('nf_preprocess:InvalidThresholdTimes', ...
        'thresholdTimes must fall inside the final epoch limits.');
end
if ~ismember(options.channelMethod, {'faster', 'cleanrawdata', 'none'})
    error('nf_preprocess:UnknownChannelMethod', ...
        'channelMethod must be faster, cleanrawdata, or none.');
end
if ~ismember(options.precleanMethod, {'gedai', 'none'})
    error('nf_preprocess:UnknownPrecleanMethod', ...
        'precleanMethod must be gedai or none.');
end
if ~ismember(options.icaAlgorithm, {'runamica15', 'runica'})
    error('nf_preprocess:UnknownICAAlgorithm', ...
        'icaAlgorithm must be runamica15 or runica.');
end
if options.cleanHighpass >= options.resample / 2
    error('nf_preprocess:InvalidCleanHighpass', ...
        'cleanHighpass must be below the post-resampling Nyquist frequency.');
end
if options.icaTrainingHighpass >= options.resample / 2
    error('nf_preprocess:InvalidICATrainingHighpass', ...
        'icaTrainingHighpass must be below the post-resampling Nyquist frequency.');
end
if strcmp(options.precleanMethod, 'gedai')
    if ~strcmp(options.gedaiOptions.outputReferenceChannel, 'AvgRef')
        error('nf_preprocess:UnsupportedGEDAIReference', ...
            ['nf_preprocess requires GEDAI outputReferenceChannel=''AvgRef''. ' ...
            'The final pipeline reference is applied only after ICA and ' ...
            'montage restoration.']);
    end
    if options.gedaiOptions.lowcutFrequency >= EEG.srate / 2 || ...
            options.gedaiOptions.lowcutFrequency >= options.resample / 2
        error('nf_preprocess:InvalidGEDAILowcut', ...
            ['gedaiOptions.lowcutFrequency must be below both input and ' ...
            'post-resampling Nyquist frequencies.']);
    end
    if ~isinf(options.gedaiOptions.enovaThresholdPerEpoch) || ...
            ~isinf(options.gedaiOptions.enovaThresholdPerChannel)
        error('nf_preprocess:UnsafeGEDAIRejection', ...
            ['nf_preprocess requires both GEDAI ENOVA rejection thresholds ' ...
            'to be Inf. GEDAI v1.7 removes samples without EEGLAB boundary ' ...
            'events; FASTER and nf_thresh own rejection in this pipeline.']);
    end
end
if strcmp(options.preset, 'child') && options.resample < 100
    error('nf_preprocess:ChildSamplingRate', ...
        ['The child adjusted_ADJUST feature extractor requires a sampling ' ...
        'rate of at least 100 Hz.']);
end
if strcmp(options.preset, 'child') && ...
        options.icaTrainingEpochLength * options.resample < 100
    error('nf_preprocess:ChildEpochLength', ...
        ['icaTrainingEpochLength must yield at least 100 samples per MADE ' ...
        'adjusted_ADJUST feature epoch.']);
end
if nf_contains_boundary_event(EEG)
    error('nf_preprocess:BoundaryUnsupported', ...
        ['GEDAI v1.7 and the common ICA-training filter do not honor EEGLAB ' ...
        'boundary events. Segment the recording before nf_preprocess.']);
end
nf_validate_faster_options(options.fasterOptions);
if strcmp(options.icaMethod, 'adjustedadjust') && options.aggressiveICA
    error('nf_preprocess:ClassifierOptionMismatch', ...
        'aggressiveICA applies only to ICLabel.');
end
if strcmp(options.icaMethod, 'adjustedadjust') && ...
        ~isempty(options.iclabelThresholds)
    error('nf_preprocess:ClassifierOptionMismatch', ...
        'iclabelThresholds apply only to ICLabel.');
end
end

function nf_preflight_dependencies(options)
required = {'eeg_checkset', 'pop_select', 'pop_eegfiltnew', ...
    'pop_resample', 'eeg_regepochs', 'pop_eegthresh', 'pop_rejspec', ...
    'pop_rejepoch', 'pop_subcomp', 'eeg_interp', ...
    'pop_reref', 'convertlocs', 'nf_filter', ...
    'nf_badchans', 'nf_cleanic', 'nf_thresh'};
if strcmp(options.icaAlgorithm, 'runica')
    required{end + 1} = 'pop_runica';
end
if ~isempty(options.events)
    required{end + 1} = 'pop_epoch';
end
if isnumeric(options.baseline) && numel(options.baseline) == 2
    required{end + 1} = 'pop_rmbase';
end
for index = 1:numel(required)
    if exist(required{index}, 'file') ~= 2
        error('nf_preprocess:MissingDependency', ...
            'Required function %s was not found on the MATLAB path.', required{index});
    end
end
if strcmp(options.channelMethod, 'faster')
    fasterDependencies = {'channel_properties', 'min_z', ...
        'distancematrix', 'hurst_exponent', ...
        'nanmean', 'nanmedian'};
    for index = 1:numel(fasterDependencies)
        if ~nf_function_available(fasterDependencies{index})
            error('nf_preprocess:MissingFASTER', ...
                ['FASTER dependency %s was not found. Install the complete ' ...
                'FASTER distribution and its MATLAB dependencies.'], ...
                fasterDependencies{index});
        end
    end
end
if strcmp(options.channelMethod, 'cleanrawdata') && ...
        exist('clean_channels', 'file') ~= 2
    error('nf_preprocess:MissingCleanRawData', ...
        'clean_channels.m from clean_rawdata is required.');
end
if strcmp(options.precleanMethod, 'gedai')
    if exist('GEDAI', 'file') ~= 2
        error('nf_preprocess:MissingGEDAI', ...
            'GEDAI.m v1.7 is required.');
    end
    try
        gedaiInputs = nargin('GEDAI');
        gedaiOutputs = nargout('GEDAI');
    catch
        gedaiInputs = NaN;
        gedaiOutputs = NaN;
    end
    if gedaiInputs ~= -12 || gedaiOutputs ~= 10
        error('nf_preprocess:GEDAIVersionMismatch', ...
            ['GEDAI.m does not expose the v1.7 contract expected by this ' ...
            'pipeline (11 supplied arguments after EEG and 10 outputs).']);
    end
    if ~nf_function_available('gpuDeviceCount')
        error('nf_preprocess:MissingGEDAIDependency', ...
            ['GEDAI v1.7 calls gpuDeviceCount even in CPU mode. The MATLAB ' ...
            'Parallel Computing Toolbox function was not found.']);
    end
    if ischar(options.gedaiOptions.referenceMatrixType) && ...
            strcmp(options.gedaiOptions.referenceMatrixType, 'precomputed')
        gedaiDirectory = fileparts(which('GEDAI'));
        databasePath = fullfile(gedaiDirectory, 'auxiliaries', ...
            'fsavLEADFIELD_4_GEDAI.mat');
        if exist(databasePath, 'file') ~= 2
            error('nf_preprocess:MissingGEDAIStandardDatabase', ...
                ['GEDAI''s standard precomputed lead-field database was not ' ...
                'found at %s. Install the complete GEDAI v1.7 distribution.'], ...
                databasePath);
        end
    end
end
if strcmp(options.icaAlgorithm, 'runamica15') && ...
        (exist('pop_runamica', 'file') ~= 2 || exist('runamica15', 'file') ~= 2)
    error('nf_preprocess:MissingAMICA', ...
        'The current AMICA plugin is required for icaAlgorithm=runamica15.');
end
if strcmp(options.icaMethod, 'iclabel') && exist('pop_iclabel', 'file') ~= 2
    error('nf_preprocess:MissingICLabel', ...
        'The ICLabel plugin is required for the adult preset.');
end
if strcmp(options.icaMethod, 'adjustedadjust')
    childDependencies = {'adjusted_ADJUST', 'compute_GD_feat', ...
        'computeSED_NOnorm', 'computeSAD', 'EM', 'trim_and_mean', ...
        'trim_and_max', 'MARA_extract_time_freq_features', ...
        'beall_horizontal', 'beall_blink_detection', 'Spatial_Info_eyes', ...
        'spectopo', 'fitlm', 'findpeaks', 'kurt'};
    for index = 1:numel(childDependencies)
        if exist(childDependencies{index}, 'file') ~= 2
            error('nf_preprocess:MissingAdjustedAdjust', ...
                ['The child preset requires %s.m from the complete MADE/' ...
                'adjusted_ADJUST dependency set.'], childDependencies{index});
        end
    end
end
if nf_save_requested(options.save) && ...
        ~nf_function_available('pop_saveset')
    error('nf_preprocess:MissingSaveFunction', ...
        'EEGLAB pop_saveset.m is required when save is requested.');
end
end

function summary = nf_dataset_summary(EEG)
summary = struct();
if isfield(EEG, 'setname')
    summary.setname = EEG.setname;
else
    summary.setname = '';
end
if isfield(EEG, 'filename')
    summary.filename = EEG.filename;
else
    summary.filename = '';
end
if isfield(EEG, 'filepath')
    summary.filepath = EEG.filepath;
else
    summary.filepath = '';
end
summary.nbchan = EEG.nbchan;
summary.channelLabels = {EEG.chanlocs.labels};
summary.pntsPerEpoch = EEG.pnts;
summary.trials = EEG.trials;
summary.srate = EEG.srate;
summary.aggregateDataSeconds = (EEG.pnts * EEG.trials) / EEG.srate;
summary.epochLimitsSeconds = [EEG.xmin EEG.xmax];
end

function software = nf_software_summary()
software = struct();
software.matlab = version;
software.eeglab = '';
if exist('eeg_getversion', 'file') == 2
    software.eeglab = eeg_getversion;
end
software.neuroFreqPreprocessSchema = '3.0.0';
end

function history = nf_compact_history(report)
history = struct();
history.schemaVersion = report.schemaVersion;
history.started = report.started;
history.finished = report.finished;
if isfield(report, 'processingFinished')
    history.processingFinished = report.processingFinished;
end
if isfield(report, 'persistenceFinished')
    history.persistenceFinished = report.persistenceFinished;
end
history.preset = report.preset;
history.classifier = report.presetDefinition.classifier;
history.input = report.input;
history.output = report.output;
history.quality = report.quality;
history.persistence = report.persistence;
end

function values = nf_normalize_legacy_events(values)
if isempty(values)
    return
end
if mod(numel(values), 2) == 1
    values = [{'events', values{1}} values(2:end)];
end
end

function row = nf_row_pair(value)
if isempty(value)
    row = [];
else
    row = reshape(value, 1, 2);
end
end

function reported = nf_report_options(options)
reported = options;
if isfield(reported, 'eventValidation')
    reported = rmfield(reported, 'eventValidation');
end
if isfield(reported, 'gedaiOptions')
    reported.gedaiOptions = nf_compact_gedai_configuration( ...
        reported.gedaiOptions);
end
end

function present = nf_quality_options_contains(options, requestedName)
present = false;
for index = 1:2:numel(options)
    value = options{index};
    if nf_is_text(value) && strcmpi(char(value), requestedName)
        present = true;
        return
    end
end
end

function EEG = nf_normalize_channel_locations(EEG)
if ~nf_function_available('convertlocs')
    error('nf_preprocess:MissingLocationConverter', ...
        'EEGLAB convertlocs.m is required to normalize channel geometry.');
end
hasCartesian = true;
hasTopographic = true;
for index = 1:numel(EEG.chanlocs)
    hasCartesian = hasCartesian && ...
        nf_location_has_xyz(EEG.chanlocs(index));
    hasTopographic = hasTopographic && ...
        nf_location_has_topography(EEG.chanlocs(index));
end
try
    if hasCartesian
        EEG.chanlocs = convertlocs(EEG.chanlocs, 'cart2all');
    elseif hasTopographic
        EEG.chanlocs = convertlocs(EEG.chanlocs, 'topo2all');
    else
        error('nf_preprocess:MissingChannelCoordinates', ...
            ['Every selected channel needs either finite nonzero X/Y/Z ' ...
            'coordinates or finite theta/radius coordinates.']);
    end
catch locationException
    if strcmp(locationException.identifier, ...
            'nf_preprocess:MissingChannelCoordinates')
        rethrow(locationException)
    end
    error('nf_preprocess:InvalidChannelCoordinates', ...
        'EEGLAB could not normalize the channel locations: %s', ...
        locationException.message);
end
nf_validate_interpolation_locations(EEG.chanlocs);
end

function valid = nf_location_has_xyz(chanloc)
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

function valid = nf_location_has_topography(chanloc)
valid = isfield(chanloc, 'theta') && isfield(chanloc, 'radius') && ...
    isnumeric(chanloc.theta) && isreal(chanloc.theta) && ...
    isscalar(chanloc.theta) && isfinite(chanloc.theta) && ...
    isnumeric(chanloc.radius) && isreal(chanloc.radius) && ...
    isscalar(chanloc.radius) && isfinite(chanloc.radius) && ...
    chanloc.radius >= 0;
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
    case 'spacetime'
        error('nf_preprocess:UnsafeInterpolationMethod', ...
            ['spacetime interpolation ignores epoch boundaries and is not ' ...
            'supported by nf_preprocess.']);
    otherwise
        error('nf_preprocess:UnknownInterpolationMethod', ...
            ['interpolationMethod must be spherical, sphericalKang, ' ...
            'sphericalCRD, invdist, or v4.']);
end
end

function nf_validate_faster_options(options)
if ~nf_is_scalar_struct(options)
    error('nf_preprocess:InvalidFASTEROptions', ...
        'fasterOptions must be one scalar structure.');
end
allowed = {'measure', 'z', 'stat'};
names = fieldnames(options);
for index = 1:numel(names)
    if ~ismember(names{index}, allowed)
        error('nf_preprocess:UnknownFASTEROption', ...
            'Unknown fasterOptions field: %s.', names{index});
    end
end
if isfield(options, 'measure')
    value = options.measure;
    if ~(isnumeric(value) || islogical(value)) || numel(value) ~= 3 || ...
            any(~isfinite(value(:))) || any(~ismember(value(:), [0 1]))
        error('nf_preprocess:InvalidFASTEROptions', ...
            'fasterOptions.measure must contain three binary values.');
    end
end
if isfield(options, 'z')
    value = options.z;
    if ~isnumeric(value) || ~isreal(value) || numel(value) ~= 3 || ...
            any(~isfinite(value(:))) || any(value(:) <= 0)
        error('nf_preprocess:InvalidFASTEROptions', ...
            'fasterOptions.z must contain three positive finite values.');
    end
end
if isfield(options, 'stat')
    value = options.stat;
    if ~nf_is_text(value) || ...
            ~ismember(lower(strtrim(char(value))), {'iqr', 'z'})
        error('nf_preprocess:InvalidFASTEROptions', ...
            'fasterOptions.stat must be ''iqr'' or ''z''.');
    end
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
    if ischar(value) && strcmpi(strtrim(value), 'boundary')
        present = true;
        return
    end
end
end

function nf_validate_covariance(value, name)
if ~isnumeric(value) || ~isreal(value) || isempty(value) || ...
        ~ismatrix(value) || size(value, 1) ~= size(value, 2) || ...
        any(~isfinite(value(:)))
    error('nf_preprocess:InvalidGEDAIOptions', ...
        '%s must be a finite real square covariance matrix.', name);
end
matrix = double(value);
scale = max(1, norm(matrix, 'fro'));
symmetryTolerance = 100 * eps(scale) * max(size(matrix));
if norm(matrix - matrix', 'fro') > symmetryTolerance
    error('nf_preprocess:InvalidGEDAIOptions', ...
        '%s must be symmetric.', name);
end
eigenvalues = eig((matrix + matrix') / 2);
positiveTolerance = 100 * eps(scale) * max(size(matrix));
if min(eigenvalues) < -positiveTolerance || trace(matrix) <= 0
    error('nf_preprocess:InvalidGEDAIOptions', ...
        '%s must be positive semidefinite with positive total variance.', name);
end
end

function textValue = nf_require_text(value, name)
if ~nf_is_text(value)
    error('nf_preprocess:InvalidOption', '%s must be scalar text.', name);
end
textValue = char(value);
end

function nf_require_positive(value, name)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value <= 0
    error('nf_preprocess:InvalidOption', ...
        '%s must be a positive finite scalar.', name);
end
end

function nf_require_nonnegative_or_inf(value, name)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    ((isfinite(value) && value >= 0) || (isinf(value) && value > 0));
if ~valid
    error('nf_preprocess:InvalidOption', ...
        '%s must be nonnegative or positive Inf.', name);
end
end

function nf_require_positive_or_inf(value, name)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    ((isfinite(value) && value > 0) || (isinf(value) && value > 0));
if ~valid
    error('nf_preprocess:InvalidOption', ...
        '%s must be positive or positive Inf.', name);
end
end

function EEG = nf_remove_duplicate_helper_ledgers(EEG)
if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc)
    return
end
fields = {'nf_filter', 'nf_badchans', 'nf_cleanic', 'nf_thresh'};
for index = 1:numel(fields)
    if isfield(EEG.etc, fields{index})
        EEG.etc = rmfield(EEG.etc, fields{index});
    end
end
historyFields = {'nf_filter_history', 'nf_cleanic_history'};
for index = 1:numel(historyFields)
    fieldName = historyFields{index};
    if ~isfield(EEG.etc, fieldName) || ~iscell(EEG.etc.(fieldName)) || ...
            isempty(EEG.etc.(fieldName))
        continue
    end
    EEG.etc.(fieldName)(end) = [];
    if isempty(EEG.etc.(fieldName))
        EEG.etc = rmfield(EEG.etc, fieldName);
    end
end
end

function available = nf_function_available(name)
available = exist(name, 'file') == 2 || exist(name, 'builtin') == 5;
end

function nf_require_unique_labels(labels)
for index = 1:numel(labels)
    if ~(ischar(labels{index}) || ...
            (isstring(labels{index}) && isscalar(labels{index}))) || ...
            isempty(strtrim(char(labels{index})))
        error('nf_preprocess:InvalidChannelLabels', ...
            'Every selected channel must have a nonempty scalar-text label.');
    end
end
normalized = lower(string(labels));
if any(strlength(normalized) == 0) || ...
        numel(unique(normalized)) ~= numel(normalized)
    error('nf_preprocess:InvalidChannelLabels', ...
        'Every selected channel must have a unique nonempty label.');
end
end

function nf_validate_interpolation_locations(chanlocs)
coordinates = zeros(numel(chanlocs), 3);
for index = 1:numel(chanlocs)
    if ~nf_location_has_xyz(chanlocs(index)) || ...
            ~nf_location_has_topography(chanlocs(index))
        error('nf_preprocess:MissingChannelCoordinates', ...
            ['Channel %s lacks finite nonzero X/Y/Z and finite theta/radius ' ...
            'coordinates required by FASTER and EEGLAB interpolation.'], ...
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
        duplicateOffset = find(distances < 1e-10, 1);
        second = first + duplicateOffset;
        error('nf_preprocess:DuplicateChannelCoordinates', ...
            'Channels %s and %s occupy the same scalp coordinate.', ...
            char(chanlocs(first).labels), char(chanlocs(second).labels));
    end
end
end

function info = nf_validate_requested_events(EEG, requested)
nf_validate_event_latencies(EEG);
requestedValues = nf_pop_epoch_events(requested);
if isempty(requestedValues)
    error('nf_preprocess:EventTypesNotFound', ...
        'At least one event type must be requested for event-locked epoching.');
end
availableValues = {EEG.event.type};
availableAreText = ischar(availableValues{1});
matchedIndices = cell(1, numel(requestedValues));
normalizedRequested = cell(1, numel(requestedValues));
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
            error('nf_preprocess:InvalidEventType', ...
                'Requested event types must be scalar text or numbers.');
        end
        normalizedRequested{requestIndex} = normalized;
        availableText = cellfun(@deblank, availableValues, ...
            'UniformOutput', false);
        matchedIndices{requestIndex} = find(strcmp(availableText, normalized));
    else
        if ischar(requestedValue) || ...
                (isstring(requestedValue) && isscalar(requestedValue))
            normalized = str2double(strtrim(requestedValue));
        else
            normalized = requestedValue;
        end
        if ~isnumeric(normalized) || ~isscalar(normalized) || ...
                ~isfinite(normalized)
            error('nf_preprocess:InvalidEventType', ...
                'A requested event type cannot be matched to numeric EEG events.');
        end
        normalizedRequested{requestIndex} = normalized;
        numericAvailable = cell2mat(availableValues);
        matchedIndices{requestIndex} = find(numericAvailable == normalized);
    end
    if isempty(matchedIndices{requestIndex})
        error('nf_preprocess:EventTypeNotFound', ...
            ['Requested event type %s does not occur with pop_epoch''s exact, ' ...
            'case-sensitive matching.'], nf_event_display(requestedValue));
    end
end
info = struct();
info.matching = 'pop_epoch exact and case-sensitive';
info.requested = normalizedRequested;
info.matchCounts = cellfun(@numel, matchedIndices);
info.matchedEventIndices = matchedIndices;
info.nUniqueMatchedEvents = numel(unique([matchedIndices{:}]));
info.popEpochEvents = normalizedRequested;
end

function displayValue = nf_event_display(value)
if isnumeric(value) && isscalar(value)
    displayValue = num2str(value, 15);
elseif ischar(value)
    displayValue = ['''' value ''''];
else
    displayValue = ['<' class(value) '>'];
end
end

function EEG = nf_normalize_event_types(EEG)
if ~isfield(EEG, 'event') || isempty(EEG.event)
    return
end
if ~isfield(EEG.event, 'type')
    error('nf_preprocess:InvalidEvents', ...
        'Every EEG.event entry must contain a type field.');
end
types = {EEG.event.type};
containsText = any(cellfun(@(value) ischar(value) || ...
    (isstring(value) && isscalar(value)), types));
for index = 1:numel(types)
    value = types{index};
    if isstring(value) && isscalar(value)
        value = char(value);
    end
    if containsText && isnumeric(value) && isscalar(value) && isfinite(value)
        value = num2str(value);
    end
    if containsText
        valid = ischar(value);
    else
        valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
            isfinite(value);
    end
    if ~valid
        error('nf_preprocess:InvalidEvents', ...
            'EEG event types must be homogeneous scalar text or finite numbers.');
    end
    EEG.event(index).type = value;
end
nf_validate_event_latencies(EEG);
end

function nf_validate_event_latencies(EEG)
if ~isfield(EEG, 'event') || isempty(EEG.event)
    return
end
if ~isfield(EEG.event, 'latency')
    error('nf_preprocess:InvalidEvents', ...
        'Every EEG.event entry must contain a latency field.');
end
for index = 1:numel(EEG.event)
    latency = EEG.event(index).latency;
    if ~isnumeric(latency) || ~isreal(latency) || ~isscalar(latency) || ...
            ~isfinite(latency) || latency < 0.5 || latency > EEG.pnts + 0.5
        error('nf_preprocess:InvalidEvents', ...
            'EEG.event(%d).latency is outside the continuous sample range.', index);
    end
    if isfield(EEG.event, 'urevent') && ~isempty(EEG.event(index).urevent)
        ureventIndex = EEG.event(index).urevent;
        if ~isnumeric(ureventIndex) || ~isscalar(ureventIndex) || ...
                ~isfinite(ureventIndex) || ureventIndex ~= round(ureventIndex) || ...
                ureventIndex < 1 || ~isfield(EEG, 'urevent') || ...
                ureventIndex > numel(EEG.urevent)
            error('nf_preprocess:InvalidEvents', ...
                'EEG.event(%d).urevent is not a valid urevent index.', index);
        end
    end
end
if isfield(EEG, 'urevent') && ~isempty(EEG.urevent)
    if ~isstruct(EEG.urevent) || ~isfield(EEG.urevent, 'latency')
        error('nf_preprocess:InvalidEvents', ...
            'EEG.urevent must contain latency fields.');
    end
    for index = 1:numel(EEG.urevent)
        latency = EEG.urevent(index).latency;
        if ~isnumeric(latency) || ~isreal(latency) || ...
                ~isscalar(latency) || ~isfinite(latency)
            error('nf_preprocess:InvalidEvents', ...
                'EEG.urevent(%d).latency must be a finite scalar.', index);
        end
    end
end
end

function nf_validate_input_eeg(EEG)
if ~isstruct(EEG) || numel(EEG) ~= 1
    error('nf_preprocess:InvalidEEG', ...
        'EEG must be one EEGLAB dataset structure.');
end
required = {'data', 'srate', 'nbchan', 'pnts', 'trials', ...
    'chanlocs', 'xmin', 'xmax'};
for index = 1:numel(required)
    if ~isfield(EEG, required{index})
        error('nf_preprocess:InvalidEEG', ...
            'EEG.%s is required.', required{index});
    end
end
dimensions = {'nbchan', 'pnts', 'trials'};
for index = 1:numel(dimensions)
    value = EEG.(dimensions{index});
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
            ~isfinite(value) || value < 1 || value ~= round(value)
        error('nf_preprocess:InvalidEEG', ...
            'EEG.%s must be a positive integer scalar.', dimensions{index});
    end
end
if ~isnumeric(EEG.srate) || ~isreal(EEG.srate) || ...
        ~isscalar(EEG.srate) || ~isfinite(EEG.srate) || EEG.srate <= 0
    error('nf_preprocess:InvalidEEG', ...
        'EEG.srate must be a positive finite scalar.');
end
if ~isnumeric(EEG.xmin) || ~isreal(EEG.xmin) || ...
        ~isscalar(EEG.xmin) || ~isfinite(EEG.xmin) || ...
        ~isnumeric(EEG.xmax) || ~isreal(EEG.xmax) || ...
        ~isscalar(EEG.xmax) || ~isfinite(EEG.xmax) || ...
        EEG.xmin >= EEG.xmax
    error('nf_preprocess:InvalidEEG', ...
        'EEG.xmin and EEG.xmax must be ordered finite scalars.');
end
if ~isnumeric(EEG.data) || ~isreal(EEG.data) || isempty(EEG.data) || ...
        ndims(EEG.data) > 3 || EEG.nbchan ~= size(EEG.data, 1) || ...
        EEG.pnts ~= size(EEG.data, 2) || ...
        EEG.trials ~= size(EEG.data, 3) || ...
        numel(EEG.chanlocs) ~= EEG.nbchan || any(~isfinite(EEG.data(:)))
    error('nf_preprocess:InvalidEEG', ...
        'EEG data, dimensions, channel locations, or sampling rate are invalid.');
end
expectedSpan = (EEG.pnts - 1) / EEG.srate;
timeTolerance = max(1e-9, 0.51 / EEG.srate);
if EEG.trials == 1 && (abs(EEG.xmin) > timeTolerance || ...
        abs((EEG.xmax - EEG.xmin) - expectedSpan) > timeTolerance)
    error('nf_preprocess:InvalidContinuousTimeAxis', ...
        ['Continuous EEG must begin at zero and have an xmax consistent ' ...
        'with pnts and srate within half a sample.']);
end
if ~isstruct(EEG.chanlocs) || ~isfield(EEG.chanlocs, 'labels')
    error('nf_preprocess:InvalidEEG', ...
        'EEG.chanlocs must contain a labels field for every channel.');
end
end

function valid = nf_is_text(value)
valid = ischar(value) || (isstring(value) && isscalar(value));
end

function valid = nf_is_save_request(value)
valid = nf_is_logical_scalar(value) || ...
    (nf_is_text(value) && ~isempty(strtrim(char(value))));
end

function requested = nf_save_requested(value)
if nf_is_text(value)
    requested = true;
else
    requested = logical(value);
end
end

function valid = nf_is_text_or_empty(value)
valid = isempty(value) || nf_is_text(value);
end

function valid = nf_is_events(value)
valid = isempty(value) || iscell(value) || isnumeric(value) || ...
    ischar(value) || isstring(value);
end

function valid = nf_is_limits_or_empty(value)
valid = isempty(value) || nf_is_increasing_pair(value);
end

function valid = nf_is_baseline(value)
valid = (islogical(value) && isscalar(value) && ~value) || ...
    (isnumeric(value) && numel(value) == 2 && nf_is_increasing_pair(value));
end

function valid = nf_is_positive_scalar(value)
valid = isnumeric(value) && isscalar(value) && isfinite(value) && value > 0;
end

function valid = nf_is_nonnegative_scalar(value)
valid = isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0;
end

function valid = nf_is_positive_integer(value)
valid = nf_is_positive_scalar(value) && value == round(value);
end

function valid = nf_is_nonnegative_integer(value)
valid = nf_is_nonnegative_scalar(value) && value == round(value);
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

function valid = nf_is_reference(value)
valid = isempty(value) || nf_is_text(value) || ...
    (isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value >= 1 && value == round(value));
end

function valid = nf_is_labels(value)
valid = isempty(value) || ischar(value) || isstring(value) || iscellstr(value);
end

function valid = nf_is_channel_specification(value)
valid = isempty(value) || islogical(value) || isnumeric(value) || ...
    ischar(value) || isstring(value) || iscellstr(value);
end

function valid = nf_is_correlation(value)
valid = isnumeric(value) && isscalar(value) && isfinite(value) && ...
    value > 0 && value <= 1;
end

function valid = nf_is_visibility(value)
valid = nf_is_text(value) && ismember(lower(char(value)), {'on', 'off'});
end

function valid = nf_is_gedai_options(value)
valid = (isstruct(value) && isscalar(value)) || iscell(value);
end

function valid = nf_is_scalar_struct(value)
valid = isstruct(value) && isscalar(value);
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
