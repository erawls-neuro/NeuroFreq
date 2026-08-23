function [EEG, report, quality, figureHandle] = nf_preprocess(EEG, varargin)
% NF_PREPROCESS  Run reproducible flexible NeuroFreq EEG preprocessing.
%
% [EEG, REPORT, QUALITY, FIGUREHANDLE] = NF_PREPROCESS(EEG, ...)
% [EEG, REPORT, QUALITY, FIGUREHANDLE] = NF_PREPROCESS(EEG, BEHAVIOR, ...)
% [EEG, REPORT, QUALITY, FIGUREHANDLE] = ...
%     NF_PREPROCESS(EEG, EVENTS, BEHAVIOR, ...)
%
% Pipeline presets:
%
%   BDC          FASTER channels + GEDAI + ICLabel + threshold/FFT rules
%   MADE         FASTER channels + adjusted_ADJUST + threshold rules
%   PREP         official PREP channel/reference code
%   FASTER       official FASTER channel, epoch, and IC classifiers
%   HAPPE+ER     installed HAPPE channel/wavelet stages + MARA and
%                HAPPE-compatible epoch rules, composed by NeuroFreq
%   cleanrawdata official clean_rawdata channel and ASR functions
%   EEGLAB       legacy EEGLAB channel and epoch-statistics functions
%
% Common name/value inputs:
%   preset                 'BDC' (default), 'MADE', 'PREP', 'FASTER',
%                          'HAPPE+ER', 'cleanrawdata', or 'EEGLAB'
%   events                 event type(s) for pop_epoch; [] makes fixed epochs
%   behavior               NeuroFreq behavior struct array with one entry per
%                          requested event/trial, ordered by increasing event
%                          latency. It is attached immediately
%                          after final epoching as EEG.etc.behav; the legacy
%                          EEG.etc.behavior alias is kept synchronized.
%                          BEHAVIOR may instead be supplied as the first
%                          positional input after EEG, or after the legacy
%                          positional EVENTS input.
%   epochLimits            [start end] seconds when events are supplied
%   continuousEpochLength  1 second when events are empty
%   epochStage             'auto' (default), 'beforeica', or 'afterica'.
%                          For BDC with exactly one requested event type,
%                          auto uses those event epochs for ICA fitting.
%                          Other presets that use ICA retain fixed training
%                          epochs unless beforeica is explicitly requested.
%   baseline               false (default) or [start end] milliseconds
%   lowpass                45 Hz; retains the complete 20-40-Hz muscle band
%   highpass               0.3 Hz
%   notch                  60 Hz; skipped when redundant below low-pass
%   resample               250 Hz
%   eegChannels            [] auto-selects channels typed EEG/empty
%   channelMethod          faster, cleanrawdata, prep, happeer, eeglab,
%                          or none
%   fasterOptions          FASTER channel measure/z settings
%   cleanrawdataOptions    clean_flatlines/clean_channels settings
%   prepOptions            PREP removeTrend/findNoisyChannels settings
%   happeerOptions         HAPPE bad-channel stage settings; acquisition
%                          density is resolved from the original EEG count
%                          unless supplied explicitly
%   eeglabOptions          legacy pop_rejchan measures/thresholds
%   maxBadChannels         floor(10%% of selected channels), used as QC limit
%   badChannelReference    [] automatic, channel label, or index in the
%                          selected EEG montage
%   cleanHighpass          1-Hz diagnostic copy for cleanrawdata channels
%   precleanMethod         gedai, asr, prep, happeer, or none
%   asrOptions             direct clean_asr settings
%   prepPipelineOptions    parameters passed to official prepPipeline
%   happeerPrecleanOptions HAPPE wavelet settings including erpMode,
%                          softThreshold, QC passband, and QC frequencies
%   gedaiOptions           GEDAI v1.7 options; default uses the standard
%                          precomputed 10-5 database. Use referenceMatrixType
%                          'interpolated' for nonstandard coordinate montages.
%   icaMethod              '' uses the preset classifier; iclabel,
%                          adjustedadjust, adjust, mara, faster, or none
%   icaAlgorithm           'runica' (default), 'picard', or 'runamica15'
%   randomSeed             1; Picard uses deterministic identity
%                          initialization and does not consume this seed
%   minimumSamplesPerRankSquared 20 for production ICA sufficiency
%   iclabelThresholds      [] uses dominant ICLabel class, or a 7x2 matrix
%   adjustOptions          ADJUST report/behavior settings
%   maraOptions            MARA rejection controls
%   fasterICAOptions       FASTER component measure/z and EOG settings;
%                          numeric EOG indices refer to the selected raw
%                          montage and are remapped by channel label
%   picardMode             'standard' (Infomax-like) or 'ortho'
%   picardMaxIterations    500
%   picardTolerance        1e-8
%   voltageThreshold       125 microvolts
%   powerThreshold         [-100 30] dB
%   muscleRange            [20 40] Hz legacy FFT band
%   epochDetectors         composable detector names: threshold, fft,
%                          peak2peak, step, gradient, flatline, clipping,
%                          faster, eeglabstats, jointprobability, or none
%   fftBands               [] uses muscleRange; otherwise N-by-2 bands,
%                          N-by-4 [low high lowerDb upperDb] rows, or
%                          structures with frequency and power-threshold fields
%   thresholdTimes         [] uses the full final epoch
%   amplitudeOptions, fftOptions, peakToPeakOptions, stepOptions,
%                          gradientOptions, flatlineOptions,
%                          clippingOptions, fasterEpochOptions,
%                          eeglabEpochOptions, and
%                          jointProbabilityOptions expose detector settings
%                          Numeric channel fields in these structures refer
%                          to the selected raw-input montage and are remapped
%                          by label after channel/ICA preparation.
%   epochRepairOptions     per-detector and fallback mark/interpolate/reject
%                          policies for composable detections
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


nf_validate_input_eeg(EEG);
EEG = nf_ensure_event_fields(EEG);
EEG = nf_normalize_event_types(EEG);

parser = inputParser;
parser.FunctionName = 'nf_preprocess';
parser.KeepUnmatched = false;
addParameter(parser, 'preset', 'BDC', @nf_is_text);
addParameter(parser, 'events', [], @nf_is_events);
addParameter(parser, 'behavior', []);
addParameter(parser, 'epochLimits', [], @nf_is_limits_or_empty);
addParameter(parser, 'continuousEpochLength', 1, @nf_is_positive_scalar);
addParameter(parser, 'epochStage', 'auto', @nf_is_text);
addParameter(parser, 'baseline', false, @nf_is_baseline);
addParameter(parser, 'lowpass', 45, @nf_is_positive_scalar);
addParameter(parser, 'highpass', 0.3, @nf_is_nonnegative_scalar);
addParameter(parser, 'notch', 60, @nf_is_nonnegative_scalar);
addParameter(parser, 'resample', 250, @nf_is_positive_scalar);
addParameter(parser, 'eegChannels', [], @nf_is_channel_specification);
addParameter(parser, 'channelMethod', '', @nf_is_text);
addParameter(parser, 'maxBadChannels', floor(EEG.nbchan / 5), ...
    @nf_is_nonnegative_integer);
addParameter(parser, 'badChannelReference', [], @nf_is_reference);
addParameter(parser, 'fasterOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'cleanCorrelation', 0.8, @nf_is_correlation);
addParameter(parser, 'cleanHighpass', 1, @nf_is_positive_scalar);
addParameter(parser, 'cleanrawdataOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'prepOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'happeerOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'eeglabOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'precleanMethod', '', @nf_is_text);
addParameter(parser, 'gedaiOptions', struct(), @nf_is_gedai_options);
addParameter(parser, 'asrOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'prepPipelineOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'happeerPrecleanOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'icaMethod', '', @nf_is_text);
addParameter(parser, 'icaAlgorithm', 'runica', @nf_is_text);
addParameter(parser, 'aggressiveICA', false, @nf_is_logical_scalar);
addParameter(parser, 'randomSeed', 1, @nf_is_nonnegative_integer);
addParameter(parser, 'iclabelThresholds', [], @nf_is_iclabel_thresholds);
addParameter(parser, 'adjustReportFile', '', @nf_is_text_or_empty);
addParameter(parser, 'adjustOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'maraOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'fasterICAOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'icaTrainingHighpass', 1, @nf_is_positive_scalar);
addParameter(parser, 'icaTrainingEpochLength', 1, @nf_is_positive_scalar);
addParameter(parser, 'icaTrainingVoltage', 1000, @nf_is_positive_scalar);
addParameter(parser, 'icaTrainingPower', [-100 30], @nf_is_increasing_pair);
addParameter(parser, 'icaTrainingFrequencies', [20 40], @nf_is_increasing_pair);
addParameter(parser, 'icaBadChannelFraction', 0.20, @nf_is_fraction);
addParameter(parser, 'minimumTrainingEpochs', 1, @nf_is_positive_integer);
addParameter(parser, 'minimumSamplesPerRankSquared', 20, ...
    @nf_is_positive_scalar);
addParameter(parser, 'amicaMaxIterations', 2000, @nf_is_positive_integer);
addParameter(parser, 'amicaThreads', 4, @nf_is_positive_integer);
addParameter(parser, 'amicaProcesses', 1, @nf_is_positive_integer);
addParameter(parser, 'picardMode', 'standard', @nf_is_text);
addParameter(parser, 'picardMaxIterations', 500, @nf_is_positive_integer);
addParameter(parser, 'picardTolerance', 1e-8, @nf_is_positive_scalar);
addParameter(parser, 'runicaStop', 1e-7, @nf_is_positive_scalar);
addParameter(parser, 'voltageThreshold', 125, @nf_is_positive_scalar);
addParameter(parser, 'powerThreshold', [-100 30], @nf_is_increasing_pair);
addParameter(parser, 'muscleRange', [20 40], @nf_is_increasing_pair);
addParameter(parser, 'epochDetectors', {}, @nf_is_method_list);
addParameter(parser, 'fftBands', [], @nf_is_frequency_bands);
addParameter(parser, 'amplitudeOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'fftOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'peakToPeakOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'stepOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'gradientOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'flatlineOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'clippingOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'fasterEpochOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'eeglabEpochOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'jointProbabilityOptions', struct(), ...
    @nf_is_scalar_struct);
addParameter(parser, 'epochRepairOptions', struct(), @nf_is_scalar_struct);
addParameter(parser, 'thresholdTimes', [], @nf_is_limits_or_empty);
addParameter(parser, 'localInterp', true, @nf_is_logical_scalar);
addParameter(parser, 'maxLocalBad', floor(EEG.nbchan / 5), ...
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
varargin = nf_normalize_positional_inputs( ...
    varargin, parser.Parameters);
parse(parser, varargin{:});
options = parser.Results;
usingDefaults = parser.UsingDefaults;
behaviorExplicitlySupplied = ...
    ~nf_was_default(usingDefaults, 'behavior');
[inputBehavior, behaviorInfo] = nf_resolve_eeg_behavior( ...
    EEG, options.behavior, behaviorExplicitlySupplied);
options.behavior = inputBehavior;
behaviorInfo.attachmentPolicy = ...
    'deferred-until-epoch-materialization';
EEG = nf_clear_eeg_behavior(EEG);
usingDefaultMaxBad = nf_was_default(usingDefaults, 'maxBadChannels');
usingDefaultMaxLocal = nf_was_default(usingDefaults, 'maxLocalBad');

options.preset = nf_normalize_preset(options.preset);
[options, presetDefinition] = nf_apply_preset_defaults( ...
    options, usingDefaults);
options.channelMethod = nf_normalize_channel_method(options.channelMethod);
options.precleanMethod = nf_normalize_preclean_method(options.precleanMethod);
options.icaAlgorithm = lower(char(options.icaAlgorithm));
options.picardMode = lower(char(options.picardMode));
options.icaMethod = nf_resolve_ica_method(options.preset, options.icaMethod);
options.epochStageRequested = nf_normalize_epoch_stage(options.epochStage);
[options.epochStage, options.epochStageResolution] = ...
    nf_resolve_epoch_stage(options.epochStageRequested, options);
options.epochDetectors = nf_normalize_epoch_detectors( ...
    options.epochDetectors);
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
if ~isempty(options.fftBands)
    options.fftBands = nf_normalize_frequency_bands( ...
        options.fftBands, options.muscleRange);
end
options.thresholdTimes = nf_row_pair(options.thresholdTimes);
if strcmp(options.precleanMethod, 'gedai')
    options.gedaiOptions = ...
        nf_resolve_gedai_configuration(options.gedaiOptions);
end

artifactPlan = nf_resolve_artifact_plan(EEG, options);
% Retain both onCleanup objects until the job reaches its finalization path.
bundleCleanup = []; %#ok
[artifactPlan, bundleCleanup] = ...
    nf_prepare_artifact_destination(artifactPlan); %#ok
diaryCleanup = []; %#ok
diaryState = struct();
if artifactPlan.log
    diaryState = nf_start_job_log(artifactPlan, options);
    diaryCleanup = onCleanup(@() nf_restore_diary(diaryState));
end
datasetCommitted = false; %#ok
qualityFigureOwned = false; %#ok
qualityPdfOwned = false; %#ok

try
    nf_validate_options(EEG, options);
    nf_check_job_log(artifactPlan, 'option validation');
    if isempty(options.events)
        options.eventValidation = struct();
    else
        options.eventValidation = nf_validate_requested_events(EEG, options.events);
    end
    if behaviorInfo.present && ~isempty(options.events)
        nf_validate_event_behavior_count( ...
            inputBehavior, options.eventValidation.nUniqueMatchedEvents);
    end
    nf_preflight_dependencies(options);

    report = struct();
    report.schemaVersion = '4.4.0';
    report.started = datestr(now, 30); %#ok<TNOW1,DATST>
    report.preset = presetDefinition.displayName;
    report.pipelineMode = 'preset-with-explicit-overrides';
    report.presetDefinition = presetDefinition;
    report.presetDefinition.resolvedChannelMethod = options.channelMethod;
    report.presetDefinition.resolvedPrecleanMethod = options.precleanMethod;
    report.presetDefinition.resolvedClassifier = options.icaMethod;
    report.presetDefinition.resolvedIcaAlgorithm = options.icaAlgorithm;
    report.presetDefinition.resolvedEpochStage = options.epochStage;
    report.presetDefinition.resolvedEpochDetectors = ...
        options.epochDetectors;
    report.options = nf_report_options(options);
    report.input = nf_dataset_summary(EEG);
    report.steps = struct();
    report.channels = struct();
    report.gedai = struct();
    report.preclean = struct();
    report.ica = struct();
    report.epochs = struct();
    report.epochs.timing = options.epochStageResolution;
    report.provenance = struct();
    report.quality = struct('computed', false, 'plotted', false, 'error', '');
    report.quality.earlierStage = ...
        'selected raw input before filtering and artifact correction';
    if strcmp(options.icaMethod, 'none')
        report.quality.laterStage = ...
            'post-precleaning data before final epoch cleaning; ICA skipped';
        report.quality.comparisonScope = ...
            ['The dashboard spans filtering, bad-channel handling, and ' ...
            'precleaning. Final epoch actions are reported separately.'];
    else
        report.quality.laterStage = ...
            'post-ICA subtraction before final epoch cleaning';
        report.quality.comparisonScope = ...
            ['The dashboard spans filtering, bad-channel handling, ' ...
            'precleaning, and ICA subtraction. Final epoch actions are ' ...
            'reported separately.'];
    end
    report.persistence = nf_initial_persistence_report(artifactPlan);

    [EEG, channelSelection] = nf_select_eeg_channels(EEG, options.eegChannels);
    if nf_pipeline_requires_geometry(options)
        EEG = nf_normalize_channel_locations(EEG);
    else
        nf_require_unique_labels({EEG.chanlocs.labels});
    end
    report.channels.selection = channelSelection;
    if strcmp(options.channelMethod, 'happeer')
        hasHappeDensity = isfield(options.happeerOptions, 'lowDensity') && ...
            ~isempty(options.happeerOptions.lowDensity);
        hasHappeParameterDensity = ...
            isfield(options.happeerOptions, 'params') && ...
            isstruct(options.happeerOptions.params) && ...
            isfield(options.happeerOptions.params, 'lowDensity') && ...
            ~isempty(options.happeerOptions.params.lowDensity);
        if ~hasHappeDensity && ~hasHappeParameterDensity
            options.happeerOptions.lowDensity = ...
                channelSelection.nAcquisitionEEGChannels <= 32;
            report.channels.happeerDensityResolution = struct( ...
                'source', 'input chanlocs.type EEG-or-empty count', ...
                'channelCount', ...
                channelSelection.nAcquisitionEEGChannels, ...
                'lowDensity', options.happeerOptions.lowDensity, ...
                'inferredByNeuroFreq', true);
        else
            report.channels.happeerDensityResolution = struct( ...
                'source', 'explicit happeerOptions', ...
                'channelCount', ...
                channelSelection.nAcquisitionEEGChannels, ...
                'lowDensity', [], ...
                'inferredByNeuroFreq', false);
        end
    end
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
    epochRepairRequirements = nf_epoch_repair_requirements(options);
    if ~strcmp(options.channelMethod, 'none') && ...
            options.maxBadChannels >= EEG.nbchan - 2
        error('nf_preprocess:InvalidBadChannelLimit', ...
            'maxBadChannels must leave at least three selected EEG channels.');
    end
    if epochRepairRequirements.channelRepairCapable && ...
            options.maxLocalBad >= EEG.nbchan
        error('nf_preprocess:InvalidLocalChannelLimit', ...
            'maxLocalBad must be below the selected EEG channel count.');
    end
    if epochRepairRequirements.requiresInterpolation && ...
            options.maxLocalBad > EEG.nbchan - 3
        error('nf_preprocess:InvalidLocalChannelLimit', ...
            'maxLocalBad must leave at least three interpolation donor channels.');
    end

    originalMontage = EEG.chanlocs;
    if ~isfield(EEG, 'etc') || isempty(EEG.etc)
        EEG.etc = struct();
    end
    EEG.etc.ogchan = originalMontage;
    if options.qualityCompute
        EEG_preclean = EEG;
    else
        EEG_preclean = [];
    end
    useUnfilteredChannelDiagnostic = ...
        nf_channel_detection_uses_unfiltered_input(options.channelMethod);
    if useUnfilteredChannelDiagnostic
        EEG_channelDiagnostic = EEG;
    else
        EEG_channelDiagnostic = [];
    end
    precleanRanEarly = false;
    precleanInfo = struct();
    gedaiInfo = struct();
    gedaiInfo.applied = false;
    if strcmp(options.precleanMethod, 'prep')
        [EEG, precleanInfo] = nf_run_prep_pipeline( ...
            EEG, options.prepPipelineOptions);
        precleanRanEarly = true;
        nf_check_job_log(artifactPlan, 'PREP pipeline');
    end

    deferHappeErpFilter = strcmp(options.precleanMethod, 'happeer') && ...
        nf_happeer_erp_mode( ...
        options.happeerPrecleanOptions, options.events);
    if deferHappeErpFilter
        preWaveletNotch = options.notch;
        skippedPreWaveletNotch = false;
        if preWaveletNotch > 0 && ...
                preWaveletNotch + 2 >= EEG.srate / 2
            preWaveletNotch = 0;
            skippedPreWaveletNotch = true;
        end
        [EEG, filterInfo] = nf_filter( ...
            EEG, 0, 0, preWaveletNotch, options.resample);
        filterInfo.requested.notchHz = options.notch;
        if skippedPreWaveletNotch
            filterInfo.skipped{end + 1} = ...
                ['Pre-wavelet notch skipped because its stop band was ' ...
                'not below the input Nyquist frequency.'];
        end
        filterInfo.deferredAnalysisPassband = true;
        filterInfo.order = ...
            ['HAPPE ERP composition: line-noise notch and resampling before ' ...
            'the HAPPE channel/wavelet stages; analysis passband afterward.'];
    else
        [EEG, filterInfo] = nf_filter( ...
            EEG, options.lowpass, options.highpass, ...
            options.notch, options.resample);
        filterInfo.deferredAnalysisPassband = false;
    end
    report.steps.filter = filterInfo;
    nf_check_job_log(artifactPlan, 'filtering');

    if useUnfilteredChannelDiagnostic
        channelInput = EEG_channelDiagnostic;
    else
        channelInput = EEG;
    end
    [channelOutput, channelInfo] = nf_badchans( ...
        channelInput, options.maxBadChannels, ...
        false, options.channelMethod, ...
        'reference', options.badChannelReference, ...
        'fasterOptions', options.fasterOptions, ...
        'cleanCorrelation', options.cleanCorrelation, ...
        'cleanHighpass', options.cleanHighpass, ...
        'cleanrawdataOptions', options.cleanrawdataOptions, ...
        'prepOptions', options.prepOptions, ...
        'happeerOptions', options.happeerOptions, ...
        'eeglabOptions', options.eeglabOptions, ...
        'interpolationMethod', options.interpolationMethod);
    if useUnfilteredChannelDiagnostic
        [EEG, channelTransfer] = nf_apply_channel_detection( ...
            EEG, channelOutput, channelInfo);
        channelInfo.inputStage = ...
            'selected raw input before analysis low-pass and resampling';
        channelInfo.transferToAnalysisData = channelTransfer;
    else
        EEG = channelOutput;
        if deferHappeErpFilter && strcmp(options.channelMethod, 'happeer')
            channelInfo.inputStage = ...
                ['pre-wavelet line-noise/resampled data before the ERP ' ...
                'analysis passband'];
        else
            channelInfo.inputStage = 'analysis-filtered working data';
        end
    end
    if strcmp(options.preset, 'faster') && ...
            strcmp(options.channelMethod, 'faster') && ...
            isfield(channelInfo, 'provenance') && ...
            isfield(channelInfo.provenance, 'contractLevel') && ...
            ~strcmp(channelInfo.provenance.contractLevel, ...
            'vendor-exact-stage')
        error('nf_preprocess:FASTERPresetReferenceRequired', ...
            ['The FASTER comparison preset requires vendor-exact ' ...
            'channel_properties. NeuroFreq could not recover a recording ' ...
            'reference and therefore refused its compatible no-reference ' ...
            'implementation. Supply badChannelReference or repair the ' ...
            'reference metadata.']);
    end
    report.channels.detection = channelInfo;
    report.steps.badChannels.applied = ~strcmp(options.channelMethod, 'none');
    report.steps.badChannels.method = options.channelMethod;
    nf_check_job_log(artifactPlan, 'bad-channel detection');
    clear channelInput channelOutput EEG_channelDiagnostic;

    if precleanRanEarly
        gedaiInfo.applied = false;
    elseif strcmp(options.precleanMethod, 'gedai')
        gedaiConfiguration = nf_subset_gedai_reference( ...
            options.gedaiOptions, originalMontage, EEG.chanlocs);
        if ischar(gedaiConfiguration.referenceMatrixType) && ...
                strcmp(gedaiConfiguration.referenceMatrixType, 'precomputed')
            nf_validate_gedai_standard_labels(EEG.chanlocs);
        end
        [EEG, gedaiInfo] = nf_run_gedai(EEG, gedaiConfiguration);
        precleanInfo = gedaiInfo;
    elseif strcmp(options.precleanMethod, 'asr')
        [EEG, precleanInfo] = nf_run_asr(EEG, options.asrOptions);
        gedaiInfo = struct();
        gedaiInfo.applied = false;
    elseif strcmp(options.precleanMethod, 'happeer')
        [EEG, precleanInfo] = nf_run_happeer_preclean( ...
            EEG, options.happeerPrecleanOptions, options);
        if deferHappeErpFilter
            [EEG, filterInfo] = nf_apply_deferred_happe_filter( ...
                EEG, filterInfo, options);
            report.steps.filter = filterInfo;
            nf_check_job_log( ...
                artifactPlan, 'post-wavelet ERP analysis filtering');
        end
        gedaiInfo = struct();
        gedaiInfo.applied = false;
    else
        gedaiInfo = struct();
        gedaiInfo.applied = false;
        precleanInfo = struct();
        precleanInfo.applied = false;
        precleanInfo.method = 'none';
    end
    EEG = eeg_checkset(EEG);
    report.gedai = gedaiInfo;
    report.preclean = precleanInfo;
    report.steps.preclean.applied = precleanInfo.applied;
    report.steps.preclean.method = options.precleanMethod;
    report.steps.gedai.applied = gedaiInfo.applied;
    nf_check_job_log(artifactPlan, 'precleaning');

    [options.fasterICAOptions, fasterIcaOptionRemapping] = ...
        nf_remap_faster_ica_options( ...
        options.fasterICAOptions, options.icaMethod, ...
        originalMontage, EEG.chanlocs);
    report.channels.fasterICAOptionRemapping = ...
        fasterIcaOptionRemapping;
    if strcmp(options.epochStage, 'beforeica')
        icaTrainingEvents = options.eventValidation.popEpochEvents;
        icaTrainingEpochLimits = options.epochLimits;
    else
        icaTrainingEvents = {};
        icaTrainingEpochLimits = [];
    end
    cleanicBehaviorArguments = {};
    if behaviorInfo.present
        cleanicBehaviorArguments = {'behavior', inputBehavior};
    end
    cleanicArguments = [ ...
        {options.icaMethod}, ...
        {options.aggressiveICA}, ...
        cleanicBehaviorArguments(:)', ...
        {'algorithm'}, {options.icaAlgorithm}, ...
        {'randomSeed'}, {options.randomSeed}, ...
        {'trainingHighpass'}, {options.icaTrainingHighpass}, ...
        {'trainingEpochLength'}, {options.icaTrainingEpochLength}, ...
        {'trainingVoltage'}, {options.icaTrainingVoltage}, ...
        {'trainingPower'}, {options.icaTrainingPower}, ...
        {'trainingFrequencies'}, {options.icaTrainingFrequencies}, ...
        {'badChannelFraction'}, {options.icaBadChannelFraction}, ...
        {'minimumTrainingEpochs'}, {options.minimumTrainingEpochs}, ...
        {'minimumSamplesPerRankSquared'}, ...
        {options.minimumSamplesPerRankSquared}, ...
        {'trainingEvents'}, {icaTrainingEvents}, ...
        {'trainingEpochLimits'}, {icaTrainingEpochLimits}, ...
        {'iclabelThresholds'}, {options.iclabelThresholds}, ...
        {'adjustReportFile'}, {options.adjustReportFile}, ...
        {'adjustOptions'}, {options.adjustOptions}, ...
        {'maraOptions'}, {options.maraOptions}, ...
        {'fasterOptions'}, {options.fasterICAOptions}, ...
        {'amicaMaxIterations'}, {options.amicaMaxIterations}, ...
        {'amicaThreads'}, {options.amicaThreads}, ...
        {'amicaProcesses'}, {options.amicaProcesses}, ...
        {'picardMode'}, {options.picardMode}, ...
        {'picardMaxIterations'}, {options.picardMaxIterations}, ...
        {'picardTolerance'}, {options.picardTolerance}, ...
        {'runicaStop'}, {options.runicaStop}];
    if options.qualityCompute
        [EEG, icaInfo, EEG_icaModel] = nf_cleanic( ...
            EEG, ...
            cleanicArguments{:});
    else
        [EEG, icaInfo] = nf_cleanic( ...
            EEG, ...
            cleanicArguments{:});
        EEG_icaModel = [];
    end
    EEG = nf_clear_eeg_behavior(EEG);
    if isfield(icaInfo, 'behavior')
        icaInfo.behavior.explicitlySuppliedToCleanic = ...
            icaInfo.behavior.explicitlySupplied;
        icaInfo.behavior.explicitlySupplied = ...
            behaviorInfo.explicitlySupplied;
        icaInfo.behavior.preprocessSource = behaviorInfo.source;
        icaInfo.behavior.forwardedByPreprocess = behaviorInfo.present;
        if isfield(EEG, 'etc') && isstruct(EEG.etc) && ...
                isfield(EEG.etc, 'nf_cleanic') && ...
                isstruct(EEG.etc.nf_cleanic)
            EEG.etc.nf_cleanic.behavior = icaInfo.behavior;
        end
    end
    EEG.icaact = [];
    nf_validate_input_eeg(EEG);
    if options.qualityCompute
        EEG_postclean = EEG;
    else
        EEG_postclean = [];
    end
    report.ica = icaInfo;
    report.steps.ica.applied = ~strcmp(options.icaMethod, 'none');
    report.steps.ica.method = options.icaMethod;
    report.steps.ica.algorithm = icaInfo.algorithm;
    report.steps.ica.epochStage = options.epochStage;
    if isfield(icaInfo.training, 'source')
        report.steps.ica.trainingEpochSource = icaInfo.training.source;
    else
        report.steps.ica.trainingEpochSource = 'none';
    end
    if isfield(icaInfo, 'classification') && ...
            isfield(icaInfo.classification, 'dataScope')
        report.steps.ica.classifierDataScope = ...
            icaInfo.classification.dataScope;
    end
    if isfield(icaInfo, 'componentSubtractionDataScope')
        report.steps.ica.componentSubtractionDataScope = ...
            icaInfo.componentSubtractionDataScope;
    end
    nf_check_job_log(artifactPlan, 'ICA cleaning');

    [EEG, epochInfo] = nf_make_final_epochs( ...
        EEG, options, behaviorInfo);
    if strcmp(options.epochStage, 'beforeica')
        trainingIndices = icaInfo.training.candidateEventIndices;
        finalIndices = epochInfo.candidateEventIndices;
        trainingLatencies = ...
            icaInfo.training.candidateEventLatenciesSamples;
        finalLatencies = epochInfo.candidateEventLatenciesSamples;
        epochInfo.matchesIcaCandidateEventSelection = ...
            isequal(trainingIndices, finalIndices) && ...
            isequal(trainingLatencies, finalLatencies);
        epochInfo.icaTrainingCandidateEventIndices = trainingIndices;
        epochInfo.icaTrainingRetainedEventIndices = ...
            icaInfo.training.retainedEventIndices;
        if ~epochInfo.matchesIcaCandidateEventSelection
            error('nf_preprocess:ICAEventEpochMismatch', ...
                ['Final epoching did not reproduce the complete candidate ' ...
                'event selection used to construct the ICA-fitting epochs. ' ...
                'No output was produced because event identity changed ' ...
                'unexpectedly between ICA preparation and final epoching.']);
        end
    else
        epochInfo.matchesIcaCandidateEventSelection = [];
        epochInfo.icaTrainingCandidateEventIndices = [];
        epochInfo.icaTrainingRetainedEventIndices = [];
    end
    epochInfo.icaTrainingEpochStage = options.epochStage;
    if strcmp(icaInfo.algorithm, 'none')
        epochInfo.finalMaterializationStage = ...
            ['after filtering and precleaning; ICA fitting and subtraction ' ...
            'were skipped'];
    else
        epochInfo.finalMaterializationStage = ...
            'after ICA component subtraction';
    end
    report.epochs.creation = epochInfo;
    report.steps.epoch = epochInfo;
    nf_check_job_log(artifactPlan, 'epoch creation');

    if isnumeric(options.baseline) && numel(options.baseline) == 2
        [baselineWindow, baselineAdjustment] = nf_fit_epoch_window( ...
            options.baseline, [EEG.xmin EEG.xmax] * 1000, ...
            1000 / EEG.srate, 'baseline');
        EEG = pop_rmbase(EEG, baselineWindow);
        [EEG, baselineDecomposition] = ...
            nf_detach_ica_for_epoching(EEG);
        EEG = eeg_checkset(EEG);
        EEG = nf_restore_ica_after_epoching( ...
            EEG, baselineDecomposition);
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
    [options, epochChannelRemapping] = ...
        nf_remap_epoch_channel_options( ...
        options, originalMontage, EEG.chanlocs);
    report.epochs.channelOptionRemapping = epochChannelRemapping;
    [EEG, ~, thresholdInfo] = nf_thresh(EEG, ...
        options.voltageThreshold, options.powerThreshold, ...
        options.muscleRange, thresholdTimes, options.localInterp, ...
        options.maxLocalBad, options.frontalChannels, ...
        'interpolationMethod', options.interpolationMethod, ...
        'detectors', options.epochDetectors, ...
        'fftBands', options.fftBands, ...
        'amplitudeOptions', options.amplitudeOptions, ...
        'fftOptions', options.fftOptions, ...
        'peak2peakOptions', options.peakToPeakOptions, ...
        'stepOptions', options.stepOptions, ...
        'gradientOptions', options.gradientOptions, ...
        'flatlineOptions', options.flatlineOptions, ...
        'clippingOptions', options.clippingOptions, ...
        'fasterOptions', options.fasterEpochOptions, ...
        'eeglabOptions', options.eeglabEpochOptions, ...
        'jointProbabilityOptions', options.jointProbabilityOptions, ...
        'repairOptions', options.epochRepairOptions);
    report.epochs.threshold = thresholdInfo;
    report.epochs.threshold.requestedTimesSeconds = options.thresholdTimes;
    report.epochs.threshold.appliedTimesSeconds = thresholdTimes;
    report.epochs.threshold.endpointAdjusted = thresholdAdjustment;
    report.steps.threshold.applied = ...
        ~all(strcmp(options.epochDetectors, 'none'));
    report.steps.threshold.nRejected = thresholdInfo.nRejected;
    report.steps.threshold.nLocallyRepaired = ...
        thresholdInfo.nLocallyRepairedEpochs;
    nf_check_job_log(artifactPlan, 'epoch thresholding');

    [EEG, interpolationInfo] = nf_finalize_montage( ...
        EEG, originalMontage, options, channelInfo);
    report.channels.finalization = interpolationInfo;
    report.steps.finalization = interpolationInfo;
    report.provenance = nf_collect_pipeline_provenance( ...
        channelInfo, precleanInfo, icaInfo, thresholdInfo);
    [EEG, finalDecomposition] = nf_detach_ica_for_epoching(EEG);
    EEG = eeg_checkset(EEG);
    EEG = nf_restore_ica_after_epoching(EEG, finalDecomposition);
    nf_check_job_log(artifactPlan, 'montage finalization');

    quality = struct();
    figureHandle = [];
    if options.qualityCompute
        try
            [quality, figureHandle] = nf_eegquality(EEG_preclean, ...
                EEG_postclean, 'final', EEG, 'report', report, ...
                'icaModel', EEG_icaModel, ...
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
        qualityFigureOwned = report.persistence.quality.fig.saved; %#ok
        qualityPdfOwned = report.persistence.quality.pdf.saved; %#ok
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
    if ~isfield(EEG.etc, 'neurofreq') || ...
            ~isstruct(EEG.etc.neurofreq)
        EEG.etc.neurofreq = struct();
    end
    EEG.etc.neurofreq.provenance = report.provenance;
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
        datasetCommitted = true; %#ok
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
    fprintf(2, '[nf_preprocess WARNING]:  %s: %s\n', ...
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
        fprintf(2, ['[nf_preprocess]: Dataset save failed; removing partial dataset files ' ...
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
fprintf(2, ['[nf_preprocess]: Preprocessing failed before dataset commit; removing QC files ' ...
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
        fprintf('\n[nf_preprocess]:  Redirecting active diary to job log: %s\n', ...
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
    reportedOptions = nf_report_options(options); %#ok <I'm gonna use it later bro I swear bro>
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
        datestr(now, 31), ... 
        plan.sourceStem, ...
        plan.outputDirectory, ...
        optionsText); %#ok<TNOW1,DATST>
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
        fprintf(['\n[nf_preprocess]:  Resumed diary after nf_preprocess job ' ...
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
        datestr(now, 31), ... 
        diagnostic, ...
        plan.logPath); %#ok<TNOW1,DATST>
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
contract = nf_vendor_function_contract('GEDAI', -12, 10);
contract.level = 'vendor-exact-stage';
contract.completePipeline = true;
contract.provider = 'GEDAI';
contract.statement = ...
    ['The uniquely resolved installed GEDAI entry point made the ' ...
    'artifact-modeling decision; no NeuroFreq substitute was used.'];
if ischar(configuration.referenceMatrixType) && ...
        strcmp(configuration.referenceMatrixType, 'precomputed')
    gedaiDirectory = fileparts(contract.path);
    databasePath = fullfile(gedaiDirectory, 'auxiliaries', ...
        'fsavLEADFIELD_4_GEDAI.mat');
    contract.referenceDatabase.path = databasePath;
    contract.referenceDatabase.sha256 = nf_file_sha256(databasePath);
    if isempty(contract.referenceDatabase.sha256)
        error('nf_preprocess:GEDAIDatabaseHashFailed', ...
            ['Could not calculate the SHA-256 identity of GEDAI''s ' ...
            'precomputed lead-field database at %s.'], databasePath);
    end
end
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
info.method = 'gedai';
info.contract = contract;
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

function [EEG, info] = nf_run_asr(EEG, supplied)
configuration = nf_resolve_asr_options(supplied);
contract = nf_vendor_function_contract('clean_asr', 11, 1);
normalizedContractPath = lower(strrep(contract.path, '\', '/'));
if ~contains(normalizedContractPath, 'clean_rawdata') && ...
        ~contains(normalizedContractPath, 'cleanrawdata')
    error('nf_preprocess:InvalidCleanRawDataProvider', ...
        ['clean_asr did not resolve inside an identifiable clean_rawdata ' ...
        'installation. Resolved path: %s'], contract.path);
end
if isfield(EEG, 'icaweights') && ~isempty(EEG.icaweights)
    error('nf_preprocess:StaleICAForASR', ...
        ['ASR changes the channel data. Remove any existing ICA ' ...
        'decomposition before nf_preprocess.']);
end
before = double(EEG.data);
EEG = clean_asr(EEG, configuration.cutoff, ...
    configuration.windowLength, configuration.stepSize, ...
    configuration.maxDimensions, ...
    configuration.referenceMaxBadChannels, ...
    configuration.referenceTolerances, ...
    configuration.referenceWindowLength, configuration.useGPU, ...
    configuration.useRiemannian, configuration.maxMemory);
EEG = eeg_checkset(EEG);
after = double(EEG.data);
if ~isequal(size(before), size(after))
    error('nf_preprocess:ASRDimensionChange', ...
        'clean_asr unexpectedly changed the data dimensions.');
end
difference = before - after;

info = struct();
info.applied = true;
info.method = 'asr';
info.configuration = configuration;
info.contract = contract;
info.contract.level = 'vendor-exact-stage';
info.contract.completePipeline = false;
info.contract.statement = ...
    ['The installed clean_rawdata clean_asr function was called directly. ' ...
    'No native substitute or silent fallback was used.'];
info.changedSampleFraction = nnz(difference ~= 0) / numel(difference);
info.removedRMSMicrovolts = sqrt(mean(difference(:) .^ 2));
end

function configuration = nf_resolve_asr_options(supplied)
configuration = struct();
configuration.cutoff = 20;
configuration.windowLength = [];
configuration.stepSize = [];
configuration.maxDimensions = 0.66;
configuration.referenceMaxBadChannels = 0.075;
configuration.referenceTolerances = [-3.5 5.5];
configuration.referenceWindowLength = 1;
configuration.useGPU = false;
configuration.useRiemannian = false;
configuration.maxMemory = 64;
configuration = nf_merge_option_struct( ...
    configuration, supplied, 'asrOptions');
nf_require_positive(configuration.cutoff, 'asrOptions.cutoff');
if ~isempty(configuration.windowLength)
    nf_require_positive( ...
        configuration.windowLength, 'asrOptions.windowLength');
end
if ~isempty(configuration.stepSize)
    nf_require_positive(configuration.stepSize, 'asrOptions.stepSize');
end
nf_require_positive( ...
    configuration.maxDimensions, 'asrOptions.maxDimensions');
if isnumeric(configuration.referenceMaxBadChannels)
    nf_require_positive( ...
        configuration.referenceMaxBadChannels, ...
        'asrOptions.referenceMaxBadChannels');
elseif nf_is_text(configuration.referenceMaxBadChannels) && ...
        strcmpi(strtrim(char(configuration.referenceMaxBadChannels)), 'off')
    configuration.referenceMaxBadChannels = 'off';
else
    error('nf_preprocess:InvalidASROptions', ...
        'asrOptions.referenceMaxBadChannels must be positive or off.');
end
if isnumeric(configuration.referenceTolerances) && ...
        (~isreal(configuration.referenceTolerances) || ...
        numel(configuration.referenceTolerances) ~= 2 || ...
        any(~isfinite(configuration.referenceTolerances)))
    error('nf_preprocess:InvalidASROptions', ...
        'asrOptions.referenceTolerances must be two finite values or off.');
elseif nf_is_text(configuration.referenceTolerances) && ...
        strcmpi(strtrim(char(configuration.referenceTolerances)), 'off')
    configuration.referenceTolerances = 'off';
elseif ~isnumeric(configuration.referenceTolerances)
    error('nf_preprocess:InvalidASROptions', ...
        'asrOptions.referenceTolerances must be two finite values or off.');
end
if isnumeric(configuration.referenceWindowLength)
    nf_require_positive( ...
        configuration.referenceWindowLength, ...
        'asrOptions.referenceWindowLength');
elseif nf_is_text(configuration.referenceWindowLength) && ...
        strcmpi(strtrim(char(configuration.referenceWindowLength)), 'off')
    configuration.referenceWindowLength = 'off';
else
    error('nf_preprocess:InvalidASROptions', ...
        'asrOptions.referenceWindowLength must be positive or off.');
end
if ~nf_is_logical_scalar(configuration.useGPU) || ...
        ~nf_is_logical_scalar(configuration.useRiemannian)
    error('nf_preprocess:InvalidASROptions', ...
        'ASR useGPU and useRiemannian settings must be logical scalars.');
end
configuration.useGPU = logical(configuration.useGPU);
configuration.useRiemannian = logical(configuration.useRiemannian);
nf_require_positive(configuration.maxMemory, 'asrOptions.maxMemory');
end

function [EEG, info] = nf_run_prep_pipeline(EEG, supplied)
pipelineContract = nf_vendor_function_contract('prepPipeline', 2, 3);
versionContract = nf_vendor_function_contract('getPrepVersion', 0, 3);
normalizedPath = lower(strrep(pipelineContract.path, '\', '/'));
if ~contains(normalizedPath, '/preppipeline/preppipeline.m')
    error('nf_preprocess:InvalidPREPProvider', ...
        ['prepPipeline did not resolve to the expected PREP directory. ' ...
        'Resolved path: %s'], pipelineContract.path);
end
prepRoot = fileparts(pipelineContract.path);
normalizedPrepRoot = lower(strrep(prepRoot, '\', '/'));
normalizedVersionPath = lower(strrep(versionContract.path, '\', '/'));
if ~startsWith(normalizedVersionPath, [normalizedPrepRoot '/'])
    error('nf_preprocess:MixedPREPInstall', ...
        ['getPrepVersion did not resolve inside the same PREP package as ' ...
        'prepPipeline. Remove mixed or shadowed PREP installations.']);
end
if EEG.trials ~= 1 || ~ismatrix(EEG.data)
    error('nf_preprocess:PREPRequiresContinuousData', ...
        'The complete PREP pipeline requires continuous two-dimensional EEG.');
end
if isfield(EEG, 'icaweights') && ~isempty(EEG.icaweights)
    error('nf_preprocess:StaleICAForPREP', ...
        ['PREP changes the channel data. Remove any existing ICA ' ...
        'decomposition before nf_preprocess.']);
end
params = supplied;
if ~isfield(params, 'name') || isempty(params.name)
    if isfield(EEG, 'setname') && ~isempty(EEG.setname)
        params.name = EEG.setname;
    elseif isfield(EEG, 'filename') && ~isempty(EEG.filename)
        params.name = EEG.filename;
    else
        params.name = 'NeuroFreq_PREP_input';
    end
end
identity = nf_dataset_identity(EEG);
warningState = warning;
pathState = path;
workingDirectory = pwd;
randomState = rng;
environmentCleanup = onCleanup(@() nf_restore_matlab_environment( ...
    warningState, pathState, workingDirectory, randomState));
[EEG, resolvedParams, computationTimes] = prepPipeline(EEG, params);
clear environmentCleanup
EEG = nf_restore_dataset_identity(EEG, identity);
if ~isfield(EEG, 'etc') || ...
        ~isfield(EEG.etc, 'noiseDetection') || ...
        ~isstruct(EEG.etc.noiseDetection) || ...
        ~isfield(EEG.etc.noiseDetection, 'errors')
    error('nf_preprocess:PREPIncomplete', ...
        ['prepPipeline returned without the required ' ...
        'EEG.etc.noiseDetection.errors report.']);
end
errors = EEG.etc.noiseDetection.errors;
status = '';
if isstruct(errors) && isfield(errors, 'status')
    status = char(string(errors.status));
end
if ~strcmpi(strtrim(status), 'good')
    diagnostic = strtrim(evalc('disp(errors)'));
    error('nf_preprocess:PREPFailed', ...
        'prepPipeline reported status %s:\n%s', status, diagnostic);
end
EEG = eeg_checkset(EEG);
try
    prepVersion = getPrepVersion();
catch
    prepVersion = '';
end

info = struct();
info.applied = true;
info.method = 'prep';
info.resolvedParameters = resolvedParams;
info.computationTimesSeconds = computationTimes;
info.nativeNoiseDetectionLocation = 'EEG.etc.noiseDetection';
info.nativeNoiseDetectionSummary = nf_prep_noise_summary( ...
    EEG.etc.noiseDetection);
info.contract = pipelineContract;
info.contract.level = 'vendor-exact-pipeline';
info.contract.completePipeline = true;
info.contract.provider = 'PREP / VisLab EEG-Clean-Tools';
info.contract.release = prepVersion;
info.contract.versionFunction = versionContract;
info.contract.packageRoot = prepRoot;
info.contract.coLocatedVersionFunctionVerified = true;
info.contract.upstreamReleaseVerified = false;
info.contract.statement = ...
    ['The uniquely resolved installed prepPipeline entry point performed PREP ' ...
    'detrending for detection, line-noise removal, robust referencing, ' ...
    'noisy-channel detection, and interpolation. NeuroFreq then continues ' ...
    'with only the explicitly resolved downstream stages.'];
end

function summary = nf_prep_noise_summary(noiseDetection)
summary = struct();
summary.version = '';
if isfield(noiseDetection, 'version')
    summary.version = noiseDetection.version;
end
summary.errors = noiseDetection.errors;
numericFields = {'interpolatedChannelNumbers', ...
    'removedChannelNumbers', 'stillNoisyChannelNumbers'};
for index = 1:numel(numericFields)
    fieldName = numericFields{index};
    if isfield(noiseDetection, fieldName)
        summary.(fieldName) = noiseDetection.(fieldName);
    else
        summary.(fieldName) = [];
    end
end
end

function [EEG, info] = nf_run_happeer_preclean( ...
    EEG, supplied, pipelineOptions)
configuration = nf_resolve_happeer_preclean_options( ...
    supplied, pipelineOptions);
contract = nf_vendor_function_contract('happe_wavThresh', 5, 3);
normalizedContractPath = lower(strrep(contract.path, '\', '/'));
if ~contains(normalizedContractPath, 'happe')
    error('nf_preprocess:InvalidHAPPEProvider', ...
        ['happe_wavThresh did not resolve inside an identifiable HAPPE ' ...
        'installation. Resolved path: %s'], contract.path);
end
dependencyNames = {'assessPipelineStep', 'calcSNR_PSNR', ...
    'wdenoise', 'mscohere', 'pop_eegfiltnew'};
dependencyContracts = cell(1, numel(dependencyNames));
for index = 1:numel(dependencyNames)
    dependencyContracts{index} = nf_vendor_function_contract( ...
        dependencyNames{index}, [], []);
end
dependencies = [dependencyContracts{:}];
happeStageRoot = fileparts(contract.path);
for index = 1:2
    if ~strcmpi(fileparts(dependencies(index).path), happeStageRoot)
        error('nf_preprocess:MixedHAPPEInstall', ...
            ['%s did not resolve beside happe_wavThresh. Remove mixed or ' ...
            'shadowed HAPPE installations.'], dependencies(index).function);
    end
end
if EEG.trials ~= 1 || ~ismatrix(EEG.data)
    error('nf_preprocess:HAPPERequiresContinuousData', ...
        ['The official HAPPE wavelet helper flattens input and does not ' ...
        'restore epoch dimensions. Continuous two-dimensional EEG is required.']);
end
if isfield(EEG, 'icaweights') && ~isempty(EEG.icaweights)
    error('nf_preprocess:StaleICAForHAPPE', ...
        ['HAPPE wavelet correction changes the channel data but does not ' ...
        'clear ICA fields. Remove the existing decomposition first.']);
end
if EEG.pnts < 1000
    error('nf_preprocess:HAPPERecordingTooShort', ...
        ['The official HAPPE QC calculation uses a 1000-sample coherence ' ...
        'window; this recording is too short.']);
end
if any(configuration.qcFrequencies >= EEG.srate / 2) || ...
        any(configuration.qcFrequencies > configuration.lowpass)
    error('nf_preprocess:HAPPEQCFrequency', ...
        ['Every HAPPE QC frequency must be below the data Nyquist ' ...
        'frequency and inside the configured retained passband.']);
end

params = struct();
params.wavelet = struct();
params.wavelet.legacy = false;
params.wavelet.softThresh = logical(configuration.softThreshold);
params.paradigm = struct();
params.paradigm.ERP = struct();
params.paradigm.ERP.on = configuration.erpMode;
params.filt = struct();
params.filt.highpass = configuration.highpass;
params.filt.lowpass = configuration.lowpass;
params.QCfreqs = configuration.qcFrequencies;
wavMeans = [];
dataQC = cell(1, 6);
currFile = 1;
identity = nf_dataset_identity(EEG);
before = double(EEG.data);
warningState = warning;
pathState = path;
workingDirectory = pwd;
randomState = rng;
environmentCleanup = onCleanup(@() nf_restore_matlab_environment( ...
    warningState, pathState, workingDirectory, randomState));
[EEG, wavMeans, dataQC] = happe_wavThresh( ...
    EEG, params, wavMeans, dataQC, currFile);
clear environmentCleanup
EEG = nf_restore_dataset_identity(EEG, identity);
EEG = eeg_checkset(EEG);
after = double(EEG.data);
if ~isequal(size(before), size(after))
    error('nf_preprocess:HAPPEDimensionChange', ...
        'happe_wavThresh unexpectedly changed the data dimensions.');
end
difference = before - after;

info = struct();
info.applied = true;
info.method = 'happeer';
info.configuration = configuration;
info.contract = contract;
info.contract.level = 'vendor-exact-stage';
info.contract.completePipeline = false;
info.contract.dependencies = dependencies;
info.contract.packageStageRoot = happeStageRoot;
info.contract.coLocatedHappeHelpersVerified = true;
info.contract.upstreamReleaseVerified = false;
info.contract.statement = ...
    ['The uniquely resolved installed happe_wavThresh function was called ' ...
    'directly with its co-located HAPPE QC helpers; exact file hashes are ' ...
    'recorded. This proves stage-code identity, not upstream release ' ...
    'authenticity or complete ordered HAPPE+ER pipeline equivalence.'];
info.waveletMetrics = wavMeans;
info.waveletMetricNames = [{'RMSE', 'MAE', 'SNR', 'peakSNR', ...
    'allDataPearsonR'} nf_happe_coherence_names( ...
    configuration.qcFrequencies)];
info.rawDataQC = dataQC;
if size(dataQC, 1) >= 1 && size(dataQC, 2) >= 6
    info.retainedVariancePercent = dataQC{1, 6};
else
    info.retainedVariancePercent = NaN;
end
info.changedSampleFraction = nnz(difference ~= 0) / numel(difference);
info.removedRMSMicrovolts = sqrt(mean(difference(:) .^ 2));
end

function configuration = nf_resolve_happeer_preclean_options( ...
    supplied, pipelineOptions)
configuration = struct();
configuration.softThreshold = false;
configuration.erpMode = [];
configuration.highpass = pipelineOptions.highpass;
configuration.lowpass = pipelineOptions.lowpass;
configuration.qcFrequencies = [0.5 1 2 5 8 12 20 30 40 45 70];
qcFrequenciesSupplied = isfield(supplied, 'qcFrequencies');
configuration = nf_merge_option_struct( ...
    configuration, supplied, 'happeerPrecleanOptions');
if ~nf_is_logical_scalar(configuration.softThreshold)
    error('nf_preprocess:InvalidHAPPEOptions', ...
        'happeerPrecleanOptions.softThreshold must be a logical scalar.');
end
configuration.softThreshold = logical(configuration.softThreshold);
if isempty(configuration.erpMode)
    configuration.erpMode = ~isempty(pipelineOptions.events);
    configuration.erpModeSource = ...
        'inferred from whether event-locked epochs were requested';
elseif nf_is_logical_scalar(configuration.erpMode)
    configuration.erpMode = logical(configuration.erpMode);
    configuration.erpModeSource = 'explicit happeerPrecleanOptions.erpMode';
else
    error('nf_preprocess:InvalidHAPPEOptions', ...
        'happeerPrecleanOptions.erpMode must be empty or binary.');
end
nf_require_positive( ...
    configuration.highpass, 'happeerPrecleanOptions.highpass');
nf_require_positive( ...
    configuration.lowpass, 'happeerPrecleanOptions.lowpass');
if configuration.highpass >= configuration.lowpass
    error('nf_preprocess:InvalidHAPPEOptions', ...
        'The HAPPE QC highpass must be below its lowpass.');
end
if ~qcFrequenciesSupplied
    qcUpper = min(configuration.lowpass, pipelineOptions.resample / 2);
    configuration.qcFrequencies = ...
        configuration.qcFrequencies( ...
        configuration.qcFrequencies < qcUpper);
end
frequencies = configuration.qcFrequencies;
if ~isnumeric(frequencies) || ~isreal(frequencies) || ...
        isempty(frequencies) || any(~isfinite(frequencies(:))) || ...
        any(frequencies(:) <= 0)
    error('nf_preprocess:InvalidHAPPEOptions', ...
        'happeerPrecleanOptions.qcFrequencies must be positive finite values.');
end
configuration.qcFrequencies = reshape(double(frequencies), 1, []);
end

function erpMode = nf_happeer_erp_mode(supplied, events)
erpMode = ~isempty(events);
if ~isfield(supplied, 'erpMode') || isempty(supplied.erpMode)
    return
end
if ~nf_is_logical_scalar(supplied.erpMode)
    error('nf_preprocess:InvalidHAPPEOptions', ...
        'happeerPrecleanOptions.erpMode must be empty or binary.');
end
erpMode = logical(supplied.erpMode);
end

function [EEG, info] = nf_apply_deferred_happe_filter( ...
        EEG, preliminaryInfo, options)
info = preliminaryInfo;
preWaveletInfo = preliminaryInfo;
postWaveletInfo = struct();
postWaveletInfo.started = datestr(now, 30); %#ok<TNOW1,DATST>
postWaveletInfo.inputSrate = EEG.srate;
postWaveletInfo.inputPnts = EEG.pnts;
postWaveletInfo.inputTrials = EEG.trials;
postWaveletInfo.applied = struct( ...
    'highpass', false, 'notch', false, ...
    'lowpass', false, 'resample', false);

if options.highpass > 0
    EEG = pop_eegfiltnew(EEG, 'locutoff', options.highpass);
    EEG = eeg_checkset(EEG);
    info.applied.highpass = true;
    postWaveletInfo.applied.highpass = true;
end
if options.lowpass > 0
    EEG = pop_eegfiltnew(EEG, 'hicutoff', options.lowpass);
    EEG = eeg_checkset(EEG);
    info.applied.lowpass = true;
    postWaveletInfo.applied.lowpass = true;
end

postWaveletInfo.finished = datestr(now, 30); %#ok<TNOW1,DATST>
postWaveletInfo.outputSrate = EEG.srate;
postWaveletInfo.outputPnts = EEG.pnts;
postWaveletInfo.outputTrials = EEG.trials;
info.requested.lowpassHz = options.lowpass;
info.requested.highpassHz = options.highpass;
info.requested.notchHz = options.notch;
info.requested.targetRateHz = options.resample;
info.finished = postWaveletInfo.finished;
info.outputSrate = EEG.srate;
info.outputPnts = EEG.pnts;
info.outputTrials = EEG.trials;
info.eeglabHistoryLength = 0;
if isfield(EEG, 'history') && ~isempty(EEG.history)
    info.eeglabHistoryLength = numel(EEG.history);
end
info.deferredAnalysisPassband = true;
info.stages = struct();
info.stages.preWaveletLineNoiseAndResample = preWaveletInfo;
info.stages.postWaveletAnalysisPassband = postWaveletInfo;

if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
EEG.etc.nf_filter = info;
if isfield(EEG.etc, 'nf_filter_history') && ...
        iscell(EEG.etc.nf_filter_history) && ...
        ~isempty(EEG.etc.nf_filter_history)
    EEG.etc.nf_filter_history{end} = info;
end
end

function names = nf_happe_coherence_names(frequencies)
names = cell(1, numel(frequencies));
for index = 1:numel(frequencies)
    names{index} = sprintf('coherence_%gHz', frequencies(index));
end
end

function identity = nf_dataset_identity(EEG)
identity = struct();
fields = {'setname', 'filename', 'filepath'};
for index = 1:numel(fields)
    fieldName = fields{index};
    if isfield(EEG, fieldName)
        identity.(fieldName) = EEG.(fieldName);
    end
end
end

function EEG = nf_restore_dataset_identity(EEG, identity)
fields = fieldnames(identity);
for index = 1:numel(fields)
    fieldName = fields{index};
    EEG.(fieldName) = identity.(fieldName);
end
end

function configuration = nf_merge_option_struct( ...
    configuration, supplied, optionName)
allowed = fieldnames(configuration);
names = fieldnames(supplied);
for index = 1:numel(names)
    name = names{index};
    if ~ismember(name, allowed)
        error('nf_preprocess:UnknownOptionField', ...
            'Unknown %s field: %s.', optionName, name);
    end
    configuration.(name) = supplied.(name);
end
end

function contract = nf_vendor_function_contract( ...
    functionName, expectedInputs, expectedOutputs)
resolved = which(functionName, '-all');
if isempty(resolved)
    error('nf_preprocess:MissingVendorFunction', ...
        'Required vendor function %s was not found.', functionName);
end
if ischar(resolved) && size(resolved, 1) > 1
    paths = cellstr(resolved);
elseif ischar(resolved)
    paths = {strtrim(resolved)};
else
    paths = resolved;
end
paths = cellfun(@strtrim, paths, 'UniformOutput', false);
paths = unique(paths, 'stable');
if numel(paths) ~= 1
    error('nf_preprocess:AmbiguousVendorFunction', ...
        ['Multiple copies of %s are on the MATLAB path. Remove path ' ...
        'ambiguity before running a named vendor method.'], functionName);
end
try
    observedInputs = nargin(functionName);
catch
    observedInputs = NaN;
end
try
    observedOutputs = nargout(functionName);
catch
    observedOutputs = NaN;
end
if ~isempty(expectedInputs) && observedInputs ~= expectedInputs
    error('nf_preprocess:VendorSignatureMismatch', ...
        '%s has %g inputs; this adapter requires %g.', ...
        functionName, observedInputs, expectedInputs);
end
if ~isempty(expectedOutputs) && observedOutputs ~= expectedOutputs
    error('nf_preprocess:VendorSignatureMismatch', ...
        '%s has %g outputs; this adapter requires %g.', ...
        functionName, observedOutputs, expectedOutputs);
end
metadata = dir(paths{1});
if isempty(metadata) || metadata.isdir
    error('nf_preprocess:UnreadableVendorFunction', ...
        '%s did not resolve to a readable source file.', functionName);
end
contract = struct();
contract.function = functionName;
contract.path = paths{1};
contract.nargin = observedInputs;
contract.nargout = observedOutputs;
contract.fileBytes = metadata.bytes;
contract.fileModified = metadata.date;
contract.sha256 = nf_file_sha256(paths{1});
if isempty(contract.sha256)
    error('nf_preprocess:VendorHashFailed', ...
        ['Could not calculate a SHA-256 identity for %s. A named vendor ' ...
        'stage will not run without a verifiable source-file identity.'], ...
        functionName);
end
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

function [EEG, info] = nf_make_final_epochs( ...
        EEG, options, behaviorInfo)
EEG = nf_ensure_event_fields(EEG);
[EEG, epochDecomposition] = nf_detach_ica_for_epoching(EEG);
if ~isempty(options.events)
    EEG = nf_sort_continuous_events(EEG);
    survivingValidation = nf_validate_requested_events(EEG, options.events);
    candidateEventIndices = ...
        sort(unique([survivingValidation.matchedEventIndices{:}]));
    candidateEventCount = numel(candidateEventIndices);
    sourceEventLatencies = double([EEG.event.latency]);
    [EEG, candidateEventPositions] = pop_epoch(EEG, ...
        survivingValidation.popEpochEvents, ...
        options.epochLimits, ...
        'eventindices', candidateEventIndices, ...
        'epochinfo', 'yes');
    candidateEventPositions = reshape(candidateEventPositions, 1, []);
    nf_validate_accepted_event_positions( ...
        candidateEventPositions, candidateEventCount, EEG.trials);
    candidateEventIndices = ...
        candidateEventIndices(candidateEventPositions);
    candidateEventLatencies = ...
        sourceEventLatencies(candidateEventIndices);
    mode = 'event-locked';
    eventTypes = options.events;
    sourceEpochIds = candidateEventPositions;
else
    survivingValidation = struct();
    candidateEventPositions = [];
    candidateEventIndices = [];
    candidateEventLatencies = [];
    EEG = eeg_regepochs(EEG, 'recurrence', options.continuousEpochLength, ...
        'limits', [0 options.continuousEpochLength], 'rmbase', NaN, ...
        'eventtype', 'nf_fixed_epoch');
    mode = 'fixed-length';
    eventTypes = {};
    candidateEventCount = EEG.trials;
    sourceEpochIds = 1:EEG.trials;
end
EEG = nf_restore_ica_after_epoching(EEG, epochDecomposition);
[EEG, behaviorAttachment] = nf_attach_final_behavior( ...
    EEG, options.behavior, behaviorInfo, mode, ...
    candidateEventCount, candidateEventPositions);
[EEG, epochCheckDecomposition] = nf_detach_ica_for_epoching(EEG);
EEG = eeg_checkset(EEG);
EEG = nf_restore_ica_after_epoching( ...
    EEG, epochCheckDecomposition);
populationDetectors = {'faster', 'jointprobability', 'eeglabstats'};
minimumEpochs = 1 + ...
    double(any(ismember(options.epochDetectors, populationDetectors)));
if EEG.trials < minimumEpochs
    error('nf_preprocess:InsufficientEpochs', ...
        ['Final epoching produced fewer than %d complete epoch(s), the ' ...
        'minimum for the selected detector configuration.'], minimumEpochs);
end
if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
EEG.etc.nf_epoch_ids = sourceEpochIds;

info = struct();
info.mode = mode;
info.events = eventTypes;
info.inputEventValidation = options.eventValidation;
info.preEpochEventValidation = survivingValidation;
info.candidateEventIndices = candidateEventIndices;
info.candidateEventPositions = candidateEventPositions;
info.candidateEventLatenciesSamples = candidateEventLatencies;
info.candidateEventLatenciesSeconds = ...
    (candidateEventLatencies - 1) / EEG.srate;
info.acceptedEventIndices = candidateEventIndices;
info.acceptedEventPositions = candidateEventPositions;
if strcmp(mode, 'event-locked')
    info.acceptedEventSemantics = ...
        ['Complete final candidate epochs accepted by pop_epoch before the ' ...
        'configured final artifact detectors and repair/rejection policies.'];
else
    info.acceptedEventSemantics = ...
        ['Not applicable to fixed-length eeg_regepochs creation; event ' ...
        'identity ledgers are empty.'];
end
info.nAcceptedEvents = numel(candidateEventIndices);
info.requestedLimitsSeconds = options.epochLimits;
info.actualLimitsSeconds = [EEG.xmin EEG.xmax];
info.nCreated = EEG.trials;
info.minimumRequired = minimumEpochs;
info.pnts = EEG.pnts;
info.sourceEpochIds = sourceEpochIds;
info.behavior = behaviorAttachment;
end

function [EEG, decomposition] = nf_detach_ica_for_epoching(EEG)
fieldNames = { ...
    'icaweights', ...
    'icasphere', ...
    'icawinv', ...
    'icachansind'};
decomposition = struct();
decomposition.gcompreject = [];

if isfield(EEG, 'reject') && isstruct(EEG.reject) && ...
        isfield(EEG.reject, 'gcompreject')
    decomposition.gcompreject = EEG.reject.gcompreject;
end

for fieldIndex = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIndex};

    if isfield(EEG, fieldName)
        decomposition.(fieldName) = EEG.(fieldName);
        EEG.(fieldName) = [];
    else
        decomposition.(fieldName) = [];
    end
end

if isfield(EEG, 'icaact')
    EEG.icaact = [];
end
end

function EEG = nf_restore_ica_after_epoching(EEG, decomposition)
fieldNames = { ...
    'icaweights', ...
    'icasphere', ...
    'icawinv', ...
    'icachansind'};

for fieldIndex = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIndex};
    EEG.(fieldName) = decomposition.(fieldName);
end

EEG.icaact = [];

componentCount = size(decomposition.icaweights, 1);

if componentCount > 0
    if ~isfield(EEG, 'reject') || ~isstruct(EEG.reject)
        EEG.reject = struct();
    end

    if numel(decomposition.gcompreject) == componentCount
        EEG.reject.gcompreject = reshape( ...
            logical(decomposition.gcompreject), ...
            1, ...
            componentCount);
    else
        EEG.reject.gcompreject = false(1, componentCount);
    end
elseif isfield(EEG, 'reject') && isstruct(EEG.reject)
    EEG.reject.gcompreject = [];
end

if isfield(EEG, 'chaninfo') && isstruct(EEG.chaninfo)
    EEG.chaninfo.icachansind = decomposition.icachansind;
end
end

function [EEG, attachment] = nf_attach_final_behavior( ...
        EEG, behavior, behaviorInfo, mode, candidateCount, ...
        acceptedPositions)
attachment = behaviorInfo;
attachment.mapping = 'not-present';
attachment.sourceEntries = numel(behavior);
attachment.outputEntries = 0;
attachment.outputField = 'EEG.etc.behav';
attachment.compatibilityField = 'EEG.etc.behavior';

if ~behaviorInfo.present
    EEG = nf_clear_eeg_behavior(EEG);
    return
end

if strcmp(mode, 'event-locked')
    if numel(behavior) ~= candidateCount
        error('nf_preprocess:BehaviorEventMismatch', ...
            ['Behavior has %d entries, but final epoching found %d ' ...
            'requested event candidates. Supply one behavior entry per ' ...
            'requested event in increasing event-latency order.'], ...
            numel(behavior), candidateCount);
    end
    behavior = behavior(acceptedPositions);
    attachment.mapping = ...
        'requested-event-order-through-pop-epoch-positions';
else
    if numel(behavior) ~= EEG.trials
        error('nf_preprocess:BehaviorFixedEpochMismatch', ...
            ['Behavior has %d entries, but fixed-length final epoching ' ...
            'created %d complete epochs.'], numel(behavior), EEG.trials);
    end
    attachment.mapping = 'fixed-epoch-order';
end

behavior = reshape(behavior, 1, []);
if numel(behavior) ~= EEG.trials
    error('nf_preprocess:BehaviorEpochMappingFailed', ...
        ['Mapped behavior has %d entries, but final epoching created %d ' ...
        'epochs.'], numel(behavior), EEG.trials);
end
EEG = nf_set_eeg_behavior(EEG, behavior);
attachment.outputEntries = numel(behavior);
end

function [options, mapping] = nf_remap_epoch_channel_options( ...
        options, sourceMontage, targetMontage)
mapping = struct();
mapping.source = 'selected raw-input montage';
mapping.target = 'post-ICA working montage';
mapping.policy = ...
    ['Explicit numeric channel fields are resolved to raw-input labels and ' ...
    'then remapped. A requested channel removed upstream is an error.'];
mapping.entries = {};
optionNames = { ...
    'amplitudeOptions', ...
    'fftOptions', ...
    'peakToPeakOptions', ...
    'stepOptions', ...
    'gradientOptions', ...
    'flatlineOptions', ...
    'clippingOptions', ...
    'fasterEpochOptions', ...
    'eeglabEpochOptions', ...
    'jointProbabilityOptions'};
detectorNames = {'threshold', 'fft', 'peak2peak', 'step', 'gradient', ...
    'flatline', 'clipping', 'faster', 'eeglabstats', ...
    'jointprobability'};
channelFields = {'channelIndices', 'localChannelIndices', 'ignoredChannels'};
sourceLabels = string({sourceMontage.labels});
targetLabels = string({targetMontage.labels});

for optionIndex = 1:numel(optionNames)
    if ~any(strcmp(options.epochDetectors, detectorNames{optionIndex}))
        continue
    end
    optionName = optionNames{optionIndex};
    settings = options.(optionName);
    for fieldIndex = 1:numel(channelFields)
        fieldName = channelFields{fieldIndex};
        if ~isfield(settings, fieldName) || isempty(settings.(fieldName))
            continue
        end
        requested = settings.(fieldName);
        if ~isnumeric(requested) || ~isreal(requested) || ...
                ~isvector(requested) || any(~isfinite(requested)) || ...
                any(requested ~= round(requested)) || any(requested < 1) || ...
                any(requested > numel(sourceLabels)) || ...
                numel(unique(requested)) ~= numel(requested)
            error('nf_preprocess:InvalidEpochChannelIndices', ...
                ['%s.%s must contain unique indices into the selected ' ...
                'raw-input montage.'], optionName, fieldName);
        end
        requested = reshape(double(requested), 1, []);
        requestedLabels = sourceLabels(requested);
        [present, remapped] = ismember( ...
            lower(requestedLabels), lower(targetLabels));
        if any(~present)
            missingLabels = strjoin(cellstr(requestedLabels(~present)), ', ');
            error('nf_preprocess:EpochChannelUnavailable', ...
                ['%s.%s requested channel(s) removed before epoch ' ...
                'detection: %s. Change the upstream channel method or ' ...
                'remove those channels from the detector selection.'], ...
                optionName, fieldName, missingLabels);
        end
        settings.(fieldName) = reshape(remapped, 1, []);
        entry = struct();
        entry.option = optionName;
        entry.field = fieldName;
        entry.requestedIndices = requested;
        entry.requestedLabels = cellstr(requestedLabels);
        entry.resolvedIndices = reshape(remapped, 1, []);
        mapping.entries{end + 1} = entry;
    end
    options.(optionName) = settings;
end
end

function [settings, mapping] = nf_remap_faster_ica_options( ...
        settings, method, sourceMontage, targetMontage)
mapping = struct();
mapping.applied = false;
mapping.source = 'selected raw-input montage';
mapping.target = 'pre-ICA working montage';
mapping.requestedIndices = [];
mapping.requestedLabels = {};
mapping.resolvedIndices = [];
if ~strcmp(method, 'faster') || ...
        ~isfield(settings, 'eogChannels') || ...
        isempty(settings.eogChannels)
    return
end

requested = settings.eogChannels;
if ~isnumeric(requested) || ~isreal(requested) || ...
        ~isvector(requested) || any(~isfinite(requested)) || ...
        any(requested ~= round(requested)) || any(requested < 1) || ...
        any(requested > numel(sourceMontage)) || ...
        numel(unique(requested)) ~= numel(requested)
    error('nf_preprocess:InvalidFASTERICAChannels', ...
        ['fasterICAOptions.eogChannels must contain unique indices into ' ...
        'the selected raw-input montage.']);
end
requested = reshape(double(requested), 1, []);
sourceLabels = string({sourceMontage.labels});
targetLabels = string({targetMontage.labels});
requestedLabels = sourceLabels(requested);
[present, remapped] = ismember(lower(requestedLabels), lower(targetLabels));
if any(~present)
    missingLabels = strjoin(cellstr(requestedLabels(~present)), ', ');
    error('nf_preprocess:FASTERICAChannelUnavailable', ...
        ['fasterICAOptions.eogChannels requested channel(s) removed before ' ...
        'ICA: %s. Change upstream channel handling or the EOG selection.'], ...
        missingLabels);
end

settings.eogChannels = reshape(double(remapped), 1, []);
mapping.applied = true;
mapping.requestedIndices = requested;
mapping.requestedLabels = cellstr(requestedLabels);
mapping.resolvedIndices = settings.eogChannels;
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

function EEG = nf_sort_continuous_events(EEG)
if ~isfield(EEG, 'event') || isempty(EEG.event)
    return
end
% nf_validate_event_latencies(EEG);
latencies = double([EEG.event.latency]);
[~, eventOrder] = sort(latencies);
EEG.event = EEG.event(eventOrder);
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
acquisitionEEG = true(1, EEG.nbchan);
typesAvailable = isfield(EEG.chanlocs, 'type');
if typesAvailable
    for index = 1:EEG.nbchan
        channelType = EEG.chanlocs(index).type;
        if isempty(channelType)
            acquisitionEEG(index) = true;
        elseif nf_is_text(channelType)
            acquisitionEEG(index) = ...
                strcmpi(strtrim(char(channelType)), 'eeg');
        else
            acquisitionEEG(index) = false;
        end
    end
    if ~any(acquisitionEEG)
        acquisitionEEG = true(1, EEG.nbchan);
    end
end

if isempty(requested)
    selected = acquisitionEEG;
    if typesAvailable
        selectionSource = 'chanlocs.type (EEG or empty)';
    else
        selectionSource = 'all input channels (no chanlocs.type field)';
    end
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
info.nAcquisitionEEGChannels = sum(acquisitionEEG);
end

function useUnfiltered = nf_channel_detection_uses_unfiltered_input(method)
useUnfiltered = ismember(method, {'cleanrawdata', 'prep'});
end

function [EEG, transfer] = nf_apply_channel_detection( ...
    EEG, detectorEEG, channelInfo)
workingLabels = {EEG.chanlocs.labels};
detectorLabels = {detectorEEG.chanlocs.labels};
nf_require_unique_labels(workingLabels);
nf_require_unique_labels(detectorLabels);
[present, keepIndices] = ismember( ...
    lower(string(detectorLabels)), lower(string(workingLabels)));
if any(~present)
    missing = strjoin(cellstr(string(detectorLabels(~present))), ', ');
    error('nf_preprocess:ChannelDetectionTransferFailed', ...
        ['The unfiltered detector returned channel labels that are absent ' ...
        'from the analysis-filtered data: %s'], missing);
end
if any(diff(keepIndices) <= 0)
    error('nf_preprocess:ChannelDetectionOrderChanged', ...
        ['The unfiltered detector changed retained channel order. ' ...
        'NeuroFreq requires label-preserving channel removal.']);
end
removedMask = true(1, numel(workingLabels));
removedMask(keepIndices) = false;
if any(removedMask)
    EEG = pop_select(EEG, 'channel', keepIndices);
    EEG = eeg_checkset(EEG);
end
if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
EEG.etc.badchans = channelInfo.artifact.nDetected;
EEG.etc.badchanindices = channelInfo.artifact.indices;
EEG.etc.badchanlabels = channelInfo.artifact.labels;
EEG.etc.nf_removed_reference = struct( ...
    'removed', channelInfo.reference.removedZeroReference, ...
    'index', channelInfo.reference.index, ...
    'label', channelInfo.reference.label);
EEG.etc.nf_badchans = channelInfo;

transfer = struct();
transfer.applied = any(removedMask);
transfer.matching = 'case-insensitive unique channel labels';
transfer.retainedIndicesInAnalysisData = keepIndices;
transfer.retainedLabels = detectorLabels;
transfer.removedIndicesInAnalysisData = find(removedMask);
transfer.removedLabels = workingLabels(removedMask);
transfer.nRemoved = sum(removedMask);
transfer.sampleDecisionsTransferred = false;
transfer.statement = ...
    ['Only channel-removal decisions were transferred. No samples or ' ...
    'filtered values from the diagnostic copy entered the analysis data.'];
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

function preset = nf_normalize_preset(value)
normalized = lower(regexprep(strtrim(char(value)), '[^a-zA-Z0-9]', ''));
switch normalized
    case {'bdc', 'adult'}
        preset = 'bdc';
    case {'made', 'child'}
        preset = 'made';
    case 'prep'
        preset = 'prep';
    case 'faster'
        preset = 'faster';
    case {'happeer', 'happeerp'}
        preset = 'happeer';
    case {'cleanrawdata', 'cleanraw'}
        preset = 'cleanrawdata';
    case {'eeglab', 'legacy', 'legacyeeglab'}
        preset = 'eeglab';
    otherwise
        error('nf_preprocess:UnknownPreset', ...
            ['preset must be BDC, MADE, PREP, FASTER, HAPPE+ER, ' ...
            'cleanrawdata, or EEGLAB.']);
end
end

function stage = nf_normalize_epoch_stage(value)
normalized = lower(regexprep(strtrim(char(value)), '[^a-zA-Z0-9]', ''));
switch normalized
    case 'auto'
        stage = 'auto';
    case {'beforeica', 'preica', 'before'}
        stage = 'beforeica';
    case {'afterica', 'postica', 'after'}
        stage = 'afterica';
    otherwise
        error('nf_preprocess:UnknownEpochStage', ...
            'epochStage must be auto, beforeica, or afterica.');
end
end

function [resolved, details] = nf_resolve_epoch_stage(requested, options)
eventTypes = nf_pop_epoch_events(options.events);
eventTypeCount = nf_unique_event_type_count(eventTypes);

details = struct();
details.requested = requested;
details.resolved = '';
details.reason = '';
details.eventTypeCount = eventTypeCount;

if strcmp(requested, 'beforeica')
    resolved = 'beforeica';
    details.reason = 'Explicit epochStage override.';
elseif strcmp(requested, 'afterica')
    resolved = 'afterica';
    details.reason = 'Explicit epochStage override.';
elseif strcmp(options.preset, 'bdc') && ...
        eventTypeCount == 1 && ~strcmp(options.icaMethod, 'none')
    resolved = 'beforeica';
    details.reason = ...
        ['BDC auto behavior with exactly one requested event type: use the ' ...
        'requested event epochs for ICA fitting.'];
else
    resolved = 'afterica';
    if isempty(eventTypes)
        details.reason = ...
            'No requested event type was available for event-scoped ICA fitting.';
    elseif strcmp(options.icaMethod, 'none')
        details.reason = ...
            'ICA is disabled, so only final analysis epoching is performed.';
    elseif ~strcmp(options.preset, 'bdc')
        details.reason = ...
            'This preset retains its existing fixed-epoch ICA preparation.';
    elseif eventTypeCount ~= 1
        details.reason = ...
            ['BDC auto behavior requires exactly one requested event type; ' ...
            'multiple types retain fixed-epoch ICA preparation.'];
    end
end

details.resolved = resolved;
details.eventScopedIcaTraining = strcmp(resolved, 'beforeica');
if strcmp(resolved, 'beforeica')
    details.scope = ...
        ['The requested event epochs define the ICA fitting copy. Filtering ' ...
        'and precleaning remain continuous; final epochs are materialized ' ...
        'after component subtraction.'];
elseif strcmp(options.icaMethod, 'none')
    details.scope = ...
        ['ICA is disabled. Requested event epochs or fixed analysis epochs ' ...
        'are materialized after filtering and precleaning.'];
else
    details.scope = ...
        ['ICA uses the existing fixed-length training preparation. Final ' ...
        'analysis epochs are materialized after component subtraction.'];
end
if strcmp(options.icaMethod, 'none')
    details.finalEpochMaterialization = ...
        ['after filtering and precleaning; ICA fitting and subtraction ' ...
        'are skipped'];
else
    details.finalEpochMaterialization = ...
        'after ICA component subtraction';
end
end

function count = nf_unique_event_type_count(events)
if isempty(events)
    count = 0;
    return
end
keys = cell(1, numel(events));
for eventIndex = 1:numel(events)
    value = events{eventIndex};
    if isstring(value) && isscalar(value)
        value = char(value);
    end
    if ischar(value)
        keys{eventIndex} = deblank(value);
    elseif isnumeric(value) && isscalar(value) && isfinite(value)
        keys{eventIndex} = num2str(value, 15);
    else
        keys{eventIndex} = sprintf( ...
            '<invalid-%s-%d>', class(value), eventIndex);
    end
end
count = numel(unique(keys));
end

function [options, definition] = nf_apply_preset_defaults( ...
    options, usingDefaults)
definition = struct();
definition.name = options.preset;
definition.fullPublishedPipelineClaim = false;
definition.contract = ...
    ['NeuroFreq composition. Each named external stage must execute the ' ...
    'installed vendor implementation and records its own provenance.'];
definition.futureComparisonNote = ...
    ['Compare resolved stages and provenance, not the preset label alone. ' ...
    'A vendor-stage call is not reported as a complete vendor pipeline.'];

switch options.preset
    case 'bdc'
        definition.displayName = 'BDC';
        definition.description = ...
            'BDC lab pipeline: FASTER channels, GEDAI, ICLabel, threshold and FFT.';
        defaults = nf_preset_values( ...
            'faster', 'gedai', 'iclabel', {'threshold', 'fft'});
    case 'made'
        definition.displayName = 'MADE';
        definition.description = ...
            ['MADE-oriented pipeline: FASTER channels and adjusted_ADJUST. ' ...
            'GEDAI is intentionally absent.'];
        defaults = nf_preset_values( ...
            'faster', 'none', 'adjustedadjust', {'threshold', 'fft'});
    case 'prep'
        definition.displayName = 'PREP';
        definition.fullPublishedPipelineClaim = true;
        definition.description = ...
            ['Official complete PREP routine, followed only by the resolved ' ...
            'NeuroFreq filter/epoch container; no ICA or artifact detector ' ...
            'is added by default.'];
        defaults = nf_preset_values( ...
            'none', 'prep', 'none', {'none'});
        defaults.globalInterpolation = false;
        defaults.rereference = false;
    case 'faster'
        definition.displayName = 'FASTER';
        definition.description = ...
            ['Vendor-exact FASTER channel, epoch, local epoch-channel, and ' ...
            'component stages embedded in NeuroFreq orchestration.'];
        defaults = nf_preset_values( ...
            'faster', 'none', 'faster', {'faster'});
        defaults.fasterEpochOptions = struct( ...
            'localChannelDetection', true, ...
            'exactVendorInterpolation', true);
    case 'happeer'
        definition.displayName = 'HAPPE+ER';
        definition.description = ...
            ['NeuroFreq composition calling installed HAPPE channel and ' ...
            'wavelet stages, MARA, amplitude, and joint-probability epoch ' ...
            'detection. It does not claim complete HAPPE+ER stage-order ' ...
            'or full-pipeline equivalence.'];
        defaults = nf_preset_values( ...
            'happeer', 'happeer', 'mara', ...
            {'threshold', 'jointprobability'});
    case 'cleanrawdata'
        definition.displayName = 'clean_rawdata';
        definition.description = ...
            ['Official clean_flatlines, clean_channels, and clean_asr ' ...
            'functions without an added ICA or epoch classifier.'];
        defaults = nf_preset_values( ...
            'cleanrawdata', 'asr', 'none', {'none'});
        defaults.globalInterpolation = false;
    case 'eeglab'
        definition.displayName = 'EEGLAB';
        definition.description = ...
            ['Legacy EEGLAB channel statistics and composable amplitude, ' ...
            'spectral, and epoch-statistics rules.'];
        defaults = nf_preset_values( ...
            'eeglab', 'none', 'none', ...
            {'threshold', 'fft', 'eeglabstats'});
end

names = fieldnames(defaults);
for index = 1:numel(names)
    name = names{index};
    if strcmp(options.preset, 'faster') && ...
            strcmp(name, 'fasterEpochOptions')
        resolvedFasterEpochOptions = defaults.fasterEpochOptions;
        suppliedFields = fieldnames(options.fasterEpochOptions);
        for fieldIndex = 1:numel(suppliedFields)
            fieldName = suppliedFields{fieldIndex};
            resolvedFasterEpochOptions.(fieldName) = ...
                options.fasterEpochOptions.(fieldName);
        end
        localDetectionSupplied = ...
            any(strcmp(suppliedFields, 'localChannelDetection'));
        exactInterpolationSupplied = ...
            any(strcmp(suppliedFields, 'exactVendorInterpolation'));
        if localDetectionSupplied && ~exactInterpolationSupplied && ...
                nf_is_logical_scalar( ...
                resolvedFasterEpochOptions.localChannelDetection) && ...
                ~logical( ...
                resolvedFasterEpochOptions.localChannelDetection)
            resolvedFasterEpochOptions.exactVendorInterpolation = false;
        end
        options.fasterEpochOptions = resolvedFasterEpochOptions;
        continue
    end
    applyDefault = nf_was_default(usingDefaults, name);
    if ismember(name, {'channelMethod', 'precleanMethod', 'icaMethod'}) && ...
            isempty(strtrim(char(options.(name))))
        applyDefault = true;
    end
    if applyDefault
        options.(name) = defaults.(name);
    end
end
definition.defaultStages = defaults;
definition.explicitOverrides = nf_explicit_option_names( ...
    options, usingDefaults);
end

function values = nf_preset_values( ...
    channelMethod, precleanMethod, icaMethod, epochDetectors)
values = struct();
values.channelMethod = channelMethod;
values.precleanMethod = precleanMethod;
values.icaMethod = icaMethod;
values.epochDetectors = epochDetectors;
end

function names = nf_explicit_option_names(options, usingDefaults)
allNames = fieldnames(options);
defaultNames = lower(string(usingDefaults));
explicitMask = ~ismember(lower(string(allNames)), defaultNames);
names = allNames(explicitMask);
end

function wasDefault = nf_was_default(usingDefaults, requestedName)
wasDefault = any(strcmpi(usingDefaults, requestedName));
end

function method = nf_normalize_channel_method(value)
normalized = lower(regexprep(strtrim(char(value)), '[^a-zA-Z0-9]', ''));
switch normalized
    case 'faster'
        method = 'faster';
    case {'cleanrawdata', 'cleanraw'}
        method = 'cleanrawdata';
    case 'prep'
        method = 'prep';
    case {'happeer', 'happeerp'}
        method = 'happeer';
    case {'eeglab', 'legacy', 'legacyeeglab'}
        method = 'eeglab';
    case 'none'
        method = 'none';
    otherwise
        error('nf_preprocess:UnknownChannelMethod', ...
            ['channelMethod must be faster, cleanrawdata, prep, happeer, ' ...
            'eeglab, or none.']);
end
end

function method = nf_normalize_preclean_method(value)
normalized = lower(regexprep(strtrim(char(value)), '[^a-zA-Z0-9]', ''));
switch normalized
    case 'gedai'
        method = 'gedai';
    case {'asr', 'cleanasr'}
        method = 'asr';
    case 'prep'
        method = 'prep';
    case {'happeer', 'happeerp', 'happewavelet'}
        method = 'happeer';
    case 'none'
        method = 'none';
    otherwise
        error('nf_preprocess:UnknownPrecleanMethod', ...
            'precleanMethod must be gedai, asr, prep, happeer, or none.');
end
end

function method = nf_resolve_ica_method(~, supplied)
normalized = lower(regexprep(strtrim(char(supplied)), '[^a-zA-Z0-9]', ''));
switch normalized
    case 'iclabel'
        method = 'iclabel';
    case {'made', 'adjustedadjust', 'adjustedadust'}
        method = 'adjustedadjust';
    case 'adjust'
        method = 'adjust';
    case 'mara'
        method = 'mara';
    case 'faster'
        method = 'faster';
    case 'none'
        method = 'none';
    otherwise
        error('nf_preprocess:UnknownICAMethod', ...
            ['icaMethod must be iclabel, adjustedadjust, adjust, mara, ' ...
            'faster, or none.']);
end
end

function detectors = nf_normalize_epoch_detectors(value)
if isempty(value)
    detectors = {'none'};
elseif ischar(value)
    detectors = {value};
elseif isstring(value)
    detectors = cellstr(value(:));
else
    detectors = value(:)';
end
normalized = cell(1, numel(detectors));
for index = 1:numel(detectors)
    token = lower(regexprep(strtrim(char(detectors{index})), ...
        '[^a-zA-Z0-9]', ''));
    switch token
        case {'threshold', 'amplitude', 'voltage'}
            normalized{index} = 'threshold';
        case {'fft', 'spectral', 'muscle'}
            normalized{index} = 'fft';
        case {'peak2peak', 'peaktopeak', 'movingwindowpeaktopeak', 'mwpp'}
            normalized{index} = 'peak2peak';
        case {'step', 'electrodepop', 'pop'}
            normalized{index} = 'step';
        case {'gradient', 'fasttransition', 'derivative'}
            normalized{index} = 'gradient';
        case 'flatline'
            normalized{index} = 'flatline';
        case {'clipping', 'clip'}
            normalized{index} = 'clipping';
        case 'faster'
            normalized{index} = 'faster';
        case {'eeglabstats', 'eeglabstatistics', 'eeglab'}
            normalized{index} = 'eeglabstats';
        case {'jointprobability', 'jointprob', 'joint'}
            normalized{index} = 'jointprobability';
        case {'none', 'off'}
            normalized{index} = 'none';
        otherwise
            error('nf_preprocess:UnknownEpochDetector', ...
                'Unknown epoch detector: %s.', char(detectors{index}));
    end
end
detectors = unique(normalized, 'stable');
if numel(detectors) > 1 && any(strcmp(detectors, 'none'))
    error('nf_preprocess:ConflictingEpochDetectors', ...
        'epochDetectors cannot combine none with active detectors.');
end
end

function bands = nf_normalize_frequency_bands(value, defaultBand)
if isempty(value)
    bands = reshape(defaultBand, 1, 2);
else
    bands = value;
    if isnumeric(bands) && isvector(bands)
        if numel(bands) == 2
            bands = reshape(bands, 1, 2);
        elseif numel(bands) == 4
            bands = reshape(bands, 1, 4);
        end
    end
end
end

function limits = nf_fft_frequency_limits(value)
if isnumeric(value)
    limits = double(value(:, 1:2));
    return
end
limits = zeros(numel(value), 2);
for index = 1:numel(value)
    if isfield(value(index), 'frequencyRangeHz')
        range = value(index).frequencyRangeHz;
    else
        range = value(index).frequencyRange;
    end
    limits(index, :) = reshape(double(range), 1, 2);
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
reservedQualityOptions = {'final', 'report', 'plot', 'visible', 'icaModel'};
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
if any(strcmp(options.epochDetectors, 'fft'))
    if ~isempty(options.fftBands) && ...
            isfield(options.fftOptions, 'bands')
        error('nf_preprocess:ConflictingFFTBands', ...
            ['Specify FFT bands with either fftBands or fftOptions.bands, ' ...
            'not both.']);
    end
    effectiveFFTBands = options.fftBands;
    if isempty(effectiveFFTBands)
        if isfield(options.fftOptions, 'bands') && ...
                ~isempty(options.fftOptions.bands)
            effectiveFFTBands = options.fftOptions.bands;
        else
            effectiveFFTBands = options.muscleRange;
        end
    end
    fftFrequencyLimits = nf_fft_frequency_limits(effectiveFFTBands);
    for bandIndex = 1:size(fftFrequencyLimits, 1)
        band = fftFrequencyLimits(bandIndex, :);
        if band(1) < 0 || band(2) > options.lowpass || ...
                band(2) >= options.resample / 2
            error('nf_preprocess:InvalidFFTBand', ...
                ['Every complete fftBands row must be retained by lowpass ' ...
                'and resampling. Invalid row: [%g %g].'], ...
                band(1), band(2));
        end
    end
end
if ~strcmp(options.icaMethod, 'none') && ...
        (options.icaTrainingFrequencies(1) < 0 || ...
        options.icaTrainingFrequencies(2) > options.lowpass || ...
        options.icaTrainingFrequencies(2) >= options.resample / 2)
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
if strcmp(options.epochStage, 'beforeica') && isempty(options.events)
    error('nf_preprocess:MissingICAEvents', ...
        ['epochStage=''beforeica'' requires events and epochLimits so the ' ...
        'requested task epochs can define the ICA fitting data.']);
end
if strcmp(options.epochStage, 'beforeica') && ...
        strcmp(options.icaMethod, 'none')
    error('nf_preprocess:ICAStageWithoutICA', ...
        'epochStage=''beforeica'' cannot be used when icaMethod=''none''.');
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
if ~ismember(options.channelMethod, ...
        {'faster', 'cleanrawdata', 'prep', 'happeer', 'eeglab', 'none'})
    error('nf_preprocess:UnknownChannelMethod', ...
        ['channelMethod must be faster, cleanrawdata, prep, happeer, ' ...
        'eeglab, or none.']);
end
if ~ismember(options.precleanMethod, ...
        {'gedai', 'asr', 'prep', 'happeer', 'none'})
    error('nf_preprocess:UnknownPrecleanMethod', ...
        'precleanMethod must be gedai, asr, prep, happeer, or none.');
end
if ~strcmp(options.icaMethod, 'none') && ...
        ~ismember(options.icaAlgorithm, {'runamica15', 'runica', 'picard'})
    error('nf_preprocess:UnknownICAAlgorithm', ...
        'icaAlgorithm must be runamica15, runica, or picard.');
end
if ~strcmp(options.icaMethod, 'none') && ...
        strcmp(options.icaAlgorithm, 'picard') && ...
        ~ismember(options.picardMode, {'standard', 'ortho'})
    error('nf_preprocess:InvalidPicardMode', ...
        'picardMode must be ''standard'' or ''ortho''.');
end
if strcmp(options.channelMethod, 'cleanrawdata') && ...
        options.cleanHighpass >= options.resample / 2
    error('nf_preprocess:InvalidCleanHighpass', ...
        'cleanHighpass must be below the post-resampling Nyquist frequency.');
end
if ~strcmp(options.icaMethod, 'none') && ...
        options.icaTrainingHighpass >= options.resample / 2
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
if strcmp(options.icaMethod, 'adjustedadjust') && options.resample < 100
    error('nf_preprocess:MADESamplingRate', ...
        ['The MADE adjusted_ADJUST feature extractor requires a sampling ' ...
        'rate of at least 100 Hz.']);
end
if strcmp(options.icaMethod, 'adjustedadjust') && ...
        options.icaTrainingEpochLength * options.resample < 100
    error('nf_preprocess:MADEEpochLength', ...
        ['icaTrainingEpochLength must yield at least 100 samples per MADE ' ...
        'adjusted_ADJUST feature epoch.']);
end
if nf_contains_boundary_event(EEG) && ...
        (ismember(options.precleanMethod, {'gedai', 'asr', 'happeer'}) || ...
        ~strcmp(options.icaMethod, 'none'))
    error('nf_preprocess:BoundaryUnsupported', ...
        ['The selected precleaner or common ICA-training filter does not ' ...
        'honor EEGLAB boundary events. Segment the recording before ' ...
        'nf_preprocess, or use PREP/none with ICA disabled.']);
end
if strcmp(options.channelMethod, 'faster')
    nf_validate_faster_options(options.fasterOptions);
end
if ~ismember(options.icaMethod, {'iclabel', 'none'}) && ...
        options.aggressiveICA
    error('nf_preprocess:ClassifierOptionMismatch', ...
        'aggressiveICA applies only to ICLabel.');
end
if ~ismember(options.icaMethod, {'iclabel', 'none'}) && ...
        ~isempty(options.iclabelThresholds)
    error('nf_preprocess:ClassifierOptionMismatch', ...
        'iclabelThresholds apply only to ICLabel.');
end
end

function nf_preflight_dependencies(options)
required = {'eeg_checkset', 'pop_select', 'pop_eegfiltnew', ...
    'pop_resample', 'nf_filter', ...
    'nf_badchans', 'nf_cleanic', 'nf_thresh'};
if nf_pipeline_requires_geometry(options)
    required{end + 1} = 'convertlocs';
end
mayRestoreChannels = options.globalInterpolation && ...
    (~strcmp(options.channelMethod, 'none') || ...
    ~strcmp(options.icaMethod, 'none'));
epochRepairRequirements = nf_epoch_repair_requirements(options);
if mayRestoreChannels || ...
        epochRepairRequirements.requiresEeglabInterpolation
    required{end + 1} = 'eeg_interp';
end
if options.rereference || strcmp(options.channelMethod, 'faster')
    required{end + 1} = 'pop_reref';
end
if ~strcmp(options.icaMethod, 'none')
    required{end + 1} = 'pop_subcomp';
end
if ~strcmp(options.icaMethod, 'none') && ...
        ismember(options.icaAlgorithm, {'runica', 'picard'})
    required{end + 1} = 'pop_runica';
end
if any(strcmp(options.epochDetectors, 'threshold'))
    required{end + 1} = 'pop_eegthresh';
end
if any(strcmp(options.epochDetectors, 'fft'))
    required{end + 1} = 'pop_rejspec';
end
if ~all(strcmp(options.epochDetectors, 'none'))
    required{end + 1} = 'pop_rejepoch';
end
if ~isempty(options.events)
    required{end + 1} = 'pop_epoch';
else
    required{end + 1} = 'eeg_regepochs';
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
    fasterChannelDependencies = {'channel_properties', 'min_z', ...
        'distancematrix', 'hurst_exponent', 'nanmean'};
    for index = 1:numel(fasterChannelDependencies)
        if ~nf_function_available(fasterChannelDependencies{index})
            error('nf_preprocess:MissingFASTER', ...
                ['FASTER channel dependency %s was not found. Install the ' ...
                'complete FASTER distribution and its MATLAB dependencies.'], ...
                fasterChannelDependencies{index});
        end
    end
end
if strcmp(options.icaMethod, 'faster')
    fasterICAdependencies = {'component_properties', 'min_z', ...
        'hurst_exponent', 'eeg_getica', 'pwelch', 'kurt', 'nanmean'};
    for index = 1:numel(fasterICAdependencies)
        if ~nf_function_available(fasterICAdependencies{index})
            error('nf_preprocess:MissingFASTER', ...
                ['FASTER ICA dependency %s was not found. Install the ' ...
                'complete FASTER distribution and its MATLAB dependencies.'], ...
                fasterICAdependencies{index});
        end
    end
end
if any(strcmp(options.epochDetectors, 'faster'))
    fasterEpochDependencies = {'epoch_properties', 'min_z'};
    if isfield(options.fasterEpochOptions, 'localChannelDetection') && ...
            nf_is_logical_scalar( ...
            options.fasterEpochOptions.localChannelDetection) && ...
            logical(options.fasterEpochOptions.localChannelDetection)
        fasterEpochDependencies{end + 1} = ...
            'single_epoch_channel_properties';
    end
    if isfield(options.fasterEpochOptions, 'exactVendorInterpolation') && ...
            nf_is_logical_scalar( ...
            options.fasterEpochOptions.exactVendorInterpolation) && ...
            logical(options.fasterEpochOptions.exactVendorInterpolation)
        fasterEpochDependencies{end + 1} = 'h_epoch_interp_spl';
    end
    for index = 1:numel(fasterEpochDependencies)
        if ~nf_function_available(fasterEpochDependencies{index})
            error('nf_preprocess:MissingFASTER', ...
                ['FASTER epoch dependency %s was not found. Install the ' ...
                'complete FASTER distribution and its MATLAB dependencies.'], ...
                fasterEpochDependencies{index});
        end
    end
end
if strcmp(options.channelMethod, 'cleanrawdata')
    cleanRawDependencies = {'clean_flatlines', 'clean_channels'};
    for index = 1:numel(cleanRawDependencies)
        if exist(cleanRawDependencies{index}, 'file') ~= 2
            error('nf_preprocess:MissingCleanRawData', ...
                ['%s.m from the official clean_rawdata distribution is ' ...
                'required.'], cleanRawDependencies{index});
        end
    end
end
if strcmp(options.channelMethod, 'prep') && ...
        (~nf_function_available('removeTrend') || ...
        ~nf_function_available('findNoisyChannels'))
    error('nf_preprocess:MissingPREP', ...
        ['removeTrend.m and findNoisyChannels.m from the PREP ' ...
        'distribution are required for channelMethod=prep.']);
end
if strcmp(options.channelMethod, 'happeer')
    happeChannelDependencies = { ...
        'happe_detectBadChans', 'pop_clean_rawdata', 'pop_rejchan'};
    for index = 1:numel(happeChannelDependencies)
        if ~nf_function_available(happeChannelDependencies{index})
            error('nf_preprocess:MissingHAPPE', ...
                ['%s is required for channelMethod=HAPPE+ER. Install the ' ...
                'complete compatible HAPPE distribution.'], ...
                happeChannelDependencies{index});
        end
    end
end
if strcmp(options.channelMethod, 'eeglab') && ...
        ~nf_function_available('pop_rejchan')
    error('nf_preprocess:MissingEEGLABChannelFunction', ...
        'pop_rejchan.m is required for channelMethod=EEGLAB.');
end
if strcmp(options.precleanMethod, 'asr') && ...
        exist('clean_asr', 'file') ~= 2
    error('nf_preprocess:MissingCleanRawData', ...
        'clean_asr.m from the official clean_rawdata distribution is required.');
end
if strcmp(options.precleanMethod, 'prep')
    prepDependencies = {'prepPipeline', 'getPrepVersion'};
    for index = 1:numel(prepDependencies)
        if exist(prepDependencies{index}, 'file') ~= 2
            error('nf_preprocess:MissingPREP', ...
                '%s.m from the installed PREP distribution is required.', ...
                prepDependencies{index});
        end
    end
end
if strcmp(options.precleanMethod, 'happeer')
    happeDependencies = {'happe_wavThresh', 'assessPipelineStep', ...
        'calcSNR_PSNR', 'wdenoise', 'mscohere'};
    for index = 1:numel(happeDependencies)
        if ~nf_function_available(happeDependencies{index})
            error('nf_preprocess:MissingHAPPE', ...
                ['%s is required for the installed HAPPE+ER ' ...
                'wavelet stage.'], happeDependencies{index});
        end
    end
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
if ~strcmp(options.icaMethod, 'none') && ...
        strcmp(options.icaAlgorithm, 'runamica15') && ...
        (exist('pop_runamica', 'file') ~= 2 || exist('runamica15', 'file') ~= 2)
    error('nf_preprocess:MissingAMICA', ...
        'The current AMICA plugin is required for icaAlgorithm=runamica15.');
end
if ~strcmp(options.icaMethod, 'none') && ...
        strcmp(options.icaAlgorithm, 'picard')
    picardDependencies = ...
        {'picard', 'picardo', 'picard_standard', 'whitening'};
    for index = 1:numel(picardDependencies)
        if exist(picardDependencies{index}, 'file') ~= 2
            error('nf_preprocess:MissingPicard', ...
                ['The complete Picard EEGLAB plugin is required for ' ...
                'icaAlgorithm=''picard''; %s.m was not found. Install ' ...
                'Picard through the EEGLAB Extension Manager.'], ...
                picardDependencies{index});
        end
    end
end
if strcmp(options.icaMethod, 'iclabel') && exist('pop_iclabel', 'file') ~= 2
    error('nf_preprocess:MissingICLabel', ...
        'The ICLabel plugin is required for icaMethod=iclabel.');
end
if strcmp(options.icaMethod, 'adjustedadjust')
    madeDependencies = {'adjusted_ADJUST', 'compute_GD_feat', ...
        'computeSED_NOnorm', 'computeSAD', 'EM', 'trim_and_mean', ...
        'trim_and_max', 'MARA_extract_time_freq_features', ...
        'beall_horizontal', 'beall_blink_detection', 'Spatial_Info_eyes', ...
        'spectopo', 'fitlm', 'findpeaks', 'kurt'};
    for index = 1:numel(madeDependencies)
        if exist(madeDependencies{index}, 'file') ~= 2
            error('nf_preprocess:MissingAdjustedAdjust', ...
                ['MADE adjusted_ADJUST requires %s.m from the complete MADE/' ...
                'adjusted_ADJUST dependency set.'], madeDependencies{index});
        end
    end
end
if nf_save_requested(options.save) && ...
        ~nf_function_available('pop_saveset')
    error('nf_preprocess:MissingSaveFunction', ...
        'EEGLAB pop_saveset.m is required when save is requested.');
end
end

function requirements = nf_epoch_repair_requirements(options)
nativeDetectors = {'threshold', 'fft', 'peak2peak', 'step', ...
    'gradient', 'flatline', 'clipping'};
nativeRepairCapable = ...
    any(ismember(options.epochDetectors, nativeDetectors));
fasterActive = any(strcmp(options.epochDetectors, 'faster'));
fasterLocalDetection = false;
fasterVendorInterpolation = false;
if fasterActive && ...
        isfield(options.fasterEpochOptions, 'localChannelDetection')
    value = options.fasterEpochOptions.localChannelDetection;
    if ~nf_is_logical_scalar(value)
        error('nf_preprocess:InvalidFasterEpochOptions', ...
            ['fasterEpochOptions.localChannelDetection must be a logical ' ...
            'scalar.']);
    end
    fasterLocalDetection = logical(value);
end
if fasterActive && ...
        isfield(options.fasterEpochOptions, 'exactVendorInterpolation')
    value = options.fasterEpochOptions.exactVendorInterpolation;
    if ~nf_is_logical_scalar(value)
        error('nf_preprocess:InvalidFasterEpochOptions', ...
            ['fasterEpochOptions.exactVendorInterpolation must be a ' ...
            'logical scalar.']);
    end
    fasterVendorInterpolation = logical(value);
end
if fasterVendorInterpolation && ~fasterLocalDetection
    error('nf_preprocess:InvalidFasterEpochOptions', ...
        ['fasterEpochOptions.exactVendorInterpolation requires ' ...
        'localChannelDetection=true.']);
end

requirements = struct();
requirements.channelRepairCapable = ...
    nativeRepairCapable || fasterLocalDetection;
requirements.requiresInterpolation = false;
requirements.requiresEeglabInterpolation = false;
requirements.requiresGeometry = false;
if ~requirements.channelRepairCapable
    return
end

if options.localInterp
    sparseAction = 'interpolate';
else
    sparseAction = 'rejectepoch';
end
frontalAction = 'rejectepoch';
excessAction = 'rejectepoch';
overrides = options.epochRepairOptions;
if isfield(overrides, 'sparseChannelAction')
    sparseAction = nf_normalize_epoch_channel_action( ...
        overrides.sparseChannelAction, ...
        'epochRepairOptions.sparseChannelAction');
end
if isfield(overrides, 'frontalAction')
    frontalAction = nf_normalize_epoch_channel_action( ...
        overrides.frontalAction, 'epochRepairOptions.frontalAction');
end
if isfield(overrides, 'excessChannelAction')
    excessAction = nf_normalize_epoch_channel_action( ...
        overrides.excessChannelAction, ...
        'epochRepairOptions.excessChannelAction');
end
actions = {sparseAction, frontalAction, excessAction};
requirements.requiresEeglabInterpolation = ...
    any(strcmp(actions, 'interpolate'));
requirements.requiresInterpolation = ...
    requirements.requiresEeglabInterpolation || fasterVendorInterpolation;
requirements.requiresGeometry = requirements.requiresInterpolation || ...
    ~strcmp(frontalAction, sparseAction);
end

function action = nf_normalize_epoch_channel_action(value, fieldName)
if ~nf_is_text(value)
    error('nf_preprocess:InvalidEpochRepairAction', ...
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
        error('nf_preprocess:InvalidEpochRepairAction', ...
            '%s must be interpolate, rejectepoch, or mark.', fieldName);
end
end

function required = nf_pipeline_requires_geometry(options)
channelGeometryMethods = {'faster', 'cleanrawdata', 'prep', 'happeer'};
precleanGeometryMethods = {'gedai', 'prep'};
mayRestoreChannels = options.globalInterpolation && ...
    ~strcmp(options.channelMethod, 'none');
epochRepairRequirements = nf_epoch_repair_requirements(options);
required = any(strcmp(options.channelMethod, channelGeometryMethods)) || ...
    any(strcmp(options.precleanMethod, precleanGeometryMethods)) || ...
    ~strcmp(options.icaMethod, 'none') || options.qualityCompute || ...
    mayRestoreChannels || epochRepairRequirements.requiresGeometry;
end

function provenance = nf_collect_pipeline_provenance( ...
    channelInfo, precleanInfo, icaInfo, thresholdInfo)
provenance = struct();
provenance.contractVocabulary = struct();
provenance.contractVocabulary.vendorExact = ...
    'The named installed vendor entry point made the scientific decision.';
provenance.contractVocabulary.vendorPrimitive = ...
    ['Named vendor primitives ran, with NeuroFreq composing decisions or ' ...
    'stage order.'];
provenance.contractVocabulary.native = ...
    'NeuroFreq implemented the explicitly described rule.';
provenance.completePipelineEquivalenceRequiresExplicitClaim = true;
provenance.channels = nf_extract_provenance(channelInfo);
if isfield(precleanInfo, 'contract')
    provenance.preclean = precleanInfo.contract;
else
    provenance.preclean = nf_extract_provenance(precleanInfo);
end
provenance.ica = nf_extract_provenance(icaInfo);
provenance.epochs = nf_extract_provenance(thresholdInfo);
end

function provenance = nf_extract_provenance(info)
if isstruct(info) && isfield(info, 'provenance')
    provenance = info.provenance;
else
    provenance = struct();
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
software.neuroFreqPreprocessSchema = '4.4.0';
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
history.classifier = report.presetDefinition.resolvedClassifier;
history.icaAlgorithm = report.steps.ica.algorithm;
history.epochTiming = report.epochs.timing;
history.input = report.input;
history.output = report.output;
history.quality = report.quality;
history.persistence = report.persistence;
end

function values = nf_normalize_positional_inputs(values, parameterNames)
if isempty(values)
    return
end

firstValue = values{1};
if mod(numel(values), 2) == 1
    if nf_is_behavior_input(firstValue)
        values = [{'behavior', firstValue} values(2:end)];
    else
        values = [{'events', firstValue} values(2:end)];
    end
    return
end

firstIsParameterName = nf_is_text(firstValue) && ...
    any(strcmpi(char(firstValue), parameterNames));
if firstIsParameterName
    return
end

if numel(values) >= 2 && nf_is_behavior_input(values{2})
    values = [{'events', firstValue, 'behavior', values{2}} ...
        values(3:end)];
end
end

function valid = nf_is_behavior_input(value)
valid = istable(value) || ...
    (isstruct(value) && isvector(value));
end

function nf_validate_accepted_event_positions( ...
        positions, candidateCount, trialCount)
if ~isnumeric(positions) || ~isreal(positions) || ...
        (~isempty(positions) && ~isvector(positions))
    error('nf_preprocess:InvalidAcceptedEventPositions', ...
        ['pop_epoch returned invalid accepted-event positions for the ' ...
        'requested candidate list.']);
end
invalidValues = any(~isfinite(positions)) || ...
    any(positions ~= round(positions)) || ...
    any(positions < 1) || ...
    any(positions > candidateCount);
invalidOrder = any(diff(positions) <= 0);
invalidCount = numel(positions) ~= trialCount;
if invalidValues || invalidOrder || invalidCount
    error('nf_preprocess:InvalidAcceptedEventPositions', ...
        ['pop_epoch returned invalid accepted-event positions for the ' ...
        'requested candidate list.']);
end
end

function row = nf_row_pair(value)
if isempty(value)
    row = [];
else
    row = reshape(value, 1, 2);
end
end

function [behavior, info] = nf_resolve_eeg_behavior( ...
        EEG, requestedBehavior, explicitlySupplied)
behavior = struct([]);
info = struct();
info.present = false;
info.explicitlySupplied = logical(explicitlySupplied);
info.source = 'none';
info.inputEntries = 0;

hasBehav = isfield(EEG, 'etc') && isstruct(EEG.etc) && ...
    isfield(EEG.etc, 'behav') && ~isempty(EEG.etc.behav);
hasBehavior = isfield(EEG, 'etc') && isstruct(EEG.etc) && ...
    isfield(EEG.etc, 'behavior') && ~isempty(EEG.etc.behavior);

behavValue = struct([]);
behaviorValue = struct([]);
if hasBehav
    behavValue = nf_normalize_behavior_array( ...
        EEG.etc.behav, 'EEG.etc.behav');
end
if hasBehavior
    behaviorValue = nf_normalize_behavior_array( ...
        EEG.etc.behavior, 'EEG.etc.behavior');
end
if hasBehav && hasBehavior && ...
        ~isequaln(orderfields(behavValue), orderfields(behaviorValue))
    error('nf_preprocess:BehaviorFieldConflict', ...
        ['EEG.etc.behav and EEG.etc.behavior contain different trial ' ...
        'metadata.']);
end

existingBehavior = struct([]);
if hasBehav
    existingBehavior = behavValue;
    info.source = 'EEG.etc.behav';
elseif hasBehavior
    existingBehavior = behaviorValue;
    info.source = 'EEG.etc.behavior';
end

requestedPresent = explicitlySupplied && ~isempty(requestedBehavior);
if requestedPresent
    requestedBehavior = nf_normalize_behavior_array( ...
        requestedBehavior, 'behavior');
    if ~isempty(existingBehavior) && ...
            ~isequaln(orderfields(requestedBehavior), ...
            orderfields(existingBehavior))
        error('nf_preprocess:BehaviorInputConflict', ...
            ['The supplied behavior differs from behavior already stored ' ...
            'in the input EEG.']);
    end
    behavior = requestedBehavior;
    info.source = 'input behavior';
elseif ~isempty(existingBehavior)
    behavior = existingBehavior;
end

info.present = ~isempty(behavior);
info.inputEntries = numel(behavior);
if info.present
    info.fieldNames = fieldnames(behavior);
else
    info.fieldNames = {};
end
end

function behavior = nf_normalize_behavior_array(behavior, fieldName)
if isempty(behavior)
    behavior = struct([]);
    return
end
if istable(behavior)
    behavior = table2struct(behavior, 'ToScalar', false);
end
if ~isstruct(behavior) || ~isvector(behavior)
    error('nf_preprocess:InvalidBehavior', ...
        ['%s must be a NeuroFreq behavior struct array with one element ' ...
        'per task trial.'], fieldName);
end
behavior = reshape(behavior, 1, []);
end

function nf_validate_event_behavior_count(behavior, eventCount)
if numel(behavior) ~= eventCount
    error('nf_preprocess:BehaviorEventMismatch', ...
        ['Behavior has %d entries, but the requested event types match %d ' ...
        'source events. Supply one behavior entry per requested event in ' ...
        'in increasing event-latency order.'], ...
        numel(behavior), eventCount);
end
end

function EEG = nf_set_eeg_behavior(EEG, behavior)
if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
behavior = reshape(behavior, 1, []);
EEG.etc.behav = behavior;
EEG.etc.behavior = behavior;
end

function EEG = nf_clear_eeg_behavior(EEG)
if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc)
    return
end
fields = {'behav', 'behavior'};
for fieldIndex = 1:numel(fields)
    if isfield(EEG.etc, fields{fieldIndex})
        EEG.etc = rmfield(EEG.etc, fields{fieldIndex});
    end
end
end

function summary = nf_behavior_option_summary(behavior)
summary = struct();
summary.provided = ~isempty(behavior);
summary.nEntries = numel(behavior);
summary.format = 'NeuroFreq struct array';
if isempty(behavior)
    summary.fieldNames = {};
else
    summary.fieldNames = fieldnames(behavior);
end
end

function reported = nf_report_options(options)
reported = options;
if isfield(reported, 'behavior')
    reported.behavior = nf_behavior_option_summary(reported.behavior);
end
if isfield(reported, 'eventValidation')
    reported = rmfield(reported, 'eventValidation');
end
if isfield(reported, 'gedaiOptions') && ...
        isstruct(reported.gedaiOptions) && ...
        isfield(reported.gedaiOptions, 'referenceMatrixType')
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
allowed = {'measure', 'z'};
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
% nf_validate_event_latencies(EEG);
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
% nf_validate_event_latencies(EEG);
end

% function nf_validate_event_latencies(EEG)
% if ~isfield(EEG, 'event') || isempty(EEG.event)
%     return
% end
% if ~isfield(EEG.event, 'latency')
%     error('nf_preprocess:InvalidEvents', ...
%         'Every EEG.event entry must contain a latency field.');
% end
% for index = 1:numel(EEG.event)
%     latency = EEG.event(index).latency;
%     if ~isnumeric(latency) || ~isreal(latency) || ~isscalar(latency) || ...
%             ~isfinite(latency) || latency < 0.5 || latency > EEG.pnts + 0.5
%         error('nf_preprocess:InvalidEvents', ...
%             'EEG.event(%d).latency is outside the continuous sample range.', index);
%     end
%     if isfield(EEG.event, 'urevent') && ~isempty(EEG.event(index).urevent)
%         ureventIndex = EEG.event(index).urevent;
%         if ~isnumeric(ureventIndex) || ~isscalar(ureventIndex) || ...
%                 ~isfinite(ureventIndex) || ureventIndex ~= round(ureventIndex) || ...
%                 ureventIndex < 1 || ~isfield(EEG, 'urevent') || ...
%                 ureventIndex > numel(EEG.urevent)
%             error('nf_preprocess:InvalidEvents', ...
%                 'EEG.event(%d).urevent is not a valid urevent index.', index);
%         end
%     end
% end
% if isfield(EEG, 'urevent') && ~isempty(EEG.urevent)
%     if ~isstruct(EEG.urevent) || ~isfield(EEG.urevent, 'latency')
%         error('nf_preprocess:InvalidEvents', ...
%             'EEG.urevent must contain latency fields.');
%     end
%     for index = 1:numel(EEG.urevent)
%         latency = EEG.urevent(index).latency;
%         if ~isnumeric(latency) || ~isreal(latency) || ...
%                 ~isscalar(latency) || ~isfinite(latency)
%             error('nf_preprocess:InvalidEvents', ...
%                 'EEG.urevent(%d).latency must be a finite scalar.', index);
%         end
%     end
% end
% end

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

function valid = nf_is_method_list(value)
valid = isempty(value) || ischar(value) || isstring(value) || ...
    iscellstr(value);
end

function valid = nf_is_frequency_bands(value)
if isempty(value)
    valid = true;
    return
end
if isnumeric(value)
    valid = isreal(value) && ...
        (size(value, 2) == 2 || size(value, 2) == 4 || ...
        (isvector(value) && ismember(numel(value), [2 4]))) && ...
        all(isfinite(value(:)));
    if ~valid
        return
    end
    if isvector(value)
        if numel(value) == 2
            value = reshape(value, 1, 2);
        else
            value = reshape(value, 1, 4);
        end
    end
    valid = all(value(:, 1) >= 0) && ...
        all(value(:, 1) < value(:, 2));
    if size(value, 2) == 4
        valid = valid && all(value(:, 3) < value(:, 4));
    end
    return
end
if ~isstruct(value) || isempty(value)
    valid = false;
    return
end
valid = true;
allowed = {'frequencyRangeHz', 'frequencyRange', ...
    'powerThresholdDb', 'powerThreshold'};
for index = 1:numel(value)
    names = fieldnames(value(index));
    if any(~ismember(names, allowed))
        valid = false;
        return
    end
    if isfield(value(index), 'frequencyRangeHz')
        frequency = value(index).frequencyRangeHz;
    elseif isfield(value(index), 'frequencyRange')
        frequency = value(index).frequencyRange;
    else
        valid = false;
        return
    end
    if ~nf_is_increasing_pair(frequency) || frequency(1) < 0
        valid = false;
        return
    end
    if isfield(value(index), 'powerThresholdDb')
        power = value(index).powerThresholdDb;
    elseif isfield(value(index), 'powerThreshold')
        power = value(index).powerThreshold;
    else
        power = [];
    end
    if ~isempty(power) && ~nf_is_increasing_pair(power)
        valid = false;
        return
    end
end
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
