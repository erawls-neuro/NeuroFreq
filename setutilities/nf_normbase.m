function TF = nf_normbase(TF, varargin)
% NF_NORMBASE  Average, decompose, normalize, and baseline TF data.
%
% GENERAL
% -------
% NF_NORMBASE averages trials within user-defined groups, optionally uses
% the concurrent phaser method (CPM) to estimate phase-locked and
% non-phase-locked ensemble power, and applies the single-trial
% normalization procedures described by Grandchamp and Delorme (2011).
% The three single-trial methods are:
%
%   'zscore'  Per-trial temporal centering and sample-SD scaling.
%   'percent' Per-trial division by temporal mean (the gain model).
%   'db'      The same gain normalization as 'percent', with 10*log10
%             applied only after trial aggregation and final baselining.
%
% Thus, 'percent' and 'db' share one pre-averaging gain operation. Trialwise
% logarithms are never averaged. By default, normalization uses the entire
% epoch and BaselineMethod='auto' applies the corresponding classical
% prestimulus correction after averaging. An explicit BaselineMethod='none'
% returns the intermediate normalized ERS.
%
% AverageMethod='cpm' runs CPM on complex linear coefficients before any
% group baseline. It returns TF.power, TF.power_phase_locked, and
% TF.power_non_phase_locked plus raw CPM parameters and diagnostics in
% TF.cpm. Power-domain trial z-normalization cannot define coefficient
% magnitudes for CPM and is rejected. Trial gain normalization ('percent'
% or 'db') is permitted as
% an experimental extension because it preserves phase, but it changes the
% estimand and may violate CPM's weighted-phase assumptions. The function
% warns unless Verbosity='silent'; full-epoch gain receives a stronger
% warning because event activity enters its own normalizer.
% Post-CPM baselines use each output map's own finite reference samples.
% The raw unbaselined CPM moments are additive when valid; no additivity
% constraint is imposed after baseline correction.
%
% Power arithmetic is performed in linear power even when the input uses
% TF.scale='log10'. If neither normalization nor baseline correction is
% requested, ordinary averaging preserves the input scale. Raw CPM output
% is linear because CPM is defined on complex linear coefficients.
% Normalized or baselined output uses TF.scale='normalized'. TF.normbase
% records the exact measure, units, input scale, windows, operation order,
% grouping, scientific cautions, and publication provenance.
%
% Every compatible top-level numeric phase field is reduced to inter-trial
% phase coherence (ITPC). Input phase and phase_* fields become output itpc
% and itpc_* fields so phase angles are never mislabeled as coherence.
% Undefined phase values stored as NaN are omitted; zero radians is valid.
% Single-trial NF_MODELERP complex signals and their now-stale model ledger
% are removed after averaging.
%
% USAGE
% -----
%   TF = nf_normbase(TF)
%   TF = nf_normbase(TF, ...
%       'SingleTrialNormalization', 'zscore', ...
%       'NormalizationTimes', 'epoch', ...
%       'BaselineMethod', 'auto', ...
%       'BaselineTimes', [-0.5 -0.1], ...
%       'TrialGroups', conditionVector, ...
%       'AverageMethod', 'cpm', ...
%       'Verbosity', 'auto')
%
% OPTIONS
% -------
% SingleTrialNormalization
%       'none' (default), 'zscore', 'percent', or 'db'. 'percent' and
%       'db' use P/mean(P) within every trial.
%
% NormalizationTimes
%       'epoch' (default) or [minimum maximum] in TF.times units. 'epoch'
%       implements full-epoch normalization. A numeric interval also
%       supports the paper's single-trial prestimulus procedures.
%
% BaselineMethod
%       'auto' (default), 'none', 'subtract', 'zscore', 'ratio', 'percent',
%       or 'db'.
%       With full-epoch normalization, auto maps zscore -> zscore,
%       percent -> ratio, db -> db, and none -> none. With a numeric
%       single-trial normalization interval, auto chooses none because the
%       interval itself is already the trialwise baseline. 'subtract'
%       subtracts the condition-level baseline mean. 'ratio' is centered at
%       one and 'percent' returns 100*(ratio-1).
%
% BaselineTimes
%       Empty (default) uses all TF.times <= 0 whenever a group baseline
%       is active. Otherwise supply inclusive [minimum maximum] limits in
%       TF.times units. No baseline samples are required for 'none'.
%
% TrialGroups
%       Empty (default) averages all trials together. Otherwise supply one
%       finite numeric or logical label per trial. Conditions are returned
%       in sorted label order and recorded in TF.normbase.conditions.
%
% AverageMethod
%       'mean' (default), 'median', or 'cpm'. CPM requires TF.phase and at
%       least two trials in every condition. It decomposes only TF.power;
%       other power fields are arithmetically averaged. CPM always uses
%       arithmetic moments and requires a verified linear complex transform
%       method or the explicit analyticCoefficientSemantics token. Known
%       bilinear/RID distributions are rejected.
%
% Verbosity
%       'auto' (default), 'verbose', or 'silent'. Auto reports warnings and
%       one concise completion summary. Verbose additionally reports the
%       resolved configuration, input dimensions, detected fields, removed
%       stale fields, and condition-level progress. Silent emits no status
%       text and suppresses function-owned warnings; errors always throw.
%
% All optional inputs are strict name-value pairs. Positional syntax is not
% supported.
%
% REFERENCE
% ---------
% Grandchamp R, Delorme A. Front Psychol. 2011;2:236.
% doi:10.3389/fpsyg.2011.00236
%
% Singhal S, Ghosh P, Kumar N, Banerjee A. J Neurophysiol.
% 2023;129(1):199-210. doi:10.1152/jn.00467.2022

if nargin < 1
    error('nf_normbase:TFRequired', ...
        'A nonempty TF structure is required.');
end

if isempty(TF) == true
    error('nf_normbase:TFRequired', ...
        'A nonempty TF structure is required.');
end

if isstruct(TF) == false || isscalar(TF) == false
    error('nf_normbase:InvalidTF', ...
        'TF must be a scalar structure.');
end

options = local_parse_inputs(TF, varargin{:});
isCPM = strcmp(options.AverageMethod, 'cpm');
options.ScientificCautions = local_cpm_scientific_cautions(options);

for cautionIndex = 1:numel(options.ScientificCautions)
    caution = options.ScientificCautions(cautionIndex);
    local_issue_warning( ...
        options.Verbosity, ...
        caution.identifier, ...
        '%s', ...
        caution.message);
end

local_report( ...
    options.Verbosity, ...
    'verbose', ...
    'Starting input validation.');

requiredFields = {'power', 'times', 'freqs', 'scale', 'ntrls', 'nsensor'};

for fieldIndex = 1:numel(requiredFields)
    fieldName = requiredFields{fieldIndex};

    if isfield(TF, fieldName) == false
        error('nf_normbase:RequiredTFField', ...
            'TF.%s is required.', fieldName);
    end
end

if isnumeric(TF.power) == false
    error('nf_normbase:InvalidPowerField', ...
        'TF.power must be a real floating-point numeric array.');
end

derivedFields = { ...
    'normbase', ...
    'cpm', ...
    'power_phase_locked', ...
    'power_non_phase_locked'};
inputFieldNames = fieldnames(TF);

for fieldIndex = 1:numel(derivedFields)
    fieldName = derivedFields{fieldIndex};
    conflictIndex = find(strcmpi(inputFieldNames, fieldName), 1);

    if isempty(conflictIndex) == false
        conflictName = inputFieldNames{conflictIndex};
        error('nf_normbase:DerivedOutputFieldConflict', ...
            ['TF.%s is an NF_NORMBASE output field. Supply unprocessed ' ...
            'single-trial TF data.'], ...
            conflictName);
    end
end

inputScale = lower(local_text_scalar(TF.scale, 'TF.scale'));

if isfield(TF, 'analyticCoefficientSemantics') == true
    inputAnalyticCoefficientSemantics = local_text_scalar( ...
        TF.analyticCoefficientSemantics, ...
        'TF.analyticCoefficientSemantics');
else
    inputAnalyticCoefficientSemantics = '';
end

if any(strcmp(inputScale, {'linear', 'log10'})) == false
    error('nf_normbase:TFScaleNotRecognized', ...
        'Input TF.scale must be ''linear'' or ''log10''.');
end

nSensors = local_positive_integer(TF.nsensor, 'TF.nsensor');
nTrials = local_positive_integer(TF.ntrls, 'TF.ntrls');

if nTrials == 1 && isfield(TF, 'conds') == true
    error('nf_normbase:TFAlreadyAveraged', ...
        'TF already contains condition averages.');
end

times = local_finite_axis(TF.times, 'TF.times');
frequencies = local_finite_axis(TF.freqs, 'TF.freqs');
nFrequencies = numel(frequencies);
nTimes = numel(times);
local_report( ...
    options.Verbosity, ...
    'verbose', ...
    ['Input: %d sensor(s), %d frequency bin(s), %d time point(s), ' ...
    '%d trial(s), scale=%s.'], ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nTrials, ...
    inputScale);

if strcmp(options.SingleTrialNormalization, 'none') == true
    normalizationIndex = [];
    normalizationTimesApplied = [];
else
    [normalizationIndex, normalizationTimesApplied] = ...
        local_normalization_index(times, options.NormalizationTimes);
end

[baselineIndex, baselineTimesApplied] = local_baseline_index( ...
    times, options.BaselineTimes, options.BaselineMethodEffective);

if strcmp(options.SingleTrialNormalization, 'zscore') == true && ...
        numel(normalizationIndex) < 2
    error('nf_normbase:InsufficientNormalizationSamples', ...
        'Single-trial z-normalization requires at least two time points.');
end

if strcmp(options.BaselineMethodEffective, 'zscore') == true && ...
        strcmp(options.SingleTrialNormalization, 'zscore') == false && ...
        numel(baselineIndex) < 2
    error('nf_normbase:InsufficientBaselineSamples', ...
        'Z-score baselining requires at least two baseline time points.');
end

trialGroups = local_trial_groups(options.TrialGroups, nTrials);
conditions = unique(trialGroups);
nConditions = numel(conditions);
local_report( ...
    options.Verbosity, ...
    'verbose', ...
    ['Configuration: normalization=%s, baseline=%s, averaging=%s, ' ...
    'conditions=%d.'], ...
    options.SingleTrialNormalization, ...
    options.BaselineMethodEffective, ...
    options.AverageMethod, ...
    nConditions);
local_report( ...
    options.Verbosity, ...
    'verbose', ...
    'Selected %d normalization sample(s) and %d baseline sample(s).', ...
    numel(normalizationIndex), ...
    numel(baselineIndex));

staleFields = {'complex_erp', 'complex_residual', 'erpmodel'};

for fieldIndex = 1:numel(staleFields)
    fieldName = staleFields{fieldIndex};

    if isfield(TF, fieldName) == true
        TF = rmfield(TF, fieldName);
        local_report( ...
            options.Verbosity, ...
            'verbose', ...
            'Removed stale field TF.%s.', ...
            fieldName);
    end
end

powerFields = nf_find_power_fields_local(TF);
phaseFields = nf_find_phase_fields_local(TF);
newPower = struct();
newITPC = struct();
sourcePower = struct();
sourcePhase = struct();
phaseOutputFields = cell(size(phaseFields));
local_report( ...
    options.Verbosity, ...
    'verbose', ...
    'Detected %d power field(s) and %d phase field(s).', ...
    numel(powerFields), ...
    numel(phaseFields));

if isCPM == true
    local_report( ...
        options.Verbosity, ...
        'verbose', ...
        ['CPM will decompose TF.power using TF.phase; ancillary power ' ...
        'fields will use arithmetic means.']);
end

for fieldIndex = 1:numel(powerFields)
    fieldName = powerFields{fieldIndex};
    values = TF.(fieldName);
    values4 = nf_as_canonical4_local( ...
        values, fieldName, nSensors, nFrequencies, nTimes, nTrials);
    local_validate_power_values(values4, inputScale, fieldName);
    sourcePower.(fieldName) = values4;
    newPower.(fieldName) = zeros( ...
        nSensors, nFrequencies, nTimes, nConditions, 'like', values);
end

for fieldIndex = 1:numel(phaseFields)
    fieldName = phaseFields{fieldIndex};
    outputFieldName = local_itpc_field_name(fieldName);
    values = TF.(fieldName);
    conflictIndex = find(strcmpi( ...
        fieldnames(TF), outputFieldName), 1);

    if isempty(conflictIndex) == false
        conflictFields = fieldnames(TF);
        conflictName = conflictFields{conflictIndex};
        error('nf_normbase:ITPCOutputFieldConflict', ...
            ['TF.%s already exists and would conflict with the ITPC ' ...
            'derived from TF.%s.'], ...
            conflictName, ...
            fieldName);
    end

    if isfloat(values) == false || isreal(values) == false
        error('nf_normbase:InvalidPhaseField', ...
            'Phase field %s must be a real floating-point array.', ...
            fieldName);
    end

    if any(isinf(values(:))) == true
        error('nf_normbase:InfinitePhase', ...
            ['Phase field %s may contain NaN for missing phase, but not ' ...
            'positive or negative infinity.'], ...
            fieldName);
    end

    if numel(values) ~= numel(TF.power)
        error('nf_normbase:PhaseSizeMismatch', ...
            ['Phase field %s must contain one value for every value in ' ...
            'TF.power.'], fieldName);
    end

    values4 = nf_as_canonical4_local( ...
        values, fieldName, nSensors, nFrequencies, nTimes, nTrials);
    sourcePhase.(fieldName) = values4;
    phaseOutputFields{fieldIndex} = outputFieldName;
    newITPC.(outputFieldName) = zeros( ...
        nSensors, nFrequencies, nTimes, nConditions, 'like', values);
    newITPC.(outputFieldName)(:) = NaN;
end

if isCPM == true
    primaryPhaseIndex = find(strcmpi(phaseFields, 'phase'), 1);

    if isempty(primaryPhaseIndex) == true
        error('nf_normbase:CPMRequiresPhase', ...
            'AverageMethod=''cpm'' requires TF.phase.');
    end

    options.CPMPhaseField = phaseFields{primaryPhaseIndex};

    options.CPMTransform = ...
        local_validate_cpm_transform_semantics(TF);

    newPower.power_phase_locked = zeros( ...
        nSensors, ...
        nFrequencies, ...
        nTimes, ...
        nConditions, ...
        'like', ...
        TF.power);
    newPower.power_non_phase_locked = zeros( ...
        nSensors, ...
        nFrequencies, ...
        nTimes, ...
        nConditions, ...
        'like', ...
        TF.power);
    cpmModels = cell(1, nConditions);
    cpmBaselineDiagnostics = repmat( ...
        local_empty_cpm_baseline_diagnostic(), ...
        1, ...
        nConditions);
end

if isfield(TF, 'behavior') == true
    if isstruct(TF.behavior) == false || numel(TF.behavior) ~= nTrials
        error('nf_normbase:BehaviorLengthMismatch', ...
            'TF.behavior must contain one structure per trial.');
    end

    newBehavior = TF.behavior(1:nConditions);
end

trialCounts = zeros(1, nConditions);

for conditionIndex = 1:nConditions
    trialIndex = find(trialGroups == conditions(conditionIndex));
    trialCounts(conditionIndex) = numel(trialIndex);
    local_report( ...
        options.Verbosity, ...
        'verbose', ...
        'Processing condition %d/%d (label=%g, trials=%d).', ...
        conditionIndex, ...
        nConditions, ...
        double(conditions(conditionIndex)), ...
        trialCounts(conditionIndex));

    if isscalar(trialIndex) == true
        if isCPM == true
            error('nf_normbase:CPMRequiresMultipleTrials', ...
                ['Condition %g contains one trial. CPM requires at ' ...
                'least two trials in every condition.'], ...
                double(conditions(conditionIndex)));
        else
            local_issue_warning( ...
                options.Verbosity, ...
                'nf_normbase:SingleTrialCondition', ...
                'Condition %g contains only one trial.', ...
                double(conditions(conditionIndex)));
        end
    end

    if isCPM == true
        conditionPower = sourcePower.power(:, :, :, trialIndex);
        primaryPhase = sourcePhase.(options.CPMPhaseField);
        conditionPhase = primaryPhase(:, :, :, trialIndex);
        [processedTotal, processedPhaseLocked, ...
            processedNonPhaseLocked, cpmModel, baselineDiagnostic] = ...
            local_process_cpm_condition( ...
            conditionPower, ...
            conditionPhase, ...
            inputScale, ...
            options, ...
            normalizationIndex, ...
            baselineIndex);
        newPower.power(:, :, :, conditionIndex) = reshape( ...
            processedTotal, ...
            [nSensors, nFrequencies, nTimes, 1]);
        newPower.power_phase_locked(:, :, :, conditionIndex) = ...
            reshape( ...
            processedPhaseLocked, ...
            [nSensors, nFrequencies, nTimes, 1]);
        newPower.power_non_phase_locked(:, :, :, conditionIndex) = ...
            reshape( ...
            processedNonPhaseLocked, ...
            [nSensors, nFrequencies, nTimes, 1]);
        cpmModels{conditionIndex} = cpmModel;
        cpmBaselineDiagnostics(conditionIndex) = baselineDiagnostic;

        ordinaryOptions = options;
        ordinaryOptions.AverageMethod = 'mean';

        for fieldIndex = 1:numel(powerFields)
            fieldName = powerFields{fieldIndex};

            if strcmp(fieldName, 'power') == true
                continue
            end

            values4 = sourcePower.(fieldName);
            conditionAncillaryPower = values4(:, :, :, trialIndex);

            if strcmp( ...
                    options.SingleTrialNormalization, ...
                    'none') == true && ...
                    strcmp( ...
                    options.BaselineMethodEffective, ...
                    'none') == true
                conditionAncillaryPower = local_cpm_linear_power( ...
                    conditionAncillaryPower, ...
                    inputScale);
                ancillaryInputScale = 'linear';
            else
                ancillaryInputScale = inputScale;
            end

            processedPower = local_process_power( ...
                conditionAncillaryPower, ...
                ancillaryInputScale, ...
                ordinaryOptions, ...
                normalizationIndex, ...
                baselineIndex, ...
                fieldName);
            newPower.(fieldName)(:, :, :, conditionIndex) = ...
                reshape( ...
                processedPower, ...
                [nSensors, nFrequencies, nTimes, 1]);
        end
    else
        for fieldIndex = 1:numel(powerFields)
            fieldName = powerFields{fieldIndex};
            values4 = sourcePower.(fieldName);
            conditionPower = values4(:, :, :, trialIndex);
            processedPower = local_process_power( ...
                conditionPower, inputScale, options, ...
                normalizationIndex, baselineIndex, fieldName);
            newPower.(fieldName)(:, :, :, conditionIndex) = ...
                reshape(processedPower, ...
                [nSensors, nFrequencies, nTimes, 1]);
        end
    end

    for fieldIndex = 1:numel(phaseFields)
        fieldName = phaseFields{fieldIndex};
        outputFieldName = phaseOutputFields{fieldIndex};
        values4 = sourcePhase.(fieldName);
        conditionPhase = values4(:, :, :, trialIndex);
        conditionITPC = nf_phase_itpc_local(conditionPhase, 4);
        newITPC.(outputFieldName)(:, :, :, conditionIndex) = ...
            reshape(conditionITPC, ...
            [nSensors, nFrequencies, nTimes, 1]);
    end

    if isfield(TF, 'behavior') == true
        newBehavior(conditionIndex) = ...
            local_average_behavior(TF.behavior(trialIndex));
    end
end

TF.ntrls = 1;
TF.conds = nConditions;
TF.trlerp = trialCounts;

transformationApplied = ...
    strcmp(options.SingleTrialNormalization, 'none') == false || ...
    strcmp(options.BaselineMethodEffective, 'none') == false;

if isCPM == true && transformationApplied == false
    TF.scale = 'linear';
elseif transformationApplied == true
    TF.scale = 'normalized';
else
    TF.scale = inputScale;
end

outputPowerFields = fieldnames(newPower);

for fieldIndex = 1:numel(outputPowerFields)
    fieldName = outputPowerFields{fieldIndex};
    TF.(fieldName) = local_output_layout( ...
        newPower.(fieldName), nSensors, nFrequencies, ...
        nTimes, nConditions);
end

for fieldIndex = 1:numel(phaseFields)
    fieldName = phaseFields{fieldIndex};
    outputFieldName = phaseOutputFields{fieldIndex};
    TF = rmfield(TF, fieldName);
    TF.(outputFieldName) = local_output_layout( ...
        newITPC.(outputFieldName), nSensors, nFrequencies, ...
        nTimes, nConditions);
end

if isfield(TF, 'behavior') == true
    TF.behavior = newBehavior;
end

if isCPM == true
    TF.cpm = local_pack_cpm_models( ...
        cpmModels, ...
        cpmBaselineDiagnostics, ...
        nSensors, ...
        nFrequencies, ...
        nTimes, ...
        nConditions, ...
        conditions, ...
        trialCounts, ...
        options);
    local_report_cpm_validity( ...
        TF.cpm, ...
        options.Verbosity);

    if TF.cpm.additivity.outputPreserved == false
        local_report( ...
            options.Verbosity, ...
            'verbose', ...
            ['The selected post-CPM transformation does not guarantee ' ...
            'total=phase-locked+non-phase-locked. The unbaselined valid ' ...
            'CPM moments in TF.cpm retain the raw identity.']);
    end
end

removableFields = { ...
    'event', ...
    'epoch', ...
    'analyticCoefficientSemantics'};

for fieldIndex = 1:numel(removableFields)
    fieldName = removableFields{fieldIndex};

    if isfield(TF, fieldName) == true
        TF = rmfield(TF, fieldName);
    end
end

TF.normbase = local_build_ledger( ...
    options, inputScale, TF.scale, normalizationIndex, ...
    normalizationTimesApplied, baselineIndex, baselineTimesApplied, ...
    conditions, trialCounts, transformationApplied);
TF.normbase.inputAnalyticCoefficientSemantics = ...
    inputAnalyticCoefficientSemantics;
TF.normbase.phaseAggregation = struct();
TF.normbase.phaseAggregation.inputFields = phaseFields;
TF.normbase.phaseAggregation.outputFields = phaseOutputFields;
TF.normbase.phaseAggregation.measure = ...
    'inter_trial_phase_coherence';
TF.normbase.phaseAggregation.formula = ...
    'abs(mean(exp(1i*phase),trials))';
TF = nf_powerfront_local(TF);
local_report( ...
    options.Verbosity, ...
    'auto', ...
    ['Completed: %d trial(s) -> %d condition average(s); ' ...
    'normalization=%s, baseline=%s, averaging=%s, scale=%s.'], ...
    nTrials, ...
    nConditions, ...
    options.SingleTrialNormalization, ...
    options.BaselineMethodEffective, ...
    options.AverageMethod, ...
    TF.scale);

end


function local_report(verbosity, level, formatSpec, varargin)

if strcmp(verbosity, 'silent') == true
    return
end

if strcmp(level, 'verbose') == true && ...
        strcmp(verbosity, 'verbose') == false
    return
end

fprintf('[nf_normbase] ');
fprintf([formatSpec, '\n'], varargin{:});

end


function local_issue_warning(verbosity, identifier, message, varargin)

if strcmp(verbosity, 'silent') == true
    return
end

warning(identifier, message, varargin{:});

end


function [processedTotal, processedPhaseLocked, ...
    processedNonPhaseLocked, model, diagnostic] = ...
    local_process_cpm_condition( ...
    power4, phase4, inputScale, options, normalizationIndex, baselineIndex)

if strcmp(options.SingleTrialNormalization, 'none') == true
    linearPower = local_cpm_linear_power(power4, inputScale);
else
    scaledTrials = local_scale_trials(power4, inputScale);

    if strcmp(inputScale, 'log10') == true && ...
            any(isfinite(power4(:)) & scaledTrials(:) == 0)
        error('nf_normbase:CPMLinearPowerRange', ...
            ['Finite log10 power underflowed while constructing ' ...
            'trial-gain-normalized CPM coefficients. Rescale the input ' ...
            'power units or reduce its dynamic range.']);
    end

    linearPower = local_single_trial_normalize( ...
        scaledTrials, ...
        normalizationIndex, ...
        options.SingleTrialNormalization, ...
        'power');
end

finitePhase = isfinite(phase4);
coefficients = sqrt(linearPower) .* exp(1i .* phase4);
zeroWithoutPhase = finitePhase == false & linearPower == 0;
coefficients(zeroWithoutPhase) = 0;
coefficients(finitePhase == false & zeroWithoutPhase == false) = NaN;
model = nf_cpm( ...
    coefficients, ...
    'TrialDimension', 4, ...
    'Verbosity', 'silent');

[processedTotal, processedPhaseLocked, ...
    processedNonPhaseLocked, diagnostic] = ...
    local_finish_cpm_components( ...
    model.total.powerMean, ...
    model.phaseLocked.powerMean, ...
    model.nonPhaseLocked.powerMean, ...
    baselineIndex, ...
    options);

processedTotal = cast(processedTotal, 'like', power4);
processedPhaseLocked = cast(processedPhaseLocked, 'like', power4);
processedNonPhaseLocked = cast(processedNonPhaseLocked, 'like', power4);

end


function linearPower = local_cpm_linear_power(power4, inputScale)

if strcmp(inputScale, 'linear') == true
    linearPower = power4;
else
    linearPower = 10 .^ power4;
    finiteLogPower = isfinite(power4);
    rangeFailure = finiteLogPower & ...
        (isfinite(linearPower) == false | linearPower == 0);

    if any(rangeFailure(:)) || any(isnan(linearPower(:)))
        error('nf_normbase:CPMLinearPowerRange', ...
            ['Finite log10 power cannot be represented as positive ' ...
            'finite linear power for CPM. Rescale the input power units ' ...
            'or reduce its dynamic range.']);
    end
end

end


function verification = local_validate_cpm_transform_semantics(TF)

requiredSemantics = ...
    'power_abs_squared_phase_angle_linear_coefficient';
hasVerifiedSemantics = false;
hasVerifiedMethod = false;
method = '';

if isfield(TF, 'analyticCoefficientSemantics') == true
    semantics = local_text_scalar( ...
        TF.analyticCoefficientSemantics, ...
        'TF.analyticCoefficientSemantics');

    if strcmpi(semantics, requiredSemantics) == false
        error('nf_normbase:UnverifiedTransform', ...
            ['TF.analyticCoefficientSemantics must equal ''%s'' for ' ...
            'CPM.'], ...
            requiredSemantics);
    end

    hasVerifiedSemantics = true;
end

if isfield(TF, 'method') == true
    method = local_text_scalar(TF.method, 'TF.method');
    methodKey = lower(method);
    methodKey = strrep(methodKey, '_', '');
    methodKey = strrep(methodKey, '-', '');
    methodKey = strrep(methodKey, ' ', '');

    unsupportedMethods = { ...
        'ridrihaczek', ...
        'rihaczek', ...
        'binomial2', ...
        'ridbinomial', ...
        'bornjordan2', ...
        'ridbornjordan'};

    if any(strcmp(methodKey, unsupportedMethods)) == true
        error('nf_normbase:NonanalyticTransform', ...
            ['TF.method ''%s'' is not a linear complex-coefficient ' ...
            'transform and cannot be decomposed with CPM.'], ...
            method);
    end

    supportedMethods = { ...
        'stft', ...
        'stransform', ...
        'demodulation', ...
        'wavelet', ...
        'dcwt', ...
        'cwt', ...
        'filterhilbert'};
    hasVerifiedMethod = any(strcmp(methodKey, supportedMethods));

    if hasVerifiedMethod == false && ...
            hasVerifiedSemantics == false
        error('nf_normbase:UnverifiedTransform', ...
            ['TF.method ''%s'' is not a verified NeuroFreq linear ' ...
            'complex transform. Supply the explicit ' ...
            'analyticCoefficientSemantics token only when TF.power and ' ...
            'TF.phase encode one linear complex coefficient.'], ...
            method);
    end
elseif hasVerifiedSemantics == false
    error('nf_normbase:UnverifiedTransform', ...
        ['CPM requires TF.method for a verified NeuroFreq linear ' ...
        'complex transform or the explicit ' ...
        'TF.analyticCoefficientSemantics token.']);
end

verification = struct();
verification.requiredSemantics = requiredSemantics;
verification.inputMethod = method;

if hasVerifiedMethod == true && hasVerifiedSemantics == true
    verification.verifiedBy = 'method_and_explicit_semantics';
elseif hasVerifiedMethod == true
    verification.verifiedBy = 'verified_neurofreq_method';
else
    verification.verifiedBy = 'explicit_semantics';
end

end


function [processedTotal, processedPhaseLocked, ...
    processedNonPhaseLocked, diagnostic] = ...
    local_finish_cpm_components( ...
    totalPower, phaseLockedPower, nonPhaseLockedPower, ...
    baselineIndex, options)

baselineMethod = options.BaselineMethodEffective;
normalizationMethod = options.SingleTrialNormalization;
diagnostic = local_empty_cpm_baseline_diagnostic();

if strcmp(baselineMethod, 'none') == true
    processedTotal = totalPower;
    processedPhaseLocked = phaseLockedPower;
    processedNonPhaseLocked = nonPhaseLockedPower;

    if strcmp(normalizationMethod, 'db') == true
        processedTotal = local_cpm_log10_power(processedTotal);
        processedPhaseLocked = ...
            local_cpm_log10_power(processedPhaseLocked);
        processedNonPhaseLocked = ...
            local_cpm_log10_power(processedNonPhaseLocked);
    end

    diagnostic.totalInvalidOutputCount = ...
        local_count_invalid_cpm_output(processedTotal);
    diagnostic.phaseLockedInvalidOutputCount = ...
        local_count_invalid_cpm_output(processedPhaseLocked);
    diagnostic.nonPhaseLockedInvalidOutputCount = ...
        local_count_invalid_cpm_output(processedNonPhaseLocked);
    return
end

totalBaseline = totalPower(:, :, baselineIndex);
phaseLockedBaseline = phaseLockedPower(:, :, baselineIndex);
nonPhaseLockedBaseline = nonPhaseLockedPower(:, :, baselineIndex);
totalValid = isfinite(totalBaseline);
phaseLockedValid = isfinite(phaseLockedBaseline);
nonPhaseLockedValid = isfinite(nonPhaseLockedBaseline);
diagnostic.totalExcludedBaselineSampleCount = ...
    sum(totalValid(:) == false);
diagnostic.phaseLockedExcludedBaselineSampleCount = ...
    sum(phaseLockedValid(:) == false);
diagnostic.nonPhaseLockedExcludedBaselineSampleCount = ...
    sum(nonPhaseLockedValid(:) == false);

if strcmp(baselineMethod, 'subtract') == true
    [totalReference, totalReferenceCount] = local_masked_mean( ...
        totalBaseline, ...
        totalValid, ...
        3);
    [phaseLockedReference, phaseLockedReferenceCount] = ...
        local_masked_mean( ...
        phaseLockedBaseline, ...
        phaseLockedValid, ...
        3);
    [nonPhaseLockedReference, nonPhaseLockedReferenceCount] = ...
        local_masked_mean( ...
        nonPhaseLockedBaseline, ...
        nonPhaseLockedValid, ...
        3);
    totalReferenceValid = totalReferenceCount > 0 & ...
        isfinite(totalReference);
    phaseLockedReferenceValid = phaseLockedReferenceCount > 0 & ...
        isfinite(phaseLockedReference);
    nonPhaseLockedReferenceValid = ...
        nonPhaseLockedReferenceCount > 0 & ...
        isfinite(nonPhaseLockedReference);
    commonReferenceValid = ...
        totalReferenceValid & ...
        phaseLockedReferenceValid & ...
        nonPhaseLockedReferenceValid;
    totalReference(totalReferenceValid == false) = NaN;
    phaseLockedReference(phaseLockedReferenceValid == false) = NaN;
    nonPhaseLockedReference( ...
        nonPhaseLockedReferenceValid == false) = NaN;
    processedTotal = bsxfun( ...
        @minus, totalPower, totalReference);
    processedPhaseLocked = bsxfun( ...
        @minus, phaseLockedPower, phaseLockedReference);
    processedNonPhaseLocked = bsxfun( ...
        @minus, nonPhaseLockedPower, nonPhaseLockedReference);
elseif strcmp(baselineMethod, 'zscore') == true
    [totalReference, totalSpread, totalReferenceCount] = ...
        local_masked_sample_stats( ...
        totalBaseline, ...
        totalValid, ...
        3);
    [phaseLockedReference, phaseLockedSpread, ...
        phaseLockedReferenceCount] = ...
        local_masked_sample_stats( ...
        phaseLockedBaseline, ...
        phaseLockedValid, ...
        3);
    [nonPhaseLockedReference, nonPhaseLockedSpread, ...
        nonPhaseLockedReferenceCount] = ...
        local_masked_sample_stats( ...
        nonPhaseLockedBaseline, ...
        nonPhaseLockedValid, ...
        3);
    totalReferenceValid = totalReferenceCount >= 2 & ...
        isfinite(totalReference) & ...
        isfinite(totalSpread) & ...
        local_cpm_spread_is_resolved( ...
        totalSpread, totalReference);
    phaseLockedReferenceValid = ...
        phaseLockedReferenceCount >= 2 & ...
        isfinite(phaseLockedReference) & ...
        isfinite(phaseLockedSpread) & ...
        local_cpm_spread_is_resolved( ...
        phaseLockedSpread, phaseLockedReference);
    nonPhaseLockedReferenceValid = ...
        nonPhaseLockedReferenceCount >= 2 & ...
        isfinite(nonPhaseLockedReference) & ...
        isfinite(nonPhaseLockedSpread) & ...
        local_cpm_spread_is_resolved( ...
        nonPhaseLockedSpread, nonPhaseLockedReference);
    commonReferenceValid = ...
        totalReferenceValid & ...
        phaseLockedReferenceValid & ...
        nonPhaseLockedReferenceValid;
    totalReference(totalReferenceValid == false) = NaN;
    phaseLockedReference(phaseLockedReferenceValid == false) = NaN;
    nonPhaseLockedReference( ...
        nonPhaseLockedReferenceValid == false) = NaN;
    totalSpread(totalReferenceValid == false) = NaN;
    phaseLockedSpread(phaseLockedReferenceValid == false) = NaN;
    nonPhaseLockedSpread( ...
        nonPhaseLockedReferenceValid == false) = NaN;
    processedTotal = bsxfun( ...
        @rdivide, ...
        bsxfun(@minus, totalPower, totalReference), ...
        totalSpread);
    processedPhaseLocked = bsxfun( ...
        @rdivide, ...
        bsxfun(@minus, phaseLockedPower, phaseLockedReference), ...
        phaseLockedSpread);
    processedNonPhaseLocked = bsxfun( ...
        @rdivide, ...
        bsxfun( ...
        @minus, nonPhaseLockedPower, nonPhaseLockedReference), ...
        nonPhaseLockedSpread);
else
    [totalReference, totalReferenceCount] = local_masked_mean( ...
        totalBaseline, ...
        totalValid, ...
        3);
    [phaseLockedReference, phaseLockedReferenceCount] = ...
        local_masked_mean( ...
        phaseLockedBaseline, ...
        phaseLockedValid, ...
        3);
    [nonPhaseLockedReference, nonPhaseLockedReferenceCount] = ...
        local_masked_mean( ...
        nonPhaseLockedBaseline, ...
        nonPhaseLockedValid, ...
        3);
    totalReferenceValid = totalReferenceCount > 0 & ...
        isfinite(totalReference) & ...
        totalReference > 0;
    phaseLockedReferenceValid = phaseLockedReferenceCount > 0 & ...
        isfinite(phaseLockedReference) & ...
        phaseLockedReference > 0;
    nonPhaseLockedReferenceValid = ...
        nonPhaseLockedReferenceCount > 0 & ...
        isfinite(nonPhaseLockedReference) & ...
        nonPhaseLockedReference > 0;
    commonReferenceValid = ...
        totalReferenceValid & ...
        phaseLockedReferenceValid & ...
        nonPhaseLockedReferenceValid;
    totalReference(totalReferenceValid == false) = NaN;
    phaseLockedReference(phaseLockedReferenceValid == false) = NaN;
    nonPhaseLockedReference( ...
        nonPhaseLockedReferenceValid == false) = NaN;
    if strcmp(baselineMethod, 'db') == true
        processedTotal = local_cpm_db_ratio( ...
            totalPower, totalReference);
        processedPhaseLocked = local_cpm_db_ratio( ...
            phaseLockedPower, phaseLockedReference);
        processedNonPhaseLocked = local_cpm_db_ratio( ...
            nonPhaseLockedPower, nonPhaseLockedReference);
    else
        totalRatio = local_cpm_ratio(totalPower, totalReference);
        phaseLockedRatio = local_cpm_ratio( ...
            phaseLockedPower, phaseLockedReference);
        nonPhaseLockedRatio = local_cpm_ratio( ...
            nonPhaseLockedPower, nonPhaseLockedReference);
        processedTotal = local_represent_cpm_ratio( ...
            totalRatio, baselineMethod, 'power');
        processedPhaseLocked = local_represent_cpm_ratio( ...
            phaseLockedRatio, ...
            baselineMethod, ...
            'power_phase_locked');
        processedNonPhaseLocked = local_represent_cpm_ratio( ...
            nonPhaseLockedRatio, ...
            baselineMethod, ...
            'power_non_phase_locked');
    end
end

processedTotal = local_sanitize_cpm_transformed_output( ...
    processedTotal, baselineMethod);
processedPhaseLocked = local_sanitize_cpm_transformed_output( ...
    processedPhaseLocked, baselineMethod);
processedNonPhaseLocked = local_sanitize_cpm_transformed_output( ...
    processedNonPhaseLocked, baselineMethod);

diagnostic.commonReferenceInvalidCount = ...
    sum(commonReferenceValid(:) == false);
diagnostic.totalInvalidReferenceCount = ...
    sum(totalReferenceValid(:) == false);
diagnostic.phaseLockedInvalidReferenceCount = ...
    sum(phaseLockedReferenceValid(:) == false);
diagnostic.nonPhaseLockedInvalidReferenceCount = ...
    sum(nonPhaseLockedReferenceValid(:) == false);
diagnostic.totalInvalidOutputCount = ...
    local_count_invalid_cpm_output(processedTotal);
diagnostic.phaseLockedInvalidOutputCount = ...
    local_count_invalid_cpm_output(processedPhaseLocked);
diagnostic.nonPhaseLockedInvalidOutputCount = ...
    local_count_invalid_cpm_output(processedNonPhaseLocked);

end


function resolved = local_cpm_spread_is_resolved(spread, reference)

referenceScale = abs(reference);
tolerance = eps(referenceScale);
resolved = spread > tolerance;

end


function ratio = local_cpm_ratio(power, reference)

ratio = bsxfun(@rdivide, power, reference);
invalid = isfinite(ratio) == false | ratio < 0;
ratio(invalid) = NaN;

end


function represented = local_cpm_db_ratio(power, reference)

validReference = isfinite(reference) & reference > 0;
positive = bsxfun( ...
    @and, ...
    isfinite(power) & power > 0, ...
    validReference);
zero = bsxfun(@and, power == 0, validReference);
logPower = power;
logPower(positive) = log10(power(positive));
logPower(positive == false) = NaN;
logReference = reference;
logReference(validReference) = log10(reference(validReference));
logReference(validReference == false) = NaN;
logRatio = bsxfun(@minus, logPower, logReference);
represented = power;
represented(positive) = 10 .* logRatio(positive);
represented(zero) = -Inf;
represented(positive == false & zero == false) = NaN;

end


function represented = local_cpm_log10_power(power)

positive = isfinite(power) & power > 0;
zero = power == 0;
represented = power;
represented(positive) = 10 .* log10(power(positive));
represented(zero) = -Inf;
represented(positive == false & zero == false) = NaN;

end


function output = local_sanitize_cpm_transformed_output(values, method)

output = values;

if strcmp(method, 'db') == true
    unintendedInfinity = output == Inf;
else
    unintendedInfinity = isinf(output);
end

output(unintendedInfinity) = NaN;

end


function count = local_count_invalid_cpm_output(values)

invalid = isnan(values) | values == Inf;
count = sum(invalid(:));

end


function represented = local_represent_cpm_ratio(ratio, method, fieldName)

switch method
    case 'ratio'
        represented = ratio;
    case 'percent'
        represented = 100 .* (ratio - 1);
    case 'db'
        represented = local_cpm_log10_power(ratio);
    otherwise
        error('nf_normbase:InternalCPMBaselineMethod', ...
            'Cannot represent CPM field %s with baseline method %s.', ...
            fieldName, ...
            method);
end

end


function [center, count] = local_masked_mean( ...
    values, valid, dimension)

working = values;
working(valid == false) = 0;
count = sum(valid, dimension);
countLike = cast(count, 'like', values);
scale = max(abs(working), [], dimension);
divisor = scale;
divisor(scale == 0) = 1;
scaled = bsxfun(@rdivide, working, divisor);
center = scale .* (sum(scaled, dimension) ./ countLike);
center(count == 0) = NaN;

end


function [center, spread, count] = ...
    local_masked_sample_stats(values, valid, dimension)

[center, count] = local_masked_mean(values, valid, dimension);
working = values;
working(valid == false) = 0;
scale = max(abs(working), [], dimension);
divisor = scale;
divisor(scale == 0) = 1;
scaled = bsxfun(@rdivide, working, divisor);
scaledCenter = center ./ divisor;
deviation = bsxfun(@minus, scaled, scaledCenter);
deviation(valid == false) = 0;
denominator = cast(count - 1, 'like', values);
scaledVariance = sum(deviation .^ 2, dimension) ./ denominator;
spread = scale .* sqrt(scaledVariance);
spread(count < 2) = NaN;

end


function diagnostic = local_empty_cpm_baseline_diagnostic()

diagnostic = struct();
diagnostic.totalExcludedBaselineSampleCount = 0;
diagnostic.phaseLockedExcludedBaselineSampleCount = 0;
diagnostic.nonPhaseLockedExcludedBaselineSampleCount = 0;
diagnostic.commonReferenceInvalidCount = 0;
diagnostic.totalInvalidReferenceCount = 0;
diagnostic.phaseLockedInvalidReferenceCount = 0;
diagnostic.nonPhaseLockedInvalidReferenceCount = 0;
diagnostic.totalInvalidOutputCount = 0;
diagnostic.phaseLockedInvalidOutputCount = 0;
diagnostic.nonPhaseLockedInvalidOutputCount = 0;

end


function processed = local_process_power( ...
    power4, inputScale, options, normalizationIndex, baselineIndex, ...
    fieldName)

if strcmp(options.SingleTrialNormalization, 'none') == true
    if strcmp(inputScale, 'linear') == true
        if strcmp(options.BaselineMethodEffective, 'db') == true
            logPower4 = local_linear_power_log10(power4);
            aggregateLog10 = local_log10_aggregate( ...
                logPower4, 4, options.AverageMethod);
            processed = local_finish_raw_log10( ...
                aggregateLog10, baselineIndex, ...
                options.BaselineMethodEffective, ...
                options.AverageMethod, fieldName);
        else
            aggregate = local_linear_aggregate( ...
                power4, 4, options.AverageMethod);
            processed = local_finish_raw_linear( ...
                aggregate, baselineIndex, ...
                options.BaselineMethodEffective, ...
                options.AverageMethod, fieldName);
        end
    else
        aggregate = local_log10_aggregate( ...
            power4, 4, options.AverageMethod);
        processed = local_finish_raw_log10( ...
            aggregate, baselineIndex, options.BaselineMethodEffective, ...
            options.AverageMethod, fieldName);
    end
else
    useLogDomainGain = ...
        strcmp(options.SingleTrialNormalization, 'db') || ...
        (strcmp(options.SingleTrialNormalization, 'percent') && ...
        strcmp(options.BaselineMethodEffective, 'db'));

    if useLogDomainGain == true
        normalizedLog10Trials = local_single_trial_gain_log10( ...
            power4, inputScale, normalizationIndex, fieldName);
        aggregateLog10 = local_log10_aggregate( ...
            normalizedLog10Trials, 4, options.AverageMethod);
        processed = local_finish_db_normalized_log10( ...
            aggregateLog10, baselineIndex, options, fieldName);
    else
        scaledTrials = local_scale_trials(power4, inputScale);
        normalizedTrials = local_single_trial_normalize( ...
            scaledTrials, normalizationIndex, ...
            options.SingleTrialNormalization, fieldName);
        aggregate = local_dimensionless_aggregate( ...
            normalizedTrials, 4, options.AverageMethod, fieldName);
        processed = local_finish_normalized( ...
            aggregate, normalizedTrials, baselineIndex, options, ...
            fieldName);
    end
end

processed = reshape(processed, ...
    [size(power4, 1), size(power4, 2), size(power4, 3)]);
processed = cast(processed, 'like', power4);

end


function local_validate_power_values(power, inputScale, fieldName)

if isfloat(power) == false || isreal(power) == false
    error('nf_normbase:InvalidPowerField', ...
        'Power field %s must be a real floating-point array.', fieldName);
end

if strcmp(inputScale, 'linear') == true
    if any(isfinite(power(:)) == false) || any(power(:) < 0)
        error('nf_normbase:InvalidLinearPower', ...
            ['Linear power field %s must contain finite, nonnegative ' ...
            'values.'], fieldName);
    end
else
    if any(isnan(power(:))) || any(power(:) == Inf)
        error('nf_normbase:InvalidLog10Power', ...
            ['Log10 power field %s may contain -Inf for exact zero, ' ...
            'but not NaN or +Inf.'], fieldName);
    end
end

end


function scaled = local_scale_trials(power4, inputScale)

if strcmp(inputScale, 'linear') == true
    maximum = max(power4, [], 3);
    divisor = maximum;
    divisor(divisor == 0) = 1;
    scaled = bsxfun(@rdivide, power4, divisor);
else
    maximumLog = max(power4, [], 3);
    shifted = bsxfun(@minus, power4, maximumLog);
    shifted(isnan(shifted)) = -Inf;
    scaled = 10 .^ shifted;
end

scaled = cast(scaled, 'like', power4);

end


function normalized = local_single_trial_normalize( ...
    scaledTrials, normalizationIndex, method, fieldName)

referenceValues = scaledTrials(:, :, normalizationIndex, :);

if strcmp(method, 'zscore') == true
    [referenceMean, referenceSpread, ~] = ...
        local_masked_sample_stats( ...
        referenceValues, isfinite(referenceValues), 3);
    local_assert_nonzero_spread( ...
        referenceSpread, referenceMean, ...
        'nf_normbase:ZeroTrialNormalizationVariance', ...
        ['Single-trial z-normalization is undefined because at least ' ...
        'one trial in power field %s has zero or numerically zero ' ...
        'variance in the normalization interval.'], fieldName);
    normalized = bsxfun(@rdivide, ...
        bsxfun(@minus, scaledTrials, referenceMean), referenceSpread);
else
    [referenceMean, ~] = local_masked_mean( ...
        referenceValues, isfinite(referenceValues), 3);
    local_assert_positive_reference( ...
        referenceMean, ...
        'nf_normbase:ZeroTrialNormalizationMean', ...
        ['Gain normalization is undefined because at least one trial ' ...
        'in power field %s has zero mean power in the normalization ' ...
        'interval.'], fieldName);
    normalized = bsxfun(@rdivide, scaledTrials, referenceMean);
end

if any(isfinite(normalized(:)) == false)
    error('nf_normbase:NonfiniteSingleTrialNormalization', ...
        'Single-trial normalization produced a nonfinite value in %s.', ...
        fieldName);
end

end


function aggregate = local_linear_aggregate(values, dimension, method)

if strcmp(method, 'mean') == true
    maximum = max(values, [], dimension);
    divisor = maximum;
    divisor(divisor == 0) = 1;
    scaled = bsxfun(@rdivide, values, divisor);
    aggregate = maximum .* mean(scaled, dimension);
else
    aggregate = median(values, dimension);
end

end


function aggregate = local_dimensionless_aggregate( ...
    values, dimension, method, fieldName)

if strcmp(method, 'mean') == true
    [aggregate, ~] = local_masked_mean( ...
        values, isfinite(values), dimension);
else
    aggregate = median(values, dimension);
end

if any(isfinite(aggregate(:)) == false)
    error('nf_normbase:NonfiniteTrialAggregate', ...
        'Trial aggregation produced a nonfinite value in %s.', fieldName);
end

end


function aggregate = local_log10_aggregate(logValues, dimension, method)

maximumLog = max(logValues, [], dimension);
shifted = bsxfun(@minus, logValues, maximumLog);
shifted(isnan(shifted)) = -Inf;
scaled = 10 .^ shifted;

if strcmp(method, 'mean') == true
    scaledAggregate = mean(scaled, dimension);
else
    aggregate = local_log10_median(logValues, dimension);
    return
end

aggregate = maximumLog + log10(scaledAggregate);
allZero = maximumLog == -Inf;
aggregate(allZero) = -Inf;

end


function aggregate = local_log10_median(logValues, dimension)

sortedValues = sort(logValues, dimension);
nValues = size(logValues, dimension);
upperIndex = floor(nValues ./ 2) + 1;
subscripts = repmat( ...
    {':'}, 1, max(ndims(logValues), dimension));
subscripts{dimension} = upperIndex;
upper = sortedValues(subscripts{:});

if mod(nValues, 2) == 1
    aggregate = upper;
    return
end

subscripts{dimension} = upperIndex - 1;
lower = sortedValues(subscripts{:});
bothZero = upper == -Inf;
logDifference = lower - upper;
logDifference(bothZero) = 0;
aggregate = upper + log10((10 .^ logDifference + 1) ./ 2);
aggregate(bothZero) = -Inf;

end


function normalizedLog10 = local_single_trial_gain_log10( ...
    power4, inputScale, normalizationIndex, fieldName)

if strcmp(inputScale, 'linear') == true
    logPower = local_linear_power_log10(power4);
else
    logPower = power4;
end

referenceLog10 = local_log10_aggregate( ...
    logPower(:, :, normalizationIndex, :), 3, 'mean');

if any(isfinite(referenceLog10(:)) == false)
    error('nf_normbase:ZeroTrialNormalizationMean', ...
        ['Gain normalization is undefined because at least one trial ' ...
        'in power field %s has zero mean power in the normalization ' ...
        'interval.'], fieldName);
end

normalizedLog10 = bsxfun(@minus, logPower, referenceLog10);
finiteInputPower = isfinite(logPower);

if any(isnan(normalizedLog10(:))) || ...
        any(finiteInputPower(:) & isfinite(normalizedLog10(:)) == false)
    error('nf_normbase:DBDynamicRange', ...
        ['Single-trial dB normalization exceeds the numeric range in ' ...
        '%s. Rescale the input representation.'], fieldName);
end

end


function processed = local_finish_db_normalized_log10( ...
    aggregateLog10, baselineIndex, options, fieldName)

if strcmp(options.BaselineMethodEffective, 'none') == true
    processed = 10 .* aggregateLog10;
else
    processed = local_finish_raw_log10( ...
        aggregateLog10, ...
        baselineIndex, ...
        options.BaselineMethodEffective, ...
        options.AverageMethod, ...
        fieldName);
end

finiteAggregate = isfinite(aggregateLog10);

if any(isnan(processed(:))) || ...
        any(finiteAggregate(:) & isfinite(processed(:)) == false)
    error('nf_normbase:DBDynamicRange', ...
        ['Decibel normalization exceeds the numeric range in %s. ' ...
        'Rescale the input representation.'], fieldName);
end

end


function processed = local_finish_raw_linear( ...
    aggregate, baselineIndex, baselineMethod, averageMethod, fieldName)

if strcmp(baselineMethod, 'none') == true
    processed = aggregate;
    return
end

if strcmp(baselineMethod, 'subtract') == true
    processed = local_baseline_subtract( ...
        aggregate, ...
        baselineIndex, ...
        averageMethod);
    return
end

if strcmp(baselineMethod, 'db') == true
    processed = local_linear_baseline_db( ...
        aggregate, baselineIndex, averageMethod, fieldName);
    return
end

scaled = local_scale_linear_over_time(aggregate);

if strcmp(baselineMethod, 'zscore') == true
    processed = local_classical_zscore( ...
        scaled, baselineIndex, averageMethod, fieldName);
    return
end

ratio = local_linear_baseline_ratio( ...
    scaled, baselineIndex, averageMethod, fieldName);
processed = local_represent_ratio(ratio, baselineMethod, fieldName);

end


function processed = local_finish_raw_log10( ...
    aggregateLog10, baselineIndex, baselineMethod, averageMethod, ...
    fieldName)

if strcmp(baselineMethod, 'none') == true
    processed = aggregateLog10;
    return
end

if strcmp(baselineMethod, 'subtract') == true
    aggregateLinear = 10 .^ aggregateLog10;
    finiteLogPower = isfinite(aggregateLog10);

    if any(isnan(aggregateLinear(:))) || ...
            any(aggregateLinear(:) == Inf) || ...
            any(finiteLogPower(:) & aggregateLinear(:) == 0)
        error('nf_normbase:LinearPowerRange', ...
            ['Finite log10 power cannot be represented as positive ' ...
            'finite linear power for subtractive baselining of %s.'], ...
            fieldName);
    end

    processed = local_baseline_subtract( ...
        aggregateLinear, ...
        baselineIndex, ...
        averageMethod);
    return
end

if strcmp(baselineMethod, 'zscore') == true
    scaled = local_scale_log10_over_time(aggregateLog10);
    processed = local_classical_zscore( ...
        scaled, baselineIndex, averageMethod, fieldName);
    return
end

baselineLog10 = local_log10_aggregate( ...
    aggregateLog10(:, :, baselineIndex), 3, averageMethod);

if any(baselineLog10(:) == -Inf)
    error('nf_normbase:ZeroGroupBaselineMean', ...
        ['Baseline correction is undefined because power field %s has ' ...
        'zero mean baseline power.'], fieldName);
end

logRatio = bsxfun(@minus, aggregateLog10, baselineLog10);

if strcmp(baselineMethod, 'db') == true
    processed = 10 .* logRatio;
    finiteSignal = isfinite(aggregateLog10);

    if any(finiteSignal(:) & isfinite(processed(:)) == false)
        error('nf_normbase:DBDynamicRange', ...
            ['Decibel baselining exceeds the numeric range in %s. ' ...
            'Rescale the input representation.'], fieldName);
    end
else
    ratio = 10 .^ logRatio;

    if any(isnan(ratio(:))) || any(ratio(:) == Inf)
        error('nf_normbase:RatioOverflow', ...
            'Baseline ratios overflowed their numeric type in %s.', ...
            fieldName);
    end

    processed = local_represent_ratio(ratio, baselineMethod, fieldName);
end

end


function processed = local_finish_normalized( ...
    aggregate, normalizedTrials, baselineIndex, options, fieldName)

normalizationMethod = options.SingleTrialNormalization;
baselineMethod = options.BaselineMethodEffective;
averageMethod = options.AverageMethod;

if strcmp(normalizationMethod, 'zscore') == true
    if strcmp(baselineMethod, 'none') == true
        processed = aggregate;
    elseif strcmp(baselineMethod, 'subtract') == true
        processed = local_baseline_subtract( ...
            aggregate, ...
            baselineIndex, ...
            averageMethod);
    else
        processed = local_pooled_normalized_zscore( ...
            aggregate, normalizedTrials, baselineIndex, ...
            averageMethod, fieldName);
    end

    return
end

if strcmp(baselineMethod, 'none') == true
    ratio = aggregate;
elseif strcmp(baselineMethod, 'subtract') == true
    processed = local_baseline_subtract( ...
        aggregate, ...
        baselineIndex, ...
        averageMethod);
    return
elseif strcmp(baselineMethod, 'zscore') == true
    processed = local_classical_zscore( ...
        aggregate, ...
        baselineIndex, ...
        averageMethod, ...
        fieldName);
    return
elseif strcmp(baselineMethod, 'db') == true
    processed = local_linear_baseline_db( ...
        aggregate, baselineIndex, averageMethod, fieldName);
    return
else
    ratio = local_linear_baseline_ratio( ...
        aggregate, baselineIndex, averageMethod, fieldName);
end

if strcmp(normalizationMethod, 'db') == true
    processed = local_represent_ratio(ratio, 'db', fieldName);
else
    processed = local_represent_ratio(ratio, baselineMethod, fieldName);
end

end


function processed = local_baseline_subtract( ...
    aggregate, baselineIndex, averageMethod)

baselineValues = aggregate(:, :, baselineIndex);

if strcmp(averageMethod, 'mean') == true
    [reference, ~] = local_masked_mean( ...
        baselineValues, isfinite(baselineValues), 3);
else
    reference = median(baselineValues, 3);
end

processed = bsxfun(@minus, aggregate, reference);

end


function processed = local_classical_zscore( ...
    aggregate, baselineIndex, averageMethod, fieldName)

baselineValues = aggregate(:, :, baselineIndex);

if strcmp(averageMethod, 'mean') == true
    [center, spread, ~] = local_masked_sample_stats( ...
        baselineValues, isfinite(baselineValues), 3);
else
    center = median(baselineValues, 3);
    spread = local_raw_mad(baselineValues, center, 3);
end

local_assert_nonzero_spread( ...
    spread, center, 'nf_normbase:ZeroGroupBaselineVariance', ...
    ['Z-score baselining is undefined because power field %s has zero ' ...
    'or numerically zero baseline dispersion.'], fieldName);
processed = bsxfun(@rdivide, ...
    bsxfun(@minus, aggregate, center), spread);

if any(isfinite(processed(:)) == false)
    error('nf_normbase:ZScoreDynamicRange', ...
        ['Z-score baselining exceeds the numeric range in %s. Rescale ' ...
        'the input power units.'], fieldName);
end

end


function processed = local_pooled_normalized_zscore( ...
    aggregate, normalizedTrials, baselineIndex, averageMethod, fieldName)

baselineTrials = normalizedTrials(:, :, baselineIndex, :);
nPooled = numel(baselineIndex) .* size(normalizedTrials, 4);

if nPooled < 2 %#ok
    error('nf_normbase:InsufficientPooledBaselineSamples', ...
        ['FullTB-z requires at least two pooled trial-by-baseline ' ...
        'observations in every condition.']);
end

pooled = reshape(baselineTrials, ...
    [size(normalizedTrials, 1), size(normalizedTrials, 2), nPooled]);

if strcmp(averageMethod, 'mean') == true
    [center, spread, ~] = local_masked_sample_stats( ...
        pooled, isfinite(pooled), 3);
else
    center = median(pooled, 3);
    spread = local_raw_mad(pooled, center, 3);
end

local_assert_nonzero_spread( ...
    spread, center, 'nf_normbase:ZeroPooledBaselineVariance', ...
    ['The pooled trial-by-baseline dispersion is zero or numerically ' ...
    'zero in power field %s. FullTB-z is undefined.'], fieldName);
processed = bsxfun(@rdivide, ...
    bsxfun(@minus, aggregate, center), spread);

if any(isfinite(processed(:)) == false)
    error('nf_normbase:ZScoreDynamicRange', ...
        ['Pooled z-score baselining exceeds the numeric range in %s. ' ...
        'Rescale the input power units.'], fieldName);
end

end


function ratio = local_linear_baseline_ratio( ...
    aggregate, baselineIndex, averageMethod, fieldName)

scaled = local_scale_linear_over_time(aggregate);
baselineValues = scaled(:, :, baselineIndex);

if strcmp(averageMethod, 'mean') == true
    [reference, ~] = local_masked_mean( ...
        baselineValues, isfinite(baselineValues), 3);
else
    reference = median(baselineValues, 3);
end

local_assert_positive_reference( ...
    reference, 'nf_normbase:ZeroGroupBaselineMean', ...
    ['Baseline division is undefined because power field %s has zero ' ...
    'baseline reference power.'], fieldName);
ratio = bsxfun(@rdivide, scaled, reference);

if any(isnan(ratio(:))) || any(ratio(:) == Inf) || any(ratio(:) < 0)
    error('nf_normbase:InvalidBaselineRatio', ...
        'Baseline division produced an invalid ratio in %s.', fieldName);
end

end


function processed = local_linear_baseline_db( ...
    values, baselineIndex, averageMethod, fieldName)

logValues = local_linear_power_log10(values);
processed = local_finish_raw_log10( ...
    logValues, baselineIndex, 'db', averageMethod, fieldName);

end


function logValues = local_linear_power_log10(values)

logValues = values;
positive = values > 0;
logValues(positive) = log10(values(positive));
logValues(positive == false) = -Inf;

end


function represented = local_represent_ratio(ratio, method, fieldName)

switch method
    case {'none', 'ratio'}
        represented = ratio;
    case 'percent'
        represented = 100 .* (ratio - 1);
    case 'db'
        represented = 10 .* log10(ratio);
    otherwise
        error('nf_normbase:InternalRatioRepresentation', ...
            'Cannot represent a ratio with method %s.', method);
end

if any(isnan(represented(:))) || any(represented(:) == Inf)
    error('nf_normbase:InvalidNormalizedOutput', ...
        'Normalized representation produced NaN or +Inf in %s.', ...
        fieldName);
end

end


function scaled = local_scale_linear_over_time(values)

maximum = max(values, [], 3);
divisor = maximum;
divisor(divisor == 0) = 1;
scaled = bsxfun(@rdivide, values, divisor);

end


function scaled = local_scale_log10_over_time(logValues)

maximumLog = max(logValues, [], 3);
shifted = bsxfun(@minus, logValues, maximumLog);
shifted(isnan(shifted)) = -Inf;
scaled = 10 .^ shifted;

end


function deviation = local_raw_mad(values, center, dimension)

absoluteDeviation = abs(bsxfun(@minus, values, center));
deviation = median(absoluteDeviation, dimension);

end


function local_assert_positive_reference( ...
    reference, identifier, message, fieldName)

if any(isfinite(reference(:)) == false) || any(reference(:) <= 0)
    error(identifier, message, fieldName);
end

end


function local_assert_nonzero_spread( ...
    spread, center, identifier, message, fieldName)

referenceScale = abs(center);
tolerance = eps(referenceScale);

if any(isfinite(spread(:)) == false) || ...
        any(spread(:) <= tolerance(:))
    error(identifier, message, fieldName);
end

end


function options = local_parse_inputs(TF, varargin)

parser = inputParser;
parser.FunctionName = mfilename;
parser.CaseSensitive = false;
parser.PartialMatching = false;
parser.KeepUnmatched = false;

addRequired(parser, 'TF', @local_valid_tf_input);
addParameter(parser, 'SingleTrialNormalization', 'none', ...
    @local_is_text_scalar);
addParameter(parser, 'NormalizationTimes', 'epoch', ...
    @local_valid_normalization_times);
addParameter(parser, 'BaselineMethod', 'auto', ...
    @local_is_text_scalar);
addParameter(parser, 'BaselineTimes', [], ...
    @local_valid_optional_interval);
addParameter(parser, 'TrialGroups', [], ...
    @local_valid_trial_groups_input);
addParameter(parser, 'AverageMethod', 'mean', ...
    @local_is_text_scalar);
addParameter(parser, 'Verbosity', 'auto', ...
    @local_is_text_scalar);
parse(parser, TF, varargin{:});

options = struct();
options.SingleTrialNormalizationRequested = lower(local_text_scalar( ...
    parser.Results.SingleTrialNormalization, ...
    'SingleTrialNormalization'));
options.SingleTrialNormalization = ...
    local_canonical_single_trial_method( ...
    parser.Results.SingleTrialNormalization);
options.NormalizationTimes = local_canonical_normalization_times( ...
    parser.Results.NormalizationTimes);
options.BaselineMethodRequestedInput = lower(local_text_scalar( ...
    parser.Results.BaselineMethod, 'BaselineMethod'));
options.BaselineMethodRequested = local_canonical_baseline_method( ...
    parser.Results.BaselineMethod);
options.BaselineMethodEffective = local_resolve_baseline_method( ...
    options.SingleTrialNormalization, options.NormalizationTimes, ...
    options.BaselineMethodRequested);
options.BaselineTimes = local_canonical_interval( ...
    parser.Results.BaselineTimes, 'BaselineTimes');
options.TrialGroups = parser.Results.TrialGroups;
options.AverageMethod = local_canonical_average_method( ...
    parser.Results.AverageMethod);
options.Verbosity = local_canonical_verbosity( ...
    parser.Results.Verbosity);
local_validate_method_combination(options);

end


function tf = local_valid_tf_input(value)

tf = isstruct(value) && isscalar(value) && isempty(value) == false;

end


function method = local_canonical_single_trial_method(value)

method = lower(local_text_scalar(value, 'SingleTrialNormalization'));

switch method
    case 'none'
        method = 'none';
    case 'zscore'
        method = 'zscore';
    case 'percent'
        method = 'percent';
    case 'db'
        method = 'db';
    otherwise
        error('nf_normbase:InvalidSingleTrialNormalization', ...
            ['SingleTrialNormalization must be none, zscore, percent, ' ...
            'or db.']);
end

end


function method = local_canonical_baseline_method(value)

method = lower(local_text_scalar(value, 'BaselineMethod'));

switch method
    case 'auto'
        method = 'auto';
    case 'none'
        method = 'none';
    case 'subtract'
        method = 'subtract';
    case 'zscore'
        method = 'zscore';
    case 'ratio'
        method = 'ratio';
    case 'percent'
        method = 'percent';
    case 'db'
        method = 'db';
    otherwise
        error('nf_normbase:InvalidBaselineMethod', ...
            ['BaselineMethod must be auto, none, subtract, zscore, ' ...
            'ratio, percent, or db.']);
end

end


function method = local_resolve_baseline_method( ...
    singleTrialMethod, normalizationTimes, requestedMethod)

if strcmp(requestedMethod, 'auto') == false
    method = requestedMethod;
    return
end

if strcmp(singleTrialMethod, 'none') == false && ...
        isnumeric(normalizationTimes) == true
    method = 'none';
    return
end

switch singleTrialMethod
    case 'none'
        method = 'none';
    case 'zscore'
        method = 'zscore';
    case 'percent'
        method = 'ratio';
    case 'db'
        method = 'db';
    otherwise
        error('nf_normbase:InternalBaselineResolution', ...
            'Could not resolve BaselineMethod=auto.');
end

end


function method = local_canonical_average_method(value)

method = lower(local_text_scalar(value, 'AverageMethod'));

if any(strcmp(method, {'mean', 'median', 'cpm'})) == false
    error('nf_normbase:InvalidAverageMethod', ...
        'AverageMethod must be ''mean'', ''median'', or ''cpm''.');
end

end


function verbosity = local_canonical_verbosity(value)

verbosity = lower(local_text_scalar(value, 'Verbosity'));

if any(strcmp(verbosity, {'verbose', 'auto', 'silent'})) == false
    error('nf_normbase:InvalidVerbosity', ...
        'Verbosity must be ''verbose'', ''auto'', or ''silent''.');
end

end


function local_validate_method_combination(options)

singleTrialMethod = options.SingleTrialNormalization;
baselineMethod = options.BaselineMethodEffective;

if strcmp(options.AverageMethod, 'cpm') == true && ...
        strcmp(singleTrialMethod, 'zscore') == true
    error('nf_normbase:CPMZScoreIncompatible', ...
        ['CPM cannot follow single-trial z-normalization because signed ' ...
        'power cannot define the complex coefficients required by CPM. ' ...
        'Use SingleTrialNormalization=''none'' and, if desired, apply ' ...
        'BaselineMethod=''zscore'' after CPM.']);
end

if strcmp(singleTrialMethod, 'zscore') == true && ...
        any(strcmp(baselineMethod, ...
        {'none', 'subtract', 'zscore'})) == false
    error('nf_normbase:IncompatibleNormalizationMethods', ...
        ['Z-normalized trials are signed. Their group baseline must be ' ...
        'none, subtract, or zscore.']);
end

if strcmp(singleTrialMethod, 'db') == true && ...
        any(strcmp(baselineMethod, {'none', 'db'})) == false
    error('nf_normbase:IncompatibleNormalizationMethods', ...
        ['SingleTrialNormalization=''db'' requires BaselineMethod none, ' ...
        'db, or auto.']);
end

end


function cautions = local_cpm_scientific_cautions(options)

cautions = repmat(local_empty_caution(), 1, 0);

if strcmp(options.AverageMethod, 'cpm') == false
    return
end

singleTrialMethod = options.SingleTrialNormalization;
baselineMethod = options.BaselineMethodEffective;

if any(strcmp(singleTrialMethod, {'percent', 'db'})) == true
    caution = local_empty_caution();
    caution.identifier = 'nf_normbase:CPMTrialGainNormalization';
    caution.category = 'experimental_combination';
    caution.message = [ ...
        'Trial gain normalization before CPM preserves coefficient phase ' ...
        'but changes the estimand to gain-normalized component power. ' ...
        'It is valid only when the trial gain is independent of the ' ...
        'non-phase-locked phase; this combination has not been validated ' ...
        'by either source publication.'];
    cautions(end + 1) = caution;

    if ischar(options.NormalizationTimes) == true
        caution = local_empty_caution();
        caution.identifier = ...
            'nf_normbase:CPMFullEpochGainNormalization';
        caution.category = 'exploratory_combination';
        caution.message = [ ...
            'Full-epoch trial gain normalization before CPM is ' ...
            'exploratory: event activity and PL/NPL interference enter ' ...
            'the normalizer and can violate CPM weighted-phase ' ...
            'assumptions. A strictly prestimulus normalization interval ' ...
            'is more defensible.'];
        cautions(end + 1) = caution;
    end
end

if any(strcmp(baselineMethod, {'ratio', 'percent', 'db'})) == true
    caution = local_empty_caution();
    caution.identifier = ...
        'nf_normbase:CPMDivisivePhaseLockedBaseline';
    caution.category = 'interpretation_caution';
    caution.message = [ ...
        'Divisive baselining of phase-locked CPM power can be unstable ' ...
        'because its prestimulus expectation is commonly near zero. ' ...
        'Use subtraction unless a genuine nonzero phase-locked baseline ' ...
        'is scientifically expected.'];
    cautions(end + 1) = caution;
end

end


function caution = local_empty_caution()

caution = struct();
caution.identifier = '';
caution.category = '';
caution.message = '';

end


function value = local_canonical_normalization_times(value)

if local_is_text_scalar(value) == true
    textValue = local_text_scalar(value, 'NormalizationTimes');

    if strcmpi(textValue, 'epoch') == false
        error('nf_normbase:InvalidNormalizationTimes', ...
            ['NormalizationTimes must be ''epoch'' or an inclusive ' ...
            '[minimum maximum] interval.']);
    end

    value = 'epoch';
else
    value = local_canonical_interval(value, 'NormalizationTimes');
end

end


function value = local_canonical_interval(value, optionName)

if isempty(value) == true
    value = [];
    return
end

if isnumeric(value) == false || isreal(value) == false || ...
        numel(value) ~= 2 || any(isfinite(value(:)) == false)
    error('nf_normbase:InvalidTimeInterval', ...
        '%s must be an empty value or two finite real limits.', ...
        optionName);
end

value = double(value(:).');

if value(1) > value(2)
    error('nf_normbase:ReversedTimeInterval', ...
        '%s limits must be ordered from minimum to maximum.', ...
        optionName);
end

end


function tf = local_valid_normalization_times(value)

if local_is_text_scalar(value) == true
    tf = true;
else
    tf = local_valid_optional_interval(value) && isempty(value) == false;
end

end


function tf = local_valid_optional_interval(value)

tf = isempty(value) || ...
    (isnumeric(value) && isreal(value) && numel(value) == 2 && ...
    all(isfinite(value(:))));

end


function tf = local_valid_trial_groups_input(value)

tf = isempty(value) || ...
    ((isnumeric(value) || islogical(value)) && isvector(value) && ...
    isreal(value));

end


function tf = local_is_text_scalar(value)

tf = ischar(value) || (isstring(value) && isscalar(value));

end


function value = local_text_scalar(value, name)

if local_is_text_scalar(value) == false
    error('nf_normbase:TextScalarRequired', ...
        '%s must be a character vector or string scalar.', name);
end

value = strtrim(char(value));

if isempty(value) == true
    error('nf_normbase:EmptyTextOption', ...
        '%s cannot be empty.', name);
end

end


function value = local_positive_integer(value, name)

if isnumeric(value) == false || isreal(value) == false || ...
        isscalar(value) == false || isfinite(value) == false || ...
        value < 1 || value ~= fix(value)
    error('nf_normbase:PositiveIntegerRequired', ...
        '%s must be a finite positive integer scalar.', name);
end

value = double(value);

end


function axisValues = local_finite_axis(value, name)

if isnumeric(value) == false || isreal(value) == false || ...
        isvector(value) == false || isempty(value) == true || ...
        any(isfinite(value(:)) == false)
    error('nf_normbase:InvalidAxis', ...
        '%s must be a nonempty finite real numeric vector.', name);
end

axisValues = double(value(:).');

end


function [index, appliedTimes] = ...
    local_normalization_index(times, requestedTimes)

if ischar(requestedTimes) == true
    index = 1:numel(times);
else
    index = find(times >= requestedTimes(1) & ...
        times <= requestedTimes(2));
end

if isempty(index) == true
    error('nf_normbase:EmptyNormalizationInterval', ...
        'NormalizationTimes selects no samples from TF.times.');
end

appliedTimes = [min(times(index)), max(times(index))];

end


function [index, appliedTimes] = ...
    local_baseline_index(times, requestedTimes, baselineMethod)

if strcmp(baselineMethod, 'none') == true
    index = [];
    appliedTimes = [];
    return
end

if isempty(requestedTimes) == true
    index = find(times <= 0);
else
    index = find(times >= requestedTimes(1) & ...
        times <= requestedTimes(2));
end

if isempty(index) == true
    error('nf_normbase:EmptyBaselineInterval', ...
        'BaselineTimes selects no samples from TF.times.');
end

appliedTimes = [min(times(index)), max(times(index))];

end


function groups = local_trial_groups(value, nTrials)

if isempty(value) == true
    groups = ones(1, nTrials);
    return
end

if numel(value) ~= nTrials
    error('nf_normbase:TrialGroupLengthMismatch', ...
        ['TrialGroups contains %d labels, but TF metadata specifies ' ...
        '%d trials.'], numel(value), nTrials);
end

groups = value(:).';

if any(isfinite(double(groups)) == false)
    error('nf_normbase:InvalidTrialGroupLabel', ...
        'TrialGroups must contain finite labels.');
end

end


function output = local_output_layout( ...
    canonical, nSensors, nFrequencies, nTimes, nConditions)

if nSensors == 1
    output = reshape(canonical, ...
        [nFrequencies, nTimes, nConditions]);
else
    output = reshape(canonical, ...
        [nSensors, nFrequencies, nTimes, nConditions]);
end

end


function cpm = local_pack_cpm_models( ...
    models, baselineDiagnostics, nSensors, nFrequencies, nTimes, ...
    nConditions, conditions, trialCounts, options)

firstModel = models{1};
cpm = struct();
cpm.schema = 'nf-cpm-tf/1.0';
cpm.citation = firstModel.citation;
cpm.doi = firstModel.doi;
cpm.estimator = firstModel.estimator;
cpm.sourceAnalyticCoefficientSemantics = ...
    options.CPMTransform.requiredSemantics;
cpm.inputTransformMethod = options.CPMTransform.inputMethod;
cpm.transformVerifiedBy = options.CPMTransform.verifiedBy;
cpm.estimand = ...
    'condition_ensemble_phase_locked_and_non_phase_locked_power';
cpm.decomposedPowerField = 'power';
cpm.pairedPhaseField = options.CPMPhaseField;
cpm.ancillaryPowerAverageMethod = 'arithmetic_mean';
cpm.normalization = options.SingleTrialNormalization;
cpm.baselineMethod = options.BaselineMethodEffective;
cpm.singleTrialGainCombinationApplied = any(strcmp( ...
    options.SingleTrialNormalization, ...
    {'percent', 'db'}));
cpm.normalizationCombinationValidatedBySources = ...
    strcmp(options.SingleTrialNormalization, 'none');

if strcmp(options.SingleTrialNormalization, 'none') == true
    cpm.rawPowerScale = 'linear_input_power';
    cpm.rawAmplitudeScale = 'sqrt_linear_input_power';
else
    cpm.rawPowerScale = 'linear_trial_gain_normalized_power';
    cpm.rawAmplitudeScale = ...
        'sqrt_linear_trial_gain_normalized_power';
end

cpm.scientificCautions = options.ScientificCautions;
cpm.conditions = struct();
cpm.conditions.values = conditions;
cpm.conditions.trialCounts = trialCounts;

cpm.total = struct();
cpm.total.powerMean = local_pack_cpm_field( ...
    models, ...
    {'total', 'powerMean'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);

cpm.phaseLocked = struct();
cpm.phaseLocked.phase = local_pack_cpm_field( ...
    models, ...
    {'phaseLocked', 'phase'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.phaseLocked.amplitudeMean = local_pack_cpm_field( ...
    models, ...
    {'phaseLocked', 'amplitudeMean'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.phaseLocked.amplitudeVariance = local_pack_cpm_field( ...
    models, ...
    {'phaseLocked', 'amplitudeVariance'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.phaseLocked.amplitudeSecondMoment = local_pack_cpm_field( ...
    models, ...
    {'phaseLocked', 'amplitudeSecondMoment'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.phaseLocked.powerMean = local_pack_cpm_field( ...
    models, ...
    {'phaseLocked', 'powerMean'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);

cpm.nonPhaseLocked = struct();
cpm.nonPhaseLocked.amplitudeMean = local_pack_cpm_field( ...
    models, ...
    {'nonPhaseLocked', 'amplitudeMean'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.nonPhaseLocked.amplitudeVariance = local_pack_cpm_field( ...
    models, ...
    {'nonPhaseLocked', 'amplitudeVariance'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.nonPhaseLocked.amplitudeSecondMoment = local_pack_cpm_field( ...
    models, ...
    {'nonPhaseLocked', 'amplitudeSecondMoment'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.nonPhaseLocked.powerMean = local_pack_cpm_field( ...
    models, ...
    {'nonPhaseLocked', 'powerMean'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.leakageRatio = local_pack_cpm_field( ...
    models, ...
    {'leakageRatio'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);

cpm.raw = struct();
cpm.raw.phaseLockedAmplitudeSecondMoment = local_pack_cpm_field( ...
    models, ...
    {'raw', 'phaseLockedAmplitudeSecondMoment'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.raw.nonPhaseLockedAmplitudeSecondMoment = local_pack_cpm_field( ...
    models, ...
    {'raw', 'nonPhaseLockedAmplitudeSecondMoment'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.raw.phaseLockedAmplitudeVariance = local_pack_cpm_field( ...
    models, ...
    {'raw', 'phaseLockedAmplitudeVariance'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.raw.nonPhaseLockedAmplitudeVariance = local_pack_cpm_field( ...
    models, ...
    {'raw', 'nonPhaseLockedAmplitudeVariance'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.raw.closureError = local_pack_cpm_field( ...
    models, ...
    {'raw', 'closureError'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);

cpm.validity = struct();
cpm.validity.validTrialCount = local_pack_cpm_field( ...
    models, ...
    {'validity', 'validTrialCount'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.validity.sufficientTrials = local_pack_cpm_field( ...
    models, ...
    {'validity', 'sufficientTrials'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.validity.phaseDefined = local_pack_cpm_field( ...
    models, ...
    {'validity', 'phaseDefined'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.validity.positiveProjectionCount = local_pack_cpm_field( ...
    models, ...
    {'validity', 'positiveProjectionCount'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.validity.conditioningDenominator = local_pack_cpm_field( ...
    models, ...
    {'validity', 'conditioningDenominator'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.validity.inverseAbsoluteDenominator = local_pack_cpm_field( ...
    models, ...
    {'validity', 'inverseAbsoluteDenominator'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.validity.nonsingular = local_pack_cpm_field( ...
    models, ...
    {'validity', 'nonsingular'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.validity.phaseLockedSecondMoment = local_pack_cpm_field( ...
    models, ...
    {'validity', 'phaseLockedSecondMoment'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.validity.nonPhaseLockedSecondMoment = local_pack_cpm_field( ...
    models, ...
    {'validity', 'nonPhaseLockedSecondMoment'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.validity.phaseLockedAmplitudeVariance = local_pack_cpm_field( ...
    models, ...
    {'validity', 'phaseLockedAmplitudeVariance'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
cpm.validity.nonPhaseLockedAmplitudeVariance = local_pack_cpm_field( ...
    models, ...
    {'validity', 'nonPhaseLockedAmplitudeVariance'}, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);
conditionSummaries = repmat( ...
    models{1}.validity.summary, ...
    1, ...
    nConditions);

for conditionIndex = 1:nConditions
    conditionSummaries(conditionIndex) = ...
        models{conditionIndex}.validity.summary;
end

cpm.validity.conditionSummaries = conditionSummaries;
cpm.validity.summary = local_sum_cpm_summaries(models);
cpm.baselineDiagnosticsByCondition = baselineDiagnostics;
cpm.baselineDiagnostics = local_sum_cpm_baseline_diagnostics( ...
    baselineDiagnostics);
cpm.additivity = struct();
cpm.additivity.rawIdentity = ...
    'total_power=phase_locked_power+non_phase_locked_power';
cpm.additivity.rawPreservedWhenMomentsValid = true;
cpm.additivity.rawIdentityDomain = ...
    'locations_where_raw_total_and_both_raw_component_moments_are_finite';
cpm.additivity.returnedCommonFiniteDomain = ...
    'locations_where_total_and_both_returned_components_are_finite';

if strcmp(options.BaselineMethodEffective, 'none') == true
    cpm.additivity.baselineReferencePolicy = 'not_applicable';
else
    cpm.additivity.baselineReferencePolicy = ...
        'component_specific_finite_baseline_samples';
end

cpm.additivity.postBaselineConstraint = 'none';
cpm.additivity.outputPreserved = ...
    strcmp(options.SingleTrialNormalization, 'db') == false && ...
    strcmp(options.BaselineMethodEffective, 'none');

end


function packed = local_pack_cpm_field( ...
    models, path, nSensors, nFrequencies, nTimes, nConditions)

firstValue = local_get_nested_field(models{1}, path);

if islogical(firstValue) == true
    canonical = false( ...
        nSensors, ...
        nFrequencies, ...
        nTimes, ...
        nConditions);
else
    canonical = zeros( ...
        nSensors, ...
        nFrequencies, ...
        nTimes, ...
        nConditions, ...
        'like', ...
        firstValue);
end

for conditionIndex = 1:nConditions
    currentValue = local_get_nested_field( ...
        models{conditionIndex}, ...
        path);
    canonical(:, :, :, conditionIndex) = reshape( ...
        currentValue, ...
        [nSensors, nFrequencies, nTimes, 1]);
end

packed = local_output_layout( ...
    canonical, ...
    nSensors, ...
    nFrequencies, ...
    nTimes, ...
    nConditions);

end


function value = local_get_nested_field(structure, path)

value = structure;

for pathIndex = 1:numel(path)
    value = value.(path{pathIndex});
end

end


function summary = local_sum_cpm_summaries(models)

summary = models{1}.validity.summary;
names = fieldnames(summary);

for nameIndex = 1:numel(names)
    fieldName = names{nameIndex};
    total = 0;

    for conditionIndex = 1:numel(models)
        total = total + ...
            models{conditionIndex}.validity.summary.(fieldName);
    end

    summary.(fieldName) = total;
end

end


function summary = local_sum_cpm_baseline_diagnostics(diagnostics)

summary = local_empty_cpm_baseline_diagnostic();
names = fieldnames(summary);

for nameIndex = 1:numel(names)
    fieldName = names{nameIndex};
    values = [diagnostics.(fieldName)];
    summary.(fieldName) = sum(values);
end

end


function local_report_cpm_validity(cpm, verbosity)

summary = cpm.validity.summary;
baseline = cpm.baselineDiagnostics;
undefinedCount = summary.insufficientTrialLocations + ...
    summary.incompleteTrialLocations + ...
    summary.undefinedPhaseLocations;
invalidMomentCount = ...
    summary.invalidPhaseLockedSecondMoments + ...
    summary.invalidNonPhaseLockedSecondMoments + ...
    summary.invalidPhaseLockedAmplitudeVariances + ...
    summary.invalidNonPhaseLockedAmplitudeVariances;
invalidReferenceCount = ...
    baseline.totalInvalidReferenceCount + ...
    baseline.phaseLockedInvalidReferenceCount + ...
    baseline.nonPhaseLockedInvalidReferenceCount;

if undefinedCount > 0
    local_issue_warning( ...
        verbosity, ...
        'nf_normbase:CPMUndefinedPhase', ...
        ['%d CPM validity issue(s) involved incomplete complex trials, ' ...
        'fewer than two finite trials, or no defined PL phase. Missing ' ...
        'trials were omitted and non-estimable outputs are NaN.'], ...
        undefinedCount);
end

if summary.singularLocations > 0
    local_issue_warning( ...
        verbosity, ...
        'nf_normbase:CPMSingularMomentEquation', ...
        ['%d CPM coefficient location(s) reached the published ' ...
        'cos(2*alpha)=0 singularity; affected component estimates are ' ...
        'NaN.'], ...
        summary.singularLocations);
end

if invalidMomentCount > 0
    local_issue_warning( ...
        verbosity, ...
        'nf_normbase:CPMInvalidMomentEstimate', ...
        ['%d negative finite-sample CPM moment estimate(s) were kept ' ...
        'in TF.cpm.raw and represented as NaN in physical outputs.'], ...
        invalidMomentCount);
end

if invalidReferenceCount > 0
    local_issue_warning( ...
        verbosity, ...
        'nf_normbase:CPMInvalidBaselineReference', ...
        ['%d CPM baseline reference(s) were undefined, nonpositive for ' ...
        'division, or had zero dispersion for z-scoring; affected ' ...
        'normalized outputs are NaN.'], ...
        invalidReferenceCount);
end

end


function ledger = local_build_ledger( ...
    options, inputScale, outputScale, normalizationIndex, ...
    normalizationTimesApplied, baselineIndex, baselineTimesApplied, ...
    conditions, trialCounts, transformationApplied)

[outputMeasure, outputUnits, percentConvention] = ...
    local_output_description(options, inputScale);

ledger = struct();
ledger.schema = 'nf-normbase/1.0';
ledger.citation = ...
    'Grandchamp R, Delorme A. Front Psychol. 2011;2:236.';
ledger.doi = '10.3389/fpsyg.2011.00236';
ledger.inputScale = inputScale;
ledger.outputScale = outputScale;
ledger.outputMeasure = outputMeasure;
ledger.outputUnits = outputUnits;
ledger.transformationApplied = transformationApplied;
ledger.inputRepresentationConvertedToLinear = ...
    strcmp(inputScale, 'log10') == true;

if transformationApplied == true
    ledger.outputKind = 'normalized_power';
else
    ledger.outputKind = 'power';
end

ledger.percentConvention = percentConvention;
ledger.averageMethod = options.AverageMethod;
ledger.scientificCautions = options.ScientificCautions;

singleTrial = struct();
singleTrial.requestedMethod = ...
    options.SingleTrialNormalizationRequested;
singleTrial.effectiveMethod = options.SingleTrialNormalization;

switch options.SingleTrialNormalization
    case 'none'
        singleTrial.internalOperator = 'none';
    case 'zscore'
        singleTrial.internalOperator = 'temporal_zscore';
    otherwise
        singleTrial.internalOperator = 'temporal_gain_ratio';
end

singleTrial.windowRequested = options.NormalizationTimes;
singleTrial.indices = normalizationIndex;
singleTrial.timesApplied = normalizationTimesApplied;

if strcmp(options.SingleTrialNormalization, 'none') == true
    singleTrial.statistic = 'none';
    singleTrial.spread = 'none';
elseif strcmp(options.SingleTrialNormalization, 'zscore') == true
    singleTrial.statistic = 'arithmetic_mean';
    singleTrial.spread = 'sample_sd_m_minus_1';
else
    singleTrial.statistic = 'arithmetic_mean';
    singleTrial.spread = 'none';
end

ledger.singleTrial = singleTrial;

baseline = struct();
baseline.requestedMethod = options.BaselineMethodRequestedInput;
baseline.effectiveMethod = options.BaselineMethodEffective;
baseline.windowRequested = options.BaselineTimes;
baseline.indices = baselineIndex;
baseline.timesApplied = baselineTimesApplied;

if strcmp(options.BaselineMethodEffective, 'zscore') == true && ...
        strcmp(options.SingleTrialNormalization, 'zscore') == true
    baseline.referenceScope = 'pooled_condition_trials_x_baseline_time';
elseif strcmp(options.BaselineMethodEffective, 'none') == true
    baseline.referenceScope = 'none';
elseif strcmp(options.AverageMethod, 'cpm') == true
    baseline.referenceScope = ...
        'each_cpm_output_x_own_finite_baseline_time';
else
    baseline.referenceScope = 'condition_aggregate_x_baseline_time';
end

if strcmp(options.AverageMethod, 'cpm') == true && ...
        strcmp(options.BaselineMethodEffective, 'none') == false
    baseline.componentReferenceCoupling = 'independent';
else
    baseline.componentReferenceCoupling = 'not_applicable';
end

if strcmp(options.BaselineMethodEffective, 'none') == true
    baseline.referenceStatistic = 'none';
elseif strcmp(options.BaselineMethodEffective, 'zscore') == true && ...
        any(strcmp(options.AverageMethod, {'mean', 'cpm'})) == true
    baseline.referenceStatistic = 'arithmetic_mean_and_sample_sd';
elseif strcmp(options.BaselineMethodEffective, 'zscore') == true
    baseline.referenceStatistic = 'median_and_raw_mad';
elseif any(strcmp(options.AverageMethod, {'mean', 'cpm'})) == true
    baseline.referenceStatistic = 'arithmetic_mean';
else
    baseline.referenceStatistic = 'median';
end

ledger.baseline = baseline;
ledger.conditions = struct();
ledger.conditions.values = conditions;
ledger.conditions.trialCounts = trialCounts;
ledger.operationOrder = local_operation_order(options);
baselineIsStrictlyPrestimulus = isempty(baselineTimesApplied) == false && ...
    all(baselineTimesApplied < 0);
ledger.baselineIsStrictlyPrestimulus = baselineIsStrictlyPrestimulus;

isFullEpoch = ischar(options.NormalizationTimes) && ...
    strcmp(options.NormalizationTimes, 'epoch');
isMean = strcmp(options.AverageMethod, 'mean');
isMatchingPaperOutput = ...
    (strcmp(options.SingleTrialNormalization, 'zscore') && ...
    strcmp(options.BaselineMethodEffective, 'zscore')) || ...
    (strcmp(options.SingleTrialNormalization, 'percent') && ...
    strcmp(options.BaselineMethodEffective, 'ratio')) || ...
    (strcmp(options.SingleTrialNormalization, 'db') && ...
    strcmp(options.BaselineMethodEffective, 'db'));
ledger.grandchampDelormeFormulaExact = ...
    isFullEpoch && isMean && isMatchingPaperOutput;
ledger.grandchampDelormeExact = ...
    ledger.grandchampDelormeFormulaExact && ...
    baselineIsStrictlyPrestimulus;
ledger.fullEpochNormalization = isFullEpoch && ...
    strcmp(options.SingleTrialNormalization, 'none') == false;
ledger.neurofreqMedianExtension = ...
    strcmp(options.AverageMethod, 'median');
ledger.neurofreqPercentChangeExtension = ...
    strcmp(options.BaselineMethodEffective, 'percent');

if strcmp(options.AverageMethod, 'cpm') == true
    cpm = struct();
    cpm.citation = [ ...
        'Singhal S, Ghosh P, Kumar N, Banerjee A. ' ...
        'J Neurophysiol. 2023;129(1):199-210.'];
    cpm.doi = '10.1152/jn.00467.2022';
    cpm.estimand = ...
        'ensemble_phase_locked_and_non_phase_locked_mean_power';
    cpm.estimator = 'published_finite_sample_moments';
    cpm.decomposedPowerField = 'power';
    cpm.pairedPhaseField = options.CPMPhaseField;
    cpm.ancillaryPowerAverageMethod = 'arithmetic_mean';
    cpm.singleTrialGainCombinationApplied = any(strcmp( ...
        options.SingleTrialNormalization, ...
        {'percent', 'db'}));
    cpm.singleTrialGainCombinationValidatedBySources = false;
    cpm.rawComponentAdditivity = ...
        'phase_locked_plus_non_phase_locked_equals_total_when_valid';
    cpm.postBaselineConstraint = 'none';

    if strcmp(options.BaselineMethodEffective, 'none') == true
        cpm.baselineReferencePolicy = 'not_applicable';
    else
        cpm.baselineReferencePolicy = ...
            'component_specific_finite_baseline_samples';
    end

    cpm.outputAdditivityPreserved = ...
        strcmp(options.SingleTrialNormalization, 'db') == false && ...
        strcmp(options.BaselineMethodEffective, 'none');
    ledger.cpm = cpm;
end

end


function order = local_operation_order(options)

order = {};

if strcmp(options.SingleTrialNormalization, 'none') == false
    order{end + 1} = 'single_trial_normalization';
end

if strcmp(options.AverageMethod, 'cpm') == true
    order{end + 1} = 'within_condition_cpm_estimation';
else
    order{end + 1} = 'within_condition_trial_aggregation';
end

if strcmp(options.BaselineMethodEffective, 'none') == false
    order{end + 1} = 'condition_level_baseline_correction';
end

if strcmp(options.SingleTrialNormalization, 'db') == true || ...
        strcmp(options.BaselineMethodEffective, 'db') == true
    order{end + 1} = 'post_aggregation_decibel_transform';
end

end


function [measure, units, percentConvention] = ...
    local_output_description(options, inputScale)

singleTrialMethod = options.SingleTrialNormalization;
baselineMethod = options.BaselineMethodEffective;
percentConvention = 'not_applicable';

if strcmp(singleTrialMethod, 'none') == true && ...
        strcmp(baselineMethod, 'none') == true
    if strcmp(options.AverageMethod, 'cpm') == true
        measure = 'cpm_ensemble_mean_power';
        units = 'linear_input_power_units';
        return
    end

    measure = 'trial_average_power';

    if strcmp(inputScale, 'linear') == true
        units = 'input_power_units';
    else
        units = 'log10_input_power_units';
    end

    return
end

if strcmp(baselineMethod, 'subtract') == true
    measure = 'baseline_subtracted_power';

    if strcmp(singleTrialMethod, 'none') == true
        units = 'linear_input_power_units';
    elseif strcmp(singleTrialMethod, 'zscore') == true
        units = 'single_trial_zscore_difference';
    else
        units = 'gain_normalized_power_difference';
    end

    return
end

if strcmp(baselineMethod, 'zscore') == true || ...
        (strcmp(singleTrialMethod, 'zscore') == true && ...
        strcmp(baselineMethod, 'none') == true)
    measure = 'zscore';
    units = 'standard_deviations';
    return
end

if strcmp(singleTrialMethod, 'db') == true || ...
        strcmp(baselineMethod, 'db') == true
    measure = 'decibel_change';
    units = 'dB';
    percentConvention = 'ratio transformed as 10*log10(ratio)';
    return
end

if strcmp(baselineMethod, 'percent') == true
    measure = 'percent_change';
    units = 'percent';
    percentConvention = '100*(ratio-1), zero denotes baseline';
    return
end

measure = 'power_ratio';
units = 'proportion_of_reference';
percentConvention = 'ratio, one denotes 100 percent of reference';

end


function powerFields = nf_find_power_fields_local(TF)

powerFields = local_find_measure_fields(TF, 'power');

end


function phaseFields = nf_find_phase_fields_local(TF)

phaseFields = local_find_measure_fields(TF, 'phase');

end


function outputFieldName = local_itpc_field_name(phaseFieldName)

if strcmpi(phaseFieldName, 'phase') == true
    outputFieldName = 'itpc';
    return
end

if strncmpi(phaseFieldName, 'phase_', 6) == false
    error('nf_normbase:InternalPhaseFieldName', ...
        'Cannot derive an ITPC field name from %s.', phaseFieldName);
end

outputFieldName = ['itpc_', phaseFieldName(7:end)];

end


function fields = local_find_measure_fields(TF, prefix)

names = fieldnames(TF);
fields = {};

for fieldIndex = 1:numel(names)
    fieldName = names{fieldIndex};
    value = TF.(fieldName);
    lowerName = lower(fieldName);
    exactName = strcmp(lowerName, prefix);
    prefixName = strncmp(lowerName, [prefix '_'], numel(prefix) + 1);

    if isnumeric(value) == true && (exactName || prefixName)
        fields{end + 1, 1} = fieldName; %#ok<AGROW>
    end
end

for fieldIndex = 1:numel(fields)
    duplicateIndex = find(strcmpi(fields, fields{fieldIndex}));

    if numel(duplicateIndex) > 1
        error('nf_normbase:AmbiguousMeasureFields', ...
            ['TF contains measure fields whose names differ only by ' ...
            'letter case (%s). Field names must be unique under ' ...
            'case-insensitive matching.'], ...
            strjoin(fields(duplicateIndex), ', '));
    end
end

primaryIndex = find(strcmpi(fields, prefix));

if isempty(primaryIndex) == false
    primaryField = fields{primaryIndex(1)};
    fields(primaryIndex) = [];
    fields = [{primaryField}; fields];
end

end


function values4 = nf_as_canonical4_local( ...
    values, fieldName, nSensors, nFrequencies, nTimes, nTrials)

logicalSize = [nSensors, nFrequencies, nTimes, nTrials];
expectedElements = prod(logicalSize);

if numel(values) ~= expectedElements
    error('nf_normbase:TensorSizeMismatch', ...
        ['Field %s contains %d values; expected %d from NeuroFreq ' ...
        'dimension metadata.'], ...
        fieldName, numel(values), expectedElements);
end

storedSize = size(values);
canonicalStoredSize = storedSize;

if numel(canonicalStoredSize) < numel(logicalSize)
    canonicalStoredSize(end + 1:numel(logicalSize)) = 1;
end

squeezedSize = logicalSize(logicalSize ~= 1);

if isempty(squeezedSize) == true
    squeezedSize = [1, 1];
elseif isscalar(squeezedSize) == true
    squeezedSize = [squeezedSize, 1];
end

if isequal(canonicalStoredSize, logicalSize) == false && ...
        isequal(storedSize, squeezedSize) == false
    error('nf_normbase:TensorDimensionOrder', ...
        ['Field %s must use sensor x frequency x time x trial ordering ' ...
        'or its singleton-squeezed representation.'], fieldName);
end

values4 = reshape(values, logicalSize);

end


function itpc = nf_phase_itpc_local(phaseValues, trialDimension)

validPhase = isfinite(phaseValues);
unitPhasors = exp(1i .* phaseValues);
unitPhasors(validPhase == false) = 0;
validCount = sum(validPhase, trialDimension);
phasorSum = sum(unitPhasors, trialDimension);
validCountLikePhase = cast(validCount, 'like', phaseValues);
itpc = abs(phasorSum) ./ validCountLikePhase;
itpc = cast(itpc, 'like', phaseValues);
itpc(validCount == 0) = NaN;

end


function TF = nf_powerfront_local(TF)

names = fieldnames(TF);
powerFields = {};
otherFields = {};

for fieldIndex = 1:numel(names)
    fieldName = names{fieldIndex};
    lowerName = lower(fieldName);
    exactName = strcmp(lowerName, 'power');
    prefixName = strncmp(lowerName, 'power_', 6);

    if isnumeric(TF.(fieldName)) == true && (exactName || prefixName)
        powerFields{end + 1, 1} = fieldName; %#ok<AGROW>
    else
        otherFields{end + 1, 1} = fieldName; %#ok<AGROW>
    end
end

primaryIndex = find(strcmpi(powerFields, 'power'));

if isempty(primaryIndex) == false
    powerFields(primaryIndex) = [];
    powerFields = [{'power'}; powerFields];
end

TF = orderfields(TF, [powerFields; otherFields]);

end


function output = local_average_behavior(structures)

names = fieldnames(structures);
output = struct();

for fieldIndex = 1:numel(names)
    fieldName = names{fieldIndex};
    values = {structures.(fieldName)};
    scalarNumeric = cellfun( ...
        @(value) (isnumeric(value) || islogical(value)) && ...
        isreal(value) && isscalar(value), values);

    if all(scalarNumeric) == true
        output.(fieldName) = nanmean([values{:}]);
        continue
    end

    firstValue = values{1};
    allIdentical = all(cellfun( ...
        @(value) isequaln(value, firstValue), values));

    if allIdentical == true
        output.(fieldName) = firstValue;
    else
        error('nf_normbase:NonconstantBehaviorMetadata', ...
            ['Nonnumeric or nonscalar behavior field %s varies within ' ...
            'a condition and cannot be averaged safely.'], fieldName);
    end
end

end


function value = nanmean(values)

if nargin < 1
    error('nf_normbase:NanmeanInput', ...
        'nanmean requires an input.');
end

if isempty(values) == true
    value = NaN;
    return
end

if isfloat(values) == false
    values = double(values);
end

try
    value = mean(values, 'omitnan');
    return
catch
end

values = values(:);
values = values(isfinite(values));

if isempty(values) == true
    value = NaN;
else
    value = mean(values);
end

end
