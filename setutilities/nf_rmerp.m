function TF = nf_rmerp(TF, trlvec, mode)
% NF_RMERP  Separate modeled ERP and residual analytic TF activity.
%
% GENERAL
% -------
% NF_RMERP reconstructs the complex analytic coefficient represented by
% TF.power and TF.phase, estimates a condition-specific ERP contribution
% for every trial, and returns the modeled ERP power/phase and the power/
% phase of the complex residual.
%
% The function operates in the complex TF domain. For a linear analytic
% transform, the trial mean of the complex coefficients is the transform
% of the trial-averaged ERP. The function must not be used with a
% bilinear/quadratic TF distribution whose power and phase fields do not
% reconstruct a linear analytic coefficient.
%
% USAGE
% -----
%   TF = nf_rmerp(TF)
%   TF = nf_rmerp(TF, trlvec)
%   TF = nf_rmerp(TF, trlvec, mode)
%
% INPUTS
% ------
%   TF      NeuroFreq single-trial TF structure. Required fields are:
%             power   : linear power or base-10 log power
%             phase   : phase in radians
%             scale   : 'linear' or 'log10'
%             nsensor : number of sensors
%             ntrls   : number of trials; must exceed one
%             freqs   : frequency vector
%             times   : time vector
%             method  : NeuroFreq linear analytic transform name, or a
%                       custom transform accompanied by
%                       analyticCoefficientSemantics =
%                       'power_abs_squared_phase_angle_linear_coefficient'
%
%           The logical data dimensions are sensor x frequency x time x
%           trial. Canonical arrays and the exact singleton-squeezed shape
%           produced from that ordering are accepted and returned in their
%           original stored shape. Merely having the expected number of
%           elements is not sufficient: axis-permuted arrays are rejected.
%
%   trlvec  Trial condition labels. Numeric, logical, string, categorical,
%           and cellstr vectors are accepted. Missing labels are rejected.
%           Empty or omitted means that every trial belongs to one
%           condition.
%
%   mode    ERP model:
%             'standard' : the condition mean is removed from every trial.
%
%             'dot'      : each trial receives a nonnegative real scale of
%                          the condition mean. Scales are the joint least-
%                          squares solution constrained to sum to the
%                          number of condition trials. Therefore their mean
%                          is exactly one and the modeled ERPs average to
%                          the condition ERP.
%
%             'dotderiv' : each trial is modeled by
%
%                              a_r * mu + b_r * dmu/dt
%
%                          using real coefficients and a numerically stable
%                          centered least-squares fit. The constraints
%                          mean(a)=1 and mean(b)=0 are enforced within every
%                          condition. The derivative is orthogonalized to
%                          mu during fitting without changing the model
%                          span. For a small temporal shift,
%                          latencyShift approximately equals -b/a in the
%                          units of TF.times.
%
% OUTPUT FIELDS
% -------------
%   TF.power_evoked   Power of the modeled per-trial ERP contribution.
%   TF.phase_evoked   Phase of the modeled per-trial ERP contribution.
%   TF.power_induced  Power of the complex residual after ERP removal.
%   TF.phase_induced  Phase of the complex residual after ERP removal.
%   TF.erprem         Model, coefficient, provenance, and invariant ledger.
%
% IMPORTANT INTERPRETATION
% ------------------------
% - 'power_induced' is residual power under the selected ERP model. It is
%   not guaranteed to be purely induced neural activity.
% - All three modes use the analyzed trials to estimate their condition
%   ERP. The fit is not cross-validated.
% - The mean modeled complex ERP is required to equal the observed
%   condition ERP within a class-aware floating-point tolerance. The mean
%   complex residual is correspondingly zero.
% - Complex signals add exactly. Their powers generally do not because
%   total power contains a modeled-ERP/residual cross term.
% - Phase is undefined where the corresponding complex coefficient is
%   exactly zero. NF_RMERP stores NaN at those phase locations.

% greetings
disp(['[nf_rmerp]: single-trial evoked potential modeling using ' mode ' mode']);

if nargin < 1
    error('nf_rmerp:MissingTF', ...
        'a nonempty NeuroFreq TF structure is required.');
end

if isempty(TF)
    error('nf_rmerp:MissingTF', ...
        'a nonempty NeuroFreq TF structure is required.');
end

if nargin < 3
    mode = 'standard';
elseif isempty(mode)
    mode = 'standard';
end

mode = local_text_scalar(mode, 'mode');
mode = lower(strtrim(mode));

allowedModes = {'standard', 'dot', 'dotderiv'};

if any(strcmp(mode, allowedModes)) == false
    error('nf_rmerp:UnknownMode', ...
        'unknown mode ''%s''. allowed modes are standard, dot, and dotderiv.', ...
        mode);
end

[TF, power4, phase4, originalSize, dimensions, scale] = ...
    local_validate_tf(TF, mode);

if nargin < 2
    trlvec = [];
end

if isempty(trlvec)
    fprintf('[nf_rmerp]: no trial vector supplied - treating all trials as one condition.\n');
    trlvec = ones(dimensions.nTrials, 1);
end

[trialLabels, conditionLabels, conditionIndices] = ...
    local_group_trials(trlvec, dimensions.nTrials);

for conditionIndex = 1:numel(conditionIndices)
    if numel(conditionIndices{conditionIndex}) < 2
        warning('nf_rmerp:ConditionTooSmall', ...
            ['condition %d contains only one trial. under the exact-mean ' ...
            'constraint, that entire trial would be classified as ERP.'], ...
            conditionIndex);
        TF = [];
        return
    end
end

if strcmp(scale, 'linear') == true
    amplitude4 = sqrt(power4);
else
    amplitude4 = 10 .^ (power4 ./ 2);
end

if strcmp(scale, 'log10') == true
    underflowMask = isfinite(power4) & amplitude4 == 0;

    if any(underflowMask(:)) == true
        error('nf_rmerp:LogPowerUnderflow', ...
            ['finite TF.power values underflowed while reconstructing the ' ...
            'analytic amplitude. rescale the coefficients or use a wider ' ...
            'numeric representation.']);
    end
end

phaseForReconstruction = phase4;
phaseForReconstruction(isnan(phaseForReconstruction)) = 0;
analytic4 = amplitude4 .* exp(1i .* phaseForReconstruction);

clear amplitude4
clear phase4
clear phaseForReconstruction

if any(~isfinite(analytic4(:))) == true
    error('nf_rmerp:NonfiniteAnalyticSignal', ...
        ['power and phase reconstruction produced NaN or Inf analytic ' ...
        'coefficients. check the power range and scale.']);
end

nChannels = dimensions.nChannels;
nFrequencies = dimensions.nFrequencies;
nTrials = dimensions.nTrials;
nConditions = numel(conditionIndices);

powerEvoked4 = zeros(size(power4), 'like', power4);
powerInduced4 = zeros(size(power4), 'like', power4);
phaseEvoked4 = zeros(size(power4), 'like', TF.phase);
phaseInduced4 = zeros(size(power4), 'like', TF.phase);

amplitudeCoefficients = ones( ...
    nChannels, nFrequencies, nTrials, 'like', power4);
derivativeCoefficients = zeros( ...
    nChannels, nFrequencies, nTrials, 'like', power4);
latencyShiftApprox = nan( ...
    nChannels, nFrequencies, nTrials, 'like', power4);

conditionCounts = zeros(nConditions, 1);
conditionMeanModelError = zeros(nConditions, 1);
conditionMeanResidual = zeros(nConditions, 1);
conditionReconstructionError = zeros(nConditions, 1);
conditionInvariantTolerance = zeros(nConditions, 1);
conditionMeanAmplitudeError = zeros(nConditions, 1);
conditionMeanDerivativeError = zeros(nConditions, 1);
conditionAmplitudeConstraintTolerance = zeros(nConditions, 1);
conditionDerivativeConstraintTolerance = zeros(nConditions, 1);
conditionMinimumAmplitude = nan(nConditions, 1);
derivativeRank = zeros(nChannels, nFrequencies, nConditions, 'uint8');
derivativeConditionNumber = nan( ...
    nChannels, nFrequencies, nConditions);
derivativeOrthogonalization = zeros( ...
    nChannels, nFrequencies, nConditions, 'like', power4);
derivativeBasisRelativeNorm = nan( ...
    nChannels, nFrequencies, nConditions);
derivativeBasisColumnNorms = nan( ...
    nChannels, nFrequencies, 2, nConditions);
templateUnidentifiableCount = zeros(nConditions, 1);
derivativeBasisUnidentifiableCount = zeros(nConditions, 1);
derivativeNegativeAmplitudeCount = zeros(nConditions, 1);
derivativeUndefinedLatencyCount = zeros(nConditions, 1);
derivativeMaximumAbsoluteLatency = nan(nConditions, 1);
dotProjectionAdjustedCount = zeros(nConditions, 1);
dotUnidentifiableCount = zeros(nConditions, 1);
undefinedEvokedPhaseCount = zeros(nConditions, 1);
undefinedInducedPhaseCount = zeros(nConditions, 1);

epsilonValue = local_epsilon(power4);

for conditionIndex = 1:nConditions

    trialIndex = conditionIndices{conditionIndex};
    conditionCounts(conditionIndex) = numel(trialIndex);

    conditionSignal = analytic4(:, :, :, trialIndex);
    conditionMean = local_stable_mean(conditionSignal, 4);

    if strcmp(mode, 'standard') == true

        modeledERP = repmat( ...
            conditionMean, [1, 1, 1, numel(trialIndex)]);
        conditionAmplitude = ones( ...
            nChannels, nFrequencies, numel(trialIndex), 'like', power4);
        conditionDerivative = zeros( ...
            nChannels, nFrequencies, numel(trialIndex), 'like', power4);
        conditionLatency = nan( ...
            nChannels, nFrequencies, numel(trialIndex), 'like', power4);

    elseif strcmp(mode, 'dot') == true

        [modeledERP, conditionAmplitude, dotDiagnostics] = ...
            local_fit_dot(conditionSignal, conditionMean);
        conditionDerivative = zeros( ...
            nChannels, nFrequencies, numel(trialIndex), 'like', power4);
        conditionLatency = nan( ...
            nChannels, nFrequencies, numel(trialIndex), 'like', power4);
        dotProjectionAdjustedCount(conditionIndex) = ...
            dotDiagnostics.projectionAdjustedCount;
        dotUnidentifiableCount(conditionIndex) = ...
            dotDiagnostics.unidentifiableCount;

    else

        [modeledERP, conditionAmplitude, conditionDerivative, ...
            conditionLatency, derivativeDiagnostics] = ...
            local_fit_dotderiv( ...
            conditionSignal, conditionMean, TF.times);
        derivativeRank(:, :, conditionIndex) = ...
            derivativeDiagnostics.rank;
        derivativeConditionNumber(:, :, conditionIndex) = ...
            derivativeDiagnostics.conditionNumber;
        derivativeOrthogonalization(:, :, conditionIndex) = ...
            derivativeDiagnostics.orthogonalization;
        derivativeBasisRelativeNorm(:, :, conditionIndex) = ...
            derivativeDiagnostics.relativeNorm;
        derivativeBasisColumnNorms(:, :, :, conditionIndex) = ...
            derivativeDiagnostics.columnNorms;
        templateUnidentifiableCount(conditionIndex) = ...
            nnz(derivativeDiagnostics.templateUnidentifiable);
        derivativeBasisUnidentifiableCount(conditionIndex) = ...
            nnz(derivativeDiagnostics.rank < 2);

    end


    if any(~isfinite(modeledERP(:))) == true
        error('nf_rmerp:NonfiniteModeledERP', ...
            ['ERP fitting produced NaN or Inf modeled coefficients in ' ...
            'condition %d.'], ...
            conditionIndex);
    end

    if any(~isfinite(conditionAmplitude(:))) == true
        error('nf_rmerp:NonfiniteModelCoefficient', ...
            ['ERP fitting produced a nonfinite amplitude coefficient in ' ...
            'condition %d.'], ...
            conditionIndex);
    end

    if any(~isfinite(conditionDerivative(:))) == true
        error('nf_rmerp:NonfiniteModelCoefficient', ...
            ['ERP fitting produced a nonfinite derivative coefficient in ' ...
            'condition %d.'], ...
            conditionIndex);
    end

    residualSignal = conditionSignal - modeledERP;

    if any(~isfinite(residualSignal(:))) == true
        error('nf_rmerp:NonfiniteResidual', ...
            ['ERP subtraction produced NaN or Inf residual coefficients in ' ...
            'condition %d.'], ...
            conditionIndex);
    end

    meanModelError = local_stable_mean(modeledERP, 4) - conditionMean;
    meanResidual = local_stable_mean(residualSignal, 4);
    reconstructionError = ...
        conditionSignal - (modeledERP + residualSignal);

    conditionMeanModelError(conditionIndex) = ...
        local_max_abs(meanModelError);
    conditionMeanResidual(conditionIndex) = ...
        local_max_abs(meanResidual);
    conditionReconstructionError(conditionIndex) = ...
        local_max_abs(reconstructionError);

    conditionScale = local_max_abs(conditionSignal);
    conditionScale = max( ...
        conditionScale, local_smallest_positive(power4));
    conditionInvariantTolerance(conditionIndex) = ...
        64 .* epsilonValue .* max(1, numel(trialIndex)) .* conditionScale;

    amplitudeMeanError = mean(conditionAmplitude, 3) - 1;
    derivativeMeanError = mean(conditionDerivative, 3);
    conditionMeanAmplitudeError(conditionIndex) = ...
        local_max_abs(amplitudeMeanError);
    conditionMeanDerivativeError(conditionIndex) = ...
        local_max_abs(derivativeMeanError);
    conditionMinimumAmplitude(conditionIndex) = ...
        double(min(conditionAmplitude(:)));

    conditionAmplitudeConstraintTolerance(conditionIndex) = ...
        64 .* epsilonValue .* max(1, numel(trialIndex)) .* ...
        max(1, local_max_abs(conditionAmplitude));
    conditionDerivativeConstraintTolerance(conditionIndex) = ...
        64 .* epsilonValue .* max(1, numel(trialIndex)) .* ...
        max(1, local_max_abs(conditionDerivative));

    if conditionMeanAmplitudeError(conditionIndex) > ...
            conditionAmplitudeConstraintTolerance(conditionIndex)
        error('nf_rmerp:AmplitudeConstraintFailed', ...
            ['mean amplitude failed its within-condition constraint in ' ...
            'condition %d.'], ...
            conditionIndex);
    end

    if conditionMeanDerivativeError(conditionIndex) > ...
            conditionDerivativeConstraintTolerance(conditionIndex)
        error('nf_rmerp:DerivativeConstraintFailed', ...
            ['mean derivative coefficient failed its within-condition ' ...
            'constraint in condition %d.'], ...
            conditionIndex);
    end

    if strcmp(mode, 'dot') == true
        if conditionMinimumAmplitude(conditionIndex) < ...
                -conditionAmplitudeConstraintTolerance(conditionIndex)
            error('nf_rmerp:NonnegativeAmplitudeConstraintFailed', ...
                ['the dot model produced a materially negative amplitude ' ...
                'in condition %d.'], ...
                conditionIndex);
        end
    end

    if strcmp(mode, 'dotderiv') == true
        derivativeNegativeAmplitudeCount(conditionIndex) = ...
            nnz(conditionAmplitude < 0);
        derivativeUndefinedLatencyCount(conditionIndex) = ...
            nnz(~isfinite(conditionLatency));
        finiteLatency = abs(conditionLatency(isfinite(conditionLatency)));

        if isempty(finiteLatency) == false
            derivativeMaximumAbsoluteLatency(conditionIndex) = ...
                double(max(finiteLatency));
        end
    end

    if conditionMeanModelError(conditionIndex) > ...
            conditionInvariantTolerance(conditionIndex)
        error('nf_rmerp:MeanModelInvariantFailed', ...
            ['the modeled ERP failed the exact-mean invariant in ' ...
            'condition %d. maximum error was %.17g; tolerance was %.17g.'], ...
            conditionIndex, ...
            conditionMeanModelError(conditionIndex), ...
            conditionInvariantTolerance(conditionIndex));
    end

    if conditionMeanResidual(conditionIndex) > ...
            conditionInvariantTolerance(conditionIndex)
        error('nf_rmerp:ResidualMeanInvariantFailed', ...
            ['the residual failed the zero-mean invariant in condition %d. ' ...
            'maximum residual mean was %.17g; tolerance was %.17g.'], ...
            conditionIndex, ...
            conditionMeanResidual(conditionIndex), ...
            conditionInvariantTolerance(conditionIndex));
    end

    if conditionReconstructionError(conditionIndex) > ...
            conditionInvariantTolerance(conditionIndex)
        error('nf_rmerp:ReconstructionInvariantFailed', ...
            ['modeled ERP plus residual did not reconstruct the analytic ' ...
            'signal in condition %d. maximum error was %.17g.'], ...
            conditionIndex, ...
            conditionReconstructionError(conditionIndex));
    end

    [conditionEvokedPower, conditionEvokedPhase, evokedZeroMask] = ...
        local_encode_signal(modeledERP, scale, TF.power, TF.phase);
    [conditionInducedPower, conditionInducedPhase, inducedZeroMask] = ...
        local_encode_signal(residualSignal, scale, TF.power, TF.phase);

    powerEvoked4(:, :, :, trialIndex) = conditionEvokedPower;
    phaseEvoked4(:, :, :, trialIndex) = conditionEvokedPhase;
    powerInduced4(:, :, :, trialIndex) = conditionInducedPower;
    phaseInduced4(:, :, :, trialIndex) = conditionInducedPhase;

    amplitudeCoefficients(:, :, trialIndex) = conditionAmplitude;
    derivativeCoefficients(:, :, trialIndex) = conditionDerivative;
    latencyShiftApprox(:, :, trialIndex) = conditionLatency;

    undefinedEvokedPhaseCount(conditionIndex) = nnz(evokedZeroMask);
    undefinedInducedPhaseCount(conditionIndex) = nnz(inducedZeroMask);

end

TF.power_evoked = reshape(powerEvoked4, originalSize);
TF.phase_evoked = reshape(phaseEvoked4, originalSize);
TF.power_induced = reshape(powerInduced4, originalSize);
TF.phase_induced = reshape(phaseInduced4, originalSize);

TF.erprem = struct();
TF.erprem.schema = 'nf_rmerp_v2';
TF.erprem.method = 'nf_rmerp';
TF.erprem.mode = mode;
TF.erprem.scale = scale;

if isfield(TF, 'method') == true
    TF.erprem.inputTransform = TF.method;
else
    TF.erprem.inputTransform = 'custom';
end

if isfield(TF, 'analyticCoefficientSemantics') == true
    TF.erprem.inputTransformVerification = ...
        TF.analyticCoefficientSemantics;
else
    TF.erprem.inputTransformVerification = 'verified NeuroFreq method';
end

TF.erprem.trlvec = trialLabels;
TF.erprem.conditionLabels = conditionLabels;
TF.erprem.conditionCounts = conditionCounts;
TF.erprem.coefficientScope = 'channel_frequency';
TF.erprem.coefficients = struct();
TF.erprem.coefficients.amplitude = amplitudeCoefficients;
TF.erprem.coefficients.derivative = derivativeCoefficients;
TF.erprem.coefficients.latencyShiftApprox = latencyShiftApprox;
TF.erprem.constraints = struct();
TF.erprem.constraints.meanAmplitude = 1;
TF.erprem.constraints.meanDerivative = 0;
TF.erprem.constraints.withinCondition = true;
TF.erprem.constraints.modelMeanEqualsObservedERP = true;
TF.erprem.diagnostics = struct();
TF.erprem.diagnostics.conditionMeanModelMaximumAbsoluteError = ...
    conditionMeanModelError;
TF.erprem.diagnostics.conditionMeanResidualMaximumAbsoluteValue = ...
    conditionMeanResidual;
TF.erprem.diagnostics.conditionReconstructionMaximumAbsoluteError = ...
    conditionReconstructionError;
TF.erprem.diagnostics.conditionInvariantTolerance = ...
    conditionInvariantTolerance;
TF.erprem.diagnostics.conditionMeanAmplitudeConstraintError = ...
    conditionMeanAmplitudeError;
TF.erprem.diagnostics.conditionMeanDerivativeConstraintError = ...
    conditionMeanDerivativeError;
TF.erprem.diagnostics.conditionAmplitudeConstraintTolerance = ...
    conditionAmplitudeConstraintTolerance;
TF.erprem.diagnostics.conditionDerivativeConstraintTolerance = ...
    conditionDerivativeConstraintTolerance;
TF.erprem.diagnostics.conditionMinimumAmplitude = ...
    conditionMinimumAmplitude;
TF.erprem.diagnostics.dotProjectionAdjustedCount = ...
    dotProjectionAdjustedCount;
TF.erprem.diagnostics.dotUnidentifiableCount = ...
    dotUnidentifiableCount;
TF.erprem.diagnostics.derivativeDesignRank = derivativeRank;
TF.erprem.diagnostics.derivativeDesignConditionNumber = ...
    derivativeConditionNumber;
TF.erprem.diagnostics.derivativeOrthogonalization = ...
    derivativeOrthogonalization;
TF.erprem.diagnostics.derivativeBasisRelativeNorm = ...
    derivativeBasisRelativeNorm;
TF.erprem.diagnostics.derivativeBasisColumnNorms = ...
    derivativeBasisColumnNorms;
TF.erprem.diagnostics.templateUnidentifiableCount = ...
    templateUnidentifiableCount;
TF.erprem.diagnostics.derivativeBasisUnidentifiableCount = ...
    derivativeBasisUnidentifiableCount;
TF.erprem.diagnostics.derivativeNegativeAmplitudeCount = ...
    derivativeNegativeAmplitudeCount;
TF.erprem.diagnostics.derivativeUndefinedLatencyCount = ...
    derivativeUndefinedLatencyCount;
TF.erprem.diagnostics.derivativeMaximumAbsoluteLatency = ...
    derivativeMaximumAbsoluteLatency;
TF.erprem.diagnostics.undefinedEvokedPhaseCount = ...
    undefinedEvokedPhaseCount;
TF.erprem.diagnostics.undefinedInducedPhaseCount = ...
    undefinedInducedPhaseCount;
TF.erprem.power_evoked_field = 'power_evoked';
TF.erprem.phase_evoked_field = 'phase_evoked';
TF.erprem.power_induced_field = 'power_induced';
TF.erprem.phase_induced_field = 'phase_induced';
if strcmp(scale, 'linear') == true
    TF.erprem.analyticReconstruction = ...
        'sqrt(power).*exp(1i.*phase), with NaN phase interpreted as zero';
else
    TF.erprem.analyticReconstruction = ...
        '10.^(power./2).*exp(1i.*phase), with NaN phase interpreted as zero';
end

TF.erprem.fitUsesAnalyzedTrialsInTemplate = true;
TF.erprem.crossValidated = false;
TF.erprem.residualInterpretation = ...
    'model-dependent residual; not guaranteed purely induced activity';
TF.erprem.powerAdditivity = ...
    ['abs(total).^2 = abs(modeled).^2 + abs(residual).^2 + ' ...
    '2*real(modeled.*conj(residual))'];
TF.erprem.zeroPowerPhasePolicy = ...
    'phase is NaN where the modeled complex coefficient is exactly zero';

if strcmp(scale, 'log10') == true
    TF.erprem.log10ZeroPowerValue = -Inf;
else
    TF.erprem.log10ZeroPowerValue = [];
end

if strcmp(mode, 'dotderiv') == true
    TF.erprem.derivative = struct();
    TF.erprem.derivative.timeAxis = TF.times(:);
    TF.erprem.derivative.timeUnits = 'units of TF.times';
    TF.erprem.derivative.coefficientsAreReal = true;
    TF.erprem.derivative.latencyApproximation = '-b/a';
    TF.erprem.derivative.fitUsesOrthogonalizedDerivative = true;
    TF.erprem.derivative.normalizedTimeBasis = true;
    TF.erprem.derivative.coefficientsAreUnbounded = true;
    TF.erprem.derivative.warning = ...
        ['Negative amplitudes and large approximate latency shifts are ' ...
        'reported, not silently clipped; inspect diagnostics before ' ...
        'interpreting coefficients physiologically.'];
end

TF = nf_powerfront_local(TF);

end


function [TF, power4, phase4, originalSize, dimensions, scale] = ...
    local_validate_tf(TF, mode)

if isstruct(TF) == false
    error('nf_rmerp:BadTFType', ...
        'TF must be a scalar NeuroFreq structure.');
end

if isscalar(TF) == false
    error('nf_rmerp:BadTFType', ...
        'TF must be a scalar NeuroFreq structure.');
end

requiredFields = { ...
    'power', ...
    'phase', ...
    'scale', ...
    'ntrls', ...
    'nsensor', ...
    'freqs', ...
    'times'};

for fieldIndex = 1:numel(requiredFields)
    fieldName = requiredFields{fieldIndex};

    if isfield(TF, fieldName) == false
        error('nf_rmerp:MissingField', ...
            'TF.%s is required.', fieldName);
    end
end

reservedFields = { ...
    'power_evoked', ...
    'phase_evoked', ...
    'power_induced', ...
    'phase_induced', ...
    'erprem'};

for fieldIndex = 1:numel(reservedFields)
    fieldName = reservedFields{fieldIndex};

    if isfield(TF, fieldName) == true
        error('nf_rmerp:ExistingERPResult', ...
            ['TF already contains reserved ERP-result field ''%s''. ' ...
            'Remove prior NF_RMERP or CPM results before fitting.'], ...
            fieldName);
    end
end

if isnumeric(TF.power) == false
    error('nf_rmerp:BadPower', ...
        'TF.power must be a real floating-point numeric array.');
end

if isfloat(TF.power) == false
    error('nf_rmerp:BadPower', ...
        'TF.power must be a real floating-point numeric array.');
end

if isreal(TF.power) == false
    error('nf_rmerp:BadPower', ...
        'TF.power must be real.');
end

if issparse(TF.power) == true
    error('nf_rmerp:BadPower', ...
        'TF.power must be a full numeric array.');
end

if isnumeric(TF.phase) == false
    error('nf_rmerp:BadPhase', ...
        'TF.phase must be a real floating-point numeric array.');
end

if isfloat(TF.phase) == false
    error('nf_rmerp:BadPhase', ...
        'TF.phase must be a real floating-point numeric array.');
end

if isreal(TF.phase) == false
    error('nf_rmerp:BadPhase', ...
        'TF.phase must be real.');
end

if issparse(TF.phase) == true
    error('nf_rmerp:BadPhase', ...
        'TF.phase must be a full numeric array.');
end

if strcmp(class(TF.power), class(TF.phase)) == false
    error('nf_rmerp:PrecisionMismatch', ...
        'TF.power and TF.phase must use the same floating-point class.');
end

if isequal(size(TF.power), size(TF.phase)) == false
    error('nf_rmerp:PowerPhaseShapeMismatch', ...
        'TF.power and TF.phase must have identical stored dimensions.');
end

if isempty(TF.power) == true
    error('nf_rmerp:EmptyPower', ...
        'TF.power and TF.phase must be nonempty.');
end

scale = local_text_scalar(TF.scale, 'TF.scale');
scale = lower(strtrim(scale));

if any(strcmp(scale, {'linear', 'log10'})) == false
    error('nf_rmerp:BadScale', ...
        'TF.scale must be ''linear'' or ''log10''.');
end

if strcmp(scale, 'linear') == true
    if any(~isfinite(TF.power(:))) == true
        error('nf_rmerp:NonfinitePower', ...
            'Linear TF.power cannot contain NaN or Inf values.');
    end

    if any(TF.power(:) < 0) == true
        error('nf_rmerp:NegativeLinearPower', ...
            ['Linear power must be nonnegative. Negative values cannot ' ...
            'reconstruct an analytic coefficient.']);
    end

    exactZeroInput = TF.power == 0;
else
    invalidLogPower = isnan(TF.power) | TF.power == Inf;

    if any(invalidLogPower(:)) == true
        error('nf_rmerp:NonfinitePower', ...
            ['Log10 TF.power may contain -Inf only to represent exact ' ...
            'zero power; NaN and +Inf are invalid.']);
    end

    exactZeroInput = TF.power == -Inf;
end


if any(isinf(TF.phase(:))) == true
    error('nf_rmerp:NonfinitePhase', ...
        'TF.phase cannot contain Inf values.');
end

invalidUndefinedPhase = isnan(TF.phase) & ~exactZeroInput;

if any(invalidUndefinedPhase(:)) == true
    error('nf_rmerp:NonfinitePhase', ...
        ['TF.phase may be NaN only where the corresponding analytic ' ...
        'coefficient has exactly zero power.']);
end

local_validate_positive_integer(TF.ntrls, 'TF.ntrls');
local_validate_positive_integer(TF.nsensor, 'TF.nsensor');

nTrials = double(TF.ntrls);
nChannels = double(TF.nsensor);

if nTrials < 2
    error('nf_rmerp:AveragedInput', ...
        'ERP removal requires at least two single trials.');
end

if isnumeric(TF.freqs) == false
    error('nf_rmerp:BadFrequencies', ...
        'TF.freqs must be a finite real numeric vector.');
end

if isvector(TF.freqs) == false
    error('nf_rmerp:BadFrequencies', ...
        'TF.freqs must be a finite real numeric vector.');
end

if isreal(TF.freqs) == false
    error('nf_rmerp:BadFrequencies', ...
        'TF.freqs must be a finite real numeric vector.');
end

if isempty(TF.freqs) == true
    error('nf_rmerp:BadFrequencies', ...
        'TF.freqs must be nonempty.');
end

if any(~isfinite(TF.freqs(:))) == true
    error('nf_rmerp:BadFrequencies', ...
        'TF.freqs must contain only finite values.');
end

if isnumeric(TF.times) == false
    error('nf_rmerp:BadTimes', ...
        'TF.times must be a finite real numeric vector.');
end

if isvector(TF.times) == false
    error('nf_rmerp:BadTimes', ...
        'TF.times must be a finite real numeric vector.');
end

if isreal(TF.times) == false
    error('nf_rmerp:BadTimes', ...
        'TF.times must be a finite real numeric vector.');
end

if isempty(TF.times) == true
    error('nf_rmerp:BadTimes', ...
        'TF.times must be nonempty.');
end

if any(~isfinite(TF.times(:))) == true
    error('nf_rmerp:BadTimes', ...
        'TF.times must contain only finite values.');
end

if strcmp(mode, 'dotderiv') == true
    if numel(TF.times) < 3
        error('nf_rmerp:TooFewDerivativeSamples', ...
            'dotderiv requires at least three time samples.');
    end

    timeDifferences = diff(double(TF.times(:)));

    if any(~isfinite(timeDifferences)) == true
        error('nf_rmerp:InvalidDerivativeTimeScale', ...
            ['Differences in TF.times overflowed the numeric range. ' ...
            'Rescale the time axis before derivative fitting.']);
    end

    if any(timeDifferences <= 0) == true
        error('nf_rmerp:NonmonotonicTimes', ...
            ['TF.times must be strictly increasing for derivative ' ...
            'ERP fitting.']);
    end

    timeSpan = double(TF.times(end)) - double(TF.times(1));

    if isfinite(timeSpan) == false || timeSpan <= 0
        error('nf_rmerp:InvalidDerivativeTimeScale', ...
            'The TF.times span must be finite and positive.');
    end
end

nFrequencies = numel(TF.freqs);
nTimes = numel(TF.times);
expectedElements = nChannels .* nFrequencies .* nTimes .* nTrials;

if numel(TF.power) ~= expectedElements    %#ok
    error('nf_rmerp:DimensionMismatch', ...
        ['TF.power/phase contain %d values, but nsensor, freqs, times, ' ...
        'and ntrls require %d (%d x %d x %d x %d).'], ...
        numel(TF.power), ...
        expectedElements, ...
        nChannels, ...
        nFrequencies, ...
        nTimes, ...
        nTrials);
end


canonicalSize = [nChannels, nFrequencies, nTimes, nTrials];
squeezedSize = local_squeezed_size(canonicalSize);
storedSize = size(TF.power);
validCanonicalShape = isequal(storedSize, canonicalSize);
validSqueezedShape = isequal(storedSize, squeezedSize);

if validCanonicalShape == false && validSqueezedShape == false
    error('nf_rmerp:DimensionOrderMismatch', ...
        ['TF.power has stored size [%s]. Expected canonical size [%s] or ' ...
        'its exact singleton-squeezed size [%s]. Axis permutation is not ' ...
        'inferred from element count.'], ...
        local_size_text(storedSize), ...
        local_size_text(canonicalSize), ...
        local_size_text(squeezedSize));
end

if isfield(TF, 'chanlocs') == true
    if isempty(TF.chanlocs) == false
        if numel(TF.chanlocs) ~= nChannels
            error('nf_rmerp:ChannelMetadataMismatch', ...
                ['TF.chanlocs contains %d channels, whereas TF.nsensor ' ...
                'equals %d.'], ...
                numel(TF.chanlocs), ...
                nChannels);
        end
    end
end

requiredSemantics = ...
    'power_abs_squared_phase_angle_linear_coefficient';
hasVerifiedSemantics = false;

if isfield(TF, 'analyticCoefficientSemantics') == true
    semantics = local_text_scalar( ...
        TF.analyticCoefficientSemantics, ...
        'TF.analyticCoefficientSemantics');

    if strcmpi(strtrim(semantics), requiredSemantics) == false
        error('nf_rmerp:UnverifiedTransform', ...
            ['TF.analyticCoefficientSemantics must equal ''%s'' for a ' ...
            'custom transform.'], ...
            requiredSemantics);
    end

    hasVerifiedSemantics = true;
end

if isfield(TF, 'method') == true
    method = local_text_scalar(TF.method, 'TF.method');
    methodKey = lower(strtrim(method));
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
        error('nf_rmerp:NonanalyticTransform', ...
            ['TF.method ''%s'' does not provide a linear analytic ' ...
            'coefficient reconstructible as sqrt(power).*exp(1i.*phase).'], ...
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

    if any(strcmp(methodKey, supportedMethods)) == false && ...
            hasVerifiedSemantics == false
        error('nf_rmerp:UnverifiedTransform', ...
            ['TF.method ''%s'' is not a verified NeuroFreq linear analytic ' ...
            'transform. Custom transforms must set ' ...
            'TF.analyticCoefficientSemantics to the documented token.'], ...
            method);
    end
elseif hasVerifiedSemantics == false
    error('nf_rmerp:UnverifiedTransform', ...
        ['TF.method or an explicit TF.analyticCoefficientSemantics token ' ...
        'is required. NF_RMERP cannot infer linear analytic semantics ' ...
        'from power and phase field names.']);
end

originalSize = size(TF.power);

power4 = reshape( ...
    TF.power, [nChannels, nFrequencies, nTimes, nTrials]);
phase4 = reshape( ...
    TF.phase, [nChannels, nFrequencies, nTimes, nTrials]);

TF.scale = scale;

dimensions = struct();
dimensions.nChannels = nChannels;
dimensions.nFrequencies = nFrequencies;
dimensions.nTimes = nTimes;
dimensions.nTrials = nTrials;

end


function local_validate_positive_integer(value, name)

if isnumeric(value) == false
    error('nf_rmerp:BadDimensionMetadata', ...
        '%s must be a positive integer scalar.', name);
end

if isreal(value) == false
    error('nf_rmerp:BadDimensionMetadata', ...
        '%s must be a positive integer scalar.', name);
end

if isscalar(value) == false
    error('nf_rmerp:BadDimensionMetadata', ...
        '%s must be a positive integer scalar.', name);
end

if isfinite(value) == false
    error('nf_rmerp:BadDimensionMetadata', ...
        '%s must be a positive integer scalar.', name);
end

if value < 1
    error('nf_rmerp:BadDimensionMetadata', ...
        '%s must be a positive integer scalar.', name);
end

if value ~= round(value)
    error('nf_rmerp:BadDimensionMetadata', ...
        '%s must be a positive integer scalar.', name);
end

end


function textValue = local_text_scalar(value, name)

if isstring(value) == true
    if isscalar(value) == false
        error('nf_rmerp:BadTextScalar', ...
            '%s must be a character vector or scalar string.', name);
    end

    if ismissing(value) == true
        error('nf_rmerp:BadTextScalar', ...
            '%s cannot be missing.', name);
    end

    textValue = char(value);
elseif ischar(value) == true
    if isrow(value) == false
        error('nf_rmerp:BadTextScalar', ...
            '%s must be a character vector or scalar string.', name);
    end

    textValue = value;
else
    error('nf_rmerp:BadTextScalar', ...
        '%s must be a character vector or scalar string.', name);
end

if isempty(strtrim(textValue)) == true
    error('nf_rmerp:BadTextScalar', ...
        '%s cannot be empty text.', name);
end

end


function [trialLabels, conditionLabels, conditionIndices] = ...
    local_group_trials(trlvec, nTrials)

if isvector(trlvec) == false
    error('nf_rmerp:BadTrialVector', ...
        'trlvec must be a vector with one nonmissing label per trial.');
end

if numel(trlvec) ~= nTrials
    error('nf_rmerp:TrialVectorLength', ...
        'trlvec has %d labels; TF contains %d trials.', ...
        numel(trlvec), ...
        nTrials);
end

if isnumeric(trlvec) == true

    if isreal(trlvec) == false
        error('nf_rmerp:BadTrialVector', ...
            'Numeric trial labels must be real and finite.');
    end

    if any(~isfinite(trlvec(:))) == true
        error('nf_rmerp:MissingTrialLabel', ...
            'Numeric trial labels cannot contain NaN or Inf.');
    end

    trialLabels = trlvec(:);

elseif islogical(trlvec) == true

    trialLabels = trlvec(:);

elseif isstring(trlvec) == true

    trialLabels = trlvec(:);

    if any(ismissing(trialLabels)) == true
        error('nf_rmerp:MissingTrialLabel', ...
            'String trial labels cannot contain missing values.');
    end

    if any(strlength(trialLabels) == 0) == true
        error('nf_rmerp:MissingTrialLabel', ...
            'String trial labels cannot contain empty labels.');
    end

elseif iscategorical(trlvec) == true

    trialLabels = trlvec(:);

    if any(isundefined(trialLabels)) == true
        error('nf_rmerp:MissingTrialLabel', ...
            'Categorical trial labels cannot contain undefined values.');
    end

elseif iscell(trlvec) == true

    if iscellstr(trlvec) == false %#ok
        error('nf_rmerp:BadTrialVector', ...
            'Cell trial labels must be a cell array of character vectors.');
    end

    trialLabels = string(trlvec(:));

    if any(strlength(trialLabels) == 0) == true
        error('nf_rmerp:MissingTrialLabel', ...
            'Cellstr trial labels cannot contain empty labels.');
    end

else

    error('nf_rmerp:BadTrialVector', ...
        ['trlvec must be numeric, logical, string, categorical, or ' ...
        'cellstr.']);

end

conditionLabels = unique(trialLabels, 'stable');
conditionIndices = cell(numel(conditionLabels), 1);

for conditionIndex = 1:numel(conditionLabels)
    conditionIndices{conditionIndex} = ...
        find(trialLabels == conditionLabels(conditionIndex));

    if isempty(conditionIndices{conditionIndex}) == true
        error('nf_rmerp:InternalConditionGrouping', ...
            'A condition label produced an empty trial set.');
    end
end

end


function [erp, amplitude, diagnostics] = local_fit_dot(X, mu)

[nChannels, nFrequencies, nTimes, nTrials] = size(X);

erp = zeros(size(X), 'like', X);
amplitude = zeros( ...
    nChannels, nFrequencies, nTrials, 'like', real(X));
diagnostics = struct();
diagnostics.projectionAdjustedCount = 0;
diagnostics.unidentifiableCount = 0;

epsilonValue = local_epsilon(real(X));
minimumValue = local_realmin(real(X));

for channelIndex = 1:nChannels
    for frequencyIndex = 1:nFrequencies

        trialMatrix = reshape( ...
            X(channelIndex, frequencyIndex, :, :), ...
            [nTimes, nTrials]);
        meanWaveform = reshape( ...
            mu(channelIndex, frequencyIndex, :), ...
            [nTimes, 1]);

        sliceScale = max(abs(trialMatrix(:)));

        if sliceScale == 0
            scaledTrialMatrix = trialMatrix;
            scaledMeanWaveform = meanWaveform;
        else
            scaledTrialMatrix = trialMatrix ./ sliceScale;
            scaledMeanWaveform = meanWaveform ./ sliceScale;
        end

        meanEnergy = real( ...
            scaledMeanWaveform' * scaledMeanWaveform);
        trialEnergy = mean(sum(abs(scaledTrialMatrix) .^ 2, 1));
        energyTolerance = sqrt(epsilonValue) .* ...
            max(double(trialEnergy), double(minimumValue));

        if double(meanEnergy) <= energyTolerance %#ok

            weights = ones(nTrials, 1, 'like', real(X));
            diagnostics.unidentifiableCount = ...
                diagnostics.unidentifiableCount + 1;

        else

            rawWeights = real( ...
                scaledMeanWaveform' * scaledTrialMatrix) ./ meanEnergy;
            rawWeights = rawWeights(:);
            weights = local_project_simplex(rawWeights, nTrials);

            adjustment = max(abs(weights - rawWeights));
            adjustmentTolerance = 32 .* epsilonValue .* ...
                max(1, double(max(abs(rawWeights))));

            if double(adjustment) > adjustmentTolerance %#ok
                diagnostics.projectionAdjustedCount = ...
                    diagnostics.projectionAdjustedCount + 1;
            end

        end

        amplitude(channelIndex, frequencyIndex, :) = ...
            reshape(weights, [1, 1, nTrials]);
        erp(channelIndex, frequencyIndex, :, :) = ...
            reshape( ...
            meanWaveform * weights.', ...
            [1, 1, nTimes, nTrials]);

    end
end

end


function projected = local_project_simplex(values, targetSum)

values = values(:);
nValues = numel(values);

if nValues == 0
    error('nf_rmerp:InternalSimplexProjection', ...
        'Simplex projection received an empty vector.');
end

if targetSum <= 0
    error('nf_rmerp:InternalSimplexProjection', ...
        'Simplex projection target must be positive.');
end

sortedValues = sort(values, 'descend');
cumulative = cumsum(sortedValues) - targetSum;
indices = (1:nValues).';
active = find(sortedValues - cumulative ./ indices > 0, 1, 'last');

if isempty(active) == true
    projected = ones(size(values), 'like', values);
    projected = projected .* (targetSum ./ nValues);
else
    threshold = cumulative(active) ./ active;
    projected = max(values - threshold, 0);
end

sumCorrection = targetSum - sum(projected);
[~, maximumIndex] = max(projected);
projected(maximumIndex) = projected(maximumIndex) + sumCorrection;

epsilonValue = local_epsilon(projected);
negativeTolerance = 64 .* epsilonValue .* max(1, targetSum);

if any(projected < -negativeTolerance) == true
    error('nf_rmerp:InternalSimplexProjection', ...
        'Simplex projection produced a materially negative coefficient.');
end

projected(projected < 0) = 0;
sumCorrection = targetSum - sum(projected);
[~, maximumIndex] = max(projected);
projected(maximumIndex) = projected(maximumIndex) + sumCorrection;

end


function [erp, amplitude, derivative, latency, diagnostics] = ...
    local_fit_dotderiv(X, mu, times)

[nChannels, nFrequencies, nTimes, nTrials] = size(X);

erp = zeros(size(X), 'like', X);
amplitude = zeros( ...
    nChannels, nFrequencies, nTrials, 'like', real(X));
derivative = zeros( ...
    nChannels, nFrequencies, nTrials, 'like', real(X));
latency = nan( ...
    nChannels, nFrequencies, nTrials, 'like', real(X));

diagnostics = struct();
diagnostics.rank = zeros(nChannels, nFrequencies, 'uint8');
diagnostics.conditionNumber = nan(nChannels, nFrequencies);
diagnostics.orthogonalization = zeros( ...
    nChannels, nFrequencies, 'like', real(X));
diagnostics.relativeNorm = nan(nChannels, nFrequencies);
diagnostics.columnNorms = nan(nChannels, nFrequencies, 2);
diagnostics.templateUnidentifiable = false(nChannels, nFrequencies);

epsilonValue = local_epsilon(real(X));
minimumValue = local_realmin(real(X));
timeVector = double(times(:));
timeSpan = timeVector(end) - timeVector(1);
normalizedTimes = (timeVector - timeVector(1)) ./ timeSpan;

for channelIndex = 1:nChannels
    for frequencyIndex = 1:nFrequencies

        trialMatrix = reshape( ...
            X(channelIndex, frequencyIndex, :, :), ...
            [nTimes, nTrials]);
        meanWaveform = reshape( ...
            mu(channelIndex, frequencyIndex, :), ...
            [nTimes, 1]);

        sliceScale = max(abs(trialMatrix(:)));

        if sliceScale == 0
            scaledTrialMatrix = trialMatrix;
            scaledMeanWaveform = meanWaveform;
        else
            scaledTrialMatrix = trialMatrix ./ sliceScale;
            scaledMeanWaveform = meanWaveform ./ sliceScale;
        end

        scaledMeanArray = reshape( ...
            scaledMeanWaveform, [1, 1, nTimes]);
        normalizedDerivativeArray = ...
            local_time_derivative(scaledMeanArray, normalizedTimes);
        normalizedDerivativeWaveform = reshape( ...
            normalizedDerivativeArray, [nTimes, 1]);

        meanEnergy = real( ...
            scaledMeanWaveform' * scaledMeanWaveform);
        trialEnergy = mean(sum(abs(scaledTrialMatrix) .^ 2, 1));
        energyTolerance = sqrt(epsilonValue) .* ...
            max(double(trialEnergy), double(minimumValue));

        if double(meanEnergy) <= energyTolerance %#ok

            amplitudeRow = ones(1, nTrials, 'like', real(X));
            derivativeRow = zeros(1, nTrials, 'like', real(X));
            latencyRow = nan(1, nTrials, 'like', real(X));
            modeledMatrix = repmat(meanWaveform, [1, nTrials]);
            designRank = 0;
            designConditionNumber = Inf;
            derivativeProjection = 0;
            relativeDerivativeNorm = NaN;
            designColumnNorms = [0, 0];
            diagnostics.templateUnidentifiable( ...
                channelIndex, frequencyIndex) = true;

        else

            derivativeProjection = ...
                real(scaledMeanWaveform' * ...
                normalizedDerivativeWaveform) ./ meanEnergy;
            orthogonalDerivative = ...
                normalizedDerivativeWaveform - ...
                derivativeProjection .* scaledMeanWaveform;
            centeredTrials = ...
                scaledTrialMatrix - ...
                repmat(scaledMeanWaveform, [1, nTrials]);
            complexDesign = [scaledMeanWaveform, orthogonalDerivative];
            realDesign = [real(complexDesign); imag(complexDesign)];
            realOutcome = [real(centeredTrials); imag(centeredTrials)];

            [designRank, designConditionNumber, designColumnNorms] = ...
                local_design_diagnostics(realDesign);
            relativeDerivativeNorm = ...
                designColumnNorms(2) ./ designColumnNorms(1);

            if designRank < 2
                orthogonalDerivative(:) = 0;
                complexDesign(:, 2) = 0; %#ok
                realDesign(:, 2) = 0;
            end

            [scaledDesign, columnScales] = ...
                local_scale_design_columns(realDesign);
            scaledCoefficients = ...
                local_real_least_squares(scaledDesign, realOutcome);
            deviationCoefficients = zeros( ...
                size(scaledCoefficients), 'like', scaledCoefficients);

            for coefficientIndex = 1:size(scaledCoefficients, 1)
                if columnScales(coefficientIndex) > 0
                    deviationCoefficients(coefficientIndex, :) = ...
                        scaledCoefficients(coefficientIndex, :) ./ ...
                        columnScales(coefficientIndex);
                end
            end

            deviationCoefficients = ...
                local_enforce_zero_row_means(deviationCoefficients);

            normalizedDerivativeRow = deviationCoefficients(2, :);
            amplitudeRow = 1 + deviationCoefficients(1, :) - ...
                derivativeProjection .* normalizedDerivativeRow;

            derivativeRow = normalizedDerivativeRow .* timeSpan;

            modeledScaledMatrix = ...
                scaledMeanWaveform * ...
                (1 + deviationCoefficients(1, :)) + ...
                orthogonalDerivative * normalizedDerivativeRow;
            modeledMatrix = modeledScaledMatrix .* sliceScale;

            latencyRow = nan(1, nTrials, 'like', real(X));

            if designRank == 2
                validAmplitude = abs(amplitudeRow) > sqrt(epsilonValue);
                latencyRow(validAmplitude) = ...
                    -derivativeRow(validAmplitude) ./ ...
                    amplitudeRow(validAmplitude);
            end

        end

        amplitude(channelIndex, frequencyIndex, :) = ...
            reshape(amplitudeRow, [1, 1, nTrials]);
        derivative(channelIndex, frequencyIndex, :) = ...
            reshape(derivativeRow, [1, 1, nTrials]);
        latency(channelIndex, frequencyIndex, :) = ...
            reshape(latencyRow, [1, 1, nTrials]);
        erp(channelIndex, frequencyIndex, :, :) = ...
            reshape(modeledMatrix, [1, 1, nTimes, nTrials]);

        diagnostics.rank(channelIndex, frequencyIndex) = uint8(designRank);
        diagnostics.conditionNumber(channelIndex, frequencyIndex) = ...
            designConditionNumber;
        diagnostics.orthogonalization(channelIndex, frequencyIndex) = ...
            derivativeProjection;
        diagnostics.relativeNorm(channelIndex, frequencyIndex) = ...
            relativeDerivativeNorm;
        diagnostics.columnNorms(channelIndex, frequencyIndex, :) = ...
            reshape(designColumnNorms, [1, 1, 2]);

    end
end

end


function [scaledDesign, columnScales] = ...
    local_scale_design_columns(design)

columnScales = sqrt(sum(design .^ 2, 1));
scaledDesign = zeros(size(design), 'like', design);

for columnIndex = 1:size(design, 2)
    if columnScales(columnIndex) > 0
        scaledDesign(:, columnIndex) = ...
            design(:, columnIndex) ./ columnScales(columnIndex);
    end
end

end


function [designRank, conditionNumber, columnNorms] = ...
    local_design_diagnostics(design)

singularValues = svd(design, 'econ');
columnNorms = sqrt(sum(design .^ 2, 1));

if isempty(singularValues) == true
    designRank = 0;
    conditionNumber = Inf;
    return
end

maximumSingularValue = max(singularValues);

if maximumSingularValue == 0
    designRank = 0;
    conditionNumber = Inf;
    return
end

epsilonValue = local_epsilon(design);
relativeTolerance = max( ...
    sqrt(epsilonValue), max(size(design)) .* epsilonValue);
tolerance = relativeTolerance .* maximumSingularValue;
designRank = sum(singularValues > tolerance);

if designRank < size(design, 2)
    conditionNumber = Inf;
else
    conditionNumber = double( ...
        maximumSingularValue ./ min(singularValues));
end

end


function coefficients = local_real_least_squares(design, outcome)

singularValues = svd(design, 'econ');

if isempty(singularValues) == true
    coefficients = zeros( ...
        size(design, 2), size(outcome, 2), 'like', outcome);
    return
end

maximumSingularValue = max(singularValues);

if maximumSingularValue == 0
    coefficients = zeros( ...
        size(design, 2), size(outcome, 2), 'like', outcome);
    return
end


epsilonValue = local_epsilon(design);
relativeTolerance = max( ...
    sqrt(epsilonValue), max(size(design)) .* epsilonValue);
tolerance = relativeTolerance .* maximumSingularValue;

coefficients = pinv(design, tolerance) * outcome;

end


function centered = local_enforce_zero_row_means(values)

centered = values - mean(values, 2);
nColumns = size(centered, 2);

if nColumns == 1
    centered(:, 1) = 0;
else
    centered(:, nColumns) = -sum(centered(:, 1:(nColumns - 1)), 2);
end

end


function derivative = local_time_derivative(values, times)

[nChannels, nFrequencies, nTimes] = size(values);
derivative = zeros(size(values), 'like', values);

if nTimes == 1
    return
end

timeVector = double(times(:));

for channelIndex = 1:nChannels
    for frequencyIndex = 1:nFrequencies

        waveform = reshape( ...
            values(channelIndex, frequencyIndex, :), [nTimes, 1]);

        if nTimes == 2

            slope = (waveform(2) - waveform(1)) ./ ...
                (timeVector(2) - timeVector(1));
            derivativeWaveform = repmat(slope, [2, 1]);

        else

            derivativeWaveform = zeros(nTimes, 1, 'like', waveform);

            for timeIndex = 2:(nTimes - 1)
                leftWidth = timeVector(timeIndex) - timeVector(timeIndex - 1);
                rightWidth = timeVector(timeIndex + 1) - timeVector(timeIndex);

                derivativeWaveform(timeIndex) = ...
                    (-rightWidth ./ (leftWidth .* ...
                    (leftWidth + rightWidth))) .* waveform(timeIndex - 1) + ...
                    ((rightWidth - leftWidth) ./ ...
                    (leftWidth .* rightWidth)) .* waveform(timeIndex) + ...
                    (leftWidth ./ (rightWidth .* ...
                    (leftWidth + rightWidth))) .* waveform(timeIndex + 1);
            end

            firstWidth = timeVector(2) - timeVector(1);
            secondWidth = timeVector(3) - timeVector(2);
            derivativeWaveform(1) = ...
                (-(2 .* firstWidth + secondWidth) ./ ...
                (firstWidth .* (firstWidth + secondWidth))) .* waveform(1) + ...
                ((firstWidth + secondWidth) ./ ...
                (firstWidth .* secondWidth)) .* waveform(2) - ...
                (firstWidth ./ ...
                (secondWidth .* (firstWidth + secondWidth))) .* waveform(3);

            previousWidth = timeVector(nTimes - 1) - ...
                timeVector(nTimes - 2);
            finalWidth = timeVector(nTimes) - timeVector(nTimes - 1);
            derivativeWaveform(nTimes) = ...
                (finalWidth ./ ...
                (previousWidth .* (previousWidth + finalWidth))) .* ...
                waveform(nTimes - 2) - ...
                ((previousWidth + finalWidth) ./ ...
                (previousWidth .* finalWidth)) .* waveform(nTimes - 1) + ...
                ((previousWidth + 2 .* finalWidth) ./ ...
                (finalWidth .* (previousWidth + finalWidth))) .* ...
                waveform(nTimes);

        end

        derivative(channelIndex, frequencyIndex, :) = ...
            reshape(derivativeWaveform, [1, 1, nTimes]);

    end
end

end


function [powerOutput, phaseOutput, zeroMask] = ...
    local_encode_signal(signal, scale, powerPrototype, phasePrototype)

signalMagnitude = abs(signal);
zeroMask = real(signal) == 0 & imag(signal) == 0;

if strcmp(scale, 'log10') == true
    powerOutput = zeros(size(signalMagnitude), 'like', powerPrototype);
    nonzeroMask = ~zeroMask;
    powerOutput(nonzeroMask) = ...
        2 .* log10(signalMagnitude(nonzeroMask));
    powerOutput(zeroMask) = -Inf;
else
    powerOutput = signalMagnitude .^ 2;

    underflowMask = ~zeroMask & powerOutput == 0;

    if any(underflowMask(:)) == true
        error('nf_rmerp:OutputPowerUnderflow', ...
            ['A nonzero modeled coefficient underflowed while encoding ' ...
            'linear power. Rescale the analytic signal or use log10 power.']);
    end

    if any(~isfinite(powerOutput(:))) == true
        error('nf_rmerp:OutputPowerOverflow', ...
            ['A finite modeled coefficient produced power outside the ' ...
            'representable numeric range. Rescale the analytic signal.']);
    end
end

phaseOutput = angle(signal);
phaseOutput = cast(phaseOutput, 'like', phasePrototype);
phaseOutput(zeroMask) = NaN;
powerOutput = cast(powerOutput, 'like', powerPrototype);

end


function squeezedSize = local_squeezed_size(canonicalSize)

nonSingleton = canonicalSize(canonicalSize ~= 1);

if isempty(nonSingleton) == true
    squeezedSize = [1, 1];
elseif isscalar(nonSingleton) == true
    squeezedSize = [nonSingleton, 1];
else
    squeezedSize = nonSingleton;
end

end


function textValue = local_size_text(sizeValue)

sizeValue = sizeValue(:).';
pieces = cell(1, numel(sizeValue));

for valueIndex = 1:numel(sizeValue)
    pieces{valueIndex} = sprintf('%d', sizeValue(valueIndex));
end

textValue = strjoin(pieces, ' x ');

end


function value = local_max_abs(array)

if isempty(array) == true
    value = 0;
else
    value = double(max(abs(array(:))));
end

end


function meanValue = local_stable_mean(array, dimension)

arrayScale = max(abs(array(:)));

if arrayScale == 0
    meanValue = mean(array, dimension);
else
    meanValue = mean(array ./ arrayScale, dimension) .* arrayScale;
end

if any(~isfinite(meanValue(:))) == true
    error('nf_rmerp:NonfiniteConditionMean', ...
        ['A condition mean overflowed or became nonfinite despite finite ' ...
        'input coefficients. Rescale the analytic coefficients.']);
end

end


function value = local_epsilon(prototype)

if isa(prototype, 'single') == true
    value = eps('single');
else
    value = eps('double');
end

end


function value = local_realmin(prototype)

if isa(prototype, 'single') == true
    value = realmin('single');
else
    value = realmin('double');
end

end


function value = local_smallest_positive(prototype)

value = local_realmin(prototype) .* local_epsilon(prototype);

end


function TF = nf_powerfront_local(TF)

fieldNames = fieldnames(TF);
powerFields = {};
otherFields = {};

for fieldIndex = 1:numel(fieldNames)

    fieldName = fieldNames{fieldIndex};
    fieldValue = TF.(fieldName);

    if isnumeric(fieldValue) == true
        if contains(lower(fieldName), 'power') == true
            powerFields{end + 1, 1} = fieldName; %#ok
        else
            otherFields{end + 1, 1} = fieldName; %#ok
        end
    else
        otherFields{end + 1, 1} = fieldName; %#ok
    end

end

primaryPowerIndex = find(strcmpi(powerFields, 'power'));

if isempty(primaryPowerIndex) == false
    powerFields(primaryPowerIndex) = [];
    powerFields = [{'power'}; powerFields];
end

fieldOrder = [powerFields; otherFields];
TF = orderfields(TF, fieldOrder);

end
