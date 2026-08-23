function model = nf_cpm(coefficients, varargin)
% NF_CPM  Estimate concurrent-phaser ensemble parameters.
%
% GENERAL
% -------
% NF_CPM is the numerical concurrent phaser method (CPM) engine used by
% NF_NORMBASE. It estimates ensemble statistics of phase-locked (PL) and
% non-phase-locked (NPL) amplitudes from complex single-trial coefficients.
% It does not average a NeuroFreq TF structure and does not apply baseline
% correction. Use NF_NORMBASE with AverageMethod='cpm' for that workflow.
%
% CPM assumes, at each coefficient location,
%
%   R_r = A_r*sin(alpha) + B_r*sin(phi_r)
%   I_r = A_r*cos(alpha) + B_r*cos(phi_r)
%
% where A_r and B_r are nonnegative, alpha is constant across trials, and
% phi_r is uniform and independent of the amplitudes. NF_CPM evaluates the
% published finite-sample moment equations on their estimable, nonsingular
% domain, expressed in the native complex-coefficient coordinates C=X+iY of
% NeuroFreq, where X=R, Y=-I, and the returned PL phase angle(mean(C))
% equals alpha-pi/2. This removes the Fourier sine/cosine convention from
% the public output without changing the estimator. In these coordinates,
%
%   E[A^2] = (E[X^2]-E[Y^2])/cos(2*angle(E[C]))
%   E[B^2] = E[abs(C)^2]-E[A^2].
%
% NF_CPM does not substitute a phase-rotated moment estimator at the
% intrinsic cos(2*alpha)=0 singularity.
%
% USAGE
% -----
%   model = nf_cpm(coefficients)
%   model = nf_cpm(coefficients, ...
%       'TrialDimension', 4, ...
%       'Verbosity', 'silent')
%
% INPUT
% -----
% coefficients
%       Floating-point complex coefficients. NaN marks a missing trial at
%       that coefficient location. Infinite values are rejected.
%
% OPTIONS
% -------
% TrialDimension
%       Dimension containing trials. The default is the last nonsingleton
%       dimension. At least two trials are required.
%
% Verbosity
%       'auto' (default), 'verbose', or 'silent'. Auto reports scientific
%       validity warnings. Verbose also reports dimensions and completion.
%       Silent suppresses function-owned output; errors always throw.
%
% OUTPUT
% ------
% model.phaseLocked and model.nonPhaseLocked contain amplitudeMean,
% amplitudeVariance, amplitudeSecondMoment, and powerMean in NeuroFreq
% coefficient units. model.raw retains signed method-of-moments estimates.
% model.validity identifies undefined, singular, or nonphysical estimates;
% negative estimates are never clipped.
%
% REFERENCE
% ---------
% Singhal S, Ghosh P, Kumar N, Banerjee A. J Neurophysiol.
% 2023;129(1):199-210. doi:10.1152/jn.00467.2022

if nargin < 1
    error('nf_cpm:CoefficientsRequired', ...
        'A complex coefficient array is required.');
end

options = local_parse_inputs(coefficients, varargin{:});
trialDimension = options.TrialDimension;
nTrials = size(coefficients, trialDimension);

if nTrials < 2
    error('nf_cpm:MultipleTrialsRequired', ...
        'CPM requires at least two trials.');
end

local_report( ...
    options.Verbosity, ...
    'verbose', ...
    'Estimating CPM across %d trials in dimension %d.', ...
    nTrials, ...
    trialDimension);

nDimensions = max(ndims(coefficients), trialDimension);
dimensionOrder = [ ...
    1:(trialDimension - 1), ...
    (trialDimension + 1):nDimensions, ...
    trialDimension];
canonical = permute(coefficients, dimensionOrder);
trialDimensionCanonical = ndims(canonical);
canonicalSize = size(canonical);
outputSize = canonicalSize(1:(end - 1));

if isempty(outputSize) == true
    outputSize = [1, 1];
elseif isscalar(outputSize) == true
    outputSize = [outputSize, 1];
end

valid = isfinite(real(canonical)) & isfinite(imag(canonical));
validCount = sum(valid, trialDimensionCanonical);
validCountLike = cast(validCount, 'like', real(canonical));
sufficientTrials = validCount >= 2;
working = canonical;
working(valid == false) = 0;

meanCoefficient = local_stable_masked_mean( ...
    working, ...
    validCountLike, ...
    trialDimensionCanonical);
phaseLockedPhase = angle(meanCoefficient);
phaseLockedAmplitudeMean = abs(meanCoefficient);
phaseDefined = phaseLockedAmplitudeMean > 0 & sufficientTrials;
phaseLockedPhase(phaseDefined == false) = NaN;
phaseLockedAmplitudeMean(phaseDefined == false) = NaN;

phaseForRotation = phaseLockedPhase;
phaseForRotation(phaseDefined == false) = 0;
rotated = bsxfun(@times, working, exp(-1i .* phaseForRotation));
positiveProjection = valid & imag(rotated) > 0;
positiveProjectionCount = sum( ...
    positiveProjection, ...
    trialDimensionCanonical);
positiveProjectionCountLike = cast( ...
    positiveProjectionCount, ...
    'like', ...
    real(canonical));
positiveProjectionValues = imag(rotated);
positiveProjectionValues(positiveProjection == false) = 0;
positiveProjectionMean = local_stable_masked_mean( ...
    positiveProjectionValues, ...
    positiveProjectionCountLike, ...
    trialDimensionCanonical);
nonPhaseLockedAmplitudeMean = ...
    cast(pi ./ 2, 'like', real(canonical)) .* ...
    positiveProjectionMean;
nonPhaseLockedAmplitudeMean(positiveProjectionCount == 0) = NaN;
nonPhaseLockedAmplitudeMean(phaseDefined == false) = NaN;

secondMomentDifference = local_stable_quadratic_difference( ...
    working, ...
    validCountLike, ...
    trialDimensionCanonical);
totalPowerMean = local_stable_square_mean( ...
    abs(working), ...
    validCountLike, ...
    trialDimensionCanonical);
conditioningDenominator = cos(2 .* phaseForRotation);
singularityTolerance = 64 .* eps(cast(1, 'like', real(canonical)));
nonsingular = phaseDefined & ...
    abs(conditioningDenominator) > singularityTolerance;
estimable = sufficientTrials & phaseDefined & nonsingular;

rawPhaseLockedSecondMoment = ...
    secondMomentDifference ./ ...
    conditioningDenominator;
rawNonPhaseLockedSecondMoment = ...
    totalPowerMean - rawPhaseLockedSecondMoment;
rawPhaseLockedSecondMoment(estimable == false) = NaN;
rawNonPhaseLockedSecondMoment(estimable == false) = NaN;

rawPhaseLockedVariance = rawPhaseLockedSecondMoment - ...
    phaseLockedAmplitudeMean .^ 2;
rawNonPhaseLockedVariance = rawNonPhaseLockedSecondMoment - ...
    nonPhaseLockedAmplitudeMean .^ 2;

phaseLockedSecondMomentValid = estimable & ...
    isfinite(rawPhaseLockedSecondMoment) & ...
    rawPhaseLockedSecondMoment >= 0;
nonPhaseLockedSecondMomentValid = estimable & ...
    isfinite(rawNonPhaseLockedSecondMoment) & ...
    rawNonPhaseLockedSecondMoment >= 0;
phaseLockedVarianceValid = phaseLockedSecondMomentValid & ...
    isfinite(rawPhaseLockedVariance) & ...
    rawPhaseLockedVariance >= 0;
nonPhaseLockedVarianceValid = nonPhaseLockedSecondMomentValid & ...
    positiveProjectionCount > 0 & ...
    isfinite(rawNonPhaseLockedVariance) & ...
    rawNonPhaseLockedVariance >= 0;

phaseLockedSecondMoment = rawPhaseLockedSecondMoment;
phaseLockedSecondMoment(phaseLockedSecondMomentValid == false) = NaN;
nonPhaseLockedSecondMoment = rawNonPhaseLockedSecondMoment;
nonPhaseLockedSecondMoment( ...
    nonPhaseLockedSecondMomentValid == false) = NaN;
phaseLockedAmplitudeVariance = rawPhaseLockedVariance;
phaseLockedAmplitudeVariance(phaseLockedVarianceValid == false) = NaN;
nonPhaseLockedAmplitudeVariance = rawNonPhaseLockedVariance;
nonPhaseLockedAmplitudeVariance( ...
    nonPhaseLockedVarianceValid == false) = NaN;
totalPowerMean(sufficientTrials == false) = NaN;

leakageRatio = rawPhaseLockedVariance ./ ...
    rawNonPhaseLockedSecondMoment;
leakageValid = phaseLockedVarianceValid & ...
    nonPhaseLockedSecondMomentValid & ...
    rawNonPhaseLockedSecondMoment > 0;
leakageRatio(leakageValid == false) = NaN;

closureError = rawPhaseLockedSecondMoment + ...
    rawNonPhaseLockedSecondMoment - totalPowerMean;
inverseAbsoluteDenominator = 1 ./ abs(conditioningDenominator);
inverseAbsoluteDenominator(nonsingular == false) = Inf;

model = struct();
model.schema = 'nf-cpm/1.0';
model.citation = [ ...
    'Singhal S, Ghosh P, Kumar N, Banerjee A. ' ...
    'J Neurophysiol. 2023;129(1):199-210.'];
model.doi = '10.1152/jn.00467.2022';
model.estimator = struct();
model.estimator.model = 'concurrent_phaser_method';
model.estimator.secondMomentEquations = ...
    'Singhal_et_al_2023_equations_23_24';
model.estimator.coefficientConvention = ...
    'power=abs(C)^2, phase=angle(C)';
model.estimator.phaseConvention = ...
    'phase_locked_phase=angle(mean(C))';
model.estimator.domain = 'estimable_nonsingular_domain';
model.estimator.singularityTolerance = double(singularityTolerance);
model.estimator.conditioningDiagnostic = ...
    'inverseAbsoluteDenominator';
model.estimator.numericalEvaluation = ...
    'max_scaled_algebraically_equivalent_moments';
model.trialDimension = trialDimension;
model.numberOfTrials = nTrials;
model.outputSize = outputSize;

model.total = struct();
model.total.powerMean = local_reshape_output(totalPowerMean, outputSize);

model.phaseLocked = struct();
model.phaseLocked.phase = local_reshape_output( ...
    phaseLockedPhase, outputSize);
model.phaseLocked.amplitudeMean = local_reshape_output( ...
    phaseLockedAmplitudeMean, outputSize);
model.phaseLocked.amplitudeVariance = local_reshape_output( ...
    phaseLockedAmplitudeVariance, outputSize);
model.phaseLocked.amplitudeSecondMoment = local_reshape_output( ...
    phaseLockedSecondMoment, outputSize);
model.phaseLocked.powerMean = local_reshape_output( ...
    phaseLockedSecondMoment, outputSize);

model.nonPhaseLocked = struct();
model.nonPhaseLocked.amplitudeMean = local_reshape_output( ...
    nonPhaseLockedAmplitudeMean, outputSize);
model.nonPhaseLocked.amplitudeVariance = local_reshape_output( ...
    nonPhaseLockedAmplitudeVariance, outputSize);
model.nonPhaseLocked.amplitudeSecondMoment = local_reshape_output( ...
    nonPhaseLockedSecondMoment, outputSize);
model.nonPhaseLocked.powerMean = local_reshape_output( ...
    nonPhaseLockedSecondMoment, outputSize);

model.leakageRatio = local_reshape_output(leakageRatio, outputSize);

model.raw = struct();
model.raw.phaseLockedAmplitudeSecondMoment = local_reshape_output( ...
    rawPhaseLockedSecondMoment, outputSize);
model.raw.nonPhaseLockedAmplitudeSecondMoment = local_reshape_output( ...
    rawNonPhaseLockedSecondMoment, outputSize);
model.raw.phaseLockedAmplitudeVariance = local_reshape_output( ...
    rawPhaseLockedVariance, outputSize);
model.raw.nonPhaseLockedAmplitudeVariance = local_reshape_output( ...
    rawNonPhaseLockedVariance, outputSize);
model.raw.closureError = local_reshape_output( ...
    closureError, outputSize);

model.validity = struct();
model.validity.validTrialCount = local_reshape_output( ...
    validCount, outputSize);
model.validity.sufficientTrials = local_reshape_output( ...
    sufficientTrials, outputSize);
model.validity.phaseDefined = local_reshape_output( ...
    phaseDefined, outputSize);
model.validity.positiveProjectionCount = local_reshape_output( ...
    positiveProjectionCount, outputSize);
model.validity.conditioningDenominator = local_reshape_output( ...
    conditioningDenominator, outputSize);
model.validity.inverseAbsoluteDenominator = local_reshape_output( ...
    inverseAbsoluteDenominator, outputSize);
model.validity.nonsingular = local_reshape_output( ...
    nonsingular, outputSize);
model.validity.phaseLockedSecondMoment = local_reshape_output( ...
    phaseLockedSecondMomentValid, outputSize);
model.validity.nonPhaseLockedSecondMoment = local_reshape_output( ...
    nonPhaseLockedSecondMomentValid, outputSize);
model.validity.phaseLockedAmplitudeVariance = local_reshape_output( ...
    phaseLockedVarianceValid, outputSize);
model.validity.nonPhaseLockedAmplitudeVariance = local_reshape_output( ...
    nonPhaseLockedVarianceValid, outputSize);

summary = struct();
summary.insufficientTrialLocations = sum(sufficientTrials(:) == false);
summary.incompleteTrialLocations = sum( ...
    sufficientTrials(:) & validCount(:) < nTrials);
summary.undefinedPhaseLocations = sum( ...
    sufficientTrials(:) & phaseDefined(:) == false);
summary.singularLocations = sum( ...
    sufficientTrials(:) & ...
    phaseDefined(:) & ...
    nonsingular(:) == false);
summary.undefinedNonPhaseLockedMeanLocations = ...
    sum(positiveProjectionCount(:) == 0); %#ok
summary.invalidPhaseLockedSecondMoments = ...
    sum(estimable(:) & phaseLockedSecondMomentValid(:) == false);
summary.invalidNonPhaseLockedSecondMoments = ...
    sum(estimable(:) & nonPhaseLockedSecondMomentValid(:) == false);
summary.invalidPhaseLockedAmplitudeVariances = ...
    sum(phaseLockedSecondMomentValid(:) & ...
    phaseLockedVarianceValid(:) == false);
summary.invalidNonPhaseLockedAmplitudeVariances = ...
    sum(nonPhaseLockedSecondMomentValid(:) & ...
    nonPhaseLockedVarianceValid(:) == false);
model.validity.summary = summary;

invalidMomentCount = ...
    summary.invalidPhaseLockedSecondMoments + ...
    summary.invalidNonPhaseLockedSecondMoments + ...
    summary.invalidPhaseLockedAmplitudeVariances + ...
    summary.invalidNonPhaseLockedAmplitudeVariances;

if summary.insufficientTrialLocations > 0
    local_issue_warning( ...
        options.Verbosity, ...
        'nf_cpm:InsufficientValidTrials', ...
        ['%d coefficient location(s) had fewer than two finite trials; ' ...
        'their CPM estimates are NaN.'], ...
        summary.insufficientTrialLocations);
end

if summary.incompleteTrialLocations > 0
    local_issue_warning( ...
        options.Verbosity, ...
        'nf_cpm:IncompleteTrialEnsemble', ...
        ['%d coefficient location(s) omitted at least one nonfinite ' ...
        'trial coefficient; CPM used the remaining finite trials.'], ...
        summary.incompleteTrialLocations);
end

if summary.undefinedPhaseLocations > 0
    local_issue_warning( ...
        options.Verbosity, ...
        'nf_cpm:UndefinedPhase', ...
        ['%d coefficient location(s) had no estimable phase-locked ' ...
        'phase; their CPM component estimates are NaN.'], ...
        summary.undefinedPhaseLocations);
end

if summary.singularLocations > 0
    local_issue_warning( ...
        options.Verbosity, ...
        'nf_cpm:SingularMomentEquation', ...
        ['%d coefficient location(s) were at the intrinsic CPM ' ...
        'cos(2*alpha)=0 singularity; their component estimates are NaN.'], ...
        summary.singularLocations);
end

if invalidMomentCount > 0
    local_issue_warning( ...
        options.Verbosity, ...
        'nf_cpm:InvalidMomentEstimate', ...
        ['%d negative finite-sample moment estimate(s) ' ...
        'were retained in model.raw and represented as NaN in physical ' ...
        'model fields.'], ...
        invalidMomentCount);
end

local_report( ...
    options.Verbosity, ...
    'verbose', ...
    'Completed CPM estimation.');

end


function options = local_parse_inputs(coefficients, varargin)

parser = inputParser;
parser.FunctionName = mfilename;
parser.CaseSensitive = false;
parser.PartialMatching = false;
parser.KeepUnmatched = false;

defaultTrialDimension = find(size(coefficients) ~= 1, 1, 'last');

if isempty(defaultTrialDimension) == true
    defaultTrialDimension = ndims(coefficients);
end

addRequired(parser, 'coefficients', @local_valid_coefficients);
addParameter(parser, 'TrialDimension', defaultTrialDimension, ...
    @local_valid_dimension);
addParameter(parser, 'Verbosity', 'auto', @local_is_text_scalar);
parse(parser, coefficients, varargin{:});

options = struct();
options.TrialDimension = double(parser.Results.TrialDimension);
options.Verbosity = lower(local_text_scalar( ...
    parser.Results.Verbosity, ...
    'Verbosity'));

if options.TrialDimension > max(ndims(coefficients), 2)
    error('nf_cpm:InvalidTrialDimension', ...
        'TrialDimension exceeds the coefficient array dimensions.');
end

if any(strcmp(options.Verbosity, ...
        {'auto', 'verbose', 'silent'})) == false
    error('nf_cpm:InvalidVerbosity', ...
        'Verbosity must be ''auto'', ''verbose'', or ''silent''.');
end

end


function valid = local_valid_coefficients(value)

valid = isnumeric(value) && isfloat(value) && isempty(value) == false;

if valid == true
    valid = any(isinf(real(value(:)))) == false && ...
        any(isinf(imag(value(:)))) == false;
end

end


function valid = local_valid_dimension(value)

valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value >= 1 && value == fix(value);

end


function valid = local_is_text_scalar(value)

valid = ischar(value) || (isstring(value) && isscalar(value));

end


function value = local_text_scalar(value, name)

if local_is_text_scalar(value) == false
    error('nf_cpm:TextScalarRequired', ...
        '%s must be a character vector or string scalar.', name);
end

value = strtrim(char(value));

if isempty(value) == true
    error('nf_cpm:EmptyTextOption', ...
        '%s cannot be empty.', name);
end

end


function meanValue = local_stable_masked_mean( ...
    values, validCount, dimension)

scale = max(abs(values), [], dimension);
divisor = scale;
divisor(scale == 0) = 1;
scaled = bsxfun(@rdivide, values, divisor);
scaledMean = sum(scaled, dimension) ./ validCount;
meanValue = scale .* scaledMean;
meanValue(validCount == 0) = NaN;

end


function meanSquare = local_stable_square_mean( ...
    values, validCount, dimension)

scale = max(values, [], dimension);
divisor = scale;
divisor(scale == 0) = 1;
scaled = bsxfun(@rdivide, values, divisor);
scaledMeanSquare = sum(scaled .^ 2, dimension) ./ validCount;
meanSquare = scale .* (scale .* scaledMeanSquare);
meanSquare(validCount == 0) = NaN;

end


function meanDifference = local_stable_quadratic_difference( ...
    values, validCount, dimension)

scale = max(abs(values), [], dimension);
divisor = scale;
divisor(scale == 0) = 1;
scaled = bsxfun(@rdivide, values, divisor);
scaledDifference = real(scaled) .^ 2 - imag(scaled) .^ 2;
scaledMeanDifference = sum(scaledDifference, dimension) ./ validCount;
meanDifference = scale .* (scale .* scaledMeanDifference);
meanDifference(validCount == 0) = NaN;

end


function output = local_reshape_output(value, outputSize)

output = reshape(value, outputSize);

end


function local_report(verbosity, level, formatSpec, varargin)

if strcmp(verbosity, 'silent') == true
    return
end

if strcmp(level, 'verbose') == true && ...
        strcmp(verbosity, 'verbose') == false
    return
end

fprintf('[nf_cpm] ');
fprintf([formatSpec, '\n'], varargin{:});

end


function local_issue_warning(verbosity, identifier, message, varargin)

if strcmp(verbosity, 'silent') == true
    return
end

warning(identifier, message, varargin{:});

end
