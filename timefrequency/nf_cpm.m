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
% NeuroFreq modification, September 2026: add the optional projection
% estimator from Supplement Eq. S1.3; retain the original estimator by default.
%
% CPM assumes, at each coefficient location,
%
%   R_r = A_r*sin(alpha) + B_r*sin(phi_r)
%   I_r = A_r*cos(alpha) + B_r*cos(phi_r)
%
% where A_r and B_r are nonnegative, alpha is constant across trials, and
% phi_r is uniform and independent of the amplitudes. NF_CPM estimates the
% published population moments using sample averages, expressed in the
% native complex-coefficient coordinates C=X+iY of
% NeuroFreq, where X=R, Y=-I, and the returned PL phase angle(mean(C))
% equals alpha-pi/2. This removes the Fourier sine/cosine convention from
% the public output without changing the estimator. In these coordinates,
%
%   E[A^2] = (E[X^2]-E[Y^2])/cos(2*angle(E[C]))
%   E[B^2] = E[abs(C)^2]-E[A^2].
%
% The optional projection estimator instead uses Z=C*exp(-i*angle(E[C])):
%
%   E[B^2] = 2*E[imag(Z)^2]
%   E[A^2] = E[abs(C)^2]-E[B^2].
%
% This avoids the original formula's cos(2*alpha) denominator. Both target
% the same population moments but differ for finite samples. Selection is
% explicit: there is no automatic fallback between estimators.
%
% USAGE
% -----
%   model = nf_cpm(coefficients)
%   model = nf_cpm(coefficients, ...
%       'TrialDimension', 4, ...
%       'Estimator', 'projection', ...
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
% Estimator
%       'original' (default): main-text Eqs. 23-24.
%       'projection': Supplement Eq. S1.3 with j=2, then total minus NPL.
%       Only second-moment estimation changes; amplitude means, variance
%       calculations, and physical-output masks retain their existing rules.
%       Projection NPL power is nonnegative; PL power and amplitude variances
%       can still have negative finite-sample estimates.
%
% ExecutionEnvironment
%       'auto' (default) follows the input. 'gpu' also accepts ordinary
%       arrays, sends each block to the GPU, and returns ordinary maps.
%       If no GPU is available, host input uses the CPU. A gpuArray input
%       retains GPU maps unless 'cpu' is selected. Summary counts and
%       estimator metadata are always available on the host.
%
% MaxBlockLocations
%       Optional positive integer limiting independent coefficient locations
%       processed together. A block retains every trial for each location.
%       GPU execution always uses a memory-based cap. On the CPU the default
%       evaluates the full array, as in previous versions.
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
% model.estimator.method records the selection. For 'projection', the
% conditioning denominator is 1 (no division by a phase-dependent quantity);
% undefined phase still makes component estimates unavailable.
%
% REFERENCE
% ---------
% Singhal S, Ghosh P, Kumar N, Banerjee A. J Neurophysiol.
% 2023;129(1):199-210. doi:10.1152/jn.00467.2022
% Supplement Eq. S1.3: doi:10.5281/zenodo.7315006

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

resident = isa(coefficients, 'gpuArray');
[useGPU, device, unavailable] = hnf_gpu_execution(coefficients,options.ExecutionEnvironment);
if ~isempty(unavailable)
    local_issue_warning(options.Verbosity, 'nf_cpm:GPUUnavailable', ...
        'GPU unavailable; evaluating CPM on the CPU. %s', unavailable);
end
if resident && ~useGPU
    coefficients = gather(coefficients);
    resident = false;
end
precision = local_precision(coefficients);

inputSize = size(coefficients);
outputSize = inputSize;
outputSize(trialDimension) = [];
if isscalar(outputSize), outputSize = [outputSize, 1]; end

if useGPU || ~isempty(options.MaxBlockLocations)
    model = local_tiled(coefficients, options, outputSize, nTrials, ...
        precision, device, useGPU, resident);
else
    if any(isinf(coefficients(:)))
        error('nf_cpm:InfiniteCoefficients','Infinite coefficients are not supported.');
    end
    % Place the full trial axis last for the original vectorized CPU kernel.
    order = [1:trialDimension-1,trialDimension+1:ndims(coefficients),trialDimension];
    model = local_kernel(permute(coefficients,order), options, ...
        outputSize, nTrials, precision);
end
if useGPU, model.executionEnvironment = 'gpu';
else, model.executionEnvironment = 'cpu'; end
model.validity.summary = local_summary(model, nTrials, resident);
local_warnings(model.validity.summary, options);
local_report(options.Verbosity,'verbose','Completed CPM estimation.');
end


function model = local_kernel(canonical, options, outputSize, nTrials, precision)
% Shared estimator for a whole array or an independent block of locations.
% The trials at a location stay together, including missing coefficients.
trialDimensionCanonical = ndims(canonical);
valid = isfinite(real(canonical)) & isfinite(imag(canonical));
validCount = sum(valid, trialDimensionCanonical);
validCountLike = cast(validCount, 'like', real(canonical(1)));
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
    real(canonical(1)));
positiveProjectionValues = imag(rotated);
positiveProjectionValues(positiveProjection == false) = 0;
positiveProjectionMean = local_stable_masked_mean( ...
    positiveProjectionValues, ...
    positiveProjectionCountLike, ...
    trialDimensionCanonical);
nonPhaseLockedAmplitudeMean = ...
    cast(pi ./ 2, 'like', real(canonical(1))) .* ...
    positiveProjectionMean;
nonPhaseLockedAmplitudeMean(positiveProjectionCount == 0) = NaN;
nonPhaseLockedAmplitudeMean(phaseDefined == false) = NaN;

totalPowerMean = local_stable_square_mean( ...
    abs(working), ...
    validCountLike, ...
    trialDimensionCanonical);

if strcmp(options.Estimator, 'original')
    secondMomentDifference = local_stable_quadratic_difference( ...
        working, ...
        validCountLike, ...
        trialDimensionCanonical);
    conditioningDenominator = cos(2 .* phaseForRotation);
    singularityTolerance = 64 .* eps(cast(1, precision));
    nonsingular = phaseDefined & ...
        abs(conditioningDenominator) > singularityTolerance;
    rawPhaseLockedSecondMoment = ...
        secondMomentDifference ./ conditioningDenominator;
    rawNonPhaseLockedSecondMoment = ...
        totalPowerMean - rawPhaseLockedSecondMoment;
else
    % Supplement Eq. S1.3, j=2: use ALL finite trials, not just the positive
    % projections used above for the NPL amplitude mean. Reuse the existing
    % rotation and scaled square-mean helper to avoid the cos(2*phase) pole.
    rawNonPhaseLockedSecondMoment = 2 .* local_stable_square_mean( ...
        abs(imag(rotated)), validCountLike, trialDimensionCanonical);
    rawPhaseLockedSecondMoment = ...
        totalPowerMean - rawNonPhaseLockedSecondMoment;
    conditioningDenominator = ones(size(phaseForRotation), ...
        'like', phaseForRotation);
    singularityTolerance = 0;
    nonsingular = phaseDefined;
end

estimable = sufficientTrials & phaseDefined & nonsingular;
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
model.estimator.method = options.Estimator;
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
if strcmp(options.Estimator, 'projection')
    model.estimator.secondMomentEquations = ...
        'Singhal_et_al_2023_supplement_S1_3_j2_and_total_minus_NPL';
    model.estimator.domain = 'sufficient_trials_and_defined_phase';
    model.estimator.conditioningDiagnostic = 'not_applicable';
end
model.estimator.numericalEvaluation = ...
    'max_scaled_algebraically_equivalent_moments';
model.trialDimension = options.TrialDimension;
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


end


function model = local_tiled(coefficients, options, outputSize, nTrials, ...
    precision, device, useGPU, resident)
% The timing-critical array reductions take place on the device. A host
% caller receives results in three bulk copies per block (floats, counts,
% masks); a resident caller retains all scientific maps on the device.
nLocations = prod(outputSize);
blockSize = nLocations;
if useGPU
    % Budget for the input, rotated/normalized complex temporaries, masks,
    % reductions and (for an internal trial axis) the linear-index matrix.
    % Leave headroom for persistent result maps, CUDA and other operations.
    bytes = 4 + 4*strcmp(precision,'double');
    scratchPerLocation = 20*bytes*nTrials + 512;
    freeBytes = double(device.AvailableMemory);
    if resident
        % The result retains 16 real maps, two double counts, seven masks
        % and transient storage for population-wide validity summaries.
        freeBytes = freeBytes - nLocations*(16*bytes + 2*8 + 7 + 8);
    end
    budget = min(192*1024^2, .20*freeBytes);
    if budget < scratchPerLocation
        error('nf_cpm:GPUBlockTooLarge', ...
            'Insufficient free GPU memory for CPM output and one complete trial ensemble.');
    end
    blockSize = min(blockSize, floor(budget/scratchPerLocation));
end
if ~isempty(options.MaxBlockLocations)
    blockSize = min(blockSize,options.MaxBlockLocations);
end
dim = options.TrialDimension;
shape = size(coefficients);
before = prod(shape(1:dim-1));
after = prod(shape(dim+1:end));
if after == 1, flat = reshape(coefficients,nLocations,nTrials); end
hasInfinity = false;
for first = 1:blockSize:nLocations
    rows = first:min(nLocations,first+blockSize-1);
    if after == 1
        tile = flat(rows,:);
    else
        % Every row is one non-trial location in MATLAB's column order.
        % Index only this block; a non-final trial axis needs no full permute.
        base = 1 + mod(rows(:)-1,before) + floor((rows(:)-1)/before)*before*nTrials;
        indices = bsxfun(@plus,base,before*(0:nTrials-1));
        tile = reshape(coefficients(indices),numel(rows),nTrials);
    end
    if useGPU
        % Sparse two-dimensional inputs can still represent an ensemble.
        % MATLAB's device elementwise kernels work on this bounded dense tile.
        if issparse(tile), tile = full(tile); end
        if ~resident, tile = gpuArray(tile); end
    end
    % One scalar check at the end avoids synchronizing a resident input
    % after each block. The kernel masks nonfinite values before reductions.
    hasInfinity = hasInfinity | any(isinf(tile(:)));
    block = local_kernel(tile,options,[numel(rows) 1],nTrials,precision);
    if first == 1
        model = block;
        [sections,fields,kinds] = local_map_layout(block);
    end
    for k = 1:numel(kinds)
        key = kinds{k};
        columns = cell(1,numel(key));
        for c = 1:numel(key)
            i = key(c);
            if isempty(sections{i}), columns{c} = block.(fields{i});
            else, columns{c} = block.(sections{i}).(fields{i}); end
        end
        if useGPU && ~resident
            packed = gather(cat(2,columns{:}));
            for c = 1:numel(columns), columns{c} = packed(:,c); end
        end
        for c = 1:numel(key)
            i = key(c);
            section = sections{i};
            field = fields{i};
            if isempty(section)
                if first == 1, model.(field) = zeros(outputSize,'like',columns{c}); end
                model.(field)(rows) = columns{c};
            else
                if first == 1, model.(section).(field) = zeros(outputSize,'like',columns{c}); end
                model.(section).(field)(rows) = columns{c};
            end
        end
    end
    clear block tile packed columns;
end
if useGPU, hasInfinity = gather(hasInfinity); end
if hasInfinity
    error('nf_cpm:InfiniteCoefficients','Infinite coefficients are not supported.');
end
% These are aliases in the published output. Share the assembled arrays.
model.phaseLocked.powerMean = model.phaseLocked.amplitudeSecondMoment;
model.nonPhaseLocked.powerMean = model.nonPhaseLocked.amplitudeSecondMoment;
model.outputSize = outputSize;
end


function [sections,fields,kinds] = local_map_layout(block)
% Discover map fields once, excluding scalar estimator metadata. Keep
% counts double and masks logical, independently of coefficient precision.
sections = {''};
fields = {'leakageRatio'};
types = {local_precision(block.leakageRatio)};
for name = {'total','phaseLocked','nonPhaseLocked','raw','validity'}
    section = name{1};
    list = fieldnames(block.(section));
    for j = 1:numel(list)
        field = list{j};
        if strcmp(field,'powerMean') && ~strcmp(section,'total'), continue; end
        sections{end+1} = section; %#ok<AGROW>
        fields{end+1} = field; %#ok<AGROW>
        types{end+1} = local_precision(block.(section).(field)); %#ok<AGROW>
    end
end
names = unique(types);
kinds = cell(size(names));
for j = 1:numel(names), kinds{j} = find(strcmp(types,names{j})); end
end


function type = local_precision(x)
if isa(x,'gpuArray'), type = classUnderlying(x); else, type = class(x); end
end


function summary = local_summary(model, nTrials, resident)
v = model.validity;
counts = [sum(~v.sufficientTrials(:)), ...
    sum(v.sufficientTrials(:) & v.validTrialCount(:)<nTrials), ...
    sum(v.sufficientTrials(:) & ~v.phaseDefined(:)), ...
    sum(v.phaseDefined(:) & ~v.nonsingular(:)), ...
    sum(v.positiveProjectionCount(:)==0), ...
    sum(v.nonsingular(:) & ~v.phaseLockedSecondMoment(:)), ...
    sum(v.nonsingular(:) & ~v.nonPhaseLockedSecondMoment(:)), ...
    sum(v.phaseLockedSecondMoment(:) & ~v.phaseLockedAmplitudeVariance(:)), ...
    sum(v.nonPhaseLockedSecondMoment(:) & ~v.nonPhaseLockedAmplitudeVariance(:))];
if resident, counts = gather(counts); end
names = {'insufficientTrialLocations','incompleteTrialLocations', ...
    'undefinedPhaseLocations','singularLocations','undefinedNonPhaseLockedMeanLocations', ...
    'invalidPhaseLockedSecondMoments','invalidNonPhaseLockedSecondMoments', ...
    'invalidPhaseLockedAmplitudeVariances','invalidNonPhaseLockedAmplitudeVariances'};
summary = cell2struct(num2cell(counts),names,2);
end


function local_warnings(summary,options)
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
        ['%d coefficient location(s) were at the original CPM formula''s ' ...
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
addParameter(parser, 'Estimator', 'original', @local_is_text_scalar);
addParameter(parser, 'Verbosity', 'auto', @local_is_text_scalar);
addParameter(parser, 'ExecutionEnvironment', 'auto', @local_is_text_scalar);
addParameter(parser, 'MaxBlockLocations', [], ...
    @(n) isempty(n) || local_valid_dimension(n));
parse(parser, coefficients, varargin{:});

options = struct();
options.ExecutionEnvironment = lower(local_text_scalar( ...
    parser.Results.ExecutionEnvironment, 'ExecutionEnvironment'));
if ~any(strcmp(options.ExecutionEnvironment,{'auto','cpu','gpu'}))
    error('nf_cpm:InvalidExecutionEnvironment','Use auto, cpu or gpu.');
end
options.MaxBlockLocations = double(parser.Results.MaxBlockLocations);
options.TrialDimension = double(parser.Results.TrialDimension);
options.Estimator = lower(local_text_scalar( ...
    parser.Results.Estimator, 'Estimator'));
options.Verbosity = lower(local_text_scalar( ...
    parser.Results.Verbosity, ...
    'Verbosity'));

if any(strcmp(options.Estimator, {'original', 'projection'})) == false
    error('nf_cpm:InvalidEstimator', ...
        'Estimator must be ''original'' or ''projection''.');
end

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

if isa(value,'gpuArray')
    valid = any(strcmp(classUnderlying(value),{'single','double'})) && ~isempty(value);
else
    valid = isnumeric(value) && isfloat(value) && ~isempty(value);
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
