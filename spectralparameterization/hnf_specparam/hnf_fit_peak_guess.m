function [peakParams, usedBounds] = hnf_fit_peak_guess( ...
    guess, freqs, flatSpec, guessWeight, ...
    widthLimits, useOptimizationToolbox, peakType, approvedPeakKernel)
% HNF_FIT_PEAK_GUESS Jointly optimize a set of candidate spectral peaks.

if nargin < 7, peakType = 'gaussian'; end
if isempty(guessWeight), guessWeight = 'weak'; end
usedBounds = false;
if isempty(guess)
    peakParams = zeros(0, 3);
    return
end

% Choose the objective once per fit. Solver callbacks can run thousands of
% times. Native and MATLAB paths use the same guesses, penalties and solver.
useNative = exist('hnf_mex_available','file') == 2 && ...
    hnf_mex_available('hnf_peak_objective_mex') && ...
    isa(guess,'double') && isa(freqs,'double') && isa(flatSpec,'double') && ...
    isreal(guess) && isreal(freqs) && isreal(flatSpec) && ...
    ~issparse(guess) && ~issparse(freqs) && ~issparse(flatSpec);
if nargin >= 8 && ischar(approvedPeakKernel)
    % An empty client path disables the worker kernel. A path mismatch can
    % occur when a worker has a differently configured filesystem or path.
    useNative = useNative && ~isempty(approvedPeakKernel) && ...
        exist('hnf_peak_objective_mex','file') == 3 && ...
        strcmp(which('hnf_peak_objective_mex'),approvedPeakKernel);
end
if useNative
    shape = double(strcmp(peakType,'cauchy'));
    penalty = 1e2 * strcmp(guessWeight,'weak') + 1e7 * strcmp(guessWeight,'strong');
    % A stale or foreign binary can be discoverable but fail to load. Probe
    % before choosing the solver so fmincon can still use MATLAB gradients.
    try
        initialError = hnf_peak_objective_mex(guess, freqs, flatSpec, guess, penalty, shape); %#ok<NASGU>
    catch err
        if ~hnf_mex_load_error(err), rethrow(err); end
        useNative = false;
    end
end
if useNative
    freeObjective = @(p) hnf_peak_objective_mex(p,freqs,flatSpec,guess,penalty,shape);
else
    freeObjective = @(p) hnf_error_model(p,freqs,flatSpec,guess,guessWeight,peakType);
end

if useOptimizationToolbox == 1
    canUseFmincon = specparam_has_fmincon_local();
    fallbackOptions = specparam_requested_fallback_options_local();

    octaveRuntime = exist('OCTAVE_VERSION','builtin') ~= 0;
    if canUseFmincon || octaveRuntime
        lowerBounds = [guess(:, 1) - guess(:, 3) .* 2, ...
            zeros(size(guess(:, 2))), ...
            ones(size(guess(:, 3))) .* widthLimits(1)];
        upperBounds = [guess(:, 1) + guess(:, 3) .* 2, ...
            inf(size(guess(:, 2))), ...
            ones(size(guess(:, 3))) .* widthLimits(2)];

        if useNative
            objective = @(params) hnf_peak_objective_mex(params, freqs, flatSpec, [], 0, shape);
        else
            objective = @(params) hnf_error_model_constr(params, freqs, flatSpec, peakType);
        end
        if canUseFmincon
            try
                [peakParams,~,exitFlag] = fmincon(objective, guess(:), ...
                    [], [], [], [], lowerBounds(:), upperBounds(:), [], ...
                    specparam_constrained_options_local());
                usedBounds = ~octaveRuntime || (exitFlag > 0 && all(isfinite(peakParams)) && ...
                    all(peakParams >= lowerBounds(:)-1e-7) && all(peakParams <= upperBounds(:)+1e-7));
            catch
                if ~octaveRuntime
                    peakParams = fminsearch(freeObjective, guess, fallbackOptions);
                end
            end
        end
        if octaveRuntime && ~usedBounds
            % Core Octave also supplies a constrained solver. Keep the
            % same peak SSE and bounds when the optim package is absent.
            peakParams = octave_peak_sqp(objective,guess(:),lowerBounds(:),upperBounds(:));
            usedBounds = true;
        end
    else
        peakParams = fminsearch(freeObjective, guess, fallbackOptions);
    end
else
    unconstrainedOptions = specparam_unconstrained_options_local();
    peakParams = fminsearch(freeObjective, guess, unconstrainedOptions);
end

peakParams = reshape(peakParams, [], 3);
% Both densities depend on scale squared, so a scale always denotes a magnitude.
peakParams(:,3) = abs(peakParams(:,3));

end

function available = specparam_has_fmincon_local()

% Paths and toolbox licenses may change between calls in one session.
available = exist('fmincon', 'file') == 2 && ...
    (exist('OCTAVE_VERSION','builtin') || license('test','Optimization_Toolbox'));
end

function p = octave_peak_sqp(objective,p,lb,ub)
% Eliminate parameters fixed by equal bounds (e.g., a prescribed width).
free = lb < ub;
p(~free) = lb(~free);
if ~any(free), return; end
% Express only finite bounds as linear inequalities. Octave 9 sqp's
% lb/ub interface mishandles mixed finite and infinite upper bounds.
lower = lb(free);
upper = ub(free);
eyeFree = eye(nnz(free));
A = [eyeFree(isfinite(lower),:); -eyeFree(isfinite(upper),:)];
b = [lower(isfinite(lower)); -upper(isfinite(upper))];
constraints = {@(x) A*x-b,@(x) A};
[x,value,info] = sqp(p(free),{@value_at,@gradient_at},[],constraints, ...
    [],[],3000,1e-7);
if ~any(info==[101 104]) || ~isfinite(value) || any(~isfinite(x)) || ...
        any(x < lb(free)-1e-7) || any(x > ub(free)+1e-7)
    error('NeuroFreq:OctavePeakFit','Bounded Octave peak fit failed (sqp status %d).',info);
end
p(free) = min(ub(free),max(lb(free),x));

    function value = value_at(x)
        trial = p;
        trial(free) = x;
        value = objective(trial);
    end
    function gradient = gradient_at(x)
        trial = p;
        trial(free) = x;
        [~,g] = objective(trial);
        gradient = g(free);
    end
end

function options = specparam_constrained_options_local()

persistent cachedOptions

if isempty(cachedOptions)
    cachedOptions = optimset('Display', 'off', 'TolX', 1e-3, ...
        'TolFun', 1e-5, 'MaxFunEvals', 3000, 'MaxIter', 3000, ...
        'GradObj', 'on');
end

options = cachedOptions;

end

function options = specparam_requested_fallback_options_local()

persistent cachedOptions

if isempty(cachedOptions)
    cachedOptions = optimset('Display', 'off', 'TolX', 1e-3, ...
        'TolFun', 1e-5, 'MaxFunEvals', 3000, 'MaxIter', 3000);
end

options = cachedOptions;

end

function options = specparam_unconstrained_options_local()

persistent cachedOptions

if isempty(cachedOptions)
    cachedOptions = optimset('Display', 'off', 'TolX', 1e-4, ...
        'TolFun', 1e-5, 'MaxFunEvals', 5000, 'MaxIter', 5000);
end

options = cachedOptions;

end
