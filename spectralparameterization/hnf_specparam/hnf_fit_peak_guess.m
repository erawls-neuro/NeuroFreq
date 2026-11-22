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
    freeObjective = @hnf_peak_objective_mex;
    freeArgs = {freqs, flatSpec, guess, penalty, shape};
else
    freeObjective = @hnf_error_model;
    freeArgs = {freqs, flatSpec, guess, guessWeight, peakType};
end

if useOptimizationToolbox == 1
    canUseFmincon = specparam_has_fmincon_local();
    fallbackOptions = specparam_requested_fallback_options_local();

    if canUseFmincon == true
        lowerBounds = [guess(:, 1) - guess(:, 3) .* 2, ...
            zeros(size(guess(:, 2))), ...
            ones(size(guess(:, 3))) .* widthLimits(1)];
        upperBounds = [guess(:, 1) + guess(:, 3) .* 2, ...
            inf(size(guess(:, 2))), ...
            ones(size(guess(:, 3))) .* widthLimits(2)];

        try
            constrainedOptions = specparam_constrained_options_local();
            if useNative
                objective = @(params) hnf_peak_objective_mex(params, freqs, flatSpec, [], 0, shape);
            else
                objective = @(params) hnf_error_model_constr(params, freqs, flatSpec, peakType);
            end
            parameterVector = fmincon(objective, guess(:), ...
                [], [], [], [], lowerBounds(:), upperBounds(:), [], ...
                constrainedOptions);
            usedBounds = true;
            peakParams = reshape(parameterVector, [], 3);
        catch
            peakParams = fminsearch(freeObjective, guess, fallbackOptions, freeArgs{:});
        end
    else
        peakParams = fminsearch(freeObjective, guess, fallbackOptions, freeArgs{:});
    end
else
    unconstrainedOptions = specparam_unconstrained_options_local();
    peakParams = fminsearch(freeObjective, guess, unconstrainedOptions, freeArgs{:});
end

peakParams = reshape(peakParams, [], 3);
% Both densities depend on scale squared, so a scale always denotes a magnitude.
peakParams(:,3) = abs(peakParams(:,3));

end

function available = specparam_has_fmincon_local()

% Paths and toolbox licenses may change between calls in one session.
available = exist('fmincon', 'file') == 2 && license('test','Optimization_Toolbox');
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
