function [peakParams, usedBounds] = hnf_fit_peak_guess( ...
    guess, freqs, flatSpec, guessWeight, ...
    widthLimits, useOptimizationToolbox, peakType)
% HNF_FIT_PEAK_GUESS Jointly optimize a set of candidate spectral peaks.

if nargin < 7, peakType = 'gaussian'; end
if isempty(guessWeight), guessWeight = 'weak'; end
usedBounds = false;
if isempty(guess)
    peakParams = zeros(0, 3);
    return
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
            objective = @(params) hnf_error_model_constr(params, freqs, flatSpec, peakType);
            parameterVector = fmincon(objective, guess(:), ...
                [], [], [], [], lowerBounds(:), upperBounds(:), [], ...
                constrainedOptions);
            usedBounds = true;
            peakParams = reshape(parameterVector, [], 3);
        catch
            peakParams = fminsearch(@hnf_error_model, guess, ...
                fallbackOptions, freqs, flatSpec, guess, guessWeight, peakType);
        end
    else
        peakParams = fminsearch(@hnf_error_model, guess, fallbackOptions, ...
            freqs, flatSpec, guess, guessWeight, peakType);
    end
else
    unconstrainedOptions = specparam_unconstrained_options_local();
    peakParams = fminsearch(@hnf_error_model, guess, unconstrainedOptions, ...
        freqs, flatSpec, guess, guessWeight, peakType);
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
