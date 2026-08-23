function peakParams = hnf_fit_peak_guess( ...
    guess, freqs, flatSpec, ~, ...
    widthLimits, useOptimizationToolbox)
% HNF_FIT_PEAK_GUESS Jointly optimize a set of candidate spectral peaks.

if isempty(guess)
    peakParams = zeros(0, 3);
    return
end

if useOptimizationToolbox == 1
    canUseFmincon = specparam_has_fmincon_local();
    constrainedOptions = specparam_constrained_options_local();
    fallbackOptions = specparam_requested_fallback_options_local();

    if canUseFmincon == true
        lowerBounds = [guess(:, 1) - guess(:, 3) .* 2, ...
            zeros(size(guess(:, 2))), ...
            ones(size(guess(:, 3))) .* widthLimits(1)];
        upperBounds = [guess(:, 1) + guess(:, 3) .* 2, ...
            inf(size(guess(:, 2))), ...
            ones(size(guess(:, 3))) .* widthLimits(2)];

        try
            parameterVector = fmincon(@error_model_constr, guess(:), ...
                [], [], [], [], lowerBounds(:), upperBounds(:), [], ...
                constrainedOptions, freqs, flatSpec);
            peakParams = reshape(parameterVector, [], 3);
        catch
            peakParams = fminsearch(@error_model, guess, ...
                fallbackOptions, freqs, flatSpec, guess, 'weak');
        end
    else
        peakParams = fminsearch(@error_model, guess, fallbackOptions, ...
            freqs, flatSpec, guess, 'weak');
    end
else
    unconstrainedOptions = specparam_unconstrained_options_local();
    peakParams = fminsearch(@error_model, guess, unconstrainedOptions, ...
        freqs, flatSpec, guess, 'weak');
end

peakParams = reshape(peakParams, [], 3);

end

function available = specparam_has_fmincon_local()

persistent cachedAvailability

if isempty(cachedAvailability)
    cachedAvailability = exist('fmincon', 'file') == 2;

    if cachedAvailability == true
        if exist('license', 'builtin') == 5 || exist('license', 'file') == 2
            try
                cachedAvailability = license('test', 'Optimization_Toolbox');
            catch
                cachedAvailability = true;
            end
        end
    end
end

available = cachedAvailability;

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
