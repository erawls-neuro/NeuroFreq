function [apFit, periodicFit, apData, periodicData, rSquare, mae, apParms, peakParms] = ...
    hnf_specparam_fit_compact( ...
    spec, freqs, options, frequencyCache, computeMetrics)
% HNF_SPECPARAM_FIT_COMPACT Numerical spectral-parameterization core.
%
% This function assumes that options have already been validated. It avoids
% inputParser, plotting, and construction of the full specparam result
% structure. The fitted model and objective are unchanged. A precomputed
% cache from specparam_frequency_cache can be supplied as the fourth input;
% repeated callers are then responsible for validating the frequency grid.
% Set the optional fifth input to false when goodness-of-fit metrics will be
% discarded. All fitted signals and parameters remain unchanged.

if nargin < 5 || isempty(computeMetrics)
    computeMetrics = true;
end

if isa(spec, 'double') && isrow(spec)
    workSpec = spec;
else
    workSpec = double(spec(:).');
end

cacheWasSupplied = nargin >= 4 && ~isempty(frequencyCache);

if cacheWasSupplied == true
    workFreqs = frequencyCache.freqs;
else
    workFreqs = double(freqs(:).');
end

if numel(workSpec) ~= numel(workFreqs)
    error('hnf_specparam_fit_compact:SizeMismatch', ...
        'spec and freqs must contain the same number of elements.');
end

if any(~isfinite(workSpec))
    error('hnf_specparam_fit_compact:NonfiniteSpectrum', ...
        'spec must contain only finite log10 power values.');
end

if cacheWasSupplied == false
    if any(~isfinite(workFreqs)) || any(workFreqs <= 0)
        error('hnf_specparam_fit_compact:InvalidFrequencies', ...
            'freqs must contain only finite positive values.');
    end

    frequencyCache = hnf_specparam_frequency_cache(workFreqs);
end

if isfield(options, 'ag') && ~isempty(options.ag)
    aperiodicGuess = options.ag;
else
    aperiodicGuess = workSpec(end) - workSpec(1) ./ ...
        log10(workFreqs(end) ./ workFreqs(1));
end

if options.maxPeaks == 0 && strcmp(options.aPeriodicMode, 'fixed')
    peakParms = zeros(0, 3);
    periodicFit = zeros(size(workSpec));
    apData = workSpec;
    apParms = hnf_simple_ap_fit(workFreqs, apData, options.aPeriodicMode, ...
        aperiodicGuess, frequencyCache);
    apFit = hnf_gen_aperiodic(workFreqs, apParms, options.aPeriodicMode, ...
        frequencyCache);
    periodicData = workSpec - apFit;

    if computeMetrics == true
        mae = mean(abs(periodicData));
        rSquare = specparam_rsquare_local(workSpec, apFit);
    else
        mae = NaN;
        rSquare = NaN;
    end
    return
end

initialApParms = hnf_robust_ap_fit(workFreqs, workSpec, ...
    options.aPeriodicMode, aperiodicGuess, frequencyCache);

initialApFit = hnf_gen_aperiodic(workFreqs, initialApParms, ...
    options.aPeriodicMode, frequencyCache);
initialPeriodicData = workSpec - initialApFit;

[peakParms, ~] = hnf_fit_peaks(workFreqs, initialPeriodicData, ...
    options.maxPeaks, options.peakThreshold, options.minPeakHeight, ...
    options.peakWidthLims ./ 2, [], options.optim);

if options.threshAfter == 1 && options.optim == 0 && ~isempty(peakParms)
    keepPeak = peakParms(:, 2) >= options.minPeakHeight;
    keepPeak = keepPeak & peakParms(:, 3) >= options.peakWidthLims(1) ./ 2;
    keepPeak = keepPeak & peakParms(:, 3) <= options.peakWidthLims(2) ./ 2;
    keepPeak = keepPeak & peakParms(:, 1) >= 0;
    peakParms = peakParms(keepPeak, :);
end

if ~isempty(peakParms)
    lowerEdgeDistance = peakParms(:, 1) - workFreqs(1);
    upperEdgeDistance = workFreqs(end) - peakParms(:, 1);
    keepPeak = lowerEdgeDistance >= frequencyCache.frequencyResolution;
    keepPeak = keepPeak & upperEdgeDistance >= frequencyCache.frequencyResolution;
    peakParms = peakParms(keepPeak, :);
end

periodicFit = hnf_specparam_peak_model(workFreqs, peakParms);
apData = workSpec - periodicFit;

apParms = hnf_simple_ap_fit(workFreqs, apData, options.aPeriodicMode, ...
    initialApParms(end), frequencyCache);
apFit = hnf_gen_aperiodic(workFreqs, apParms, options.aPeriodicMode, ...
    frequencyCache);
periodicData = workSpec - apFit;

if computeMetrics == true
    modelFit = apFit + periodicFit;
    residual = workSpec - modelFit;
    mae = mean(abs(residual));
    rSquare = specparam_rsquare_local(workSpec, modelFit);
else
    mae = NaN;
    rSquare = NaN;
end

end

function rSquare = specparam_rsquare_local(observed, fitted)

if any(~isfinite(observed)) || any(~isfinite(fitted))
    rSquare = NaN;
    return
end

observedCentered = observed - mean(observed);
fittedCentered = fitted - mean(fitted);

observedSumSquares = sum(observedCentered .* observedCentered);
fittedSumSquares = sum(fittedCentered .* fittedCentered);
denominator = observedSumSquares .* fittedSumSquares;

if denominator > 0 %#ok
    correlation = sum(observedCentered .* fittedCentered) ./ sqrt(denominator);
    rSquare = correlation .* correlation;
else
    rSquare = NaN;
end


end
