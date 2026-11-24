function aperiodicParams = hnf_robust_ap_fit(freqs, powerSpectrum, aperiodicMode, aperiodicGuess, frequencyCache)
% HNF_ROBUST_AP_FIT Fit the aperiodic component while suppressing positive peaks.

freqs = double(freqs(:).');
powerSpectrum = double(powerSpectrum(:).');
if nargin < 4 || isempty(aperiodicGuess)
    aperiodicGuess = hnf_specparam_slope_guess(freqs,powerSpectrum);
end

if nargin < 5 || isempty(frequencyCache)
    frequencyCache = hnf_specparam_frequency_cache(freqs);
end

initialParams = hnf_simple_ap_fit(freqs, powerSpectrum, aperiodicMode, ...
    aperiodicGuess, frequencyCache);
initialFit = hnf_gen_aperiodic(freqs, initialParams, aperiodicMode, ...
    frequencyCache);

flatSpectrum = powerSpectrum - initialFit;
flatSpectrum(flatSpectrum < 0) = 0;

if numel(flatSpectrum) <= 2000 && all(~isnan(flatSpectrum))
    % Under MATLAB's midpoint percentile definition, the minimum occupies
    % percentile 50/N. Therefore the 0.025th percentile is exactly the
    % minimum whenever N is no greater than 2000.
    percentileThreshold = min(flatSpectrum);
else
    percentileThreshold = prctile(flatSpectrum, 0.025);
end

percentileMask = flatSpectrum <= percentileThreshold;
fitFreqs = freqs(percentileMask);
fitSpectrum = powerSpectrum(percentileMask);

if nnz(percentileMask) < numel(initialParams)
    aperiodicParams = initialParams;
    return
end

switch aperiodicMode
    case 'fixed'
        fitLog10Freqs = frequencyCache.log10Freqs(percentileMask);
        [aperiodicParams, fullRank] = hnf_specparam_fixed_ap_fit( ...
            fitLog10Freqs, fitSpectrum);

        if fullRank == false
            searchOptions = specparam_ap_search_options_local();
            loss = @(p) hnf_error_expo_nk_function(p,fitFreqs,fitSpectrum,fitLog10Freqs);
            aperiodicParams = fminsearch(loss,initialParams,searchOptions);
        end

    case 'knee'
        searchOptions = specparam_ap_search_options_local();
        fitLogFreqs = frequencyCache.logFreqs(percentileMask);
        loss = hnf_knee_loss(fitFreqs,fitSpectrum,fitLogFreqs,initialParams,frequencyCache);
        aperiodicParams = fminsearch(loss,initialParams,searchOptions);
end

end

function options = specparam_ap_search_options_local()

persistent cachedOptions

if isempty(cachedOptions)
    cachedOptions = optimset('Display', 'off', 'TolX', 1e-4, ...
        'TolFun', 1e-6, 'MaxFunEvals', 5000, 'MaxIter', 5000);
end

options = cachedOptions;

end
