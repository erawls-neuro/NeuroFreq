function aperiodicParams = hnf_simple_ap_fit(freqs, powerSpectrum, aperiodicMode, aperiodicGuess, frequencyCache)
% SIMPLE_AP_FIT Fit the aperiodic component of a log10 power spectrum.

if nargin < 4 || isempty(aperiodicGuess)
    aperiodicGuess = -(powerSpectrum(end) - powerSpectrum(1)) ./ ...
        log10(freqs(end) ./ freqs(1));
end

if nargin < 5 || isempty(frequencyCache)
    frequencyCache = hnf_specparam_frequency_cache(freqs);
end

switch aperiodicMode
    case 'fixed'
        [aperiodicParams, fullRank] = hnf_specparam_fixed_ap_fit( ...
            frequencyCache.log10Freqs, powerSpectrum, frequencyCache);

        if fullRank == false
            searchOptions = specparam_ap_search_options_local();
            guessVector = [powerSpectrum(1) aperiodicGuess];
            aperiodicParams = fminsearch(@error_expo_nk_function, ...
                guessVector, searchOptions, freqs, powerSpectrum, ...
                frequencyCache.log10Freqs);
        end

    case 'knee'
        searchOptions = specparam_ap_search_options_local();
        guessVector = [powerSpectrum(1) 0 aperiodicGuess];
        aperiodicParams = fminsearch(@error_expo_function, guessVector, ...
            searchOptions, freqs, powerSpectrum, frequencyCache.logFreqs);
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
