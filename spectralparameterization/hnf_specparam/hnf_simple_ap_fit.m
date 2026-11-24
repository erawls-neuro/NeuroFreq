function aperiodicParams = hnf_simple_ap_fit(freqs, powerSpectrum, aperiodicMode, aperiodicGuess, frequencyCache)
% SIMPLE_AP_FIT Fit the aperiodic component of a log10 power spectrum.

freqs = double(freqs(:).');
powerSpectrum = double(powerSpectrum(:).');
if nargin < 4 || isempty(aperiodicGuess)
    aperiodicGuess = hnf_specparam_slope_guess(freqs,powerSpectrum);
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
            loss = @(p) hnf_error_expo_nk_function(p,freqs,powerSpectrum,frequencyCache.log10Freqs);
            aperiodicParams = fminsearch(loss,guessVector,searchOptions);
        end

    case 'knee'
        searchOptions = specparam_ap_search_options_local();
        guessVector = [powerSpectrum(1) 0 aperiodicGuess];
        loss = hnf_knee_loss(freqs,powerSpectrum,frequencyCache.logFreqs,guessVector,frequencyCache);
        aperiodicParams = fminsearch(loss,guessVector,searchOptions);
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
