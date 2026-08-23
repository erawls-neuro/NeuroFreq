function ap_vals = hnf_gen_aperiodic(freqs, aperiodic_params, aperiodic_mode, frequencyCache)
%       Generate aperiodic values, from parameter definition.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency vector to create aperiodic component for.
%       aperiodic_params : 1x3 array
%           Parameters that define the aperiodic component.
%       aperiodic_mode : {'fixed', 'knee'}
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       ap_vals : 1d array
%           Generated aperiodic values.

    if nargin < 4 || isempty(frequencyCache)
        frequencyCache = hnf_specparam_frequency_cache(freqs);
    end

    switch aperiodic_mode
        case 'fixed'
            ap_vals = aperiodic_params(1) - ...
                aperiodic_params(2) .* frequencyCache.log10Freqs;
        case 'knee'
            ap_vals = hnf_expo_function(freqs, aperiodic_params, ...
                frequencyCache.logFreqs);
    end
end
