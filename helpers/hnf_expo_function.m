function ys = hnf_expo_function(freqs, params, logFreqs)
%       Exponential function to use for fitting 1/f, with a 'knee' (maximum at low frequencies).
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Input x-axis values.
%       params : 1x3 array (offset, knee, exp)
%           Parameters (offset, knee, exp) that define Lorentzian function:
%           y = 10^offset * (1/(knee + x^exp))
%
%       Returns
%       -------
%       ys :    1xn array
%           Output values for exponential function.

    if nargin < 3
        logFreqs = log(freqs);
    end

    frequencyPower = exp(params(3) .* logFreqs);
    ys = params(1) - log10(abs(params(2)) + frequencyPower);

end
