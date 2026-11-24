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

    % Evaluate log(knee + f^exponent) directly on the log scale.
    % For a zero knee this is exactly the fixed power law.
    logPower = params(3).*logFreqs;
    if params(2) == 0
        logSum = logPower;
    else
        logKnee = log(abs(params(2)));
        logSum = max(logKnee,logPower) + log1p(exp(-abs(logKnee-logPower)));
    end
    ys = params(1) - logSum./log(10);

end
