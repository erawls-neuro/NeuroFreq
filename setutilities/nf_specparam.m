function parm_spec = nf_specparam(spec, freqs, varargin)
% NF_SPECPARAM Parameterize one log10 power spectrum.
%
% parm_spec = nf_specparam(spec, freqs, Name, Value, ...)
%
% This MATLAB implementation was adapted from the implementation for
% Brainstorm. If used, please retain the citations below.
%
% This is the validating entry point. The numerical work is delegated
% to specparam_fit so repeated calls can bypass inputParser after their
% options have been validated once.
%
% INPUTS
% ------
% spec  : F-element row or column log10 power spectrum.
% freqs : F-element row or column frequency vector in Hz.
%
% NAME-VALUE OPTIONS
% ------------------
% peakWidthLims : [minimum maximum] peak bandwidth in Hz. Default [0.5 12].
% maxPeaks      : Maximum number of candidate peaks. Default numel(freqs).
% minPeakHeight : Minimum candidate peak height. Default 0.
% aPeriodicMode : 'fixed' or 'knee'. Default 'fixed'.
% peakThreshold : Candidate threshold in residual standard deviations.
% threshAfter   : Apply intended post-fit limits for unconstrained fits.
% optim         : Use fmincon when Optimization Toolbox is available.
% ag            : Initial aperiodic-exponent guess.
%
% REFERENCES
% ----------
% Wilson, L. E., da Silva Castanheira, J., and Baillet, S. (2022).
% Time-resolved parameterization of aperiodic and periodic brain activity.
% eLife, 11, e77348.
%
% Donoghue, T. et al. (2020). Parameterizing neural power spectra into
% periodic and aperiodic components. Nature Neuroscience, 23, 1655-1665.
%
% Tadel, F., Baillet, S., Mosher, J. C., Pantazis, D., and Leahy, R. M.
% (2011). Brainstorm: A user-friendly application for MEG/EEG analysis.
% Computational Intelligence and Neuroscience, 2011, 879716.
% https://doi.org/10.1155/2011/879716

p = inputParser;

validBinary = @(x) isscalar(x) && (x == 0 || x == 1);
validScalar = @(x) isscalar(x);
validTwoElementNumericVector = @(x) numel(x) == 2 && isnumeric(x) && isvector(x);
expectedAperiodicMode = @(x) any(validatestring(x, {'knee', 'fixed'}));

defaultPeakWidthLimits = [0.5 12];
defaultMaxPeaks = numel(freqs);
defaultMinPeakHeight = 0;
defaultAperiodicMode = 'fixed';
defaultPeakThreshold = 2.0;
defaultThresholdAfter = 1;
defaultOptimization = 1;
defaultAperiodicGuess = spec(end) - spec(1) ./ log10(freqs(end) ./ freqs(1));

addRequired(p, 'spec');
addRequired(p, 'f');
addParameter(p, 'peakWidthLims', defaultPeakWidthLimits, validTwoElementNumericVector);
addParameter(p, 'maxPeaks', defaultMaxPeaks, validScalar);
addParameter(p, 'minPeakHeight', defaultMinPeakHeight, validScalar);
addParameter(p, 'aPeriodicMode', defaultAperiodicMode, expectedAperiodicMode);
addParameter(p, 'peakThreshold', defaultPeakThreshold, validScalar);
addParameter(p, 'threshAfter', defaultThresholdAfter, validBinary);
addParameter(p, 'optim', defaultOptimization, validBinary);
addParameter(p, 'ag', defaultAperiodicGuess, validScalar);

parse(p, spec, freqs, varargin{:});

parm_spec = hnf_specparam_fit(spec, freqs, p.Results);

end
