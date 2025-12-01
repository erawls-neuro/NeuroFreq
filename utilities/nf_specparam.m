function parm_spec = nf_specparam( spec, freqs, varargin )
% NF_SPECPARAM    MATLAB implementation of specparam.
%
% GENERAL
% -------
% MATLAB implementation of specparam.
% Parameterizes power spectra into separate 
% oscillatory and aperiodic components.
% Code is adapted as originally from implementation for BrainStorm.
%
% If code is used please cite:
% Wilson, L. E., da Silva Castanheira, J., & Baillet, S. (2022). 
% Time-resolved parameterization of aperiodic and periodic brain activity. 
% Elife, 11, e77348.
%
% Please also cite BrainStorm software:
% François Tadel, Sylvain Baillet, John C. Mosher, Dimitrios Pantazis, 
% Richard M. Leahy, "Brainstorm: A User-Friendly Application for MEG/EEG 
% Analysis", Computational Intelligence and Neuroscience, vol. 2011, 
% Article ID 879716, 13 pages, 2011. https://doi.org/10.1155/2011/879716
%
% The spectral parameterization method should be cited as 
%
% Donoghue, T., Haller, M., Peterson, E. J., Varma, P., Sebastian, P., 
% Gao, R., ... & Voytek, B. (2020). Parameterizing neural power spectra 
% into periodic and aperiodic components. Nature neuroscience, 23(12), 
% 1655-1665.
%
%
% OUTPUT
% ------
%
% parm_spec - contains results of spectral parameterization analysis
%
% INPUT
% -----
% 
% spec - power spectrum, must be log10!
% f - vector of frequencies in spec
% various keyword-arg pairs: specparam settings detailed below
% 
% -----
% E. Rawls, erawls89@gmail.com, rawls017@umn.edu. 
% July 2023
% Copyright (c) 2023 by E. Rawls.
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA


%make input parser
p                      = inputParser;
%define valid inputs
validBinary            = @(x) x==0 || x==1;
validScalar            = @(x) length(x)==1;
valid2Scalar           = @(x) numel(x)==2 && isnumeric(x) && isvector(x);
expectedPeakType       = @(x) any(validatestring(x,{'gaussian' 'cauchy'}));
expectedAperiodicModes = @(x) any(validatestring(x,{'knee' 'fixed'}));
% defaults for specparam
defaultPeakWidthLimits = [.5 12]; %bandwidth limits for peaks
defaultMaxPeaks        = numel(freqs); %max n peaks
defaultMinPeakHeight   = 0; %min peak height in power
defaultAperiodicMode   = 'fixed'; %aperiodic can be 'fixed' 'knee'
defaultPeakThreshold   = 2.0; %noise SD for Gaussian detection
defaultPeakType        = 'gaussian'; %can also be cauchy, but untested
defaultThreshAfter     = 1; %threshold after fitting? only used if no-optim
defaultOptim           = 1; %use optimization toolbox? mimics python
defaultPlot            = 1; %make plot?
defaultAG              = spec(end)-spec(1)./log10(freqs(end)./freqs(1)); %ap-exp guess
%add REQUIRED parameters to parser
addRequired(p,'spec'); %power spectrum
addRequired(p,'f'); %frequency vector
%add OPTIONAL parameters to parser
addParameter(p,'peakWidthLims',defaultPeakWidthLimits,valid2Scalar);
addParameter(p,'maxPeaks',defaultMaxPeaks,validScalar);
addParameter(p,'minPeakHeight',defaultMinPeakHeight,validScalar);
addParameter(p,'aPeriodicMode',defaultAperiodicMode,expectedAperiodicModes);
addParameter(p,'peakThreshold',defaultPeakThreshold,validScalar);
addParameter(p,'peakType',defaultPeakType,expectedPeakType);
addParameter(p,'threshAfter',defaultThreshAfter,validBinary);
addParameter(p,'optim',defaultOptim,validBinary);
addParameter(p,'ag',defaultAG,validScalar);
addParameter(p,'plt',defaultPlot,validBinary);
%parse it
parse(p,spec,freqs,varargin{:});
%done parsing!










%aperiodic parameters
aperiodic_pars=robust_ap_fit(freqs,spec, p.Results.aPeriodicMode, ...
    p.Results.ag);

%remove aperiodic
flat_spec = flatten_spectrum(freqs, spec, aperiodic_pars, ...
    p.Results.aPeriodicMode);

% Fit peaks
[peak_pars, peak_function] = fit_peaks(double(freqs), double(flat_spec), ...
    p.Results.maxPeaks, p.Results.peakThreshold, p.Results.minPeakHeight, ...
    p.Results.peakWidthLims/2, ...
    p.Results.peakType, [], p.Results.optim);
if p.Results.threshAfter && ~p.Results.optim  % Check thresholding requirements are met for unbounded optimization
    peak_pars(peak_pars(:,2) < opt.min_peak_height,:)     = []; % remove peaks shorter than limit
    peak_pars(peak_pars(:,3) < opt.peak_width_limits(1)/2,:)  = []; % remove peaks narrower than limit
    peak_pars(peak_pars(:,3) > opt.peak_width_limits(2)/2,:)  = []; % remove peaks broader than limit
    peak_pars(peak_pars(:,1) < 0,:) = []; % remove peaks with a centre frequency less than zero
end
peak_pars( find( (peak_pars(:,1)-freqs(1)) < mode(diff(freqs)) ),: ) = []; %sometimes the function fits peaks that are impossible?
peak_pars( find( freqs(end)-peak_pars(:,1) < mode(diff(freqs)) ),: ) = []; %sometimes the function fits peaks that are impossible?

% Refit aperiodic
aperiodic = spec;
for peak = 1:size(peak_pars,1)
    aperiodic = aperiodic - peak_function(freqs,peak_pars(peak,1), peak_pars(peak,2), peak_pars(peak,3));
end
aperiodic_pars = simple_ap_fit(freqs, aperiodic, p.Results.aPeriodicMode, aperiodic_pars(end));

%re-calculate a better flattened spectrum
flat_spec = flatten_spectrum(freqs, spec, aperiodic_pars, ...
    p.Results.aPeriodicMode);

% Generate model fit
ap_fit = gen_aperiodic(freqs, aperiodic_pars, p.Results.aPeriodicMode);
model_fit = ap_fit;
peak_fit = zeros(1,numel(ap_fit));
for peak = 1:size(peak_pars,1)
    peak_fit = peak_fit + peak_function(freqs,peak_pars(peak,1),...
        peak_pars(peak,2),peak_pars(peak,3));
    model_fit = model_fit + peak_function(freqs,peak_pars(peak,1),...
        peak_pars(peak,2),peak_pars(peak,3));
end

% Calculate model error
MAE = sum(abs(spec-model_fit))/length(model_fit);
rsq_tmp = corrcoef(spec,model_fit).^2;

%generate output
parm_spec.options        = p.Results;
parm_spec.aperiodicParms = aperiodic_pars;
parm_spec.peakParms      = peak_pars;
parm_spec.f              = freqs;
parm_spec.data           = spec;
parm_spec.modelFit       = model_fit;
parm_spec.aPeriodicData  = aperiodic;
parm_spec.aPeriodicFit   = ap_fit;
parm_spec.PeriodicData   = flat_spec;
parm_spec.PeriodicFit    = peak_fit;
parm_spec.MAE            = MAE;
parm_spec.RSquare        = rsq_tmp(2);

%make a nice plot
if p.Results.plt  
    plotMatSpecParam( parm_spec ); 
end
    
end