function TF = nf_tfspecparam(TF, varargin)
% NF_TFSPECPARAM    Computes a time-resolved spectral parameterization
%
% GENERAL
% -------
% SPRiNT (Wilson et al., 2022) rewritten for use with TF structures output
% by tfUtility/tf_fun. Computes a time-resolved spectral parameterization
% using matspecparam.m.
%
% Accepts TF transforms calculated using any method. Automatically detects
% if power values are on a linear scale and converts to log scale.
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
% TF.SPRiNT fields (all log10 scale)
%   Aperiodic (all log10):
%     TF.SPRiNT.ap_power    : C×F×T×R (fit or data per 'FoD')
%     TF.SPRiNT.osc_power   : C×F×T×R (fit or data per 'FoD')
%     TF.SPRiNT.rsquare     : C×T×R
%     TF.SPRiNT.mae         : C×T×R
%     TF.SPRiNT.ap_offset   : C×T×R
%     TF.SPRiNT.ap_knee     : C×T×R (NaN if 'fixed')
%     TF.SPRiNT.ap_exponent : C×T×R
%   Peaks (ALL peaks per slice, parfor-safe numeric arrays; padded with NaN):
%     TF.SPRiNT.pk_n        : C×T×R    (#peaks detected)
%     TF.SPRiNT.pk_f        : C×K×T×R  (Hz)        K = maxPeaks
%     TF.SPRiNT.pk_aLog10   : C×K×T×R  (log10 power height)
%     TF.SPRiNT.pk_sd       : C×K×T×R  (Hz)
%
% INPUT
% -----
%   TF.power : C×F×T[×R] power
%   TF.scale : 'linear'|'log10'  (if 'linear', converted to log10)
%   TF.freqs : F×1 (Hz)
%
% Name-Value
% ----------
%   'peakWidthLims' : [min max] Hz (default [0.5 12])
%   'maxPeaks'      : max #peaks per slice (default numel(TF.freqs))
%   'minPeakHeight' : minimum peak height (log10 power units if scale->log10) (default 0)
%   'aPeriodicMode' : 'fixed'|'knee' (default 'fixed')
%   'peakThreshold' : noise SD threshold for detection (default 2.0)
%   'peakType'      : 'gaussian'|'cauchy' (default 'gaussian')
%   'threshAfter'   : re-threshold when no optim (default 1)
%   'optim'         : use Optimization Toolbox if available (default 1)
%   'FoD'           : 'fit' returns model fit; 'data' returns modeled data (default 'fit')
%


% Parse options
p = inputParser;
valid2Scalar = @(x) numel(x)==2 && isnumeric(x) && isvector(x);
expectedPeakType = @(x) any(validatestring(x,{'gaussian','cauchy'}));
expectedAPModes  = @(x) any(validatestring(x,{'knee','fixed'}));
expectedFoD      = @(x) any(validatestring(x,{'fit','data'}));

addRequired(p,'TF');
addParameter(p,'peakWidthLims',[0.5 12],valid2Scalar);
addParameter(p,'maxPeaks',numel(TF.freqs),@(x)isscalar(x)&&x>=0);
addParameter(p,'minPeakHeight',0,@isscalar);
addParameter(p,'aPeriodicMode','fixed',expectedAPModes);
addParameter(p,'peakThreshold',2.0,@isscalar);
addParameter(p,'peakType','gaussian',expectedPeakType);
addParameter(p,'threshAfter',1,@(x) isscalar(x)&&ismember(x,[0 1]));
addParameter(p,'optim',1,@(x) isscalar(x)&&ismember(x,[0 1]));
addParameter(p,'FoD','fit',expectedFoD);
parse(p,TF,varargin{:});
opt = p.Results;

% Sanity & shape
if ~isfield(TF,'power') || ~isfield(TF,'freqs') || ~isfield(TF,'scale')
    error('TF.power, TF.freqs, TF.scale are required.');
end
freqs = TF.freqs(:);
surf  = TF.power;             % C×F×T[×R]
if ndims(surf)==3
    surf = reshape(surf, size(surf,1), size(surf,2), size(surf,3), 1);
end
[C,F,T,R] = size4d(surf);
if F ~= numel(freqs)
    error('Second dimension of TF.power (F=%d) differs from numel(TF.freqs) (%d).', F, numel(freqs));
end

% Convert to log10 if needed (we keep TF.power in log10 afterwards)
if ~isfield(TF,'scale') || ~(strcmpi(TF.scale,'linear') || strcmpi(TF.scale,'log10'))
    error('TF.scale must be ''linear'' or ''log10''.');
end
if strcmpi(TF.scale,'linear')
    surf(surf<=0) = eps;
    surf = log10(surf);
    TF.power = surf;
    TF.scale = 'log10';
end

% Allocate outputs
aperiodic_log = zeros(C,F,T,R,'like',surf);
osc_log       = zeros(C,F,T,R,'like',surf);
rsq           = zeros(C,T,R);
mae           = zeros(C,T,R);

ap_offset     = nan(C,T,R);
ap_knee       = nan(C,T,R);
ap_exponent   = nan(C,T,R);

K             = max(1, round(opt.maxPeaks));   % cap
pk_n          = zeros(C,T,R);                  % #peaks detected (untruncated)
pk_f          = nan(C,K,T,R);                  % Hz
pk_aLog10     = nan(C,K,T,R);                  % log10(power) height
pk_sd         = nan(C,K,T,R);                  % Hz

% Main loop
fprintf(1,'SPRiNT (nf_tfspecparam) progress:   0%%');
parfor ch = 1:C
    % local slices (parfor-friendly)
    aperiodic_log_ch = zeros(F,T,R);
    osc_log_ch       = zeros(F,T,R);
    rsq_ch           = zeros(T,R);
    mae_ch           = zeros(T,R);
    ap_offset_ch     = nan(T,R);
    ap_knee_ch       = nan(T,R);
    ap_exponent_ch   = nan(T,R);
    pk_n_ch          = zeros(T,R);
    pk_f_ch          = nan(K,T,R);
    pk_aLog10_ch     = nan(K,T,R);
    pk_sd_ch         = nan(K,T,R);
    
    for r = 1:R
        % extract channel x trial slab → T×F (each row one time-slice)
        spec_log10_TR = squeeze(surf(ch,:,:,r)).'; % T×F
        for t = 1:T
            dataX = spec_log10_TR(t,:);   % 1×F log10 power
            parm = nf_specparam( dataX, freqs, ...
                'aPeriodicMode', opt.aPeriodicMode, ...  
                'minPeakHeight', opt.minPeakHeight, ...
                'peakThreshold', opt.peakThreshold, ...
                'peakType',      opt.peakType, ...
                'peakWidthLims', opt.peakWidthLims, ...
                'threshAfter',   opt.threshAfter, ...
                'optim',         opt.optim, ...
                'plt',           0 );                %#ok
            
            % spectra (log10)
            if strcmpi(opt.FoD,'fit')
                aperiodic_log_ch(:,t,r) = parm.aPeriodicFit(:);
                osc_log_ch(:,t,r)       = parm.PeriodicFit(:);
            else
                aperiodic_log_ch(:,t,r) = parm.aPeriodicData(:);
                osc_log_ch(:,t,r)       = parm.PeriodicData(:);
            end
            
            % fit quality
            rsq_ch(t,r) = parm.RSquare;
            mae_ch(t,r) = parm.MAE;
            
            % aperiodic params → [offset, knee(=NaN if fixed), exponent]
            apv = parm.aperiodicParms(:).';
            if strcmpi(parm.options.aPeriodicMode,'fixed')
                ap_offset_ch(t,r)   = apv(1);
                ap_knee_ch(t,r)     = NaN;
                ap_exponent_ch(t,r) = apv(2);
            else
                ap_offset_ch(t,r)   = apv(1);
                ap_knee_ch(t,r)     = apv(2);
                ap_exponent_ch(t,r) = apv(3);
            end
            
            % ALL peaks (Nx3), truncate/pad to K
            pk = parm.peakParms;           % [f0, ampLog10, sd]
            n  = size(pk,1);
            pk_n_ch(t,r) = n;
            if n > 0
                nuse = min(n, K);
                pk_f_ch(1:nuse,t,r)      = pk(1:nuse,1);
                pk_aLog10_ch(1:nuse,t,r) = pk(1:nuse,2);
                pk_sd_ch(1:nuse,t,r)     = pk(1:nuse,3);
            end
        end
    end
    
    % write back this channel’s slab
    aperiodic_log(ch,:,:,:) = aperiodic_log_ch;
    osc_log(ch,:,:,:)       = osc_log_ch;
    rsq(ch,:,:)             = rsq_ch;
    mae(ch,:,:)             = mae_ch;
    
    ap_offset(ch,:,:)       = ap_offset_ch;
    ap_knee(ch,:,:)         = ap_knee_ch;
    ap_exponent(ch,:,:)     = ap_exponent_ch;
    
    pk_n(ch,:,:)            = pk_n_ch;
    pk_f(ch,:,:,:)          = pk_f_ch;
    pk_aLog10(ch,:,:,:)     = pk_aLog10_ch;
    pk_sd(ch,:,:,:)         = pk_sd_ch;
    
    % light progress (best-effort; harmless in parfor)
    fprintf(1,'\b\b\b\b%3d%%', round(100*ch/C));
end
fprintf(1,'\n');

% Package outputs
TF.SPRiNT = struct();
TF.SPRiNT.ap_power    = aperiodic_log;      % C×F×T×R (log10)
TF.SPRiNT.osc_power   = osc_log;            % C×F×T×R (log10)
TF.SPRiNT.rsquare     = rsq;                % C×T×R
TF.SPRiNT.mae         = mae;                % C×T×R

TF.SPRiNT.ap_offset   = ap_offset;          % C×T×R
TF.SPRiNT.ap_knee     = ap_knee;            % C×T×R
TF.SPRiNT.ap_exponent = ap_exponent;        % C×T×R

% ALL peaks
TF.SPRiNT.pk_n        = pk_n;               % C×T×R (#peaks detected)
TF.SPRiNT.pk_f        = pk_f;               % C×K×T×R (Hz)
TF.SPRiNT.pk_aLog10   = pk_aLog10;          % C×K×T×R (log10 power height)
TF.SPRiNT.pk_sd       = pk_sd;              % C×K×T×R (Hz)

end % function


% helpers
function [C,F,T,R] = size4d(X)
sz = size(X);
C = sz(1);
F = sz(2);
T = 1; R = 1;
if numel(sz)>=3, T = sz(3); end
if numel(sz)>=4, R = sz(4); end
end








