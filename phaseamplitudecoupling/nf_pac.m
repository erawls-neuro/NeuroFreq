function pac = nf_pac(TF, foiLow, foiHi, varargin)
% NF_PAC    Phase–amplitude coupling (PAC) using NeuroFreq TF output.
%
% Phase–amplitude coupling (PAC) using NeuroFreq TF output. Supports:
%   - MVL-based PAC (Canolty 2006, Cohen 2014 variant)
%   - MI-based PAC (Tort et al., 2008)
%   - ERPAC (Voytek et al., 2013)
%
% and:
%   - within-channel PAC (default)
%   - cross-channel PAC between specified channel sets
%   - optional MEX acceleration (MI / MVL / circ corr / time shift)
%
% GENERAL BEHAVIOR
% ----------------
% - Uses TF.phase and TF.power from nf_tftransform.m
% - Selects frequency bands via foiLow / foiHi ranges
% - Selects time-of-interest as the contiguous span between the first and
%   last element of 'toi'
% - Permutation count nperm controls z-scoring:
%       nperm == 0  -> raw metric (MVL / MI / ERPAC)
%       nperm  > 0  -> z-scored metric via null permutations/surrogates
%
% CROSS-CHANNEL MODE
% ------------------
% When mode is 'between' or 'cross', PAC is computed between:
%   - phase channels: phaseCh
%   - amplitude channels: ampCh
%
% Output dimensions:
%   - MVL / MI:
%       within: [chan x loFreq x hiFreq x trial]
%       cross : [nPhaseChan x nAmpChan x loFreq x hiFreq x trial]
%
%   - ERPAC:
%       within: [chan x loFreq x hiFreq x time]
%       cross : [nPhaseChan x nAmpChan x loFreq x hiFreq x time]
%
% USAGE
% -----
% pac = nf_pac(TF, foiLow, foiHi, ...
%              'toi',    [tStart tEnd], ...
%              'method', 'MVL' | 'MI' | 'ERPAC', ...
%              'mode',   'within' | 'between' | 'cross', ...
%              'nperm',  0 or positive integer, ...
%              'nbin',   positive integer (MI only), ...
%              'phaseCh', [phase channel indices], ...
%              'ampCh',   [amp channel indices], ...
%              'use_mex', true / false);
%
% REQUIRED INPUTS
%   TF      : struct from nf_tftransform.m with fields:
%             .phase [chan x freq x time x trial]
%             .power [chan x freq x time x trial]
%             .freqs [1 x nFreq]
%             .times [1 x nTimes]
%             .Fs
%             .chanlocs
%   foiLow  : [fMin fMax] low-frequency range (phase)
%   foiHi   : [fMin fMax] high-frequency range (amplitude)
%
% OPTIONAL PARAMETERS (Name–Value)
%   'toi'      : time range [tStart tEnd], default = TF.times
%   'method'   : 'MVL' (default), 'MI', or 'ERPAC'
%   'mode'     : 'within' (default), or 'between'/'cross'
%   'nperm'    : permutations for z-scoring (default 100)
%   'nbin'     : phase bins for MI (default 18)
%   'phaseCh'  : phase channel indices (cross mode; default = all)
%   'ampCh'    : amplitude channel indices (cross mode; default = all)
%   'use_mex'  : logical, enable MEX acceleration (default true)
%

% Input parsing
validNumVector = @(x) isnumeric(x) && isvector(x);
validScalarNum = @(x) isnumeric(x) && isscalar(x);
expectedMethods = {'MVL','MI','ERPAC'};
expectedModes   = {'within','between','cross'};
defaultMethod  = 'MVL';
defaultTois    = TF.times;
defaultNperm   = 100;
defaultMode    = 'within';
defaultNbin    = 18;
defaultUseMex  = true;
defaultPhaseCh = [];
defaultAmpCh   = [];
p = inputParser;
p.FunctionName = 'nf_pac';
addRequired(p, 'TF');
addRequired(p, 'foiLow', validNumVector);
addRequired(p, 'foiHi',  validNumVector);
addParameter(p, 'toi',     defaultTois, validNumVector);
addParameter(p, 'method',  defaultMethod, ...
    @(x) any(strcmpi(x, expectedMethods)));
addParameter(p, 'mode',    defaultMode, ...
    @(x) any(strcmpi(x, expectedModes)));
addParameter(p, 'nperm',   defaultNperm, validScalarNum);
addParameter(p, 'nbin',    defaultNbin,  validScalarNum);
addParameter(p, 'phaseCh', defaultPhaseCh, validNumVector);
addParameter(p, 'ampCh',   defaultAmpCh,   validNumVector);
addParameter(p, 'use_mex', defaultUseMex, ...
    @(x) (islogical(x) && isscalar(x)) || (isnumeric(x) && isscalar(x)));
parse(p, TF, foiLow, foiHi, varargin{:});

% Extract parsed values
toi     = p.Results.toi;
method  = lower(p.Results.method);
modeStr = lower(p.Results.mode);
nperm   = p.Results.nperm;
nbin    = p.Results.nbin;
phaseCh = p.Results.phaseCh;
ampCh   = p.Results.ampCh;
use_mex = logical(p.Results.use_mex);

% greetings
disp(['[nf_pac]: phase-amplitude coupling using ' method]);
% Sanity checks on scalar parameters
if ~isscalar(nperm) || nperm < 0 || nperm ~= round(nperm)
    error('nf_pac:BadNPerm', 'nperm must be a non-negative integer scalar.');
end

if ~isscalar(nbin) || nbin < 1 || nbin ~= round(nbin)
    error('nf_pac:BadNBin', 'nbin must be a positive integer scalar.');
end

isCross = strcmp(modeStr, 'between') || strcmp(modeStr, 'cross');

% Frequency and time selection
freqs = TF.freqs;
times = TF.times;

if numel(foiLow) < 2
    error('nf_pac:BadFOILow', 'foiLow must specify at least [fMin fMax].');
end

if numel(foiHi) < 2
    error('nf_pac:BadFOIHi', 'foiHi must specify at least [fMin fMax].');
end

fLowMin = min(foiLow);
fLowMax = max(foiLow);
fHiMin  = min(foiHi);
fHiMax  = max(foiHi);

loFreqs = find(freqs >= fLowMin & freqs <= fLowMax);
hiFreqs = find(freqs >= fHiMin  & freqs <= fHiMax);

if isempty(loFreqs)
    error('nf_pac:NoLow', 'no low frequencies found in requested range [%g %g].', fLowMin, fLowMax);
end

if isempty(hiFreqs)
    error('nf_pac:NoHigh', 'no high frequencies found in requested range [%g %g].', fHiMin, fHiMax);
end

tMin = min(toi);
tMax = max(toi);

timeInd = find(times >= tMin & times <= tMax);

if isempty(timeInd)
    error('nf_pac:NoTimes', 'no time points found in requested range [%g %g].', tMin, tMax);
end

% Extract TF data and basic dimension checks
dataLow = TF.phase(:, loFreqs, timeInd, :);
dataHi  = TF.power(:, hiFreqs, timeInd, :);

[nChanL, ~, nTimes, nTrialL]     = size(dataLow);
[nChanH, ~, nTimesH, nTrialH]    = size(dataHi);

if nChanL ~= nChanH
    error('nf_pac:PhasePowerDiffChan', 'phase and power must have the same number of channels.');
end

if nTimesH ~= nTimes
    error('nf_pac:PhasePowerDiffTimes', 'phase and power must have the same number of time points.');
end

if nTrialH ~= nTrialL
    error('nf_pac:PhasePowerDiffTrials', 'phase and power must have the same number of trials.');
end

nChan = nChanL;

% Channel selection for cross mode
if isCross
    if isempty(phaseCh)
        phaseCh = 1:nChan;
    end

    if isempty(ampCh)
        ampCh = 1:nChan;
    end

    phaseCh = phaseCh(:)';
    ampCh   = ampCh(:)';

    if any(phaseCh < 1) || any(phaseCh > nChan) || any(phaseCh ~= round(phaseCh))
        error('nf_pac:PhaseChBadInd', 'phaseCh must contain valid integer channel indices within [1 %d].', nChan);
    end

    if any(ampCh < 1) || any(ampCh > nChan) || any(ampCh ~= round(ampCh))
        error('nf_pac:AmpChBadInd', 'ampCh must contain valid integer channel indices within [1 %d].', nChan);
    end

    nPhaseChan = numel(phaseCh);
    nAmpChan   = numel(ampCh);
    nPairs     = nPhaseChan * nAmpChan;

    if nPairs > 1024
        warning('nf_pac:largeCrossGrid', ...
            ['Cross-channel PAC requested for %d phase × %d amp channels (%d pairs).\n' ...
             'This can be memory- and time-intensive, especially with nperm = %d.'], ...
             nPhaseChan, nAmpChan, nPairs, nperm);
    end
else
    if ~isempty(phaseCh) || ~isempty(ampCh)
        warning('nf_pac:ignoredChannels', ...
            ['phaseCh / ampCh provided but mode = ''within''; ', ...
             'ignoring channel subsets and using all channels within-channel.']);
    end
end

% PAC computation (delegated to metric-specific functions)
switch method
    case 'mvl'
        if isCross
            pacRes = nf_mvlpac_cross(dataLow, dataHi, nperm, phaseCh, ampCh, use_mex);
        else
            pacRes = nf_mvlpac_within(dataLow, dataHi, nperm, use_mex);
        end

    case 'mi'
        if isCross
            pacRes = nf_mipac_cross(dataLow, dataHi, nperm, nbin, phaseCh, ampCh, use_mex);
        else
            pacRes = nf_mipac_within(dataLow, dataHi, nperm, nbin, use_mex);
        end

    case 'erpac'
        if isCross
            pacRes = nf_erpac_cross(dataLow, dataHi, nperm, phaseCh, ampCh, use_mex);
        else
            pacRes = nf_erpac_within(dataLow, dataHi, nperm, use_mex);
        end

    otherwise
        error('nf_pac: unknown method "%s".', method);
end

% Assemble output struct
pac.PAC      = pacRes;
pac.method   = method;
pac.mode     = modeStr;
pac.nperm    = nperm;
pac.use_mex  = use_mex;
pac.nbin     = nbin;
pac.loFreqs  = freqs(loFreqs);
pac.hiFreqs  = freqs(hiFreqs);
pac.times    = times(timeInd);
pac.Fs       = TF.Fs;

% Channel metadata
pac.chanlocs = TF.chanlocs;

if isCross
    pac.phaseCh        = phaseCh;
    pac.ampCh          = ampCh;
    pac.chanlocs_phase = TF.chanlocs(phaseCh);
    pac.chanlocs_amp   = TF.chanlocs(ampCh);
end

% Method-specific extras
if strcmp(method, 'erpac')
    pac.type = 'erpac';
else
    if isfield(TF, 'event')
        pac.event = TF.event;
    end

    if isfield(TF, 'epoch')
        pac.epoch = TF.epoch;
    end

    pac.type = method;
end

end








