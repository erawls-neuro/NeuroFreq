function [tfRes, p] = nf_tftransform(EEG, varargin)
% NF_TFTRANSFORM    Quickly and flexibly return time-frequency transformations for input EEG data.
%
% GENERAL
% -------
% The main TF transformation function of NeuroFreq. Quickly and flexibly
% return time-frequency transformations for input EEG data. Inputs must be
% EEGLAB .set structures; for raw data matrices use the transform directly
% rather than nf_tftransform. Uses various methods, see the readme in the
% transforms directory for details.
%
% OUTPUT:
% -------
% 1) TF - a structure with information about the TF analysis, including
%   a) Times - times returned
%   b) Chanlocs - information about included channels
%   c) Freqs - vector of frequencies
%   d) Power - TF power
%   e) Phase - TF phase
%   f) Event - set of time-locking events (one per trial) from EEGLAB
%   g) Behavior - EEG set of trial-by-trial behaviors (if present)
%
% INPUT:
% ------
% 1) EEG - an eeglab structure, required
%

% parse input
p = inputParser;
validNumVector = @(x) isnumeric(x) && isvector(x);
validMonotonicVector = @(x) isnumeric(x) && isvector(x) && sum(x < 0) == 0 ...
    && sum(diff(x) < 0) == 0;
validScalar = @(x) isscalar(x);

%set legal methods
expectedMethods = {'stft' ...
                   'filterhilbert' ...
                   'demodulation' ...
                   'dcwt' ...
                   'cwt' ...
                   'stransform' ...
                   'ridbinomial', ...
                   'ridbornjordan', ...
                   'ridrihaczek'};

%DEFAULTS
defaultMethod = 'dcwt';
defaultFreqs = [];
defaultTimes = [EEG.xmin EEG.xmax];
defaultCycles = [];
defaultBandWidth = [];
defaultOrder = [];
defaultWindow = [];
defaultOverlap = [];
defaultlowpassF = [];
defaultMaxLag = [];
defaultCWKernel = [];
defaultfRes = [];
defaultMakePos = [];
defaultVerbose = [];

%INPUTS
addRequired(p, 'EEG');
addParameter(p, 'method', defaultMethod, @(x) any(validatestring(x, expectedMethods)));
addParameter(p, 'freqs', defaultFreqs, validMonotonicVector);
addParameter(p, 'times', defaultTimes, validNumVector);
addParameter(p, 'window', defaultWindow, validScalar);
addParameter(p, 'overlap', defaultOverlap, validScalar);
addParameter(p, 'fBandwidth', defaultBandWidth);
addParameter(p, 'order', defaultOrder, validScalar);
addParameter(p, 'lowpassF', defaultlowpassF, validScalar);
addParameter(p, 'cycles', defaultCycles);
addParameter(p, 'fRes', defaultfRes, validScalar);
addParameter(p, 'maxLags', defaultMaxLag);
addParameter(p, 'cwkernel', defaultCWKernel);
addParameter(p, 'makePos', defaultMakePos, validScalar);
addParameter(p, 'verbose', defaultVerbose, validScalar)

%parse
parse(p, EEG, varargin{:});

%get all inputs for ease
method = p.Results.method;
freqs = p.Results.freqs;
window = p.Results.window;
overlap = p.Results.overlap;
fBandwidth = p.Results.fBandwidth;
order = p.Results.order;
lowpassF = p.Results.lowpassF;
cycles = p.Results.cycles;
fRes = p.Results.fRes;
maxLags = p.Results.maxLags;
cwkernel = p.Results.cwkernel;
makePos = p.Results.makePos;
reqT = p.Results.times;
verbose = p.Results.verbose;
Fs = EEG.srate;
data = EEG.data;
times = EEG.times / 1000;

% GREETINGS 
disp(['[nf_tftransform]: beginning time-frequency transformation using ' method]);

% OPTIONS
method = lower(method);
switch method
    case 'stft'
        tfRes = nf_stft(data, Fs, ...
                window, ...
                overlap, ...
                fRes, ...
                verbose);
        times = times(1) + tfRes.times;
    case 'filterhilbert'
        tfRes = nf_filterhilbert(data, Fs, ...
                freqs, ...
                fBandwidth, ...
                order, ...
                verbose);
    case 'demodulation'
        tfRes = nf_demodulation(data, Fs, ...
                freqs, ...
                lowpassF, ...
                order, ...
                verbose);
    case {'dcwt', 'wavelet'}
        tfRes = nf_dcwt(data, Fs, ...
                freqs, ...
                cycles, ...
                verbose);
    case 'cwt'
        tfRes = nf_cwt(data, Fs, verbose);
    case 'stransform'
        tfRes = nf_stransform(data, Fs, verbose);
    case 'ridbinomial'
        tfRes = nf_ridbinomial(data, Fs, ...
                fRes, ...
                maxLags, ...
                makePos, ...
                verbose);
    case 'ridbornjordan'
        tfRes = nf_ridbornjordan(data, Fs, ...
                fRes, ...
                maxLags, ...
                makePos, ...
                verbose);
    case 'ridrihaczek'
        tfRes = nf_ridrihaczek(data, Fs, ...
                fRes, ...
                cwkernel, ...
                makePos, ...
                verbose);
end
%END OPTIONS
tfRes.times = times;

%trim freqs & times
if isempty(freqs)
    freqs = tfRes.freqs;
end
fIndices = find((tfRes.freqs >= freqs(1)) + (tfRes.freqs <= freqs(end)) == 2);
tIndices = find((tfRes.times >= reqT(1)) + (tfRes.times <= reqT(end)) == 2);

if EEG.nbchan > 1
    tfRes.power = tfRes.power(:, fIndices, tIndices, :);
    if isfield(tfRes, 'phase')
        tfRes.phase = tfRes.phase(:, fIndices, tIndices, :);
    end
else
    tfRes.power = tfRes.power(fIndices, tIndices, :);
    if isfield(tfRes, 'phase')
        tfRes.phase = tfRes.phase(fIndices, tIndices, :);
    end
end
tfRes.freqs = tfRes.freqs(fIndices);
tfRes.times = tfRes.times(tIndices);

%get channel info
tfRes.chanlocs = EEG.chanlocs;

%get events
if ~isempty(EEG.event)
    tfRes.event = EEG.event;
end
if ~isempty(EEG.epoch)
    tfRes.epoch = EEG.epoch;
end
if isfield(EEG.etc, 'behavior')
    tfRes.behavior = EEG.etc.behavior;
end

% Compatibility fields required by downstream NF functions
if ~isfield(tfRes, 'nsensor')
    tfRes.nsensor = EEG.nbchan;
end

if ~isfield(tfRes, 'ntrls')
    tfRes.ntrls = size(tfRes.power, ndims(tfRes.power));
end

if ~isfield(tfRes, 'Fs')
    tfRes.Fs = Fs;
end

if ~isfield(tfRes, 'scale') || isempty(tfRes.scale)
    tfRes.scale = 'linear';
end

disp(['[nf_tftransform]: time-frequency transformation using ' method ' complete']);

end














