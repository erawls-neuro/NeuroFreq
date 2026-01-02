function [tfRes,p] = nf_tftransform(EEG,varargin)
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
% Name-Value Arguments:
% METHOD-GENERAL:
%
%   'method': method for tf calculation. Can be any supported - new methods
%       should be made in conformance with the formatting of the transforms.
%       Find supported methods below and in the transforms folder.
%
%   'freqs': REQUIRED for wavelet,filterhilbert,demodulation. Must be a 
%                       vector of numbers (i.e. 1:30). Only some methods 
%                       return exact freq vectors - other methods will 
%                       return native freqs within the range specified.
%
%   'times': will trim times off the ends to accomodate the requested
%                       range. Can be a vector of times [0:4:1000] or a 
%                       max-min pair, but only the ends are used.
%
% METHOD-SPECIFIC
%
%   wavelet: 'cycles': Must be either a vector of length freqs, OR 
%                       a max-min pair i.e. [3 8], OR a single number i.e. 
%                       4. If a max-min pair is provided a linear space 
%                       over the range will be used. 
%   
%   filterhilbert: 'fBandWidth': bandwidth of frequency filters.
%                       may be a single element i.e. 1 (1-Hz bandwidth at 
%                       all frequencies) or a pair of values i.e. [1 8]. 
%                       If a max-min pair is provided a linear space 
%                       over the range will be used. 
%   
%   filterhilbert,demodulation: 'order': order of filters (e.g. 3).
%   
%   stft: 'window': window in seconds (e.g. 0.5) for temporal 
%                       smoothing, 'overlap': STFT overlap in percent.
%
%   demodulation: 'lowpassF': low-pass F of butterworth
%
%   stft,binomial2,bornjordan2,ridrihaczek: 'fRes': frequency resolution 
%                       i.e. 0.5 Hz, 'makePos': make surfaces positive?
%
%   binomial2,bornjordan2: 'maxlags': length of autocorrelation window
%
%   ridrihaczek: 'cwkernel': Choi-Williams kernel size
%
%



% parse input
p = inputParser;
validNumVector = @(x) isnumeric(x) && isvector(x);
validMonotonicVector = @(x) isnumeric(x) && isvector(x) && sum(x<0)==0 ...
    && sum(diff(x)<0)==0;
validScalar = @(x) length(x)==1;

%set legal methods
expectedMethods = {'stft' ... %spectrogram
                   'filterhilbert' ... %filter and Hilbert transform
                   'demodulation'... %complex demodulation
                   'dcwt'... %Custom wavelet transform
                   'cwt'... %continuous wavelet transform
                   'stransform'... %Stockwell transform
                   'ridbinomial',... %type-II binomial RID
                   'ridbornjordan',... %type-II Born-Jordan RID
                   'ridrihaczek'}; %RID-Rihaczek
                   %any new methods added here!

%DEFAULTS                   
%all methods
defaultMethod = 'dcwt';
defaultFreqs = []; %1:1:floor(EEG.srate/4);
defaultTimes = [EEG.xmin EEG.xmax];
%pass on defaults
defaultCycles = [];%3;
defaultBandWidth = [];%1;
defaultOrder = [];%3;
defaultWindow = [];%.5;
defaultOverlap = [];%80;
defaultlowpassF = [];%3;
defaultMaxLag = [];%2*EEG.pnts;
defaultCWKernel = [];%0.001;
defaultfRes = [];%0.25;
defaultMakePos = [];%0;

%INPUTS
addRequired(p,'EEG'); %EEGLAB set
%all methods
addParameter(p,'method',defaultMethod,@(x) any(validatestring(x,expectedMethods))); %basically required
addParameter(p,'freqs',defaultFreqs,validMonotonicVector); %freq vector
addParameter(p,'times',defaultTimes,validNumVector); %time vector
%only some methods
addParameter(p,'window',defaultWindow,validScalar); %stft
addParameter(p,'overlap',defaultOverlap,validScalar); %stft
addParameter(p,'fBandwidth',defaultBandWidth); %filter-hilbert
addParameter(p,'order',defaultOrder,validScalar); %filter-hilbert,demodulation
addParameter(p,'lowpassF',defaultlowpassF,validScalar); %demodulation
addParameter(p,'cycles',defaultCycles); %wavelet
addParameter(p,'fRes',defaultfRes,validScalar); %binomial2,bornjordan2,ridrihaczek
addParameter(p,'maxLags',defaultMaxLag); %binomial2,bornjordan2
addParameter(p,'cwkernel',defaultCWKernel); %ridRihaczek
addParameter(p,'makePos',defaultMakePos,validScalar); %binomial2,bornjordan2,ridrihaczek

%parse
parse(p,EEG,varargin{:});

%get all inputs for ease
method=p.Results.method;
freqs=p.Results.freqs;
window=p.Results.window;
overlap=p.Results.overlap;
fBandwidth=p.Results.fBandwidth;
order=p.Results.order;
lowpassF=p.Results.lowpassF;
cycles=p.Results.cycles;
fRes=p.Results.fRes;
maxLags=p.Results.maxLags;
cwkernel=p.Results.cwkernel;
makePos=p.Results.makePos;
reqT=p.Results.times;
Fs=EEG.srate;
data=EEG.data;
times = EEG.times/1000;

% OPTIONS
method = lower(method);
switch method
    %linear operations
    case 'stft'
        tfRes = nf_stft(data,Fs,...
                window,...
                overlap,...
                fRes);
        times = times(1)+tfRes.times; %refigure times
    case 'filterhilbert'
        tfRes = nf_filterhilbert(data,Fs,...
                freqs,...
                fBandwidth,...
                order);
    case 'demodulation'
        tfRes = nf_demodulation(data,Fs,...
                freqs,...
                lowpassF,...
                order);
    case {'dcwt','wavelet'}
        tfRes = nf_dcwt(data,Fs,...
                freqs,...
                cycles);
    case 'cwt'
            tfRes = nf_cwt(data,Fs);
    case 'stransform'
        tfRes = nf_stransform(data,Fs);            
    %quadratic tfds
    case 'ridbinomial'
        tfRes = nf_ridbinomial(data,Fs,...
                fRes,...
                maxLags,...
                makePos);
    case 'ridbornjordan'
        tfRes = nf_ridbornjordan(data,Fs,...
                fRes,...
                maxLags,...
                makePos);
    case 'ridrihaczek'
        tfRes = nf_ridrihaczek(data,Fs,...
                fRes,...
                cwkernel,...
                makePos);

end
%END OPTIONS
tfRes.times = times;

%trim freqs & times
if isempty(freqs)
    freqs = tfRes.freqs;
end
fIndices = find( (tfRes.freqs>=freqs(1)) + (tfRes.freqs<=freqs(end)) == 2);
tIndices = find( (tfRes.times>=reqT(1)) + (tfRes.times<=reqT(end)) == 2);

if EEG.nbchan>1
    tfRes.power=tfRes.power(:,fIndices,tIndices,:);
    if isfield(tfRes,'phase') %some methods do not return phases
        tfRes.phase=tfRes.phase(:,fIndices,tIndices,:);
    end
else
    tfRes.power=tfRes.power(fIndices,tIndices,:);
    if isfield(tfRes,'phase') %some methods do not return phases
        tfRes.phase=tfRes.phase(fIndices,tIndices,:);
    end
end
tfRes.freqs = tfRes.freqs(fIndices);
tfRes.times = tfRes.times(tIndices);

%get channel info
tfRes.chanlocs = EEG.chanlocs;

%get events
if ~isempty(EEG.event) %if there are events
    tfRes.event = EEG.event;
end
if ~isempty(EEG.epoch) %if there are trials
    tfRes.epoch = EEG.epoch;
end
if isfield(EEG.etc,'behavior')
    tfRes.behavior=EEG.etc.behavior; %if behavioral structure exists
end


end