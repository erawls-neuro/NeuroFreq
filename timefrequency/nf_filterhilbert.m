function tfRes = nf_filterhilbert(data, Fs, freqs, fBandWidth, order, verbose)
% NF_FILTERHILBERT    Calculates time-frequency of an input dataset (1/2/3D) using filter-Hilbert.
%
% GENERAL
% -------
% Calculates time-frequency of an input dataset (1/2/3D) using 
% filter-Hilbert.
% 
% Expected dimensions: channel (1), time samples (2), trials/segments (3;
% optional). 
%
% Allows filter bandwidth to change by frequency, mimicking 
% wavelet smoothing that changes by frequency. Uses Butterworth filters of
% specified order.
%
% OUTPUT
% ------
% tfRes: structure with fields 
%   1) power
%   2) phase
%   3) frequencies
%   4) times
%   5) sampling rate
%   6) bandwidth
%   7) order
%
% INPUT
% -----
% 1) data: 1D(time), 2D(channelXtime) or 3D(channelXtimesample) (REQUIRED)
% 2) Fs: sampling rate of signal, in Hz (REQUIRED)
% 3) freqs: requested frequencies in Hz, defaults to 1:1:Fs/2
% 4) fBandWidth: bandwidth of Butterworth filters, formatted 1 or [1 8] 
%       i.e. start-and stop-bandwidth, defaults to 1
% 5) order: order of Butterworth filters, defaults to 3
% 6) verbose: print verbose output? 0 or 1, defaults to 0
%

% defaults
if nargin < 6 || isempty(verbose)
    verbose = 0;
end
if nargin < 5 || isempty(order)
    order = 3;
    if verbose == 1 
        disp('[nf_filterhilbert]: setting filter order = 3 (default)');
    end
end
if nargin < 4 || isempty(fBandWidth)
    fBandWidth = 1;
    if verbose == 1 
        disp('[nf_filterhilbert]: setting filter bandwidth to 1 Hz (default)');
    end
end
if nargin < 3 || isempty(freqs)
    freqs = 1:1:floor((Fs / 2) - 1);
    if verbose == 1 
        disp('[nf_filterhilbert]: setting freqs to 1:1:Fs/2 (default)');
    end
end
if nargin < 2 || isempty(Fs) || isempty(data)
    error('nf_filterhilbert:InputsRequired',...
        'at least a signal and sampling rate are required');
end

% sizes
nChan  = size(data, 1);
nTimes = size(data, 2);

if ndims(data) == 3
    nTrls = size(data, 3);
else
    nTrls = 1;
end

% concatenate time across trials
data = reshape(data, nChan, nTimes * nTrls);
nTot  = size(data, 2);

nFreq = numel(freqs);

% parse bandwidth specification
fLap = fBandWidth;

if isscalar(fLap)
    if verbose == 1 
        disp(['[nf_filterhilbert]: using ' num2str(fLap) ' Hz bandwidth for all filters']);
    end
    bwidths = repmat(fLap, 1, nFreq);
elseif numel(fLap) == 2
    if verbose == 1 
        disp(['[nf_filterhilbert]: filter bandwidth will linearly scale from ' num2str(fLap(1)) ...
          ' to ' num2str(fLap(2))]);
    end
    bwidths = linspace(fLap(1), fLap(2), nFreq);
else
    error('nf_filterhilbert:FBandwidthBad',...
        'fBandWidth must be scalar or [start end]');
end

% preallocate
filtPow  = zeros(nChan, nFreq, nTot);
filtPhas = zeros(nChan, nFreq, nTot);

% progress
disp('[nf_filterhilbert]: beginning filter-hilbert');

% loop over frequencies (filters designed per band)
for j = 1:nFreq
    
    fc = freqs(j);
    bw = bwidths(j);
    hipF  = fc - (bw / 2);
    lowpF = fc + (bw / 2);
   
    % clamp to valid passband
    if hipF <= 0
        hipF = 0.01;
    end
    nyq = Fs / 2;
    if lowpF >= nyq
        lowpF = nyq - 0.01;
    end
    if lowpF <= hipF
        warning('nf_filterhilbert:invalidBand', ...
                'Invalid band for freq %.2f Hz (hipF=%.2f, lowpF=%.2f). Skipping.', ...
                fc, hipF, lowpF);
        continue
    end
    Wn = [hipF lowpF] ./ nyq;
  
    % design bandpass Butterworth (correct b,a order)
    [b, a] = butter(order, Wn, 'bandpass');
    
    % filter all channels at once (time along rows)
    band = filtfilt(b, a, double(data.')).';
    
    % analytic signal and power/phase
    if isreal(data)
        
        % hilbert operates along columns, so transpose twice
        analytic = hilbert(band.').';
        mag2     = abs(analytic).^2;
        phs      = angle(analytic);
        
    else
        
        % input already complex/analytic
        mag2 = abs(band).^2;
        phs  = angle(band);
        
    end
    
    % store (explicit reshape to avoid surprises)
    filtPow(:, j, :)  = reshape(mag2, nChan, 1, nTot);
    filtPhas(:, j, :) = reshape(phs,   nChan, 1, nTot);
    
end

% reshape back to [chan x freq x time x trial]
tfRes.power  = squeeze(reshape(filtPow,  nChan, nFreq, nTimes, nTrls));
tfRes.phase  = squeeze(reshape(filtPhas, nChan, nFreq, nTimes, nTrls));
tfRes.freqs  = freqs;
tfRes.times  = 0:1 / Fs:((1 / Fs) * nTimes) - (1 / Fs);
tfRes.nsensor = nChan;
tfRes.ntrls   = nTrls;
tfRes.Fs      = Fs;
tfRes.bandwidth = bwidths;
tfRes.order     = order;
tfRes.method    = 'filter-hilbert';
tfRes.scale     = 'linear';

disp('[nf_filterhilbert]: filter-hilbert complete');

end









