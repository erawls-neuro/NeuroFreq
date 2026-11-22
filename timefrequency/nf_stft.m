function tfRes = nf_stft(data,Fs,window,overlap,fRes,verbose)
% NF_STFT    Calculates time-frequency of an input dataset (1/2/3D) using short-time Fourier transform.
%
% GENERAL
% -------
% Calculates time-frequency of an input dataset (1/2/3D) using short-time
% Fourier transform with symmetric Hamming windows and a frame-wise FFT.
% If 100% overlap is specified the transform
% will have the same sampling rate as the original data, although it will
% be missing 1/2 the window width from the beginning and the end.
%
% Expected dimensions: channel (1), time samples (2), trials/segments (3;
% optional).
%
%
% OUTPUT
% ------
% tfRes: structure with fields
%   1) power
%   2) phase
%   3) frequencies
%   4) times
%   5) sampling rate
%
% INPUT
% -----
% 1) data: 1D(time), 2D(channelXtime) or 3D(channelXtimesample) (REQUIRED)
% 2) Fs: sampling rate of signal, in Hz (REQUIRED)
% 3) window: window width, in s (Hamming), defaults to 0.5
% 4) overlap: overlap of STFTs in percent, defaults to 80%
% 5) fRes: frequency resolution, defaults to 1 Hz
% 6) verbose: print verbose output? 0 or 1, defaults to 0
%


narginchk(2,6);
data = hnf_tf_input(data, Fs);
if nargin < 3 || isempty(window), window = 0.5; end
if nargin < 4 || isempty(overlap), overlap = 80; end
if nargin < 5 || isempty(fRes), fRes = 1; end
if nargin < 6 || isempty(verbose), verbose = 0; end
validateattributes(window, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(overlap, {'numeric'}, {'real','scalar','finite','>=',0,'<=',100});
validateattributes(fRes, {'numeric'}, {'real','scalar','finite','positive'});
[nChan,nTimes,nTrls] = size(data);
window = round(window*Fs);
if window < 2 || window > nTimes
    error('nf_stft:WindowLength','The window must contain 2 to %d samples.',nTimes);
end
nfft = round(Fs/fRes);
if nfft < window
    warning('nf_stft:fResTooCoarse','The requested window requires at least %d FFT samples.',window);
    nfft = window;
end
hop = max(1,window-round(window*overlap/100));
starts = 1:hop:nTimes-window+1;
% Symmetric Hamming window, as in hamming(window). Use the explicit frame
% DFT to keep the same scaling and bin order for real and complex signals.
taper = 0.54 - 0.46*cos(2*pi*(0:window-1)/(window-1));
nFreq = floor(nfft/2)+1;
gain = ones(1,nFreq)/sum(taper);
if isreal(data)
    % Peak-amplitude convention for real sinusoids. Only interior bins have
    % a separate negative-frequency partner to combine with the positive one.
    gain(2:end) = 2*gain(2:end);
    if mod(nfft,2) == 0, gain(end) = gain(end)/2; end
end
power = zeros(nChan,nFreq,numel(starts),nTrls);
phase = power;
for frame = 1:numel(starts)
    indices = starts(frame)+(0:window-1);
    spectrum = fft(bsxfun(@times, data(:,indices,:), taper), nfft, 2);
    spectrum = bsxfun(@times, spectrum(:,1:nFreq,:), gain);
    power(:,:,frame,:) = reshape(abs(spectrum).^2,nChan,nFreq,1,nTrls);
    phase(:,:,frame,:) = reshape(angle(spectrum),nChan,nFreq,1,nTrls);
end
tfRes.power = hnf_tf_output(power);
tfRes.phase = hnf_tf_output(phase);
tfRes.freqs = (0:nFreq-1)*Fs/nfft;
% MATLAB spectrogram timestamps each frame half a window from its start.
tfRes.times = (starts-1+window/2)/Fs;
tfRes.nsensor = nChan;
tfRes.ntrls = nTrls;
tfRes.Fs = Fs/hop;
tfRes.window = window;
tfRes.overlap = overlap;
tfRes.method = 'stft';
% A frame DFT is referenced to the beginning of its local window.
tfRes.phaseReference = 'frame-start-carrier';
tfRes.phaseReferenceSeconds = window/(2*Fs);
tfRes.scale = 'linear';
if verbose, disp('[nf_stft]: short time Fourier transform complete'); end
end
