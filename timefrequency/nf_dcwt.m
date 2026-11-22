function tfRes = nf_dcwt(data,Fs,freqs,cycles,verbose)
% NF_DCWT    Calculates time-frequency of an input dataset (1/2/3D) using Morlet wavelets
%
% GENERAL
% -------
% Calculates time-frequency of an input dataset (1/2/3D) using discretized
% continuous wavelet transform, among the most common TF applications for
% EEG data. Wavelets are Gaussian-windowed complex sine waves (Morlet 
% wavelets), unit max-normalized in frequency domain.
%
% Expected dimensions: channel (1), time samples (2), trials/segments 
% (3; optional). 
%
% Allows either constant wavelets (e.g. 3 cycles) or wavelets that change 
% by frequency (e.g. [3 8] is 3 cycles at lowest frequency, increasing 
% linearly to 8 cycles at highest frequency).
%
% OUTPUT
% ------
% tfRes: structure with fields 
%   1) power
%   2) phase
%   3) frequencies
%   4) times
%   5) sampling rate
%   6) cycles
%
% INPUT
% -----
% 1) data: 1D(time), 2D(channelXtime) or 3D(channelXtimesample) (REQUIRED)
% 2) Fs: sampling rate of signal, in Hz (REQUIRED)
% 3) freqs: requested frequencies in Hz, defaults to 1:1:Fs/2
% 4) cycles in wavelet, can be constant (3) or changing with frequency
%       [3,8] for example, defaults to 3
% 5) verbose: print verbose output? 0 or 1, defaults to 0
%

narginchk(2,5);
data = hnf_tf_input(data, Fs);
if nargin < 3 || isempty(freqs), freqs = 1:floor(Fs/2); end
if nargin < 4 || isempty(cycles), cycles = [3 8]; end
if nargin < 5 || isempty(verbose), verbose = 0; end
validateattributes(freqs, {'numeric'}, {'real','vector','nonempty','finite','positive','<=',Fs/2});
validateattributes(cycles, {'numeric'}, {'real','vector','nonempty','finite','positive'});
freqs = freqs(:).';
cycles = cycles(:).';
if isscalar(cycles) || numel(freqs) == 1
    cycles = repmat(cycles(1), 1, numel(freqs));
elseif numel(cycles) == 2
    cycles = linspace(cycles(1), cycles(2), numel(freqs));
elseif numel(cycles) ~= numel(freqs)
    error('nf_dcwt:Cycles', 'Provide one cycle count, two endpoints, or one value per frequency.');
end
[nChan,nTimes,nTrls] = size(data);
% Symmetric, odd-length support puts time zero on an actual sample. Each
% epoch is independently zero-padded for linear (not circular) convolution.
half = ceil(Fs*max(cycles./freqs));
lag = (-half:half)/Fs;
nfft = 2^nextpow2(nTimes+2*half);
X = fft(data, nfft, 2);
power = zeros(nChan,numel(freqs),nTimes,nTrls);
phase = power;
for k = 1:numel(freqs)
    sigma = cycles(k)/(2*pi*freqs(k));
    wavelet = exp(2i*pi*freqs(k)*lag) .* exp(-0.5*(lag/sigma).^2);
    H = fft(wavelet,nfft);
    % Unit peak gain is a real scaling; dividing by a complex FFT bin
    % would rotate every estimated phase by an arbitrary constant.
    H = H/max(abs(H));
    coefficients = ifft(bsxfun(@times,X,H), [], 2);
    coefficients = coefficients(:,half+(1:nTimes),:);
    power(:,k,:,:) = reshape(abs(coefficients).^2,nChan,1,nTimes,nTrls);
    phase(:,k,:,:) = reshape(angle(coefficients),nChan,1,nTimes,nTrls);
end
tfRes.power = hnf_tf_output(power);
tfRes.phase = hnf_tf_output(phase);
tfRes.freqs = freqs;
tfRes.times = (0:nTimes-1)/Fs;
tfRes.nsensor = nChan;
tfRes.ntrls = nTrls;
tfRes.Fs = Fs;
tfRes.cycles = cycles;
tfRes.method = 'wavelet';
tfRes.phaseReference = 'instantaneous';
tfRes.phaseReferenceSeconds = 0;
tfRes.scale = 'linear';
if verbose, disp('[nf_dcwt]: complex wavelet convolution complete'); end
end
