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

%defaults
if nargin<5 || isempty(verbose)
    verbose=0;
end
if nargin<4 || isempty(cycles)
    cycles = [3 8];
    if verbose == 1
        disp('[nf_dcwt]: setting cycles = [3 8] (default)');
    end
end
if nargin<3 || isempty(freqs) 
    freqs = 1:1:floor(Fs/2);
    if verbose == 1
        disp('[nf_dcwt]: setting freqs to 1:1:Fs/2 (default)');
    end
end
if nargin<2 || isempty(Fs) || isempty(data)
    if verbose == 1
        error('nf_dcwt:InputsRequired','at least a signal and sampling rate are required');
    end
end

%get all sizes
nChan = size(data,1);
nTimes = size(data,2);
%multiple trials?
if ndims(data)==3
    nTrls = size(data,3);
else
    nTrls = 1;
end

%concatenate times
data=reshape(data,nChan,nTimes*nTrls); 

%preallocate
wavDat = zeros(nChan,numel(freqs),nTimes*nTrls);

% parse cycles
nCycles = cycles;
if isscalar(nCycles)
    if verbose == 1
        disp(['[nf_dcwt]: using ' num2str(nCycles) ' cycles at all frequencies']);
    end
    cycles = repmat(nCycles,1,numel(freqs));
elseif length(nCycles)==2
    if verbose == 1
        disp(['[nf_dcwt]: cycles will linearly scale from ' num2str(nCycles(1)) ...
        ' to ' num2str(nCycles(2))]);
    end
    cycles = nCycles(1):((nCycles(2)-nCycles(1))/(numel(freqs)-1)):nCycles(2);
end

%for each frequency, find the length of time required
reqTimes = 2./freqs.*cycles;
[longTime] = max(reqTimes); %this ensures that we can get the requested resolution
time = -longTime/2:1/Fs:longTime/2; %time might be odd here, account for this in later calculations
half_of_wavelet_size = round((length(time)-1)/2);
if mod(numel(time),2)==0
    t_even=1;
else
    t_even=0;
end
Ly=length(data(1,:))+length(time)-1; %find sizes
Ly2=pow2(nextpow2(Ly));     % Find smallest power of 2 that is > Ly

% make family of wavelets and get their power spectra
H = zeros(numel(freqs),Ly2); %preallocate for wavelet FFT
for fi=1:numel(freqs) %frequency loop
    w = (pi*freqs(fi)*sqrt(pi))^-.5 * ...
        exp(2*1i*pi*freqs(fi).*time) .* ...
        exp(-time.^2./(2*( cycles(fi) / ...
        (2*pi*freqs(fi)))^2))/freqs(fi); %morlet wavelet equation
    tFFT = fft(w, Ly2);    % Fast Fourier transform
    H(fi,:) = tFFT./max(tFFT); %normalize power
end

%do fft-convolution
disp('[nf_dcwt]: beginning complex wavelet convolution');
%sensor loop
for elec = 1:nChan  
    dataY = squeeze(data(elec,:)); % one sensor of data
    X=fft(dataY, Ly2);              % Fast Fourier transform
    Y=X.*H;                         % multiply power spectra (frequency-domain convolution)
    convDat=ifft(Y, Ly2, 2);        % Inverse fast Fourier transform
    convDat=convDat(:,1:1:Ly);      % Take just the first N elements
    if t_even==1
        wavDat(elec,:,:) = convDat(:,half_of_wavelet_size:end-half_of_wavelet_size); %remove half of wavelet from either side
    else
        wavDat(elec,:,:) = convDat(:,half_of_wavelet_size+1:end-half_of_wavelet_size); %remove half of wavelet from either side
    end
end

%format
tfRes.power = squeeze(reshape(abs(wavDat).^2,nChan,numel(freqs),nTimes,nTrls)); %power, reshape back
tfRes.phase = squeeze(reshape(angle(wavDat),nChan,numel(freqs),nTimes,nTrls)); %phase, reshape back
tfRes.freqs=freqs;
tfRes.times=0:1/Fs:((1/Fs)*nTimes)-(1/Fs);
tfRes.nsensor=nChan;
tfRes.ntrls=nTrls;
tfRes.Fs=Fs;
tfRes.cycles = cycles;
tfRes.method='wavelet';
tfRes.scale = 'linear';

disp('[nf_dcwt]: complex wavelet convolution complete');

end









