function tfRes = nf_demodulation(data,Fs,freqs,lowpassF,order,plt)
% NF_DEMODULATION    Calculates time-frequency of an input dataset (1/2/3D) using complex demodulation.
%
% GENERAL
% -------
% Calculates time-frequency of an input dataset (1/2/3D) using complex
% demodulation.
%
% Expected dimensions: channel (1), time samples (2), trials/segments (3;
% optional).
%
% Multiplies data with a carrier frequency at frequencies
% 'freqs', then low-pass filters the signal with Butterworth filter
% frequency 'bF' (frequency spread) and order 'order'.
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
% 4) lowpassF: frequency of low-pass Butterworth filter (controls spread),
%   defaults to 2 Hz
% 5) order: order of butterworth filter, defaults to 3
% 6) plt: plot result? 0 or 1, defaults to 0
%

%no plot default
if nargin<6 || isempty(plt)
    plt=0;
end
if nargin<5 || isempty(order)
    order = 3;
    disp('setting filter order = 3 (default)');
end
if nargin<4 || isempty(lowpassF)
    lowpassF = 2;
    disp('setting low-pass filter at 2 Hz (default)');
end
if nargin<3 || isempty(freqs)
    freqs = 1:1:floor(Fs/2);
    disp('setting freqs to 1:1:Fs/2 (default)');
end
if nargin<2 || isempty(Fs) || isempty(data)
    error('at least a signal and sampling rate are required');
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

%make data long in time
data=reshape(data,nChan,nTimes*nTrls); %concatenate times
demPow = zeros(nChan,numel(freqs),nTimes*nTrls); %preallocate
demPhas = zeros(nChan,numel(freqs),nTimes*nTrls); %preallocate

%make long time
ti=0:1/Fs:((1/Fs)*nTimes*nTrls)-(1/Fs);

%make bank of demodulating signals
H = zeros(numel(freqs),nTimes*nTrls); %preallocate
for fi=1:numel(freqs) %frequency loop
    H(fi,:)=exp(-2*1i*pi*freqs(fi).*ti);
end

%make Butterworth low-pass filter
[b1,a1] = butter(order, lowpassF/(Fs/2), 'low');

% make bank of demodulating signals
H = zeros(numel(freqs), nTimes * nTrls);
for fi = 1:numel(freqs)
    H(fi, :) = exp(-2 * 1i * pi * freqs(fi) .* ti);
end
Hconj = conj(H);  % <-- cache once

%demodulate data
prog=1;
fprintf(1,'complex demodulation progress: %3d%%\n',prog);

%sensor loop
for elec=1:nChan
    %one sensor of data
    dataY=data(elec,:);
    %demodulate by multiplying with complex oscillation
    demD = dataY.*H;
    %filter
    fD = filtfilt(b1,a1,double(demD)')';
    %multiply by oscillation's complex conjugate
    fD = fD.*Hconj;
    %extract
    demPow(elec,:,:) = abs(fD).^2;
    demPhas(elec,:,:) = angle(fD); %get phase angle
    %update progress bar
    prog = 100 * (elec / nChan);
    fprintf(1, '\b\b\b\b%3.0f%%', prog);
end
fprintf(1,'\n');

%format
tfRes.power = squeeze(reshape(demPow,nChan,numel(freqs),nTimes,nTrls)); %power, reshape back
tfRes.phase = squeeze(reshape(demPhas,nChan,numel(freqs),nTimes,nTrls)); %phase, reshape back
tfRes.freqs=freqs;
tfRes.times=0:1/Fs:((1/Fs)*nTimes)-(1/Fs);
tfRes.nsensor=nChan;
tfRes.ntrls=nTrls;
tfRes.Fs=Fs;
tfRes.lowpassF = lowpassF;
tfRes.order = order;
tfRes.method='demodulation';
tfRes.scale = 'linear';

% plot results
if plt == 1
    nf_tfplot(tfRes);
end

end



