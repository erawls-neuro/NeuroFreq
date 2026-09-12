function tfRes = nf_cwt(data,Fs,verbose)
% NF_CWT    Calculates time-frequency of an input dataset (1/2/3D) using matlab built-in function cwt.m.
%
% GENERAL
% -------
% Calculates time-frequency of an input dataset (1/2/3D) using matlab 
% built-in function cwt.m.
% 
% Expected dimensions: channel (1), time samples (2), trials/segments (3;
% optional). 
%
% Uses Morlet wavelets in all cases.
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
% 3) verbose: print verbose output? 0 or 1, defaults to 0
%

%defaults
if nargin<3 || isempty(verbose)
    %verbose=0;
end
if nargin<2 || isempty(Fs) || isempty(data)
    error('nf_cwt:SignalAndFsRequired',...
        'at least a signal and sampling rate are required');
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

%make data long over trials
data=reshape(data,nChan,nTimes*nTrls); %concatenate times

%test run - figure out returned frequencies for preallocation
[~,f] = cwt(single(data(1,:)),Fs,'amor'); %test frequencies
cwtPow = zeros(nChan,numel(f),nTimes*nTrls); %preallocate
cwtPhas = zeros(nChan,numel(f),nTimes*nTrls); %preallocate

%progress
disp('[nf_cwt]: beginning continuous wavelet transform');
%sensor loop
for i=1:nChan
    %one sensor of data
    dataY=data(i,:);
    %continuous wavelet
    convDat = cwt(single(dataY),Fs,'amor');
    if isreal(dataY)
        cwtPow(i,:,:) = flipud(abs(convDat).^2);
        cwtPhas(i,:,:) = flipud(angle(convDat));
    else
        cwtPow(i,:,:) = flipud(abs(convDat(:,:,1)).^2);
        cwtPhas(i,:,:) = flipud(angle(convDat(:,:,1)));
    end
end

%format
tfRes.power = squeeze(reshape(cwtPow,nChan,numel(f),nTimes,nTrls)); %power, reshape back
tfRes.phase = squeeze(reshape(cwtPhas,nChan,numel(f),nTimes,nTrls)); %phase, reshape back
tfRes.freqs=flipud(f)';
tfRes.times=0:1/Fs:((1/Fs)*nTimes)-(1/Fs);
tfRes.nsensor=nChan;
tfRes.ntrls=nTrls;
tfRes.Fs=Fs;
tfRes.method='cwt';
tfRes.scale = 'linear';

disp('[nf_cwt]: continuous wavelet transform complete');

end







