function tfRes = nf_ridbornjordan(data,Fs,fRes,maxlags,makePos,verbose)
% NF_RIDBORNJORDAN    Calculates time-frequency of an input dataset (1/2/3D) using Cohen's class RID with Born-Jordan kernel.
%
% GENERAL
% -------
% Calculates time-frequency of an input dataset (1/2/3D) using Cohen's
% class RID with Born-Jordan kernel. Includes inline functions written
% by Jeff O'Neill for RID calculation.
%
% Expected dimensions: channel (1), time samples (2), trials/segments (3;
% optional).
%
% Born-Jordan distribution does not provide phase estimates.
%
% OUTPUT
% ------
% tfRes: structure with fields
%   1) power
%   2) frequencies
%   3) times
%   4) sampling rate
%
% INPUT
% -----
% 1) data: 1D(time), 2D(channelXtime) or 3D(channelXtimesample) (REQUIRED)
% 2) Fs: sampling rate of signal, in Hz (REQUIRED)
% 3) fRes: frequency spacing of output, defaults to 2*N(times) freqs
% 4) maxlags: autocorrelation window length in s, defaults to n(times)
% 5) makePos: make result positive? 0 or 1, defaults to 0
% 6) verbose: print verbose output? 0 or 1, defaults to 0
%

%defaults
nTimes = size(data,2);
if nargin<6 || isempty(verbose)
    verbose=0;
end
if nargin<5 || isempty(makePos)
    makePos=0;
    if verbose == 1
        disp('[rid_bornjordan]: returning both positive and negative values (default)');
    end
end
if nargin<4 || isempty(maxlags)
    maxlags=2*nTimes;
    if verbose == 1
        disp('[rid_bornjordan]: setting window length to 2*N (default)');
    end
else
    maxlags = round(maxlags*Fs);
end
if nargin<3 || isempty(fRes)
    fout = linspace(0,Fs/2,nTimes);
    if verbose == 1
        disp('[rid_bornjordan]: setting frequencies output to 2*N (default)');
    end
else
    fout = 0:fRes:Fs/2; %user-specified resolution
end
if nargin<2 || isempty(Fs) || isempty(data)
    error('rid_bornjordan:InputsRequired',...
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

%preallocate
ridPowDat = zeros(nChan,numel(fout),nTimes,nTrls); %preallocate

%make window
maxlags=min(2*numel(fout),maxlags);
if mod(maxlags,2)==1
    maxlags=maxlags-1;
end

%progress
disp('[rid_bornjordan]: beginning born-jordan rid');

%sensor loop
for eloc = 1:nChan
    %trial loop
    for trl = 1:nTrls
        %get one sensor of data
        dataY = squeeze(data(eloc,:,trl));
        %Jeff O'Neills toolbox - Born-Jordan RID
        tmpPow = born_jordan2(dataY,Fs,2*numel(fout),maxlags);
        ridPowDat(eloc,:,:,trl) = tmpPow(((length(tmpPow(:,1))/2)+1):end,:);
    end
end

if makePos==1
    %set all of surface to minimal non-zero value
    tfThresh = squeeze(ridPowDat(:));
    tfThresh(tfThresh<=0)=[];
    ridPowDat(ridPowDat<min(tfThresh)) = min(tfThresh);
end
%format
tfRes.power = squeeze(ridPowDat); %power
tfRes.freqs = fout;
tfRes.times=0:1/Fs:((1/Fs)*nTimes)-(1/Fs);
tfRes.nsensor=nChan;
tfRes.ntrls=nTrls;
tfRes.Fs=Fs;
tfRes.kernel=maxlags;
tfRes.makePos=makePos;
tfRes.method='bornjordan2';
tfRes.scale = 'linear';

disp('[rid_bornjordan]: born-jordan rid complete');

end















