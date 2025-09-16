function tfRes = nf_stft(data,Fs,window,overlap,fRes,plt)
% GENERAL
% -------
% Calculates time-frequency of an input dataset (1/2/3D) using short-time
% Fourier transform. For matlab >R2019a, uses matlab built-in function 
% 'stft'. Otherwise, uses the slower but otherwise mathematically 
% equivalent 'spectrogram'. If 100% overlap is specified then spectrogram 
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
% 6) plt: plot result? 0 or 1, defaults to 0
%
% ------
% E. Rawls, erawls89@gmail.com, rawls017@umn.edu. 
% July 2023
% Copyright (c) 2023 by E. Rawls.
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%
%
%
% Change Log
% ------------
% 2/10/24 ER: made compatible with analytic signals
% 4/5/24 ER: fix window normalization
% 9/16/25 ER: bug fix for complex signals (cap at Fs/2)

%defaults
if nargin<6 || isempty(plt)
    plt=0;
end
if nargin<5 || isempty(fRes)
    fRes=1;
end
if nargin<4 || isempty(overlap)
    overlap = 80;
    disp('setting overlap to 80% (default)');
end
if nargin<3 || isempty(window)
    window = .5;
    disp('setting window to 500 ms (default)');
end
if nargin<2 || isempty(Fs) || isempty(data)
    error('at least a signal and sampling rate are required inputs');
end

%get all sizes
nChan = size(data,1);

%multiple trials?
if ndims(data)==3
    nTrls = size(data,3);
else
    nTrls = 1;
end

%make window in s
window = round(window*Fs);

%check fRes/window compatibility
if round(Fs/fRes)<window
    disp('fRes too coarse - setting fRes to allow requested window (check frequency output!)');
    fRes = Fs/window;
end

%ensure we do not split samples
if mod(window,2)==1 && mod(Fs/2,2)==1
    window = window+1;
end
w = hamming(window);

%calculate overlap
o = round(length(w)*(overlap/100));
if o>=length(w)
    disp('Setting overlap to n-1 samples');
    o = length(w)-1;
end

%is signal complex?
if ~isreal(data)
    fra = 'twosided';
    scalef = 1; %power scaling - for complex we do not need to multiply
    lengf  = 0.5; %freq cutting - for complex negative frequencies need removed
else
    fra='onesided';
    scalef = 2; %power scaling - for one-sided we need to multiply by two
    lengf  = 1; %freq cutting - for real no negative frequencies
end

%find current version of matlab - if >R2019a stft, else spectrogram
if ~verLessThan('matlab','9.6') %stft was not introduced until 2019a
    %test run - figure out returned times/frequencies for preallocation
    [~,f,t] = stft(single(squeeze(data(1,:,1))),Fs,...
        'Window',w,'OverlapLength',o,'FrequencyRange',fra,...
        'FFTLength',round(Fs/fRes)); %test
    if strcmp(fra,'twosided')
        f = f(1:(length(f)*lengf)+1); %cut frequencies
    else
        f = f(1:(length(f)*lengf)); %do not cut frequencies
    end
    specPow = zeros(nChan,numel(f),numel(t),nTrls); %preallocate
    specPhas = zeros(nChan,numel(f),numel(t),nTrls); %preallocate
    %progress
    prog=1;
    fprintf(1,'stft progress: %3d%%\n',prog);
    %sensor loop
    for i = 1:nChan
        dataY=squeeze(data(i,:,:)); %one channel of data
        specDat = stft(single(dataY),Fs,...
            'Window',w,'OverlapLength',o,'FrequencyRange',fra,...
            'FFTLength',round(Fs/fRes)); %stft
        specDat = specDat( 1:length(f), :, : ); %drop negative freqs (if complex)
        specDat = (scalef*specDat)./sum(w); %normalize and scale
        specPow(i,:,:,:) = abs(specDat).^2; %square magnitude
        specPhas(i,:,:,:) = angle(specDat); %argument
        %update progress
        prog=100*(i/nChan);
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    fprintf(1,'\n');
else
    %test run - figure out returned times/frequencies for preallocation
    [~,f,t] = spectrogram(single(squeeze(data(1,:,1))),w,o,round(Fs/fRes),Fs); %test
    if strcmp(fra,'twosided')
        f = f(1:(length(f)*lengf)+1); %cut frequencies
    else
        f = f(1:(length(f)*lengf)); %do not cut frequencies
    end
    specPow = zeros(nChan,numel(f),numel(t),nTrls); %preallocate
    specPhas = zeros(nChan,numel(f),numel(t),nTrls); %preallocate
    %progress
    prog=1;
    fprintf(1,'spectrogram progress: %3d%%\n',prog);
    %sensor loop
    for i = 1:nChan
        %trial loop
        for trl = 1:nTrls
            dataY=squeeze(data(i,:,trl)); %one trial of data
            specDat = spectrogram(single(dataY),w,o,round(Fs/fRes),Fs); %spectrogram
            specDat = specDat( 1:length(f), :, : ); %drop negative freqs (if complex)
            specDat = (scalef*specDat)./sum(w); %normalize
            specPow(i,:,:,trl) = abs(specDat).^2; %square magnitude
            specPhas(i,:,:,trl) = angle(specDat); %argument
            %update progress
            prog=100*(i/nChan);
            fprintf(1,'\b\b\b\b%3.0f%%',prog);
        end
    end
    fprintf(1,'\n');
end

%package output
tfRes.power = squeeze(specPow); %power
tfRes.phase = squeeze(specPhas); %phase
tfRes.freqs=f';
tfRes.times=t';
tfRes.nsensor=nChan;
tfRes.ntrls=nTrls;
tfRes.Fs=1/mean(diff(t));
tfRes.window=window;
tfRes.overlap=overlap;
tfRes.method='stft';
tfRes.scale = 'linear';

% plot results
if plt == 1
    nf_tfplot(tfRes);
end


end


