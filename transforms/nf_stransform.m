function tfRes = nf_stransform(data, Fs, factor, plt)
% NF_STRANSFORM    Calculates time-frequency of an input dataset (1/2/3D) using the Stockwell transform (S-transform).
%
% GENERAL
% -------
% Calculates time-frequency of an input dataset (1/2/3D) using the
% Stockwell transform (S-transform).
%
% S-transform is a middle ground between the CWT and STFT - it windows FFT
% components by frequency-domain Gaussians, and recovers time-domain data
% following windowing via inverse-FFT.
%
% Expected dimensions: channel (1), time samples (2), trials/segments (3;
% optional).
%
% OUTPUT
% ------
% tfRes: structure with fields
%   1) power
%   2) phase
%   3) frequencies
%   4) times
%   5) sampling rate
%   6) window
%
% INPUT
% -----
% 1) data: 1D(time), 2D(channelXtime) or 3D(channelXtimesample) (REQUIRED)
% 2) Fs: sampling rate of signal, in Hz (REQUIRED)
% 3) factor: Gaussian windowing factor. Usually 1, 3 for high resolution.
% 4) plt: plot result? 0 or 1, defaults to 0
%

if nargin < 4 || isempty(plt)
    plt = 0;
end

if nargin < 3 || isempty(factor)
    factor = 1;
    disp('setting factor to 1 (default)');
end

if nargin < 2 || isempty(Fs) || isempty(data)
    error('at least a signal and sampling rate are required inputs');
end

nChan  = size(data, 1);
nTimes = size(data, 2);

if ndims(data) == 3
    nTrls = size(data, 3);
else
    nTrls = 1;
end

% frequency vector from original implementation
fout = 0:floor(Fs / 2) / fix(size(data, 2) / 2):floor(Fs / 2);

stransPow  = zeros(nChan, numel(fout), nTimes, nTrls);
stransPhas = zeros(nChan, numel(fout), nTimes, nTrls);

prog = 1;
fprintf(1, 'S-transform progress: %3d%%\n', prog);

for i = 1:nChan
    
    for j = 1:nTrls
        
        dataY = squeeze(data(i, :, j));
        
        % optimized Stockwell S-transform
        tmpTF = stransf_optimized(dataY, factor);
        
        stransPow(i, :, :, j)  = abs(tmpTF).^2;
        stransPhas(i, :, :, j) = angle(tmpTF);
        
    end
    
    prog = 100 * (i / nChan);
    fprintf(1, '\b\b\b\b%3.0f%%', prog);
    
end

fprintf(1, '\n');

tfRes.power   = squeeze(reshape(stransPow,  nChan, numel(fout), nTimes, nTrls));
tfRes.phase   = squeeze(reshape(stransPhas, nChan, numel(fout), nTimes, nTrls));
tfRes.freqs   = fout;
tfRes.times   = 0:1 / Fs:((1 / Fs) * nTimes) - (1 / Fs);
tfRes.nsensor = nChan;
tfRes.ntrls   = nTrls;
tfRes.Fs      = Fs;
tfRes.factor  = factor;
tfRes.method  = 'stransform';
tfRes.scale   = 'linear';

if plt == 1
    nf_tfplot(tfRes);
end

end



function stockTF = stransf_optimized(data, factor)
%STRANSF_OPTIMIZED  Stockwell S-transform (optimized Stockwell 1996 code)
%
%   stockTF = stransf_optimized(data, factor)
%
%   data   : column or row vector (real or analytic)
%   factor : Gaussian windowing factor (usually 1, 3 for higher resolution)
%
%   stockTF: matrix of size (fix(n/2)+1) x n
%            row 1    : DC (mean)
%            rows 2.. : positive-frequency S-transform slices
if ~iscolumn(data)
    data = data(:);
end
n = length(data);
% FFT and wrap to length 2n as in original code
vector_fft = fft(data);
vector_fft = [vector_fft, vector_fft];
% output: (fix(n/2)+1) x n
stockTF = zeros(fix(n / 2) + 1, n);
% DC row: mean across time
stockTF(1, :) = mean(data) * ones(1, n);
% persistent index caches for speed across channels/trials
persistent idx1_cache
persistent idx2_cache
persistent n_cache
if isempty(n_cache) || n_cache ~= n
    idx1_cache = (0:(n - 1)).^2;
    idx2_cache = (-n:-1).^2;
    n_cache    = n;
end
idx1 = idx1_cache;
idx2 = idx2_cache;
% base constant for Gaussian exponent
baseConst = -factor * 2 * pi^2;
halfN = fix(n / 2);
for banana = 1:halfN
    % scale for this frequency bin
    scale = baseConst / (banana^2);
    % Gaussian window across time indices
    gauss = exp(scale * idx1) + exp(scale * idx2);
    % frequency-domain slice for this bin
    fftSlice = vector_fft(banana + 1:banana + n);
    % inverse FFT to recover time-domain localized component
    stockTF(banana + 1, :) = ifft(fftSlice .* gauss);
end
end






