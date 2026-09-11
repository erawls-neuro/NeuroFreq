function tfRes = nf_ridrihaczek(data, Fs, fRes, cwkernel, makePos, verbose)
% NF_RIDRIHACZEK    Calculates time-frequency of an input dataset (1/2/3D) using Cohen's class RID Rihaczek Distribution.
%
% GENERAL
% -------
% Calculates time-frequency of an input dataset (1/2/3D) using Cohen's
% class RID Rihaczek Distribution. Includes inline functions written
% originally by Selin Aviyente for RID-Rihaczek calculation.
%
% The real part of the Rihaczek RID corresponds to the Margenau-Hill
% distribution, which fulfills the marginals and is therefore used as
% the power estimate. The complex phase of the RID is retained.
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
%
% INPUT
% -----
% 1) data: 1D(time), 2D(channelXtime) or 3D(channelXtimesample) (REQUIRED)
% 2) Fs: sampling rate of signal, in Hz (REQUIRED)
% 3) fRes: frequency resolution of output in Hz, defaults to N(times) freqs
% 4) cwkernel: Choi-Williams parameter, defaults to 0.001
% 5) makePos: make distribution positive? 0 or 1, defaults to 0 
% 6) verbose: print verbose output? 0 or 1, defaults to 0
%

% defaults
if nargin < 6 || isempty(verbose)
    verbose = 0;
end

if nargin < 5 || isempty(makePos)
    makePos = 0;
    if verbose == 1
        disp('[nf_ridrihaczek]: returning both positive and negative values');
    end
end

if nargin < 4 || isempty(cwkernel)
    cwkernel = 0.001;
    if verbose == 1
        disp('[nf_ridrihaczek]: setting Choi-Williams kernel to 0.001 (default)');
    end
end

if nargin < 3 || isempty(fRes)
    nTimes = size(data, 2);
    fout   = linspace(0, Fs / 2, nTimes);
    if verbose == 1
        disp('[nf_ridrihaczek]: setting frequencies output to 2*N (default)');
    end
else
    fout = 0 : fRes : Fs / 2;
end

if nargin < 2 || isempty(Fs) || isempty(data)
    error('nf_ridrihaczek:InputsRequired',...
        'at least a signal and sampling rate are required inputs');
end

% sizes
nChan  = size(data, 1);
nTimes = size(data, 2);

if ndims(data) == 3
    nTrls = size(data, 3);
else
    nTrls = 1;
end

% preallocate
ridPowDat  = zeros(nChan, numel(fout), nTimes, nTrls);
ridPhasDat = zeros(nChan, numel(fout), nTimes, nTrls);

% fix odd times, RID-Rihaczek requires even length
if mod(nTimes, 2) == 1
    data    = cat(2, data, data(:, end, :));
    padFlag = 1;
else
    padFlag = 0;
end

% progress
disp('[nf_ridrihaczek]: beginning rihaczek rid');

% sensor loop
for eloc = 1:nChan
    
    % trial loop
    for trl = 1:nTrls
        
        % one sensor / trial
        dataY = squeeze(data(eloc, :, trl));
        
        % optimized Rihaczek RID with CW kernel
        RID_Rih = rid_rihaczek_optimized(dataY, 2 * numel(fout), cwkernel);
        
        % if we padded an extra sample, drop the trailing time bin
        if padFlag == 1
            RID_Rih(:, end) = [];
        end
        
        ridPowDat(eloc, :, :, trl)  = real(RID_Rih(1:numel(fout), :));
        ridPhasDat(eloc, :, :, trl) = angle(RID_Rih(1:numel(fout), :));
        
    end
    
end

times = 0 : 1 / Fs : ((1 / Fs) * nTimes) - (1 / Fs);

if makePos == 1
    tfThresh = squeeze(ridPowDat(:));
    tfThresh(tfThresh <= 0) = [];
    ridPowDat(ridPowDat < min(tfThresh)) = min(tfThresh);
end

% format output
tfRes.power    = squeeze(ridPowDat);
tfRes.phase    = squeeze(ridPhasDat);
tfRes.freqs    = fout;
tfRes.times    = times;
tfRes.nsensor  = nChan;
tfRes.ntrls    = nTrls;
tfRes.Fs       = Fs;
tfRes.cwkernel = cwkernel;
tfRes.makePos  = makePos;
tfRes.method   = 'RID-Rihaczek';
tfRes.scale    = 'linear';

disp('[nf_ridrihaczek]: rihaczek rid complete');

end



function tfd = rid_rihaczek_optimized(x, fbins, cwkernel)
% rid_rihaczek_optimized
% ----------------------
% Vectorized implementation of the reduced-interference Rihaczek
% distribution with Choi-Williams kernel in the ambiguity domain.
%
% This is mathematically equivalent to the classic implementation:
%   - builds an ambiguity-like surface via circular products
%   - applies CW kernel in (lag, Doppler)
%   - wraps Doppler axis to desired fbins
%   - applies 2D FFT to obtain RID-Rihaczek TFD

x = x(:);
tbins = length(x);

% build ambiguity surface via vectorized circular products
% amb(tau, k) = ifft_n{ conj(x(n)) * x((n+tau) mod N) } over n
N = tbins;

rowIdx = (0:(N - 1)).';
colIdx = 0:(N - 1);

shiftIdx = bsxfun(@plus, rowIdx, colIdx);
shiftIdx = mod(shiftIdx, N) + 1;

xConjRow = conj(x(:)).';
xShift   = x(shiftIdx);

prodMat = bsxfun(@times, xConjRow, xShift);

amb = ifft(prodMat, [], 2);

% center lag and Doppler axes as in the original code
ambTemp = [amb(:, tbins / 2 + 1 : tbins), amb(:, 1 : tbins / 2)];
amb1    = [ambTemp(tbins / 2 + 1 : tbins, :); ambTemp(1 : tbins / 2, :)];

% build or reuse Choi-Williams kernel
persistent K_cache
persistent tbins_cache
persistent cwkernel_cache

if isempty(K_cache)
    
    D = (-1 : 2 / (tbins - 1) : 1).';
    D = D * D.';
    L = D;
    
    K_cache        = chwi_krn(D, L, cwkernel);
    tbins_cache    = tbins;
    cwkernel_cache = cwkernel;
    
else
    
    if tbins_cache ~= tbins || cwkernel_cache ~= cwkernel
        
        D = (-1 : 2 / (tbins - 1) : 1).';
        D = D * D.';
        L = D;
        
        K_cache        = chwi_krn(D, L, cwkernel);
        tbins_cache    = tbins;
        cwkernel_cache = cwkernel;
        
    end
    
end

K = K_cache;

[s, d] = size(amb1);
df     = K(1:s, 1:d);

ambf = amb1 .* df;

% wrap Doppler axis to desired fbins
if tbins ~= fbins
    
    A = zeros(fbins, tbins);
    
    for tt = 1:tbins
        A(:, tt) = datawrap(ambf(:, tt), fbins);
    end
    
else
    
    A = ambf;
    
end

% final 2D FFT to obtain TFD
tfd = fft2(A);

end



function K = chwi_krn(D, L, A)
%CHWI_KRN Choi-Williams kernel function.
%   K = CHWI_KRN(D,L,A) returns the values K of the Choi-Williams kernel
%   function evaluated at the doppler-values in matrix D and the lag-
%   values in matrix L. Matrices D and L must have the same size. The
%   values in D should be in the range between -1 and +1 (with +1 being
%   the Nyquist frequency). The parameter A is optional and controls the
%   "diagonal bandwidth" of the kernel. Matrix K is of the same size as
%   the matrices D and L. Parameter A defaults to 10 if omitted.

if nargin < 3
    A = [];
end

if isempty(A)
    A = 10;
end

K = exp((-1 / (A * A)) * (D .* D .* L .* L));

end