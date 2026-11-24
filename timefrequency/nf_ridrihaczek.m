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

narginchk(2,6);
data = hnf_tf_input(data, Fs);
if nargin < 3, fRes = []; end
if nargin < 4 || isempty(cwkernel), cwkernel = 0.001; end
if nargin < 5 || isempty(makePos), makePos = 0; end
if nargin < 6 || isempty(verbose), verbose = 0; end
validateattributes(cwkernel, {'numeric'}, {'scalar','finite','real','positive'});
[nChan,nTimes,nTrls] = size(data);
[nfft,freqs] = hnf_rid_grid(nTimes,Fs,fRes);
power = zeros(nChan,numel(freqs),nTimes,nTrls);
phase = power;

% Signed DFT indices in FFT storage order (zero lag/Doppler first).
% Normalizing by Nyquist retains the original dimensionless bandwidth.
% The Choi-Williams kernel is exactly one on both axes, so on the native
% DFT grid it preserves the time and frequency marginals.
bins = ifftshift(-floor(nTimes/2):ceil(nTimes/2)-1);
normalized = 2*bins/nTimes;
kernel = exp(-((normalized.' * normalized)/cwkernel).^2);
shift = mod(bsxfun(@plus,bins.',0:nTimes-1),nTimes)+1;
useNativeFold = nfft ~= nTimes && isa(kernel,'double') && exist('hnf_rid_fold_mex','file')==3 && ...
        exist('hnf_mex_available','file')==2 && hnf_mex_available('hnf_rid_fold_mex');
for sensor = 1:nChan
    for trial = 1:nTrls
        x = reshape(double(data(sensor,:,trial)),1,[]);
        % NF's phase convention uses conj(x(t))*x(t+lag). With an open
        % kernel the native result is X(f)*conj(x(t))*exp(+i*2*pi*f*t).
        products = bsxfun(@times,x(shift),conj(x));
        ambiguity = ifft(products,[],2);
        if useNativeFold
            try
                folded = hnf_rid_fold_mex(ambiguity,kernel,double(nfft));
            catch err
                if ~hnf_mex_load_error(err), rethrow(err); end
                useNativeFold = false;
            end
        end
        if ~useNativeFold
            if nfft == nTimes
                folded = ambiguity.*kernel;
            else
                % Accumulate signed lags modulo the requested FFT length.
                % This works for coarser as well as finer frequency grids.
                ambiguity = ambiguity.*kernel;
                folded = complex(zeros(nfft,nTimes));
                for lag = 1:nTimes
                    row = mod(bins(lag),nfft)+1;
                    folded(row,:) = folded(row,:)+ambiguity(lag,:);
                end
            end
        end
        distribution = fft2(folded);
        kept = distribution(1:numel(freqs),:);
        power(sensor,:,:,trial) = real(kept);
        phase(sensor,:,:,trial) = angle(kept);
    end
end
if makePos, power = hnf_positive_floor(power); end
tfRes.power = hnf_tf_output(power);
tfRes.phase = hnf_tf_output(phase);
tfRes.freqs = freqs;
tfRes.times = (0:nTimes-1)/Fs;
tfRes.nsensor = nChan;
tfRes.ntrls = nTrls;
tfRes.Fs = Fs;
tfRes.cwkernel = cwkernel;
tfRes.makePos = makePos;
tfRes.method = 'RID-Rihaczek';
% A bilinear distribution does not supply instantaneous oscillator phase.
tfRes.phaseReference = 'rihaczek';
tfRes.phaseReferenceSeconds = [];
tfRes.scale = 'linear';
if verbose, disp('[nf_ridrihaczek]: transform complete'); end
end
