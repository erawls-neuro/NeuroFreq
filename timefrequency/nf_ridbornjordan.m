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
% 3) fRes: frequency spacing of output, default grid has n(times) bins
% 4) maxlags: autocorrelation window length in s, capped by the lag and frequency grid
% 5) makePos: make result positive? 0 or 1, defaults to 0
% 6) verbose: print verbose output? 0 or 1, defaults to 0
%

narginchk(2,6);
data = hnf_tf_input(data,Fs);
if nargin < 3, fRes = []; end
if nargin < 4, maxlags = []; end
if nargin < 5 || isempty(makePos), makePos = 0; end
if nargin < 6 || isempty(verbose), verbose = 0; end
[nChan,nTimes,nTrls] = size(data);
[nfft,freqs] = hnf_rid_grid(nTimes,Fs,fRes);
if isempty(maxlags)
    maxlags = 2*nTimes;
else
    validateattributes(maxlags, {'numeric'}, {'real','finite','scalar','positive'});
    maxlags = round(maxlags*Fs);
end
% Type II kernels use an even window and at most N positive lags.
maxlags = 2*floor(min([nfft,maxlags,2*nTimes])/2);
if maxlags < 2
    error('NeuroFreq:RIDLagWindow','The effective lag window must contain at least two samples.');
end
power = zeros(nChan,numel(freqs),nTimes,nTrls);
for sensor = 1:nChan
    for trial = 1:nTrls
        signed = hnf_born_jordan2(data(sensor,:,trial),Fs,nfft,maxlags);
        % Native DFT order puts DC in row 1, and, for even lengths, the
        % unique Nyquist coefficient in row nfft/2+1.
        native = ifftshift(signed,1);
        power(sensor,:,:,trial) = native(1:numel(freqs),:);
    end
end
if makePos, power = hnf_positive_floor(power); end
tfRes.power = hnf_tf_output(power);
tfRes.freqs = freqs;
tfRes.times = (0:nTimes-1)/Fs;
tfRes.nsensor = nChan;
tfRes.ntrls = nTrls;
tfRes.Fs = Fs;
tfRes.kernel = maxlags;
tfRes.makePos = makePos;
tfRes.method = 'bornjordan2';
tfRes.scale = 'linear';
if verbose, disp('[nf_ridbornjordan]: transform complete'); end
end
