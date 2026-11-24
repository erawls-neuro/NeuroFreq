function tfRes = nf_demodulation(data,Fs,freqs,lowpassF,order,verbose)
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
% 6) verbose: print verbose output? 0 or 1, defaults to 0
%

narginchk(2,6);
data = hnf_tf_input(data, Fs);
if nargin < 3 || isempty(freqs), freqs = 1:floor(Fs/2); end
if nargin < 4 || isempty(lowpassF), lowpassF = 2; end
if nargin < 5 || isempty(order), order = 3; end
if nargin < 6 || isempty(verbose), verbose = 0; end
validateattributes(freqs, {'numeric'}, {'real','vector','nonempty','finite','nonnegative','<=',Fs/2});
validateattributes(lowpassF, {'numeric'}, {'real','scalar','positive','finite','<',Fs/2});
validateattributes(order, {'numeric'}, {'real','scalar','integer','positive','finite'});
[nChan,nTimes,nTrls] = size(data);
times = (0:nTimes-1)/Fs;
freqs = freqs(:).';
[z,p,gain] = butter(order,lowpassF/(Fs/2),'low');
[sos,gain] = zp2sos(z,p,gain);
% Each column is one continuous trace. Filter in time, independently for
% every sensor/epoch, so the first sample cannot depend on the prior trial.
traces = reshape(permute(double(data),[2 1 3]),nTimes,[]);
power = zeros(nChan,numel(freqs),nTimes,nTrls);
phase = power;
for k = 1:numel(freqs)
    carrier = exp(-2i*pi*freqs(k)*times.');
    baseband = hnf_zero_phase(sos,gain,3*numel(p),bsxfun(@times,traces,carrier));
    dem = bsxfun(@times, baseband, conj(carrier));
    dem = permute(reshape(dem,nTimes,nChan,nTrls),[2 4 1 3]);
    [power(:,k,:,:),phase(:,k,:,:)] = hnf_powphase(dem);
end
tfRes.power = hnf_tf_output(power);
tfRes.phase = hnf_tf_output(phase);
tfRes.freqs = freqs;
tfRes.times = times;
tfRes.nsensor = nChan;
tfRes.ntrls = nTrls;
tfRes.Fs = Fs;
tfRes.lowpassF = lowpassF;
tfRes.order = order;
tfRes.method = 'demodulation';
tfRes.phaseReference = 'instantaneous';
tfRes.phaseReferenceSeconds = 0;
tfRes.scale = 'linear';
if verbose, disp('[nf_demodulation]: complex demodulation complete'); end
end
