function tfRes = nf_filterhilbert(data, Fs, freqs, fBandWidth, order, verbose)
% NF_FILTERHILBERT    Calculates time-frequency of an input dataset (1/2/3D) using filter-Hilbert.
%
% GENERAL
% -------
% Calculates time-frequency of an input dataset (1/2/3D) using 
% filter-Hilbert.
% 
% Expected dimensions: channel (1), time samples (2), trials/segments (3;
% optional). 
%
% Allows filter bandwidth to change by frequency, mimicking 
% wavelet smoothing that changes by frequency. Uses Butterworth filters of
% specified order.
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
% 4) fBandWidth: bandwidth of Butterworth filters, formatted 1 or [1 8] 
%       i.e. start-and stop-bandwidth, defaults to 1
% 5) order: order of Butterworth filters, defaults to 3
% 6) verbose: print verbose output? 0 or 1, defaults to 0
%

narginchk(2,6);
data = hnf_tf_input(data, Fs);
if nargin < 3 || isempty(freqs), freqs = 1:floor(Fs/2-1); end
if nargin < 4 || isempty(fBandWidth), fBandWidth = 1; end
if nargin < 5 || isempty(order), order = 3; end
if nargin < 6 || isempty(verbose), verbose = 0; end
validateattributes(freqs, {'numeric'}, {'real','vector','nonempty','finite','nonnegative','<=',Fs/2});
validateattributes(fBandWidth, {'numeric'}, {'real','vector','nonempty','finite','positive'});
validateattributes(order, {'numeric'}, {'real','scalar','integer','finite','positive'});
freqs = freqs(:).';
if isscalar(fBandWidth)
    widths = repmat(fBandWidth,1,numel(freqs));
elseif numel(fBandWidth) == 2
    widths = linspace(fBandWidth(1),fBandWidth(2),numel(freqs));
else
    error('nf_filterhilbert:FBandwidthBad','Use one bandwidth or a pair of endpoints.');
end
[nChan,nTimes,nTrls] = size(data);
% Separate columns prevent zero-phase filtering and Hilbert transforms from
% leaking one epoch into the next. The Hilbert transform acts along time.
traces = reshape(permute(double(data),[2 1 3]),nTimes,[]);
power = zeros(nChan,numel(freqs),nTimes,nTrls);
phase = power;
for k = 1:numel(freqs)
    nyquist = Fs/2;
    band = freqs(k)+[-0.5 0.5]*widths(k);
    if band(1) <= 0 && band(2) >= nyquist
        analytic = traces;
    else
        if band(1) <= 0
            [z,p,gain] = butter(order,band(2)/nyquist,'low');
        elseif band(2) >= nyquist
            [z,p,gain] = butter(order,band(1)/nyquist,'high');
        else
            [z,p,gain] = butter(order,band/nyquist,'bandpass');
        end
        [sos,gain] = zp2sos(z,p,gain);
        analytic = hnf_zero_phase(sos,gain,3*numel(p),traces);
    end
    if isreal(data), analytic = hilbert(analytic); end
    analytic = permute(reshape(analytic,nTimes,nChan,nTrls),[2 4 1 3]);
    [power(:,k,:,:),phase(:,k,:,:)] = hnf_powphase(analytic);
end
tfRes.power = hnf_tf_output(power);
tfRes.phase = hnf_tf_output(phase);
tfRes.freqs = freqs;
tfRes.times = (0:nTimes-1)/Fs;
tfRes.nsensor = nChan;
tfRes.ntrls = nTrls;
tfRes.Fs = Fs;
tfRes.bandwidth = widths;
tfRes.order = order;
tfRes.method = 'filter-hilbert';
tfRes.phaseReference = 'instantaneous';
tfRes.phaseReferenceSeconds = 0;
tfRes.scale = 'linear';
if verbose, disp('[nf_filterhilbert]: filter-hilbert complete'); end
end
