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

narginchk(2,3);
data = hnf_tf_input(data, Fs);
if nargin < 3 || isempty(verbose), verbose = 0; end
[nChan,nTimes,nTrls] = size(data);
for trial = 1:nTrls
    for sensor = 1:nChan
        % CWT chooses scales and reflects boundaries for an individual
        % epoch. Concatenating trials changes both the scales and the data
        % available at either edge. Page 1 is positive rotation for complex x.
        [coefficients, f] = cwt(single(data(sensor,:,trial)), 'amor', Fs);
        if trial == 1 && sensor == 1
            freqs = flipud(f(:)).';
            power = zeros(nChan,numel(f),nTimes,nTrls);
            phase = power;
        end
        coefficients = flipud(coefficients(:,:,1));
        power(sensor,:,:,trial) = abs(coefficients).^2;
        phase(sensor,:,:,trial) = angle(coefficients);
    end
end
tfRes.power = hnf_tf_output(power);
tfRes.phase = hnf_tf_output(phase);
tfRes.freqs = freqs;
tfRes.times = (0:nTimes-1)/Fs;
tfRes.nsensor = nChan;
tfRes.ntrls = nTrls;
tfRes.Fs = Fs;
tfRes.method = 'cwt';
tfRes.phaseReference = 'instantaneous';
tfRes.phaseReferenceSeconds = 0;
tfRes.scale = 'linear';
if verbose, disp('[nf_cwt]: continuous wavelet transform complete'); end
end
