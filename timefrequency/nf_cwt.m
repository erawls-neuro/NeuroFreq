function tfRes = nf_cwt(data,Fs,verbose)
% NF_CWT  Analytic Morlet time-frequency transform of a 1D, 2D or 3D dataset.
%
% GENERAL
% -------
% MATLAB uses Wavelet Toolbox cwt(...,'amor',Fs). Octave uses the
% NeuroFreq analytic Morlet bank, with ten voices per octave and symmetric
% extension of each epoch. Inspect freqs and cwtBackend when combining runs.
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
if exist('OCTAVE_VERSION','builtin')
    transform = @hnf_octave_amor;
    tfRes.cwtBackend = 'NeuroFreq analytic Morlet';
else
    transform = @(x,rate) cwt(x,'amor',rate);
    tfRes.cwtBackend = 'MATLAB Wavelet Toolbox';
end
for trial = 1:nTrls
    for sensor = 1:nChan
        % CWT chooses scales and reflects boundaries for an individual
        % epoch. Concatenating trials changes both the scales and the data
        % available at either edge. Page 1 is positive rotation for complex x.
        epoch = data(sensor,:,trial);
        % Ordinary EEG uses the established single-precision transform.
        % Preserve finite double range for traces outside single precision.
        maxValue = max(max(abs(real(epoch))),max(abs(imag(epoch))));
        if maxValue > realmax('single')/(2*nTimes) || (maxValue > 0 && maxValue < realmin('single'))
            epoch = double(epoch);
        else
            epoch = single(epoch);
        end
        [coefficients, f] = transform(epoch, Fs);
        if trial == 1 && sensor == 1
            freqs = flipud(f(:)).';
            power = zeros(nChan,numel(f),nTimes,nTrls);
            phase = power;
        end
        coefficients = flipud(coefficients(:,:,1));
        [power(sensor,:,:,trial),phase(sensor,:,:,trial)] = hnf_powphase(coefficients);
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
