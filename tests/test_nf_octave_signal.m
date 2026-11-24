function test_nf_octave_signal()
% Standalone checks of Octave-specific transform and spectral fitting paths.
if ~exist('OCTAVE_VERSION','builtin'), return; end
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'timefrequency'),fullfile(root,'spectralparameterization'), ...
    fullfile(root,'spectralparameterization','hnf_specparam'));

% Half-sample-shifted cosines extend to an exact Fourier bin on the
% double-length reflection grid. Check every sample, including the edges.
for N = [255 256]
    Fs = 128;
    k = 31;
    toneFrequency = k*Fs/(2*N);
    phase = pi*k*((0:N-1)+0.5)/N;
    x = cos(phase);
    input = cat(3,[x; 2*x],zeros(2,N));
    tf = nf_cwt(input,Fs);
    assert(strcmp(tf.cwtBackend,'NeuroFreq analytic Morlet'));
    assert(isequal(size(tf.power),[2 numel(tf.freqs) N 2]));
    assert(all(diff(tf.freqs)>0) && max(tf.freqs)<Fs/2);
    assert(max(abs(diff(log2(tf.freqs))-.1))<1e-12);
    expectedGain = exp(-.5*(6*(toneFrequency./tf.freqs(:)-1)).^2);
    expected = expectedGain*exp(1i*phase);
    z = squeeze(sqrt(tf.power(1,:,:,1)).*exp(1i*tf.phase(1,:,:,1)));
    assert(max(abs(z(:)-expected(:)))<1e-6);
    assert(all(reshape(tf.power(:,:,:,2),[],1)==0));
    assert(max(abs(reshape(tf.power(2,:,:,1)-4*tf.power(1,:,:,1),[],1)))<1e-5);
end

% A positive complex carrier retains its amplitude on interior samples.
% Its clockwise conjugate has negligible power in the same passband.
N = 513; Fs = 128; times = (0:N-1)/Fs;
x = exp(2i*pi*16*times+.43i);
tf = nf_cwt([x; conj(x)],Fs);
[~,idx] = min(abs(tf.freqs-16));
z = sqrt(squeeze(tf.power(1,idx,:,1))).'.*exp(1i*squeeze(tf.phase(1,idx,:,1)).');
gain = exp(-.5*(6*(16/tf.freqs(idx)-1))^2);
inner = 100:400;
assert(max(abs(z(inner)-gain*x(inner)))<1e-6);
assert(max(reshape(tf.power(2,idx,inner),[],1))<1e-11);

% K nonuniform spectra: absent 'contains' or 'istable' used to stop the
% public wrapper before fitting or when aligning rejected-trial metadata.
f = [2 3 4 7 11 17 23];
logPower = 1.8-1.2*log10(f);
s = struct('freqs',f,'nsensor',1,'ntrls',3,'scale','log10', ...
    'power',reshape(repmat(logPower(:),1,3),numel(f),1,3),'trial_ids',[5 6 7]);
s.power(3,1,2) = NaN;
fit = nf_tfspecparam(s,'maxPeaks',0);
assert(isequal(fit.trial_ids,[5 7]) && fit.ntrls==2);
assert(max(abs(fit.specparam.ap_exponent(:)-1.2))<1e-12);
assert(max(abs(fit.specparam.ap_offset(:)-1.8))<1e-12);

% With optim loaded this uses its fmincon front end. Without that package
% the same check exercises core sqp, including an exactly fixed width.
f = 1:.25:40;
for type = {'gaussian','cauchy'}
    kind = type{1};
    truth = [12 .7 1.5; 30 1.2 1.5];
    y = hnf_specparam_peak_model(f,truth,kind);
    for limits = {[.4 3],[1.5 1.5]}
        start = truth+[.3 .12 0; -.35 -.09 0];
        [p,bounded] = hnf_fit_peak_guess(start,f,y,[],limits{1},true,kind);
        assert(bounded && max(abs(p(:)-truth(:)))<.01);
        assert(all(p(:,3)>=limits{1}(1)-1e-8 & p(:,3)<=limits{1}(2)+1e-8));
    end
end
fprintf('Octave analytic Morlet, bounded peak fits and TF metadata passed.\n');
end
