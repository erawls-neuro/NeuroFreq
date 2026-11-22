function test_nf_timefrequency()
% Scientific regression checks using synthetic signals with known answers.
% Add NeuroFreq and its source folders, then call test_nf_timefrequency.
% Signal Processing Toolbox (or Octave signal) is needed for the filter
% checks. MATLAB Wavelet Toolbox and a compiled MEX are checked when present.
here = fileparts(fileparts(mfilename('fullpath')));
savePath = path;
restore = onCleanup(@() path(savePath)); %#ok<NASGU>
addpath(here, fullfile(here,'timefrequency'), fullfile(here,'timefrequency','hnf_timefrequency'), ...
    fullfile(here,'setutilities','hnf_setutilities'));

% All three linear methods must act independently on each trial. An impulse
% at the last sample of trial 1 must never appear in a silent trial 2.
Fs = 128;
x = zeros(1,128,2);
x(1,end,1) = 1;
methods = {@(s) nf_dcwt(s,Fs,[8 12],3), ...
           @(s) nf_filterhilbert(s,Fs,[8 12],4,2), ...
           @(s) nf_demodulation(s,Fs,[8 12],2,2)};
for method = methods
    multi = method{1}(x);
    isolated = method{1}(x(:,:,1));
    silent = multi.power(:,:,2);
    assert(all(silent(:)==0));
    close_to(multi.power(:,:,1), isolated.power, 1e-11);
end

% Verify a discretized Morlet coefficient against independent linear conv.
Fs = 64;
t = (0:239)/Fs;
tone = exp(2i*pi*8*t + 0.73i);
wave = nf_dcwt(tone,Fs,8,3);
half = ceil(Fs*3/8);
lag = (-half:half)/Fs;
w = exp(2i*pi*8*lag).*exp(-lag.^2/(2*(3/(2*pi*8))^2));
response = conv(tone,w,'same') / max(abs(fft(w,2^nextpow2(numel(t)+2*half))));
close_to(sqrt(wave.power).*exp(1i*wave.phase),response,1e-11);
assert(isequal(size(wave.power),[1 numel(t)]));
equalCycles = nf_dcwt(tone.',Fs,[8,10],[3,3]);
assert(isequal(equalCycles.cycles,[3,3]));
close_to(equalCycles.power(1,:),wave.power,1e-12);

% Explicit windowed sums at each STFT bin, including the DC/Nyquist gains.
% Check an odd FFT and a single frame, for which Fs is still defined.
for lengthFFT = [15 16]
    Fs = lengthFFT;
    n = 0:Fs-1;
    w = .54-.46*cos(2*pi*n/(Fs-1));
    for complexInput = [false true]
        z = 0.7+1.3*cos(2*pi*3*n/Fs+0.2)+0.2*(-1).^n;
        if complexInput, z = z+1i*cos(2*pi*5*n/Fs); end
        result = nf_stft(z,Fs,1,0,1);
        bins = (0:floor(Fs/2))';
        reference = (exp(-2i*pi*bins*n/Fs)*(z.*w).')/sum(w);
        if ~complexInput
            interior = bins > 0 & bins < lengthFFT/2;
            reference(interior) = 2*reference(interior);
        end
        close_to(sqrt(result.power).*exp(1i*result.phase),reference,1e-12);
        assert(result.Fs == 1 && result.times == .5);
    end
end
dc = nf_stft(ones(1,32),32,1,0,1);
ny = nf_stft((-1).^(0:31),32,1,0,1);
close_to([dc.power(1) ny.power(end)],[1 1],1e-13);

% Stockwell slices on both parities and noninteger physical sample rate.
% For a positive complex bin tone, its own slice is constant in time.
for N = [31 32]
    Fs = 127.5;
    x = exp(2i*pi*5*(0:N-1)/N + 0.4i);
    r = nf_stransform(x,Fs,1);
    close_to(r.freqs,(0:floor(N/2))*Fs/N,1e-12);
    close_to(r.power(6,:),ones(1,N),1e-10);
    close_to(exp(1i*r.phase(6,:)),exp(0.4i)*ones(1,N),1e-11);
    EEG = eeg(x,Fs);
    wrapped = nf_tftransform(EEG,'method','stransform','verbose',0);
    close_to(wrapped.power,r.power,1e-12);
    assert(wrapped.factor == 1);
end

% Recover physical oscillator phase from both frame- and epoch-referenced
% complex coefficients. Cropping a transform cannot re-zero its carrier.
Fs = 128; N = 128;
t = (0:N-1)/Fs;
x = exp(2i*pi*8*t + .37i);
for method = {'stft','stransform'}
    EEG = eeg(x,Fs);
    wrapped = nf_tftransform(EEG,'method',method{1},'freqs',8,'times',[-.2 .2], ...
        'window',.25,'overlap',75,'fRes',1);
    result = hnf_phase_at_sample(wrapped);
    close_to(exp(1i*result),exp(1i*(2*pi*8*(wrapped.times-EEG.xmin)+.37)),1e-10);
end
% Crop a continuous grid and a strided STFT. Latency starts at 1 for each
% epoch before it is flattened again onto the new, shorter time axis.
EEG = eeg(repmat(x,[2 1 2]),Fs);
EEG.event = struct('latency',{65,193,1}, 'epoch',{1,2,1},'type',{'a','b','before'});
EEG.epoch = struct('event',{[1 3],2},'eventtype',{{'a','before'},{'b'}});
wrapped = nf_tftransform(EEG,'method','stft','window',.25,'overlap',75, ...
    'freqs',[4 16],'times',[-.2 .2]);
ix = find(wrapped.times==0);
assert(isequal([wrapped.event.latency],[ix ix+numel(wrapped.times)]));
assert(isequal(wrapped.epoch(1).event,1) && isequal(wrapped.epoch(2).event,2));
% Transfer function at a unit-amplitude tone: reference is baseband power
% 1/4 for real cosine demodulation and analytic power 1 for filtering.
Fs = 128;
t = (0:Fs*8-1)/Fs;
s = cos(2*pi*8*t+0.1);
d = nf_demodulation(s,Fs,8,1.5,3);
h = nf_filterhilbert(s,Fs,8,4,2);
inner = (3*Fs):(5*Fs);
close_to(mean(d.power(inner)),.25,1e-5);
close_to(mean(h.power(inner)),1,.01);
close_to(mean(exp(1i*(d.phase(inner)-2*pi*8*t(inner)))), exp(.1i), 1e-4);
% A low, narrow order-5 filter in polynomial form yielded NaNs.
Fs = 1024;
t = (0:Fs*20-1)/Fs;
h = nf_filterhilbert(cos(2*pi*t),Fs,1,1,5);
assert(all(isfinite(h.power)));
close_to(mean(h.power(8*Fs:12*Fs)),1,.02);

% Both quadratic kernels: check each surface against a direct sum over time
% and signed lags, independent of the implementation's FFT convolutions.
x = [1+0.3i 2-0.4i -0.5+1.2i 0.8 1.7i 1-0.1i -2 0.35i];
for Nfft = [9 12]
    for fred = [false true]
        if fred, kernel = @hnf_binomial2; else, kernel = @hnf_born_jordan2; end
        [z,tt,ff] = kernel(x,18,Nfft,6);
        direct = direct_distribution(x,Nfft,3,fred);
        close_to(ifftshift(z,1), direct, 1e-12);
        close_to(sum(z,1),abs(x).^2,1e-12);
        close_to(tt,(0:7)/18,1e-14);
        close_to(ff,(-floor(Nfft/2):ceil(Nfft/2)-1)*(18/Nfft),1e-14);
    end
end
for name = {'nf_ridbinomial','nf_ridbornjordan','nf_ridrihaczek'}
    zero = feval(name{1},zeros(1,12),12,1,[],true);
    assert(all(zero.power(:) == 0));
    assert(isequal(zero.freqs,0:6));
    EEG = eeg(repmat(x,[2 1 2]),16);
    wrapped = nf_tftransform(EEG,'method',name{1}(4:end),'fRes',2);
    assert(isequal(size(wrapped.power),[2 5 numel(x) 2]));
    assert(wrapped.nsensor == 2 && wrapped.ntrls == 2);
end

% Rihaczek's open-kernel limit has a closed form using NF's lag sign.
for N = [7 8]
    x = x(1:min(numel(x),N));
    if numel(x) < N, x(end+1) = .35i; end
    r = nf_ridrihaczek(x,N,1,1e9);
    bins = (0:floor(N/2))';
    X = fft(x);
    expected = bsxfun(@times,X(bins+1).',conj(x)).*exp(2i*pi*bins*(0:N-1)/N);
    close_to(r.power,real(expected),1e-12);
    close_to(exp(1i*r.phase),exp(1i*angle(expected)),1e-12);
end
% Closed kernel preserves a bin-centered complex tone on the native grid.
N = 32;
r = nf_ridrihaczek(exp(2i*pi*4*(0:N-1)/N),N,1,.001);
expected = zeros(N/2+1,N); expected(5,:) = N;
close_to(r.power,expected,1e-12);

% TrialDimension and missing-trial handling in both CPM estimators.
c = [2+1i 2-1i NaN 3+2i 3-2i];
for estimator = {'original','projection'}
    m = nf_cpm(c, 'Estimator',estimator{1},'Verbosity','silent');
    assert(m.validity.validTrialCount==4 && m.validity.phaseDefined);
    close_to(m.total.powerMean,9,1e-12);
    close_to(m.raw.phaseLockedAmplitudeSecondMoment,4,1e-12);
    close_to(m.raw.nonPhaseLockedAmplitudeSecondMoment,5,1e-12);
    close_to(m.raw.closureError,0,1e-12);
    cm = repmat(c(:),1,3);
    mm = nf_cpm(cm,'TrialDimension',1,'Estimator',estimator{1},'Verbosity','silent');
    assert(isequal(mm.total.powerMean, repmat(m.total.powerMean,3,1)));
end

% One-frequency/one-frame selections retain the correct axes with >1 sensor.
EEG = eeg(repmat(1:16,[2 1 3]),16);
r = nf_tftransform(EEG,'method','stft','window',1,'freqs',2,'times',0);
assert(isequal(size(r.power),[2 1 1 3]) && r.freqs == 2 && r.times == 0);
r = nf_tftransform(EEG,'method','wavelet','freqs',2,'cycles',[3 3]);
assert(isequal(size(r.power),[2 1 16 3]));
if exist('cwt','file') && ~exist('OCTAVE_VERSION','builtin')
    % MATLAB's licensed reference determines supported scales.
    Fs = 128; v = cos(2*pi*8*(0:255)/Fs);
    [obs,f] = cwt(single(v),'amor',Fs);
    epochs = cat(3,v,zeros(size(v)));
    r = nf_cwt(epochs,Fs);
    close_to(r.freqs,flipud(f(:)).',1e-6);
    close_to(r.power(:,:,1),flipud(abs(obs).^2),1e-5);
    assert(all(reshape(r.power(:,:,2),1,[])==0));
    fprintf('CWT vendor comparison passed.\n');
else
    fprintf('Skipped licensed MathWorks CWT reference.\n');
end
if exist('hnf_lacf2_mex','file') == 3
    x = [1+2i -4+0.5i 7 1i];
    reference = zeros(3,4);
    for t = 1:4
        valid = 1:min(3,5-t);
        reference(valid,t) = conj(x(t))*x(t+valid-1).';
    end
    close_to(hnf_lacf2_mex(x,3),reference,1e-13);
    close_to(hnf_lacf2(x,3),reference,1e-13);
    for invalid = {0,1.5,NaN,Inf,5}
        try
            hnf_lacf2_mex(x,invalid{1});
            error('NeuroFreq:TestMiss','Accepted invalid lag count');
        catch err
            assert(~strcmp(err.identifier,'NeuroFreq:TestMiss'));
        end
    end
    fprintf('Compiled MEX comparison passed.\n');
end
fprintf('All available time-frequency reference checks passed.\n');
end

function result = direct_distribution(x,nfft,w,binomial)
N = numel(x);
result = zeros(nfft,N);
for t = 1:N
    result(:,t) = abs(x(t))^2;
    for lag = 1:w-1
        smoothed = 0;
        for m = 0:lag
            j = mod(t-1-m,N)+1;
            if j+lag > N, continue, end
            if binomial, weight = nchoosek(lag,m)/2^lag; else, weight = 1/(lag+1); end
            smoothed = smoothed + weight*conj(x(j))*x(j+lag);
        end
        result(:,t) = result(:,t)+2*real(smoothed*exp(-2i*pi*(0:nfft-1)'*lag/nfft));
    end
end
result = result/nfft;
end

function EEG = eeg(x,Fs)
EEG.data = x;
EEG.nbchan = size(x,1);
EEG.srate = Fs;
EEG.xmin = -.5;
EEG.times = 1000*(EEG.xmin+(0:size(x,2)-1)/Fs);
EEG.xmax = EEG.times(end)/1000;
EEG.chanlocs = [];
end

function close_to(actual,expected,tolerance)
assert(isequal(size(actual),size(expected)), 'Reference comparison shape mismatch');
difference = max(abs(actual(:)-expected(:)));
assert(difference <= tolerance, 'Reference discrepancy %.10g exceeds %.10g',difference,tolerance);
end
