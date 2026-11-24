function [coefficients, freqs] = hnf_octave_amor(x, Fs)
% Positive-rotation analytic Morlet for one epoch. Return descending center
% frequencies, as expected by nf_cwt. No Octave wavelet package is needed.
N = numel(x);
if N < 4
    error('nf_cwt:ShortEpoch','The Morlet CWT requires at least four samples.');
end
motherCenter = 6; % radians per sample at unit scale
smallest = (motherCenter + sqrt(2*log(10)))/pi; % 10% amplitude at Nyquist
% For exp(-t^2/2) the energy density has standard deviation 1/sqrt(2).
% Use its analytic spread to set the long-period cutoff. MATLAB estimates
% the mother-wavelet spread numerically, so compare the returned Hz axes.
largest = N/sqrt(2);
lastVoice = max(0,floor(10*log2(largest/smallest)));
scales = smallest * 2.^((0:lastVoice)/10);
freqs = (Fs*(motherCenter./(2*pi*scales))).';

% Joining the epoch to its time reversal gives a period of 2N, with
% reflection at both half-sample boundaries and no data from adjacent trials.
spectrum = fft([x,fliplr(x)]).';
omega = pi*(0:N)'/N;
coefficients = complex(zeros(numel(scales),N,'like',x));
% A real unit cosine has half its Fourier amplitude on the positive axis;
% for complex inputs retain only their original positive-rotation amplitude.
if isreal(x), spectrum = 2*spectrum; end
spectrum(1) = 0;
% The frequency responses depend only on epoch length. Reuse a bank across
% sensors and trials if it fits in 8 MiB. Batch inverse FFTs along columns.
persistent cachedN bank
if isempty(cachedN) || cachedN ~= N
    bank = [];
    if 8*(N+1)*numel(scales) <= 8*1024^2
        bank = exp(-0.5*(omega*scales-motherCenter).^2);
    end
    cachedN = N;
end
step = max(1,min(numel(scales),floor(4*1024^2/(32*N))));
for first = 1:step:numel(scales)
    cols = first:min(numel(scales),first+step-1);
    if isempty(bank)
        response = exp(-0.5*(omega*scales(cols)-motherCenter).^2);
    else
        response = bank(:,cols);
    end
    filtered = zeros(2*N,numel(cols),'like',spectrum);
    filtered(1:N+1,:) = spectrum(1:N+1).*response;
    values = ifft(filtered,[],1);
    coefficients(cols,:) = values(1:N,:).';
end
end
