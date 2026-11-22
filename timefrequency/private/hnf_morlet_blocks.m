function [power,phase] = hnf_morlet_blocks(data, frequencies, cycles, lag, half, ...
    nfft, scratchBytes, useGPU)
% Independent linear convolutions in blocks of channels, trials and frequencies.
% A full trial and the symmetric Morlet kernel always share the same FFT size.
% The CPU mode runs the same blocked algebra for reference checks.
[nChan,nTimes,nTrials] = size(data);
nFreq = numel(frequencies);
resident = isa(data,'gpuArray');
% Allow for products, complex inverse FFT, magnitude/phase, temporary
% permutations and cuFFT workspace. Each batch is bounded on every axis.
capacity = floor(scratchBytes/(96*nfft));
if capacity < 1
    error('nf_dcwt:GPUBlockTooLarge','Insufficient working memory for one Morlet FFT.');
end
if resident
    power = gpuArray.zeros(nChan,nFreq,nTimes,nTrials,'double');
    phase = gpuArray.zeros(nChan,nFreq,nTimes,nTrials,'double');
else
    power = zeros(nChan,nFreq,nTimes,nTrials);
    phase = zeros(nChan,nFreq,nTimes,nTrials);
end
freqStep = min(nFreq,max(1,floor(sqrt(capacity))));
chanStep = min(nChan,floor(capacity/freqStep));
trialStep = min(nTrials,floor(capacity/(freqStep*chanStep)));
if useGPU
    lag = gpuArray(lag);
    frequencies = gpuArray(frequencies);
    cycles = gpuArray(cycles);
end
cached = [];
% Reuse normalized wavelet responses whenever their entire bank fits.
if 16*nfft*nFreq <= min(scratchBytes/8,16*1024^2)
    for first = 1:freqStep:nFreq
        bins = first:min(first+freqStep-1,nFreq);
        response = responses(frequencies(bins),cycles(bins),lag,nfft);
        if first==1, cached = zeros(nFreq,nfft,'like',response); end
        cached(bins,:) = response;
    end
end
for ch = 1:chanStep:nChan
    chan = ch:min(nChan,ch+chanStep-1);
    for tr = 1:trialStep:nTrials
        trials = tr:min(nTrials,tr+trialStep-1);
        signal = data(chan,:,trials);
        if useGPU && ~resident, signal = gpuArray(signal); end
        X = fft(signal,nfft,2);
        for first = 1:freqStep:nFreq
            bins = first:min(first+freqStep-1,nFreq);
            if isempty(cached), H = responses(frequencies(bins),cycles(bins),lag,nfft);
            else, H = cached(bins,:); end
            % H uses the supplied frequency/cycle precision, typically
            % double. MATLAB single precision dominates mixed multiplication.
            if useGPU, xType=classUnderlying(X); hType=classUnderlying(H);
            else, xType=class(X); hType=class(H); end
            if ~strcmp(xType,hType)
                X = single(X);
                H = single(H);
            end
            weights = permute(H,[3 2 4 1]);
            coef = ifft(bsxfun(@times,X,weights),[],2);
            coef = permute(coef(:,half+(1:nTimes),:,:),[1 4 2 3]);
            p = abs(coef).^2; a = angle(coef);
            if useGPU && ~resident, p=gather(p); a=gather(a); end
            power(chan,bins,:,trials) = p;
            phase(chan,bins,:,trials) = a;
        end
    end
end
end

function H = responses(f,c,lag,nfft)
f = f(:);
sigma = c(:)./(2*pi*f);
carrier = exp(bsxfun(@times,2i*pi*f,lag));
taper = exp(-.5*bsxfun(@rdivide,lag,sigma).^2);
H = fft(carrier.*taper,nfft,2);
H = bsxfun(@rdivide,H,max(abs(H),[],2));
end
