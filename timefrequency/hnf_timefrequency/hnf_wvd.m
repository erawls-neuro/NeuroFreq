function [distribution,freq,time] = hnf_wvd(signal,Fs)
% HNF_WVD  Unsmoothed one-sided discrete Wigner–Ville map of one analytic trace.
%   [D,f,t] = hnf_wvd(x,Fs) implements the two-input call used in the demo.
% Output axes are Hz (rows) and original sample times in seconds (columns).
% Integer-sample symmetric lags use a 2/Fs lag step. Therefore one DFT
% period spans [0,Fs/2). Real inputs are first converted to analytic signals.
% Values are signed energy per bin: sum(D,1) is abs(analytic(x)).^2.
% MATLAB wvd uses an oversampled time grid and a different normalization.
validateattributes(signal,{'numeric'},{'vector','finite','nonempty'});
validateattributes(Fs,{'numeric'},{'scalar','real','finite','positive'});
x = double(signal(:));
% Scaling before the lag FFT avoids an overflowing intermediate sum even
% when each normalized distribution bin has a representable magnitude.
scale = max(max(abs(real(x))),max(abs(imag(x))));
if scale>0, x = x./scale; end
if isreal(x) && numel(x)>1, x = hilbert(x); end
n = numel(x);
distribution = zeros(n,n);
for t = 1:n
    lag = 1:min(t-1,n-t);
    product = x(t+lag).*conj(x(t-lag));
    ac = zeros(n,1);
    ac(1) = abs(x(t))^2;
    ac(lag+1) = product;
    ac(n-lag+1) = conj(product);
    distribution(:,t) = ((real(fft(ac))/n).*scale).*scale;
end
freq = (0:n-1).'*(Fs/(2*n));
time = (0:n-1)/Fs;
end
