function [nfft, freqs] = hnf_rid_grid(N, Fs, fRes)
% The full signed spectrum has nfft equally spaced bins. Keep DC through the
% last nonnegative bin, including Nyquist only when nfft is even.
if isempty(fRes)
    nfft = 2*N-2; % retain the N reported bins from DC through Nyquist
else
    validateattributes(fRes, {'numeric'}, {'scalar','real','finite','positive','<=',Fs/2});
    nfft = round(Fs/fRes);
end
freqs = (0:floor(nfft/2)) .* (Fs/nfft);
end
