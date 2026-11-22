function values = hnf_cauchy(freqs, center, height, halfWidth)
% HNF_CAUCHY  Lorentzian peak of a given height and half-width at half-height.
values = height ./ (1 + ((freqs - center) ./ halfWidth).^2);
end
