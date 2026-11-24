function values = hnf_cauchy(freqs, center, height, halfWidth)
% HNF_CAUCHY  Lorentzian peak of a given height and half-width at half-height.
if isscalar(center) && isscalar(height) && isscalar(halfWidth)
    values = reshape(hnf_peak_values(freqs(:).',center,height,halfWidth,'cauchy'),size(freqs));
else
    values = height./(1+((freqs-center)./halfWidth).^2);
end
end
