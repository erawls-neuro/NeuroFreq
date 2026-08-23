function model = hnf_specparam_peak_model(freqs, peakParams)
% HNF_SPECPARAM_PEAK_MODEL Evaluate all fitted peaks in one vectorized pass.

freqRow = double(freqs(:).');

if isempty(peakParams)
    model = zeros(size(freqRow));
    return
end

peakParams = double(peakParams);

if size(peakParams, 2) ~= 3
    peakParams = reshape(peakParams, [], 3);
end

centers = peakParams(:, 1);
heights = peakParams(:, 2);
widths = peakParams(:, 3);

distance = bsxfun(@minus, freqRow, centers);
scaledDistance = bsxfun(@rdivide, distance, widths);

basis = exp(-0.5 .* scaledDistance .* scaledDistance);

weightedBasis = bsxfun(@times, basis, heights);
model = sum(weightedBasis, 1);

end
