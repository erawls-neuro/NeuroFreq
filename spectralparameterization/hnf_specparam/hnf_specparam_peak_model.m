function model = hnf_specparam_peak_model(freqs, peakParams, peakType)
% HNF_SPECPARAM_PEAK_MODEL Evaluate all fitted peaks in one vectorized pass.

if nargin < 3, peakType = 'gaussian'; end
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

weightedBasis = hnf_peak_values(freqRow,centers,heights,widths,peakType);
model = sum(weightedBasis, 1);

end
