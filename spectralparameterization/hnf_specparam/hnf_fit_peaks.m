function [modelParams, peakFunction] = hnf_fit_peaks( ...
    freqs, flatIter, maxNPeaks, peakThreshold, minPeakHeight, ...
    widthLimits, guessWeight, useOptimizationToolbox)
% HNF_FIT_PEAKS Detect candidate peaks and jointly optimize their parameters.

peakFunction = @gaussian;
    
if maxNPeaks == 0
    modelParams = zeros(0, 3);
    return
end

guessParams = zeros(maxNPeaks, 3);
flatSpec = flatIter;
frequencyStep = freqs(2) - freqs(1);

for guessIndex = 1:maxNPeaks
    maxHeight = max(flatIter);
    maxIndices = find(flatIter == maxHeight);

    if numel(maxIndices) > 1
        middleIndex = round(numel(maxIndices) ./ 2);
        maxIndex = maxIndices(middleIndex);
    else
        maxIndex = maxIndices;
    end

    residualStandardDeviation = std(flatIter);

    if maxHeight <= peakThreshold .* residualStandardDeviation
        break
    end

    guessFrequency = freqs(maxIndex);
    guessHeight = maxHeight;

    if guessHeight <= minPeakHeight
        break
    end

    halfHeight = 0.5 .* maxHeight;
    leftIndex = sum(flatIter(1:maxIndex) <= halfHeight);

    rightIndex = numel(flatIter) - ...
        sum(flatIter(maxIndex:end) <= halfHeight) + 1;

    shortSide = min(abs([leftIndex rightIndex] - maxIndex));
    fullWidthHalfMaximum = shortSide .* 2 .* frequencyStep;

    guessWidth = fullWidthHalfMaximum ./ (2 .* sqrt(2 .* log(2)));

    if guessWidth < widthLimits(1) %#ok
        guessWidth = widthLimits(1);
    end

    if guessWidth > widthLimits(2)
        guessWidth = widthLimits(2);
    end

    guessParams(guessIndex, :) = [guessFrequency guessHeight guessWidth];

    candidatePeak = hnf_gaussian(freqs, guessFrequency, guessHeight, guessWidth);

    flatIter = flatIter - candidatePeak;
end

guessParams = guessParams(guessParams(:, 1) ~= 0, :);

if isempty(guessParams)
    modelParams = zeros(0, 3);
else
    modelParams = hnf_fit_peak_guess(guessParams, freqs, flatSpec, ...
        guessWeight, widthLimits, useOptimizationToolbox);
end

end
