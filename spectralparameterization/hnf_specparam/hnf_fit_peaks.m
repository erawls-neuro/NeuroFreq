function [modelParams, peakFunction, usedBounds] = hnf_fit_peaks( ...
    freqs, flatIter, maxNPeaks, peakThreshold, minPeakHeight, ...
    widthLimits, guessWeight, useOptimizationToolbox, peakType)
% HNF_FIT_PEAKS Detect candidate peaks and jointly optimize their parameters.

if nargin < 9, peakType = 'gaussian'; end
if strcmp(peakType,'cauchy')
    peakFunction = @hnf_cauchy;
    halfHeightToScale = 1;             % Lorentzian HWHM
else
    peakFunction = @hnf_gaussian;
    halfHeightToScale = sqrt(2*log(2)); % Gaussian HWHM = SD * sqrt(2*log(2))
end
usedBounds = false;
freqs = freqs(:).';
flatIter = flatIter(:).';
    
if maxNPeaks == 0
    modelParams = zeros(0, 3);
    return
end

guessParams = zeros(maxNPeaks, 3);
flatSpec = flatIter;


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
    % Use the nearest crossing on each side. Distant peaks need not
    % lie below this peak's half-height. Measure distances in Hz so the
    % same calculation works on log-spaced and linearly spaced grids.
    left = find(flatIter(1:maxIndex) <= halfHeight, 1, 'last');
    right = find(flatIter(maxIndex:end) <= halfHeight, 1, 'first') + maxIndex - 1;
    sides = [];
    if ~isempty(left)
        cross = freqs(left) + (halfHeight-flatIter(left)) * ...
            (freqs(left+1)-freqs(left)) / (flatIter(left+1)-flatIter(left));
        sides(end+1) = guessFrequency - cross;
    end
    if ~isempty(right)
        cross = freqs(right-1) + (halfHeight-flatIter(right-1)) * ...
            (freqs(right)-freqs(right-1)) / (flatIter(right)-flatIter(right-1));
        sides(end+1) = cross - guessFrequency;
    end
    if isempty(sides)
        guessWidth = widthLimits(2);
    else
        guessWidth = min(sides) / halfHeightToScale;
    end

    if guessWidth < widthLimits(1) %#ok
        guessWidth = widthLimits(1);
    end

    if guessWidth > widthLimits(2)
        guessWidth = widthLimits(2);
    end

    guessParams(guessIndex, :) = [guessFrequency guessHeight guessWidth];

    candidatePeak = peakFunction(freqs, guessFrequency, guessHeight, guessWidth);

    flatIter = flatIter - candidatePeak;
end

guessParams = guessParams(guessParams(:, 1) ~= 0, :);

if isempty(guessParams)
    modelParams = zeros(0, 3);
else
    [modelParams, usedBounds] = hnf_fit_peak_guess(guessParams, freqs, flatSpec, ...
        guessWeight, widthLimits, useOptimizationToolbox, peakType);
end

end
