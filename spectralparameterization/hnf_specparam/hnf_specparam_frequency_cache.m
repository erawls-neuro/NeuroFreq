function cache = hnf_specparam_frequency_cache(freqs)
% HNF_SPECPARAM_FREQUENCY_CACHE Cache invariants for a repeated frequency grid.

persistent cachedFreqs
persistent cachedValues

freqRow = double(freqs(:).');

if ~isempty(cachedFreqs)
    if isequal(freqRow, cachedFreqs)
        cache = cachedValues;
        return
    end
end

cache = struct();
cache.freqs = freqRow;
cache.log10Freqs = log10(freqRow);
cache.logFreqs = log(freqRow);
cache.fixedPredictor = -cache.log10Freqs(:);
cache.fixedPredictorMean = mean(cache.fixedPredictor);
cache.fixedPredictorCentered = ...
    cache.fixedPredictor - cache.fixedPredictorMean;
cache.fixedPredictorSumSquares = sum( ...
    cache.fixedPredictorCentered .* cache.fixedPredictorCentered);

if numel(freqRow) > 1
    cache.frequencyResolution = mode(diff(freqRow));
else
    cache.frequencyResolution = NaN;
end

cachedFreqs = freqRow;
cachedValues = cache;

end
