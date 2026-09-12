function [aperiodicParams, fullRank] = hnf_specparam_fixed_ap_fit( ...
    log10Freqs, powerSpectrum, frequencyCache)
% HNF_SPECPARAM_FIXED_AP_FIT Exact least-squares fixed aperiodic fit.
%
% The fixed model is y = offset - exponent * log10(f). Its two parameters
% therefore have an exact linear least-squares solution.

response = double(powerSpectrum(:));

fullRank = false;
aperiodicParams = [NaN NaN];

if numel(log10Freqs) < 2
    return
end

if nargin >= 3 && ~isempty(frequencyCache)
    predictor = frequencyCache.fixedPredictor;
    predictorMean = frequencyCache.fixedPredictorMean;
    predictorCentered = frequencyCache.fixedPredictorCentered;
    predictorSumSquares = frequencyCache.fixedPredictorSumSquares;
else
    predictor = -double(log10Freqs(:));
    predictorMean = mean(predictor);
    predictorCentered = predictor - predictorMean;
    predictorSumSquares = sum(predictorCentered .* predictorCentered);
end

if any(~isfinite(predictor)) || any(~isfinite(response))
    return
end

responseMean = mean(response);
responseCentered = response - responseMean;

if predictorSumSquares <= 0
    return
end

exponent = sum(predictorCentered .* responseCentered) ./ predictorSumSquares;
offset = responseMean - exponent .* predictorMean;

aperiodicParams = [offset exponent];
fullRank = true;

end
