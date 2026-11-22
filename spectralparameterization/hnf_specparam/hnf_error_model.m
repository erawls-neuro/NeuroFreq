function err = hnf_error_model(params, xVals, yVals, guess, guessWeight, peakType)
% HNF_ERROR_MODEL Peak-model SSE with optional initial-guess penalties.

if nargin < 6, peakType = 'gaussian'; end
params = reshape(params, [], 3);
fittedVals = hnf_specparam_peak_model(xVals, params, peakType);
residual = yVals(:).' - fittedVals;
err = sum(residual .* residual);

weakWeight = 1e2;
strongWeight = 1e7;

if strcmp(guessWeight, 'weak')
    penaltyWeight = weakWeight;
elseif strcmp(guessWeight, 'strong')
    penaltyWeight = strongWeight;
else
    penaltyWeight = 0;
end

if penaltyWeight > 0
    centerDifference = params(:, 1) - guess(:, 1);
    heightDifference = params(:, 2) - guess(:, 2);
    err = err + penaltyWeight .* sum(centerDifference .* centerDifference);
    err = err + penaltyWeight .* sum(heightDifference .* heightDifference);
end

end
