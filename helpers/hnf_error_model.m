function err = hnf_error_model(params, xVals, yVals, guess, guessWeight)
% HNF_ERROR_MODEL Peak-model SSE with optional initial-guess penalties.

if size(params, 2) ~= 3
    params = reshape(params, [], 3);
end

if isa(params, 'double') == false
    params = double(params);
end

if isempty(params) == false
    if isa(xVals, 'double') && isrow(xVals)
        frequencyValues = xVals;
    else
        frequencyValues = double(xVals(:).');
    end

    if size(params, 1) == 1
        scaledDistance = ...
            (frequencyValues - params(1, 1)) ./ params(1, 3);

            fittedVals = params(1, 2) .* ...
                exp(-0.5 .* scaledDistance .* scaledDistance);
    else
        centers = params(:, 1);
        heights = params(:, 2);
        widths = params(:, 3);
        distance = bsxfun(@minus, frequencyValues, centers);
        scaledDistance = bsxfun(@rdivide, distance, widths);
        basis = exp(-0.5 .* scaledDistance .* scaledDistance);
        fittedVals = sum(bsxfun(@times, basis, heights), 1);
    end
else
    fittedVals = hnf_specparam_peak_model(xVals, params);
end

residual = yVals - fittedVals;
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
