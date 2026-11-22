function [err, gradient] = hnf_error_model_constr(params, xVals, yVals, peakType)
% HNF_ERROR_MODEL_CONSTR Peak-model SSE and exact analytic gradient.

if nargin < 4, peakType = 'gaussian'; end
yVals = yVals(:).';
if isa(xVals, 'double') && isrow(xVals)
    freqRow = xVals;
else
    freqRow = double(xVals(:).');
end
originalParameterSize = size(params);

if size(params, 2) ~= 3
    params = reshape(params, [], 3);
end

centers = params(:, 1);
heights = params(:, 2);
widths = params(:, 3);

distance = bsxfun(@minus, freqRow, centers);
scaledDistance = bsxfun(@rdivide, distance, widths);

if strcmp(peakType,'cauchy')
    unitBasis = 1 ./ (1 + scaledDistance.^2);
else
    unitBasis = exp(-0.5 .* scaledDistance .* scaledDistance);
end

peakValues = bsxfun(@times, unitBasis, heights);
fittedVals = sum(peakValues, 1);
residual = yVals - fittedVals;
err = sum(residual .* residual);

if nargout > 1
    inverseWidthSquared = 1 ./ (widths .* widths);
    derivativeCenter = bsxfun(@times, peakValues .* distance, ...
        inverseWidthSquared);
    derivativeHeight = unitBasis;
    inverseWidthCubed = 1 ./ (widths .* widths .* widths);
    derivativeWidth = bsxfun(@times, ...
        peakValues .* distance .* distance, inverseWidthCubed);
    if strcmp(peakType,'cauchy')
        % d[1/(1+u^2)] has one more factor of the Lorentzian basis.
        derivativeCenter = 2 .* unitBasis .* derivativeCenter;
        derivativeWidth = 2 .* unitBasis .* derivativeWidth;
    end
    gradientMatrix = zeros(size(params));
    gradientMatrix(:, 1) = -2 .* ...
        sum(bsxfun(@times, derivativeCenter, residual), 2);
    gradientMatrix(:, 2) = -2 .* ...
        sum(bsxfun(@times, derivativeHeight, residual), 2);
    gradientMatrix(:, 3) = -2 .* ...
        sum(bsxfun(@times, derivativeWidth, residual), 2);
    gradient = reshape(gradientMatrix, originalParameterSize);
end

end
