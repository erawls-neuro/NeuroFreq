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

[peakValues,unitBasis,scaledDistance] = hnf_peak_values(freqRow,centers,heights,widths,peakType);
fittedVals = sum(peakValues, 1);
residual = yVals - fittedVals;
err = sum(residual .* residual);

if nargout > 1
    u = scaledDistance;
    if strcmp(peakType,'cauchy')
        % u/(1+u²) and u²/(1+u²) are bounded. For the tail,
        % algebra with 1/u also covers displacement ratios above realmax.
        q = u .* unitBasis;
        w = u .* q;
        tail = abs(u)>1;
        inv = 1./u(tail);
        q(tail) = inv./(1+inv.^2);
        w(tail) = 1./(1+inv.^2);
        centerNumerator = 2*(peakValues.*q);
        widthNumerator = 2*(peakValues.*w);
    else
        % On a Gaussian tail the polynomial-weighted exponential tends
        % to zero. Avoid zero * Inf for an overflowing displacement.
        u(~isfinite(u) & unitBasis==0) = 0;
        centerNumerator = peakValues.*u;
        widthNumerator = centerNumerator.*u;
    end
    derivativeCenter = bsxfun(@rdivide,centerNumerator,widths);
    derivativeWidth = bsxfun(@rdivide,widthNumerator,widths);
    derivativeHeight = unitBasis;
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
