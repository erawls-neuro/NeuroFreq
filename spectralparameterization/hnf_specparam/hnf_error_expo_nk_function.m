function err = hnf_error_expo_nk_function(params, xs, ys, log10Xs)

if nargin < 4
    log10Xs = log10(xs);
end

model = params(1) - params(2) .* log10Xs;
residual = ys - model;
err = sum(residual .* residual);

end
