function err = hnf_error_expo_function(params, xs, ys, logXs)

if nargin < 4
    logXs = log(xs);
end

model = hnf_expo_function(xs, params, logXs);
residual = ys - model;
err = sum(residual .* residual);

end
