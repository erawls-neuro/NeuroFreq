function hnf_la_scalar(v, name, lower, upper)
% Internal: finite real scalar control validation.
if ~isnumeric(v) || ~isscalar(v) || ~isreal(v) || ~isfinite(v) ...
        || v < lower || v > upper || isa(v,'gpuArray')
    error('hnf:InvalidOption','%s must be a finite real scalar in [%g,%g].',name,lower,upper);
end
end
