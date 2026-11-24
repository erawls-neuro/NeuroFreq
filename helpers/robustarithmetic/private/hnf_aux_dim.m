function dim = hnf_aux_dim(X, dim)
if isempty(dim)
    dim = find(size(X) ~= 1, 1);
    if isempty(dim), dim = 1; end
end
if ~isnumeric(dim) || ~isreal(dim) || ~isscalar(dim) || ...
        ~isfinite(dim) || dim < 1 || dim ~= fix(dim)
    error('hnf:Dimension', 'Dimension must be a positive integer scalar.');
end
end
