function dim = elem_dim(x, dim, id)
if isempty(dim)
    dim = find(size(x) ~= 1,1);
    if isempty(dim), dim = 1; end
end
if ~isnumeric(dim) || ~isscalar(dim) || ~isreal(dim) || ...
        ~isfinite(dim) || dim < 1 || dim ~= fix(dim)
    error([id ':Dimension'],'dim must be a positive integer scalar.');
end
dim = double(dim);
end
