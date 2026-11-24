function [a, valid, count, scale, origin, shape, order, invalid, raw] = stat_prepare(x, dim, nonfinite, allowComplex)
% Real arrays: translate before scaling to preserve small within-slice spreads.
if nargin < 4, allowComplex = false; end
if ~any(strcmp(class(x),{'double','single'})) || (~allowComplex && ~isreal(x)) || issparse(x)
    error('hnf:statistics:Input', 'X must be a dense CPU single/double array; ordered estimators require real X.');
end
if isempty(dim)
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end
stat_number(dim,1,flintmax,'dim',true);
shape = size(x);
shape(end+1:dim) = 1;
order = [dim, 1:dim-1, dim+1:numel(shape)];
n = shape(dim);
m = prod(shape(order(2:end)));
x = reshape(permute(x,order),n,m);
valid = isfinite(x);
invalid = any(~valid,1);
if strcmp(nonfinite,'error') && any(invalid)
    error('hnf:statistics:Nonfinite','X contains NaN or Inf.');
end
count = sum(valid,1);
x(~valid) = 0;
raw = x;
if n == 0
    scale = ones(1,m,'like',real(x)); origin = zeros(1,m,'like',x); a = x;
    return
end
scale = max(max(abs(real(x)),[],1),max(abs(imag(x)),[],1));
scale(scale == 0) = 1;
[~,first] = max(valid,[],1);
origin = x(first + (0:m-1)*n);
a = bsxfun(@minus,x,origin);
bad = ~isfinite(a);
a = bsxfun(@rdivide,a,scale);
if any(bad(:))
    fallback = bsxfun(@minus,bsxfun(@rdivide,x,scale),origin./scale);
    a(bad) = fallback(bad);
end
a(~valid) = 0;
end
