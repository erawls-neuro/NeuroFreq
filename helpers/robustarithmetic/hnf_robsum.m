function [y, info] = hnf_robsum(x, dim, opts)
% HNF_ROBSUM  Sum with selective scaled, compensated accumulation.
%   Y = HNF_ROBSUM(X,DIM,OPTS), DIM defaults to first nonsingleton axis.
%   Options: method='auto' (default)|'compensated'|'native',
%   nanPolicy='propagate'|'omit' (default 'propagate'),
%   inputPolicy='propagate'|'fill'|'error', fillValue=0,
%   overflow='error' (default)|'saturate'|'allow'.
%   Auto starts with native SUM, then uses scaled Neumaier accumulation
%   for finite slices with nonfinite sums or strong cancellation
%   (|sum| <= sqrt(eps)*sum(|real(X)|+|imag(X)|)). Complex components are
%   compensated separately. Compensated mode always takes that path.
%   When scale separation risks losing small observations, selected finite
%   slices use exponent bins and one final rounding per real component.
%   Empty/all-omitted slices sum to zero. No outlier deletion is performed.
id = 'hnf_robsum';
if nargin < 2, dim = []; end
if nargin < 3, opts = []; end
elem_float(x,id); dim = elem_dim(x,dim,id);
o = elem_options(struct('method','auto','nanPolicy','propagate', ...
    'inputPolicy','propagate','fillValue',0,'overflow','error', ...
    'underflow','allow'),opts,id);
o.method = elem_choice(o.method,{'auto','compensated','native'},id);
o.nanPolicy = elem_choice(o.nanPolicy,{'propagate','omit'},id);
if ~strcmp(o.underflow,'allow')
    error([id ':Options'],'Sum underflow detection requires exact arithmetic and is unsupported.');
end
if nargout > 1, inputCount = nnz(~isfinite(x)); end
x = elem_input(x,o,id);
if strcmp(o.nanPolicy,'omit'), x(isnan(x)) = 0; end
y = sum(x,dim); finiteSlice = all(isfinite(x),dim);
repair = false(size(y));
if ~strcmp(o.method,'native') && size(x,dim) > 1
    if strcmp(o.method,'compensated')
        repair = finiteSlice;
    else
        bound = sum(abs(real(x)),dim) + sum(abs(imag(x)),dim);
        repair = finiteSlice & (~isfinite(y) | ...
            (bound > 0 & abs(y) <= sqrt(eps(class(x))).*bound));
    end
    if any(repair(:)), y(repair) = elem_compsum(x,dim,repair); end
end
y = elem_range(y,finiteSlice,false(size(y)),o,id);
if nargout > 1
    info = struct('method',o.method,'dimension',dim,'compensatedSlices',nnz(repair), ...
        'nonfiniteInput',inputCount,'nonfiniteOutput',nnz(~isfinite(y)));
end
end
