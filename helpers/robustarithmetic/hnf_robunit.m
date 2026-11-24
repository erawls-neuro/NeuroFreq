function [y, info] = hnf_robunit(x, dim, opts)
% HNF_ROBUNIT  Stable Euclidean normalization without forming large norms.
%   Y = HNF_ROBUNIT(X,DIM,OPTS); DIM defaults to first nonsingleton axis.
%   Options: minimumNorm=0 (absolute units of X),
%   zeroPolicy='nan' (default)|'zero'|'error',
%   inputPolicy='propagate' (default)|'fill'|'error', fillValue=0.
%   A finite slice whose norm is <=minimumNorm is handled by zeroPolicy.
%   NaN/Inf makes its entire slice NaN by default; unit vectors for input
%   infinities do not have a unique limit. Finite enormous or subnormal
%   vectors are normalized by dividing twice, so their norm need not be
%   representable. Complex phase is preserved, including huge components.
id = 'hnf_robunit';
if nargin < 2, dim = []; end
if nargin < 3, opts = []; end
elem_float(x,id); dim = elem_dim(x,dim,id);
o = elem_options(struct('minimumNorm',0,'zeroPolicy','nan', ...
    'inputPolicy','propagate','fillValue',0),opts,id);
o.zeroPolicy = elem_choice(o.zeroPolicy,{'nan','zero','error'},id);
t = o.minimumNorm;
if ~isnumeric(t) || ~isscalar(t) || ~isreal(t) || ~isfinite(t) || t < 0
    error([id ':Options'],'minimumNorm must be a finite nonnegative scalar.');
end
if nargout > 1, inputCount = nnz(~isfinite(x)); end
x = elem_input(x,o,id);
if size(x,dim) == 0
    y = x; smallCount = 0;
else
    [scale,root,nanSlice,infSlice] = elem_normparts(x,dim);
    finiteSlice = ~nanSlice & ~infSlice;
    small = finiteSlice & (scale == 0 | scale <= t./root);
    if strcmp(o.zeroPolicy,'error') && any(small(:))
        error([id ':ZeroNorm'],'A norm is at or below minimumNorm.');
    end
    safeScale = scale; safeScale(safeScale == 0 | ~finiteSlice) = 1;
    safeRoot = root; safeRoot(safeRoot == 0 | ~finiteSlice) = 1;
    y = (x./safeScale)./safeRoot;
    invalid = (small | ~finiteSlice) | false(size(y));
    y(invalid) = NaN;
    if strcmp(o.zeroPolicy,'zero')
        zero = small | false(size(y)); y(zero) = 0;
    end
    if nargout > 1, smallCount = nnz(small); end
end
if nargout > 1
    info = struct('method','scaled-unit','dimension',dim, ...
        'smallNormSlices',smallCount,'nonfiniteInput',inputCount, ...
        'nonfiniteOutput',nnz(~isfinite(y)));
end
end
