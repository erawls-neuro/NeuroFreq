function [P, info] = hnf_robprob(X, dim, opts)
% HNF_ROBPROB  Scale-safe normalization of nonnegative weights.
%   P = hnf_robprob(X,dim,opts) normalizes along dim (first nonsingleton).
%   opts.pseudocount = 0: additive count in the SAME UNITS as X.
%   opts.negative = 'error' | 'clip': negative weights require explicit repair.
%   opts.zero = 'error' | 'uniform' | 'nan': behavior for zero-mass slices.
%   Inputs must be finite real full CPU single/double. No implicit omission.
%   Scaling avoids sum overflow even for weights near realmax. Positive
%   probabilities below the representable range can underflow to zero.
if nargin < 2, dim = []; end
if nargin < 3, opts = []; end
hnf_aux_real(X); dim = hnf_aux_dim(X,dim);
o = hnf_aux_options(opts,struct('pseudocount',0,'negative','error','zero','error'));
if ~isscalar(o.pseudocount) || ~isnumeric(o.pseudocount) || ...
        ~isreal(o.pseudocount) || ~isfinite(o.pseudocount) || o.pseudocount < 0 || ...
        o.pseudocount > realmax(class(X)) || ...
        ~any(strcmp(o.negative,{'error','clip'})) || ~any(strcmp(o.zero,{'error','uniform','nan'}))
    error('hnf:Options','Invalid probability options.');
end
if any(~isfinite(X(:))), error('hnf:Nonfinite','Probability weights must be finite.'); end
negative = X < 0;
if any(negative(:))
    if strcmp(o.negative,'error'), error('hnf:Domain','Probability weights must be nonnegative.'); end
    X(negative) = 0;
end
if size(X,dim) == 0, error('hnf:Empty','Cannot normalize an empty probability axis.'); end
pc = cast(o.pseudocount,'like',X);
scale = max(max(X,[],dim),pc); scale(scale == 0) = 1;
W = X ./ scale + pc ./ scale;
mass = sum(W,dim); zero = mass == 0;
if any(zero(:)) && strcmp(o.zero,'error')
    error('hnf:ZeroMass','A probability slice has zero total mass.');
end
den = mass; den(zero) = 1; P = W ./ den;
if any(zero(:))
    if strcmp(o.zero,'uniform'), fill = cast(1/size(X,dim),'like',X); else, fill = cast(NaN,'like',X); end
    mask = repmat_mask(zero,size(P),dim);
    P(mask) = fill;
end
if nargout > 1
    info = struct('zeroMass',zero,'clippedCount',nnz(negative), ...
        'pseudocount',o.pseudocount,'sum',sum(P,dim));
end
end

function mask = repmat_mask(mask, target, dim)
rep = ones(1,max(numel(target),dim)); rep(dim) = size_along(target,dim);
mask = repmat(mask,rep);
end
function n = size_along(sz,dim)
if dim > numel(sz), n = 1; else, n = sz(dim); end
end
