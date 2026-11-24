function [y, info] = hnf_roblogsumexp(x, dim, opts)
% HNF_ROBLOGSUMEXP  Stable log(sum(exp(X),DIM)) for real logits.
%   Options: nanPolicy='propagate' (default)|'omit',
%   inputPolicy='propagate' (default)|'fill'|'error', fillValue=0.
%   Maximum subtraction avoids exponential overflow. +Inf yields +Inf,
%   all -Inf or empty/omitted slices yield -Inf; NaN dominates unless
%   explicitly omitted. Extremely small relative exponential weights may
%   round to zero, without affecting a representable log sum materially.
%   Single/double class and reduction shape are preserved.
id = 'hnf_roblogsumexp';
if nargin < 2, dim = []; end
if nargin < 3, opts = []; end
elem_float(x,id); dim = elem_dim(x,dim,id);
o = elem_options(struct('nanPolicy','propagate','inputPolicy','propagate', ...
    'fillValue',0),opts,id);
o.nanPolicy = elem_choice(o.nanPolicy,{'propagate','omit'},id);
if nargout > 1, inputCount = nnz(~isfinite(x)); end
x = elem_input(x,o,id);
if ~isreal(x), error([id ':Domain'],'Log-sum-exp requires real input.'); end
if strcmp(o.nanPolicy,'omit'), x(isnan(x)) = -Inf; end
if size(x,dim) == 0
    y = sum(x,dim); y(:) = -Inf;
else
    [m,where] = max(x,[],dim); safe = m; safe(~isfinite(safe)) = 0;
    weights = exp(x-safe);
    % Remove one known unit weight before summation; LOG1P preserves the
    % contribution of tiny weights when the largest logit is near zero.
    sz = size(x); lead = prod(sz(1:min(dim-1,numel(sz))));
    outer = (0:numel(where)-1)';
    atMax = floor(outer./lead).*(lead.*size(x,dim)) + ...
        mod(outer,lead) + (double(where(:))-1).*lead + 1;
    weights(atMax) = 0;
    y = safe + log1p(sum(weights,dim));
    y(isinf(m) & m > 0) = Inf;
    y(all(x == -Inf,dim)) = -Inf;
    y(any(isnan(x),dim)) = NaN;
end
if nargout > 1
    info = struct('method','maximum-shift','dimension',dim, ...
        'nonfiniteInput',inputCount,'nonfiniteOutput',nnz(~isfinite(y)));
end
end
