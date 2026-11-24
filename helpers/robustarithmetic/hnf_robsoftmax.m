function [y, info] = hnf_robsoftmax(x, dim, opts)
% HNF_ROBSOFTMAX  Stable softmax, with explicit infinite-logit conventions.
%   Options: temperature=1 (>0),
%   infinityPolicy='nan' (default)|'equal'|'error',
%   degeneratePolicy='nan' (default)|'uniform'|'error',
%   inputPolicy='propagate' (default)|'fill'|'error', fillValue=0,
%   underflow='allow' (default)|'error'.
%   NaN makes its entire slice NaN. 'equal' explicitly assigns equal mass
%   among +Inf logits; 'uniform' explicitly assigns equal mass when every
%   logit is -Inf. Neither convention is a uniquely defined mathematical
%   limit. Finite logits use max subtraction BEFORE temperature division.
%   A temperature rescales the model and is not numerical regularization.
id = 'hnf_robsoftmax';
if nargin < 2, dim = []; end
if nargin < 3, opts = []; end
elem_float(x,id); dim = elem_dim(x,dim,id);
o = elem_options(struct('temperature',1,'infinityPolicy','nan', ...
    'degeneratePolicy','nan','inputPolicy','propagate','fillValue',0, ...
    'underflow','allow'),opts,id);
o.infinityPolicy = elem_choice(o.infinityPolicy,{'nan','equal','error'},id);
o.degeneratePolicy = elem_choice(o.degeneratePolicy,{'nan','uniform','error'},id);
o.underflow = elem_choice(o.underflow,{'allow','error'},id);
t = o.temperature;
if ~isnumeric(t) || ~isscalar(t) || ~isreal(t) || ~isfinite(t) || t <= 0
    error([id ':Options'],'temperature must be a finite strictly positive scalar.');
end
t = cast(t,'like',real(x));
if ~isfinite(t) || t <= 0
    error([id ':Options'],'temperature is not representable in the input class.');
end
if nargout > 1, inputCount = nnz(~isfinite(x)); end
x = elem_input(x,o,id);
if ~isreal(x), error([id ':Domain'],'Softmax requires real input.'); end
if size(x,dim) == 0
    y = x; posCount = 0; degenerateCount = 0;
else
    m = max(x,[],dim); nanSlice = any(isnan(x),dim);
    pos = any(x == Inf,dim) & ~nanSlice;
    degenerate = all(x == -Inf,dim);
    if strcmp(o.infinityPolicy,'error') && any(pos(:))
        error([id ':InfiniteLogit'],'A slice contains +Inf logits.');
    end
    if strcmp(o.degeneratePolicy,'error') && any(degenerate(:))
        error([id ':Degenerate'],'A slice contains only -Inf logits.');
    end
    safe = m; safe(~isfinite(safe)) = 0;
    shift = x-safe;
    rescue = isinf(shift) & isfinite(x) & isfinite(safe);
    shift = shift./t;
    if any(rescue(:))
        scaled = ((x./2-safe./2)./t).*2;
        shift(rescue) = scaled(rescue);
    end
    e = exp(shift); y = e./sum(e,dim);
    invalid = (nanSlice | pos | degenerate) | false(size(y)); y(invalid) = NaN;
    if strcmp(o.infinityPolicy,'equal') && any(pos(:))
        atInf = x == Inf; q = cast(atInf,'like',x)./sum(atInf,dim);
        mask = pos | false(size(y)); y(mask) = q(mask);
    end
    if strcmp(o.degeneratePolicy,'uniform')
        mask = degenerate | false(size(y)); y(mask) = 1./size(x,dim);
    end
    if strcmp(o.underflow,'error')
        under = isfinite(x) & ~invalid & y == 0;
        if any(under(:)), error([id ':Underflow'],'A positive probability rounded to zero.'); end
    end
    if nargout > 1, posCount = nnz(pos); degenerateCount = nnz(degenerate); end
end
if nargout > 1
    info = struct('method','maximum-shift','dimension',dim, ...
        'infiniteSlices',posCount,'degenerateSlices',degenerateCount, ...
        'nonfiniteInput',inputCount,'nonfiniteOutput',nnz(~isfinite(y)));
end
end
