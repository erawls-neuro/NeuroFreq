function [y, info] = hnf_roblog(x, opts)
% HNF_ROBLOG  Logarithm with explicit real-domain and floor policies.
%   Options: domainPolicy='nan'|'floor'|'complex'|'error' (default 'nan'),
%   floorValue=realmin(class(X)), base=exp(1), inputPolicy='propagate'|
%   'fill'|'error', fillValue=0, overflow='error'|'saturate'|'allow',
%   underflow='allow'|'error'. Complex X requires domainPolicy='complex'.
%   The default preserves log(0)=-Inf and makes negative real values NaN.
%   Floor mode replaces finite X below floorValue, including negatives;
%   this changes the quantity being estimated and is never automatic.
%   For log(1+small X), use MATLAB LOG1P directly to retain small offsets.
id = 'hnf_roblog';
if nargin < 2, opts = []; end
elem_float(x,id);
o = elem_options(struct('domainPolicy','nan','floorValue',realmin(class(x)), ...
    'base',exp(1),'inputPolicy','propagate','fillValue',0, ...
    'overflow','error','underflow','allow'),opts,id);
o.domainPolicy = elem_choice(o.domainPolicy,{'nan','floor','complex','error'},id);
o.underflow = elem_choice(o.underflow,{'allow','error'},id);
if ~isnumeric(o.base) || ~isscalar(o.base) || ~isreal(o.base) || ...
        ~isfinite(o.base) || o.base <= 0 || o.base == 1
    error([id ':Options'],'base must be finite, positive, and unequal to one.');
end
if ~isnumeric(o.floorValue) || ~isscalar(o.floorValue) || ~isreal(o.floorValue) || ...
        ~isfinite(o.floorValue) || o.floorValue <= 0
    error([id ':Options'],'floorValue must be finite and strictly positive.');
end
f = cast(o.floorValue,'like',real(x));
if ~isfinite(f) || f <= 0
    error([id ':Options'],'floorValue is not representable in the input class.');
end
if nargout > 1, inputCount = nnz(~isfinite(x)); end
x = elem_input(x,o,id);
if ~isreal(x) && ~strcmp(o.domainPolicy,'complex')
    error([id ':Domain'],'Complex input requires domainPolicy=''complex''.');
end
if strcmp(o.domainPolicy,'complex')
    bad = false(size(x));
    y = log(complex(x));
    s = max(abs(real(x)),abs(imag(x)));
    valid = isfinite(x) & s > 0;
    if any(valid(:))
        z = x(valid)./s(valid);
        y(valid) = complex(log(s(valid)) + log(hypot(real(z),imag(z))), ...
            atan2(imag(x(valid)),real(x(valid))));
    end
else
    bad = x < 0;
    if strcmp(o.domainPolicy,'error') && any(bad(:))
        error([id ':Domain'],'Logarithm received negative real values.');
    elseif strcmp(o.domainPolicy,'floor')
        bad = isfinite(x) & x < f; x(bad) = f;
    else
        x(bad) = NaN;
    end
    y = log(x);
end
if o.base == 10
    y = y./cast(log(10),'like',real(y));
elseif o.base == 2
    y = y./cast(log(2),'like',real(y));
elseif o.base ~= exp(1)
    y = y./cast(log(o.base),'like',real(y));
end
% A zero argument has an infinite logarithm by definition, not overflow.
y = elem_range(y,isfinite(x) & x ~= 0,x ~= 1,o,id);
if nargout > 1
    info = struct('modifiedDomain',nnz(bad),'nonfiniteInput',inputCount, ...
        'nonfiniteOutput',nnz(~isfinite(y)),'domainPolicy',o.domainPolicy);
end
end
