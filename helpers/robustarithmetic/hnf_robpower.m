function [y, info] = hnf_robpower(x, opts)
% HNF_ROBPOWER  Magnitude squared, or its logarithm, with safe scaling.
%   Y = HNF_ROBPOWER(X) computes |X|^2, NOT arbitrary elementwise X.^P.
%   Options: output='linear' (default)|'log'|'log10'|'db',
%   floorValue=0 (absolute power floor), overflow='error'|'saturate'|'allow'
%   (default 'error'), underflow='allow'|'floor'|'error' (default 'allow'),
%   inputPolicy='propagate'|'fill'|'error', fillValue=0.
%   Log outputs are computed before squaring, remaining finite even where
%   linear power would overflow or underflow. Zero power gives -Inf unless
%   a positive floorValue is requested. dB output is 10*log10(power).
%   A floor modifies the estimand; it is not inferred from machine epsilon.
id = 'hnf_robpower';
if nargin < 2, opts = []; end
elem_float(x,id);
o = elem_options(struct('output','linear','floorValue',0, ...
    'inputPolicy','propagate','fillValue',0,'overflow','error', ...
    'underflow','allow'),opts,id);
o.output = elem_choice(o.output,{'linear','log','log10','db'},id);
f = o.floorValue;
if ~isnumeric(f) || ~isscalar(f) || ~isreal(f) || ~isfinite(f) || f < 0
    error([id ':Options'],'floorValue must be a finite nonnegative scalar.');
end
f = cast(f,'like',real(x));
if ~isfinite(f) || (o.floorValue > 0 && f == 0)
    error([id ':Options'],'floorValue is not representable in the input class.');
end
if nargout > 1, inputCount = nnz(~isfinite(x)); end
x = elem_input(x,o,id);
finiteSource = isfinite(x);
if strcmp(o.output,'linear')
    s = max(abs(real(x)),abs(imag(x)));
    safe = s; safe(safe == 0 | ~isfinite(safe)) = 1;
    q = (real(x)./safe).^2 + (imag(x)./safe).^2;
    % Multiply the bounded factor first so two sub-half-subnormal squares
    % can still combine into a representable subnormal power.
    y = (s.*q).*s;
    y(isnan(x)) = NaN;
    y = elem_range(y,finiteSource,x ~= 0,o,id);
    if f > 0, y(isfinite(y) & y < f) = f; end
else
    s = max(abs(real(x)),abs(imag(x)));
    safe = s; safe(safe == 0 | ~isfinite(safe)) = 1;
    r = real(x)./safe; i = imag(x)./safe;
    y = 2.*log(safe) + log(r.^2 + i.^2);
    y(s == 0) = -Inf;
    y(isinf(s) & ~isnan(x)) = Inf;
    y(isnan(x)) = NaN;
    if f > 0
        below = y < log(f); y(below) = log(f);
    end
    if strcmp(o.output,'log10') || strcmp(o.output,'db')
        y = y./cast(log(10),'like',y);
        if strcmp(o.output,'db'), y = 10.*y; end
    end
end
if nargout > 1
    info = struct('output',o.output,'nonfiniteInput',inputCount, ...
        'nonfiniteOutput',nnz(~isfinite(y)),'floorValue',f);
end
end
