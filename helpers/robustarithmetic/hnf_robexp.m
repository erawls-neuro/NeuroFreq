function [y, info] = hnf_robexp(x, opts)
% HNF_ROBEXP  Exponential with explicit range and missing-input policies.
%   Options: overflow='error' (default)|'saturate'|'allow',
%   underflow='allow' (default)|'error'|'floor',
%   inputPolicy='propagate' (default)|'fill'|'error', fillValue=0.
%   Saturation is componentwise at realmax. Floor replaces a real positive
%   result rounded to zero by the smallest subnormal; complex floor mode
%   is disallowed because a scalar floor cannot preserve arbitrary phase.
%   Complex finite inputs use a half-exponential rescue near overflow.
%   Underflow policy concerns total zero, not small subnormal outputs.
%   For exp(X)-1 near zero use MATLAB EXPM1 directly.
id = 'hnf_robexp';
if nargin < 2, opts = []; end
elem_float(x,id);
o = elem_options(struct('inputPolicy','propagate','fillValue',0, ...
    'overflow','error','underflow','allow'),opts,id);
if nargout > 1, inputCount = nnz(~isfinite(x)); end
x = elem_input(x,o,id);
if ~isreal(x) && strcmp(o.underflow,'floor')
    error([id ':Options'],'underflow=''floor'' supports real input only.');
end
y = exp(x);
finiteSource = isfinite(x);
if ~isreal(x)
    % exp(real(X)/2) stays representable whenever either output component
    % can be represented. Multiplying the phase first avoids Inf*0.
    edge = finiteSource & real(x) > log(realmax(class(x))) - 2;
    if any(edge(:))
        r = real(x(edge)); p = imag(x(edge)); h = exp(r./2);
        yr = (h.*cos(p)).*h; yi = (h.*sin(p)).*h;
        % h can itself overflow only where at least one component must
        % overflow; preserve exactly zero trigonometric components.
        yr(cos(p) == 0) = 0; yi(sin(p) == 0) = 0;
        y(edge) = complex(yr,yi);
    end
end
y = elem_range(y,finiteSource,true(size(x)),o,id);
if nargout > 1
    info = struct('nonfiniteInput',inputCount,'nonfiniteOutput',nnz(~isfinite(y)), ...
        'zeroOutput',nnz(finiteSource & y == 0),'overflowPolicy',o.overflow);
end
end
