function [y, info] = hnf_robdiv(a, b, opts)
% HNF_ROBDIV  Safe elementwise division, including implicit expansion.
%   Y = HNF_ROBDIV(A,B,OPTS) supports real/complex full CPU float arrays.
%   Two inputs preserve NeuroFreq's original contract: invalid operands or
%   a zero denominator return zero. Finite real overflow returns signed Inf.
%   Pass an options struct (including struct()) to select the policies below.
%   Options (scalar struct):
%     denominatorPolicy : 'nan' (default), 'floor', 'native', or 'error'
%     denominatorFloor  : absolute denominator threshold, default 0
%     inputPolicy       : 'propagate' (default), 'fill', or 'error'
%     fillValue         : replacement for nonfinite input in fill mode (0)
%     overflow          : 'error' (default), 'saturate', or 'allow'
%     underflow         : 'allow' (default) or 'error'
%   Thresholds use B's units. Floor mode requires a strictly positive
%   threshold and preserves nonzero complex phase; exact zero maps to the
%   positive real floor. Native mode deliberately uses IEEE zero division.
%   NaN inputs propagate by default. Ordinary complex quotients use native
%   division; extreme components or collapsed native results trigger an
%   exponent-scaled rescue. This is ./, not a linear solve. An optional
%   second output includes the number of rescued complex quotients.
id = 'hnf_robdiv';
if nargin < 3 || isempty(opts)
    if isreal(a) && isreal(b) || ~(isfloat(a) && isfloat(b)) || issparse(a) || issparse(b) || ...
            isa(a,'gpuArray') || isa(b,'gpuArray') || ~isa(b,class(a))
        y = classic(a,b);
        if nargout>1
            info = struct('method','native-masked','rescuedComplexQuotients',0, ...
                'nonfiniteInput',nnz(~isfinite(a))+nnz(~isfinite(b)), ...
                'smallDenominator',nnz((b==0) | false(size(y))), ...
                'nonfiniteOutput',nnz(~isfinite(y)),'overflowPolicy','allow');
        end
        return
    end
    % Apply the same masks after rescuing ill-scaled finite complex quotients.
    historical = true;
    opts = struct('denominatorPolicy','native','overflow','allow');
else
    historical = false;
end
[a,b] = elem_pair(a,b,id);
o = elem_options(struct('denominatorPolicy','nan','denominatorFloor',0, ...
    'inputPolicy','propagate','fillValue',0,'overflow','error', ...
    'underflow','allow'),opts,id);
o.denominatorPolicy = elem_choice(o.denominatorPolicy,{'nan','floor','native','error'},id);
o.underflow = elem_choice(o.underflow,{'allow','error'},id);
f = o.denominatorFloor;
if ~isnumeric(f) || ~isscalar(f) || ~isreal(f) || ~isfinite(f) || f < 0
    error([id ':Options'],'denominatorFloor must be a finite nonnegative scalar.');
end
f = cast(f,'like',real(b));
if ~isfinite(f) || (strcmp(o.denominatorPolicy,'floor') && f <= 0)
    error([id ':Options'],'Floor mode requires a representable, positive denominatorFloor.');
end
if nargout > 1, inputCount = nnz(~isfinite(a)) + nnz(~isfinite(b)); end
a = elem_input(a,o,id); b = elem_input(b,o,id);
bad = isfinite(b) & abs(b) <= f;
if any(bad(:))
    if strcmp(o.denominatorPolicy,'error')
        error([id ':Denominator'],'A denominator is at or below denominatorFloor.');
    elseif strcmp(o.denominatorPolicy,'floor')
        small = b(bad); s = max(abs(real(small)),abs(imag(small)));
        zero = s == 0; s(zero) = 1;
        p = small./s; p = p./hypot(real(p),imag(p)); p(zero) = 1;
        b(bad) = p.*f;
    end
end
y = a./b;
finiteSource = isfinite(a) & isfinite(b) & b ~= 0;
rescued = 0;
if ~isreal(a) || ~isreal(b)
    % Native complex arithmetic is the hot path. Conservative component
    % bands flag possible squared-product overflow/underflow; zero
    % components themselves are harmless. Inspect the output as well.
    lo = 2.*sqrt(realmin(class(y))); hi = sqrt(realmax(class(y)))./2;
    risk = finiteSource & (elem_extreme(a,lo,hi) | elem_extreme(b,lo,hi) | ...
        ~isfinite(y) | (y == 0 & a ~= 0));
    if any(risk(:))
        % Allocate expanded operands only when the fallback is needed.
        aa = a + zeros(size(y),'like',a); bb = b + zeros(size(y),'like',b);
        y(risk) = elem_cdivfinite(aa(risk),bb(risk));
        if nargout > 1, rescued = nnz(risk); end
    end
end
if strcmp(o.denominatorPolicy,'nan')
    mask = bad | false(size(y)); y(mask) = NaN;
    finiteSource = finiteSource & ~mask;
end
y = elem_range(y,finiteSource,a ~= 0,o,id);
if historical, y(~isfinite(a) | ~isfinite(b) | b==0) = 0; end
if nargout > 1
    info = struct('method','native-with-scaled-complex-rescue', ...
        'rescuedComplexQuotients',rescued, ...
        'nonfiniteInput',inputCount,'smallDenominator',nnz(bad | false(size(y))), ...
        'nonfiniteOutput',nnz(~isfinite(y)),'overflowPolicy',o.overflow);
end
end

function risk = elem_extreme(x,lo,hi)
r = abs(real(x)); i = abs(imag(x));
risk = (r > 0 & r < lo) | r > hi | (i > 0 & i < lo) | i > hi;
end

function Z = classic(A, B)
% HNF_ROBDIV  Element-wise A./B with zero for invalid operands or zero B.
% The native loop fuses division and masking for large full arrays. MATLAB
% handles other types and broadcasting along individual array dimensions.
if max(numel(A),numel(B)) >= 65536 && (isa(A,'single') || isa(A,'double')) && ...
        isa(B,class(A)) && isreal(A) && isreal(B) && ~issparse(A) && ~issparse(B) && ...
        (isscalar(A) || isscalar(B) || isequal(size(A),size(B))) && ...
        exist('hnf_mex_available','file') == 2 && hnf_mex_available('hnf_robdiv_mex')
    try
        Z = hnf_robdiv_mex(A,B);
        return
    catch err
        if ~hnf_mex_load_error(err), rethrow(err); end
    end
end
Z = A ./ B;
invalid = ~isfinite(A) | ~isfinite(B) | B == 0;
Z(invalid) = 0;
end
