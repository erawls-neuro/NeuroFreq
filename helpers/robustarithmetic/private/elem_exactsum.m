function total = elem_exactsum(x)
% Finite real IEEE slice with widely separated exponents. Count signed
% mantissa bits in bins, carry in base 2, then round the final significand
% once. Intermediate sums cannot overflow or erase a tiny surviving term.
% This slow path handles slices where scaling could underflow an addend.
precision = 53; minExponent = -1074; maxExponent = 1023;
if isa(x,'single'), precision = 24; minExponent = -149; maxExponent = 127; end
v = double(x(x~=0));
total = zeros(1,'like',x);
if isempty(v), return; end
[fraction, exponent] = log2(abs(v));
mantissa = uint64(fraction*2^precision);
signs = sign(v);
bins = zeros(maxExponent-minExponent+ceil(log2(numel(v)))+3,1);
for bit = 1:precision
    take = logical(bitget(mantissa,bit));
    if any(take)
        idx = exponent(take)-precision+(bit-1)-minExponent+1;
        bins = bins + accumarray(idx(:),signs(take),size(bins));
    end
end
digits = carry(bins);
resultSign = 1;
if digits(end)<0
    digits = carry(-bins);
    resultSign = -1;
end
hi = find(digits,1,'last');
if isempty(hi), return; end
lo = max(1,hi-precision+1);
sig = sum(digits(lo:hi) .* (2.^(0:hi-lo))');
if lo>1 && digits(lo-1) && (mod(sig,2) || any(digits(1:lo-2)))
    sig = sig+1;  % nearest, ties to even
end
total = cast(resultSign*pow2(sig,minExponent+lo-1),'like',x);
end

function d = carry(d)
for i=1:numel(d)-1
    c = floor(d(i)/2);
    d(i) = d(i)-2*c;
    d(i+1) = d(i+1)+c;
end
end
