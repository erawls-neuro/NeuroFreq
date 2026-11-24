function X = hnf_la_rescale(X, numerator, denominator)
% Internal: X*numerator/denominator without overflowing the scalar ratio.
[fn,en] = log2(double(numerator));
[fd,ed] = log2(double(denominator));
if numerator == denominator
    if any(~isfinite(X(:))), error('hnf:UnrepresentableOutput','Result exceeds the output floating-point range.'); end
    return;
end
% Keep the mantissa multiplier <= 1 before applying a power of two.
f = fn/fd;
if f > 1, f = f/2; en = en+1; end
if isreal(X)
    [fx,ex] = log2(X);
    X = pow2(fx.*cast(f,'like',X),ex+cast(en-ed,'like',ex));
else
    [fr,er] = log2(real(X)); [fi,ei] = log2(imag(X));
    X = complex(pow2(fr.*cast(f,'like',fr),er+cast(en-ed,'like',er)), ...
        pow2(fi.*cast(f,'like',fi),ei+cast(en-ed,'like',ei)));
end
if any(~isfinite(X(:)))
    error('hnf:UnrepresentableOutput','Result exceeds the output floating-point range.');
end
end
