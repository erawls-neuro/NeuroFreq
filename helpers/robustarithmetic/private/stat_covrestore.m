function c = stat_covrestore(g, factor, exponent)
mantissaScale = factor.'*factor;
power = bsxfun(@plus,exponent.',exponent);
if isreal(g), c = restoreComponent(g,mantissaScale,power);
else
    c = complex(restoreComponent(real(g),mantissaScale,power), ...
        restoreComponent(imag(g),mantissaScale,power));
end
end

function y = restoreComponent(x, scale, power)
% Decompose even subnormal Gram entries before multiplying by scale.
[fraction, exponent] = log2(abs(x));
y = pow2(sign(x).*fraction.*scale,power+double(exponent));
end
