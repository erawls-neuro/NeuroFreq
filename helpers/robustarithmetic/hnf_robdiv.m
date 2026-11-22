function Z = hnf_robdiv(A, B)
% HNF_ROBDIV  Element-wise A./B with zero for invalid operands or zero B.
% Ordinary division supplies the correct scalar/array expansion and class.
Z = A ./ B;
invalid = ~isfinite(A) | ~isfinite(B) | B == 0;
Z(invalid) = 0;
end
