function Z = hnf_robdiv(A,B)
% HNF_ROBDIV  element‑wise A./B but returns 0 where B==0 and NaN‑safe
Z = zeros(size(A),class(A));
mask = (B ~= 0) & isfinite(A) & isfinite(B);
Z(mask) = A(mask) ./ B(mask);
% leave the rest at 0 – statistic = 0 → p‑value = 1 (non‑sig)
end





