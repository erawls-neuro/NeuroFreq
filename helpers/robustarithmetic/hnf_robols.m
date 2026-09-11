function [B, reciprocalCondition, usedPseudoinverse] = hnf_robols(Z, Y)
% HNF_ROBOLS  Use QR for OLS and a tolerance-truncated pseudoinverse if needed.
relativeTolerance = 1e-9;
[Q, R] = qr(Z, 0);
projectedY = Q.' * Y;
if size(R, 1) == size(R, 2)
    reciprocalCondition = rcond(R);
else
    % More coefficients than observations makes the full vector nonunique.
    reciprocalCondition = 0;
end
if reciprocalCondition > relativeTolerance
    B = R \ projectedY;
    usedPseudoinverse = false;
else
    B = pinv(R, relativeTolerance .* norm(R, 2)) * projectedY;
    usedPseudoinverse = true;
end
end
