function C = hnf_robpinv(A, tol)
% HNF_ROBPINV  Robust Moore–Penrose pseudoinverse with early‑exit fast paths.
%
%   C = hnf_robpinv(A)          – automatic tolerance (pinv default)
%   C = hnf_robpinv(A, tol)     – user‑supplied tolerance
%
% Fast paths
%   • If ‖A‖ is finite and  cond(A) ≤ 1e8 use inv(A) (square) 
%   • If rectangular but well‑conditioned use A\I or I/A (backslash)
%   • Else fall back to PINV or ridge‑regularised SVD (always finite).
%
% The threshold 1e8 is a conservative rule‑of‑thumb; change it if you
% routinely work with double‑precision matrices that have larger dynamic
% ranges.
%
if nargin < 2, tol = []; end
A(~isfinite(A)) = 0;                     % zero‑out NaN / Inf
[m,n]    = size(A);
isSquare = (m == n);
try
    % condition number (use 2‑norm estimate, cheap for small / sparse A)
    kappa = condest(A);
    if kappa < 1e8   % **well conditioned**
        if isSquare
            % plain inverse is cheapest & most accurate in this regime
            C = inv(A);
        else
            % choose the cheaper direction:  A\I  or  I/A
            if m > n                     % tall  (full column rank likely)
                C = A \ eye(m);
            else                         % wide (full row rank likely)
                C = eye(n) / A;
            end
        end
        return
    end
catch
    % condest can fail for rank‑deficient sparse A → fall through
end
try
    if isempty(tol)
        C = pinv(A);  % MATLAB’s built‑in adaptive tolerance
    else
        C = pinv(A, tol);
    end
    if all(isfinite(C(:))), return, end
catch
    % fall through to SVD rescue
end
if isempty(tol)
    tol = max(size(A))*eps(norm(A,2));   % default numeric tolerance
end
lambda = 1e-6 * trace(abs(A)) / max(size(A));   % tiny ridge
[U,S,V] = svd(A + lambda*eye(m,n),'econ');
s      = diag(S);
sInv   = zeros(size(s));
sInv(s > tol) = 1 ./ s(s > tol);
C = V * diag(sInv) * U';
end





