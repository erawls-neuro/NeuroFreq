function L = hnf_robchol(A, tol)
% HNF_ROBCHOL  NaN/Inf & near‑singular robust Cholesky decomposition.
%
%   L = hnf_robchol(A)          – automatic ridge λ·I
%   L = hnf_robchol(A, tol)     – user‑supplied tolerance for direct chol
%
% Strategy
%   1) zero out non‑finite entries
%   2) symmetrize A ← (A + A')/2
%   3) try chol(A,'lower') or chol(A+tol·I,'lower')
%   4) if that fails, add λ·I with λ = 1e‑8·trace(A)/n and retry
%
% Output
%   L  – lower‑triangular Cholesky factor (so that L*L' ≈ A)
    if nargin<2, tol = []; end
    % 1) scrub NaN / Inf
    A(~isfinite(A)) = 0;
    % 2) symmetrize
    A = (A + A')/2;
    % 3) fast path
    try
        if isempty(tol)
            L = chol(A, 'lower');
        else
            L = chol(A + tol*eye(size(A)), 'lower');
        end
        return
    catch
        % fall through to ridge‐regularized path
    end
    % 4) ridge regularization
    n    = size(A,1);
    lam  = 1e-8 * trace(A) / n;
    Areg = A + lam*eye(n);
    % ensure symmetry again
    Areg = (Areg + Areg')/2;
    L    = chol(Areg, 'lower');
end