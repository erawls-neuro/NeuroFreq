function L = hnf_robchol(A, tol)
% HNF_ROBCHOL  Lower Cholesky factor for a possibly singular covariance.
% Nonfinite entries are set to zero and the matrix is symmetrized. A
% supplied nonnegative tol adds tol*I before factorization. For singular
% positive-semidefinite input, a small scale-relative ridge permits CHOL.
% A materially indefinite covariance still raises a factorization error.
if nargin < 2 || isempty(tol), tol = 0; end
validateattributes(tol, {'numeric'}, {'scalar','real','finite','nonnegative'});
A(~isfinite(A)) = 0;
A = (A + A') / 2;
I = eye(size(A),class(A));
A = A + tol*I;
[L, failure] = chol(A,'lower');
if failure
    scale = mean(abs(diag(A)));
    ridge = max(1e-8, 10*eps(class(A))) * scale;
    if ridge == 0, ridge = eps(class(A)); end
    L = chol(A + ridge*I,'lower');
end
end
