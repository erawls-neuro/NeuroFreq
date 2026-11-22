function C = hnf_robpinv(A, tol)
% HNF_ROBPINV  Moore-Penrose pseudoinverse after replacing NaN/Inf by zero.
% Pass tol to discard singular values at or below that absolute threshold.
% PINV's default threshold depends on the input precision and dimensions.
A(~isfinite(A)) = 0;
% MATLAB's dense SVD supplies the same minimum-norm solution for tall,
% square and wide matrices, including rank-deficient inputs.
if nargin < 2 || isempty(tol)
    % Preserve a sparse fast path only when every singular direction is
    % safely above the default relative SVD tolerance.
    if issparse(A) && size(A,1)==size(A,2) && ~isempty(A) && ...
            condest(A) < 1/max(1e-8, numel(A)*eps(class(A)))
        C = A \ speye(size(A,1));
        return
    end
    C = pinv(full(A));
else
    C = pinv(full(A), tol);
end
end
