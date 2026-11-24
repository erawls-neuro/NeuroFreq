function [X, reciprocalCondition, usedSVD, rankA, method] = hnf_la_ls(A, B, o)
% Internal: normalized weighted LS, QR then economical SVD (no A'*A).
[m,n] = size(A);
reciprocalCondition = 0; usedSVD = true; rankA = 0; method = 'svd';
if m == 0 || n == 0
    X = zeros(n,size(B,2),'like',A); return;
end
if m >= n && o.ridge == 0
    [Q,R] = qr(A,0);
    reciprocalCondition = rcond(R);
    if reciprocalCondition > max(o.rcondThreshold,o.relativeTolerance)
        [X,~] = linsolve(R,Q'*B,struct('UT',true));
        if all(isfinite(X(:)))
            usedSVD = false; rankA = n; method = 'qr'; return;
        end
    end
    % SVD of the small R saves work on tall matrices.
    [U,S,V] = svd(R,'econ'); B = Q'*B;
else
    [U,S,V] = svd(A,'econ');
end
s = diag(S); k = min(m,n); s = s(1:k); U = U(:,1:k); V = V(:,1:k);
if isempty(s) || s(1) == 0
    X = zeros(n,size(B,2),'like',A); return;
end
keep = s > cast(o.relativeTolerance,'like',s)*s(1);
rankA = nnz(keep);
if m >= n, reciprocalCondition = double(s(end)/s(1)); end
if o.ridge > 0
    % Ridge is defined on the normalized design; it retains all modes.
    g = s./(s.^2 + cast(o.ridge,'like',s));
    X = V * (g.*(U'*B)); method = 'svd-ridge'; usedSVD = false;
else
    X = V(:,keep) * ((U(:,keep)'*B)./s(keep));
end
if any(~isfinite(X(:)))
    error('hnf:UnrepresentableOutput','Normalized solve could not produce finite coefficients.');
end
end
