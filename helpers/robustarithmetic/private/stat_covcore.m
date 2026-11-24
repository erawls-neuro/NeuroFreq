function [g, factor, exponent, mu, n, omitted, propagated] = stat_covcore(x, opts)
% Complete observations only: Gram matrices preserve Hermitian PSD structure.
if ~any(strcmp(class(x),{'double','single'})) || ~ismatrix(x) || issparse(x)
    error('hnf:statistics:Input','X must be a dense CPU single/double observation-by-variable matrix.');
end
if ~strcmp(opts.method,'classical') && ~isreal(x)
    error('hnf:statistics:ComplexOrder','Winsor covariance requires real X.');
end
badrows = any(~isfinite(x),2); omitted = nnz(badrows);
propagated = false;
if omitted > 0
    switch opts.nonfinite
        case 'error'
            error('hnf:statistics:Nonfinite','X contains NaN or Inf.');
        case 'propagate'
            propagated = true;
        case 'omit'
            x(badrows,:) = [];
    end
end
n = size(x,1); p = size(x,2);
if propagated, n = n-omitted; end
if propagated || n == 0
    g = nan(p,p,'like',x); factor = ones(1,p,'like',real(x));
    exponent = zeros(1,p); mu = nan(1,p,'like',x); return
end
[a,valid,~,scale,origin,~,~,~,raw] = stat_prepare(x,1,'error',true);
if strcmp(opts.method,'winsor')
    raw = stat_transform(raw,valid,'winsor',opts.trim);
    [a,~,~,scale,origin] = stat_prepare(raw,1,'omit');
end
mu = hnf_robmean(raw,1);
a = bsxfun(@minus,raw,mu);
overflowDifference = isinf(a);
a = bsxfun(@rdivide,a,scale);
if any(overflowDifference(:))
    fallback = bsxfun(@minus,bsxfun(@rdivide,raw,scale),mu./scale);
    a(overflowDifference) = fallback(overflowDifference);
end
peak = max(abs(a),[],1); peak(peak == 0) = 1;
a = bsxfun(@rdivide,a,peak);
if opts.normalization == 0, denominator = max(n-1,1); else, denominator = n; end
g = (a'*a)/denominator;
g = g/2 + g'/2;
g(1:p+1:end) = real(diag(g));
[sf,se] = log2(scale); [pf,pe] = log2(peak);
factor = sf.*pf; exponent = double(se)+double(pe);
end
