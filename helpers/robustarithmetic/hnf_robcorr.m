function [r, info] = hnf_robcorr(x, opts)
% HNF_ROBCORR  Pearson or marginally winsorized correlation without cov overflow.
%   R = HNF_ROBCORR(X,OPTS); rows=observations, columns=variables.
%   Dense CPU real/complex single/double; complex classical uses X'*X.
%   opts.method='classical' (default) or 'winsor' (real X only), trim=.1.
%   opts.nonfinite='error'|'omit'|'propagate'; omission uses complete rows;
%     propagation returns an entirely NaN matrix if any input is nonfinite.
%   opts.constant='nan' (default), 'zero', 'error' for zero-variance columns.
%     'zero' makes every entry touching a constant column zero, including
%     its diagonal. All-missing/empty input remains undefined (NaN).
%   Correlation is computed from scaled centered observations directly, so
%   it stays usable when covariance or original-unit SD would overflow.
%   No pairwise missingness, significance testing, or hidden ridge.
if nargin < 2, opts = struct; end
opts = stat_options(opts,struct('method','classical','trim',.1, ...
    'nonfinite','error','constant','nan'));
opts.method = stat_choice(opts.method,{'classical','winsor'},'method');
opts.constant = stat_choice(opts.constant,{'nan','zero','error'},'constant');
stat_number(opts.trim,0,.5-eps,'trim'); opts.normalization = 1;
[g,~,~,~,n,omitted,propagated] = stat_covcore(x,opts);
p = size(g,1); sd = sqrt(max(real(diag(g)),0));
constant = sd == 0;
if strcmp(opts.constant,'error') && any(constant) && n > 0 && ~propagated
    error('hnf:statistics:ZeroScale','At least one variable has zero variance.');
end
den = sd; den(constant) = 1;
r = bsxfun(@rdivide,bsxfun(@rdivide,g,den),den.');
% Only roundoff excursions beyond the unit disk are projected back.
magnitude = abs(r); outside = magnitude > 1;
r(outside) = r(outside)./magnitude(outside);
r = r/2 + r'/2;
if n > 0 && ~propagated
    r(1:p+1:end) = 1;
    if strcmp(opts.constant,'zero')
        r(constant,:) = 0; r(:,constant) = 0;
    else
        r(constant,:) = NaN; r(:,constant) = NaN;
    end
end
if nargout > 1
    info = struct('method',opts.method,'nonfinite',opts.nonfinite, ...
        'count',n,'invalidRowCount',omitted,'constant',constant.', ...
        'roundoffClamped',nnz(outside),'convention','conjugate first variable');
end
end
