function [c, mu, info] = hnf_robcov(x, opts)
% HNF_ROBCOV  Stable complete-row covariance, optionally regularized.
%   [C,MU,INFO] = HNF_ROBCOV(X,OPTS). Rows=observations; columns=variables.
%   Dense CPU real/complex single/double. Complex convention is (X-MU)'*
%   (X-MU)/denom, matching MATLAB COV's conjugation of the first variable.
%   opts.method='classical' (default) or 'winsor' (real X only).
%   opts.trim=.1 clips EACH marginal tail at type-7 quantiles for winsor.
%   opts.normalization=0: max(n-1,1); 1: n. Singleton covariance is zero.
%   opts.nonfinite='error'|'omit'|'propagate'. Omit removes ENTIRE rows with
%     any NaN/Inf; propagate returns an entirely NaN matrix if any are present.
%   opts.shrinkage=0 in [0,1]: (1-alpha)*C+alpha*diag(diag(C)).
%   opts.loading=0: add this nonnegative absolute value to each diagonal.
%   opts.psdrepair=false: optional eigenvalue clipping of the scaled matrix.
%   opts.eigenfloor=0 in [0,1]: floor relative to largest nonnegative
%     eigenvalue; nonzero value also enables eigenvalue projection.
%   opts.overflow='error'|'allow'|'saturate' (default error). Saturation is
%     elementwise and loses a PSD guarantee. Singular covariance is permitted.
%   Returned MU is mean of the selected/clipped observations. No pairwise
%   deletion, automatic shrinkage, or hidden ridge is applied.
if nargin < 2, opts = struct; end
opts = stat_options(opts,struct('method','classical','trim',.1, ...
    'normalization',0,'nonfinite','error','shrinkage',0,'loading',0, ...
    'psdrepair',false,'eigenfloor',0,'overflow','error'));
opts.method = stat_choice(opts.method,{'classical','winsor'},'method');
stat_number(opts.trim,0,.5-eps,'trim'); stat_number(opts.normalization,0,1,'normalization',true);
stat_number(opts.shrinkage,0,1,'shrinkage'); stat_number(opts.loading,0,Inf,'loading');
stat_number(opts.eigenfloor,0,1,'eigenfloor');
if ~(isscalar(opts.psdrepair) && (islogical(opts.psdrepair) || ...
        (isnumeric(opts.psdrepair) && any(opts.psdrepair == [0 1]))))
    error('hnf:statistics:OptionValue','psdrepair must be a logical scalar.');
end
[g,f,e,mu,n,omitted,propagated] = stat_covcore(x,opts);
p = size(g,1); projected = false;
if opts.shrinkage > 0
    diagonal = diag(g);
    g = (1-opts.shrinkage)*g;
    g(1:p+1:end) = diagonal;
end
c = stat_covrestore(g,f,e);
if n > 0 && ~propagated
    c(1:p+1:end) = real(diag(c)) + opts.loading;
end
[c,overflow] = stat_overflow(c,opts.overflow);
if (opts.psdrepair || opts.eigenfloor > 0) && ~isempty(c) && all(isfinite(c(:)))
    cscale = max(max(abs(real(c(:)))),max(abs(imag(c(:)))));
    if cscale > 0
        cn = c/cscale; cn = cn/2 + cn'/2;
        [v,d] = eig(cn); lambda = real(diag(d));
        floorvalue = opts.eigenfloor*max(max(lambda),0);
        lambda = max(lambda,floorvalue);
        cn = bsxfun(@times,v,lambda.')*v';
        cn = cn/2 + cn'/2;
        c = cn*cscale;
        c(1:p+1:end) = real(diag(c));
        [c,secondOverflow] = stat_overflow(c,opts.overflow);
        overflow = overflow | secondOverflow;
    end
    projected = true;
end
if nargout > 2
    info = struct('method',opts.method,'normalization',opts.normalization, ...
        'nonfinite',opts.nonfinite,'count',n,'invalidRowCount',omitted, ...
        'shrinkage',opts.shrinkage,'loading',opts.loading,'projected',projected, ...
        'eigenfloor',opts.eigenfloor,'overflow',overflow, ...
        'psdByConstruction',n > 0 && ~propagated && ~any(overflow(:)), ...
        'convention','conjugate first variable: centered X'' * centered X');
end
end
