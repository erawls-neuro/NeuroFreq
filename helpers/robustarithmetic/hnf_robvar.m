function [v, info] = hnf_robvar(x, dim, opts)
% HNF_ROBVAR  Variance or squared robust scale with explicit overflow policy.
%   V = HNF_ROBVAR(X,DIM,OPTS) accepts HNF_ROBSTD options (same defaults).
%   method='std' is classical variance (complex uses absolute squares).
%   method='mad' gives squared normal-consistent MAD, not sample variance.
%   method='winsor' gives variance after quantile clipping. An SD can be
%   representable even when its square cannot: opts.overflow controls this.
if nargin < 2, dim = []; end
if nargin < 3, opts = struct; end
opts = stat_options(opts,struct('method','std','normalization',0, ...
    'consistency',1.482602218505602,'trim',.1,'nonfinite','error','overflow','error'));
if nargout > 1, [s,info] = hnf_robstd(x,dim,opts);
else, s = hnf_robstd(x,dim,opts); end
[v,overflow] = stat_overflow(s.*s,opts.overflow);
if nargout > 1
    info.overflow = info.overflow | overflow;
    info.quantity = 'squared scale';
    if strcmp(opts.method,'std'), info.quantity = 'variance'; end
end
end
