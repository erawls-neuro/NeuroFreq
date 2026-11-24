function [s, info] = hnf_robmad(x, dim, opts)
% HNF_ROBMAD  Median absolute deviation, raw by default (no toolbox needed).
%   S = HNF_ROBMAD(X,DIM,OPTS); real dense single/double; DIM=[] automatic.
%   opts.consistency=1 (raw MAD); choose 1.482602218505602 for asymptotic
%   normal-model scale. opts.nonfinite='error'|'omit'|'propagate' and
%   opts.overflow='error'|'allow'|'saturate'. See HNF_ROBSTD for degenerate cases.
if nargin < 2, dim = []; end
if nargin < 3, opts = struct; end
opts = stat_options(opts,struct('consistency',1,'nonfinite','error','overflow','error'));
opts.method = 'mad';
if nargout > 1, [s,info] = hnf_robstd(x,dim,opts);
else, s = hnf_robstd(x,dim,opts); end
end
