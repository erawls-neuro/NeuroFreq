function [q, info] = hnf_robquantile(x, p, dim, opts)
% HNF_ROBQUANTILE  Safe real quantiles with an explicit interpolation rule.
%   Q = HNF_ROBQUANTILE(X,P,DIM,OPTS). P is a nonempty vector in [0,1].
%   DIM=[] selects the first nonsingleton dimension. Q replaces that
%   dimension with numel(P), preserving P's supplied order.
%   Hyndman-Fan type 7: h=1+(n-1)*p, linear interpolation between order
%   statistics. This explicit rule need not match every MATLAB release's
%   default quantile interpolation. Endpoints are the minimum and maximum.
%   opts.nonfinite='error' (default), 'omit' (NaN AND Inf), 'propagate'.
%   Empty/all-omitted slices return NaN; no Statistics Toolbox is used.
if nargin < 3, dim = []; end
if nargin < 4, opts = struct; end
opts = stat_options(opts,struct('nonfinite','error'));
if ~isnumeric(p) || ~isreal(p) || ~isvector(p) || isempty(p) || ...
        any(~isfinite(p(:)) | p(:) < 0 | p(:) > 1)
    error('hnf:statistics:Probability','P must be a nonempty finite vector in [0,1].');
end
[~,v,n,~,~,shape,order,bad,raw] = stat_prepare(x,dim,opts.nonfinite);
q = stat_quantile(raw,v,double(p(:)));
if strcmp(opts.nonfinite,'propagate'), q(:,bad) = NaN; end
q = stat_restore(q,shape,order);
if nargout > 1
    info = struct('method','Hyndman-Fan type 7','probabilities',p(:), ...
        'count',stat_restore(n,shape,order),'invalid',stat_restore(bad,shape,order), ...
        'nonfinite',opts.nonfinite);
end
end
