function [y, info] = hnf_robabs2(x, opts)
% HNF_ROBABS2  Explicitly named alias for HNF_ROBPOWER (magnitude squared).
% See HNF_ROBPOWER for options, including log-power output without overflow.
if nargin < 2, opts = []; end
if nargout > 1
    [y,info] = hnf_robpower(x,opts);
else
    y = hnf_robpower(x,opts);
end
end
