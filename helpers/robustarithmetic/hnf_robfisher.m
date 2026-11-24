function [Z, info] = hnf_robfisher(R, opts)
% HNF_ROBFISHER  Stable real Fisher atanh transform with explicit boundary policy.
%   opts.boundary='error' (default), 'clip', or 'allow' (+/-Inf at +/-1).
%   opts.margin=[]: eps(class(R)) when clipping; strictly in (0,1).
%   opts.tolerance=[]: 8*eps for roundoff overshoots; larger violations error.
%   Clipping is an estimator choice and never activates without 'clip'.
if nargin < 2, opts = []; end
hnf_aux_real(R);
o = hnf_aux_options(opts,struct('boundary','error','margin',[],'tolerance',[]));
if isempty(o.margin), o.margin = eps(class(R)); end
if isempty(o.tolerance), o.tolerance = 8*eps(class(R)); end
if ~any(strcmp(o.boundary,{'error','clip','allow'})) || ...
        ~isnumeric(o.margin) || ~isreal(o.margin) || ~isscalar(o.margin) || ...
        ~isfinite(o.margin) || o.margin<=0 || o.margin>=1 || ...
        ~isnumeric(o.tolerance) || ~isreal(o.tolerance) || ~isscalar(o.tolerance) || ...
        ~isfinite(o.tolerance) || o.tolerance<0
    error('hnf:Options','Invalid Fisher-transform options.');
end
if any(~isfinite(R(:))), error('hnf:Nonfinite','Correlations must be finite.'); end
if any(abs(R(:))>1+o.tolerance), error('hnf:Domain','Correlation outside [-1,1].'); end
boundary = abs(R)>=1;
if strcmp(o.boundary,'error') && any(boundary(:))
    error('hnf:Boundary','Finite Fisher transforms require abs(r)<1.');
end
if strcmp(o.boundary,'clip')
    cap = cast(1-o.margin,'like',R);
    if cap >= 1, error('hnf:Options','margin is too small in the input precision.'); end
    changed = abs(R)>cap; R = min(max(R,-cap),cap);
elseif strcmp(o.boundary,'allow')
    changed = abs(R)>1; R = min(max(R,-1),1);
else
    changed = false(size(R));
end
Z = 0.5 .* (log1p(R)-log1p(-R));
if nargout > 1, info = struct('clippedCount',nnz(changed),'boundaryCount',nnz(boundary)); end
end
