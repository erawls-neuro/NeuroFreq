function [H, info] = hnf_robentropy(P, dim, opts)
% HNF_ROBENTROPY  Shannon entropy with the continuous convention 0*log(0)=0.
%   H = hnf_robentropy(P,dim,opts); real finite nonnegative probabilities.
%   opts.normalize=false requires each slice to sum to one within tolerance.
%   opts.normalize=true first calls hnf_robprob (zero mass still errors).
%   opts.base=exp(1) returns nats; base=2 returns bits. No clipping/imputation.
%   opts.tolerance=[] uses 32*eps(class(P)); tiny sum errors are NOT renormalized.
if nargin < 2, dim = []; end
if nargin < 3, opts = []; end
hnf_aux_real(P); dim = hnf_aux_dim(P,dim);
o = hnf_aux_options(opts,struct('normalize',false,'base',exp(1),'tolerance',[]));
if ~(islogical(o.normalize) && isscalar(o.normalize)) || ~isnumeric(o.base) || ...
        ~isreal(o.base) || ~isscalar(o.base) || ~isfinite(o.base) || o.base <= 1
    error('hnf:Options','normalize must be logical; entropy base must exceed one.');
end
if isempty(o.tolerance), o.tolerance = 32*eps(class(P)); end
if ~isnumeric(o.tolerance) || ~isreal(o.tolerance) || ~isscalar(o.tolerance) || ...
        ~isfinite(o.tolerance) || o.tolerance < 0
    error('hnf:Options','Invalid probability-sum tolerance.');
end
if o.normalize
    P = hnf_robprob(P,dim);
else
    if any(~isfinite(P(:))) || any(P(:)<0) || any(P(:)>1)
        error('hnf:Domain','Expected finite probabilities in [0,1].');
    end
    mass = sum(P,dim);
    if any(abs(mass(:)-1)>o.tolerance) || size(P,dim)==0
        error('hnf:ProbabilitySum','Probability slices must sum to one.');
    end
end
T = zeros(size(P),'like',P); positive = P > 0;
T(positive) = -P(positive).*log(P(positive));
H = sum(T,dim) ./ log(o.base);
if nargout > 1, info = struct('normalized',o.normalize,'base',o.base,'sum',sum(P,dim)); end
end
