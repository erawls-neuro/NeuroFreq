function [r, info] = hnf_robrank(A, opts)
% HNF_ROBRANK  Scale-safe SVD rank without forming a covariance/Gram matrix.
%   [r,info] = hnf_robrank(A,opts)
%   opts.relativeTolerance=[] => max(size(A))*eps(class(A));
%   opts.absoluteTolerance=0 is an additional singular-value cutoff in A units.
%   A must be full finite CPU single/double, real or complex, 2-D.
%   INFO exposes normalized singular values and matrixScale separately so
%   diagnostics themselves need not overflow. Rank depends on tolerance.
if nargin<2, opts=[]; end
hnf_la_check(A);
o = hnf_aux_options(opts,struct('relativeTolerance',[],'absoluteTolerance',0));
if isempty(o.relativeTolerance), o.relativeTolerance=min(1,max(size(A))*eps(class(A))); end
hnf_la_scalar(o.relativeTolerance,'relativeTolerance',0,1);
hnf_la_scalar(o.absoluteTolerance,'absoluteTolerance',0,realmax);
scale = hnf_la_scale(A);
if isempty(A), s=zeros(0,1,'like',real(A)); else, s=svd(A./scale,'econ'); end
if isempty(s), threshold=o.absoluteTolerance/scale;
else, threshold=max(o.relativeTolerance*double(s(1)),o.absoluteTolerance/scale); end
r = nnz(s>threshold);
if nargout>1
    info=struct('rank',r,'singularValuesScaled',s,'matrixScale',scale, ...
        'thresholdScaled',threshold,'relativeTolerance',o.relativeTolerance, ...
        'absoluteTolerance',o.absoluteTolerance);
end
end
