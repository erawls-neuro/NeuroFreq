function [V, d, info] = hnf_robeig(A, opts)
% HNF_ROBEIG  Scaled Hermitian eigensystem with opt-in spectral controls.
%   [V,D,INFO] = hnf_robeig(A,OPTS), where D is an eigenvalue VECTOR.
%   Options: symmetry='check'/'symmetrize'; symmetryTolerance=32*n*eps;
%   order='ascend'/'descend'; psd='none'/'clip'/'error'; relativeFloor=0;
%   negativeTolerance=32*n*eps; ridge=0.
%   ridge adds ridge*componentScale(A) to each eigenvalue. relativeFloor
%   clips eigenvalues below relativeFloor*max(abs(eigenvalues)), only when
%   psd='clip'. 'error' rejects materially negative eigenvalues but retains
%   tiny negative eigenvalues rather than silently altering them.
%   General nonsymmetric eigenproblems are intentionally excluded.
if nargin<2, opts=struct; end
hnf_la_check(A); n=size(A,1); ep=eps(class(A));
o=hnf_la_opts(opts,struct('symmetry','check','symmetryTolerance',min(1,32*max(1,n)*ep), ...
    'order','ascend','psd','none','relativeFloor',0,'negativeTolerance',min(1,32*max(1,n)*ep),'ridge',0));
hnf_la_scalar(o.relativeFloor,'relativeFloor',0,1);
hnf_la_scalar(o.negativeTolerance,'negativeTolerance',0,1);
hnf_la_scalar(o.ridge,'ridge',0,realmax(class(A)));
if o.ridge>0 && cast(o.ridge,'like',A)==0, error('hnf:InvalidOption','ridge is below the input precision range.'); end
if ~any(strcmp(o.order,{'ascend','descend'})) || ~any(strcmp(o.psd,{'none','clip','error'}))
    error('hnf:InvalidOption','Invalid order or psd policy.');
end
if o.relativeFloor>0 && ~strcmp(o.psd,'clip'), error('hnf:InvalidOption','relativeFloor requires psd=clip.'); end
[H,scale,defect]=hnf_la_hermitian(A,o.symmetry,o.symmetryTolerance);
[V,D]=eig(H); raw=real(diag(D)); shifted=raw+cast(o.ridge,'like',raw);
if any(~isfinite(shifted)), error('hnf:UnrepresentableOutput','Regularized spectrum exceeds floating-point range.'); end
if isempty(raw), spectralScale=1; else, spectralScale=max(abs(shifted)); end
negative=shifted < -o.negativeTolerance*spectralScale;
if strcmp(o.psd,'error') && any(negative), error('hnf:NotPositiveSemidefinite','Matrix has materially negative eigenvalues.'); end
dn=shifted;
if strcmp(o.psd,'clip'), dn=max(dn,cast(o.relativeFloor,'like',dn)*spectralScale); end
[dn,order]=sort(dn,o.order); V=V(:,order);
d=hnf_la_rescale(dn,scale,1);
if nargout>2
    info=struct('scale',scale,'symmetryDefect',defect,'ridge',o.ridge, ...
        'psd',o.psd,'relativeFloor',o.relativeFloor,'modifiedEigenvalues',nnz(dn~=raw(order)), ...
        'rawEigenvaluesNormalized',raw(order));
end
end
