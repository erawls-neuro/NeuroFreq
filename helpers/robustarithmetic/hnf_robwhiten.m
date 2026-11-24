function [W, info, Cused] = hnf_robwhiten(C, opts)
% HNF_ROBWHITEN  Covariance whitening with explicit rank/regularization policy.
%   [W,INFO,CUSED]=hnf_robwhiten(C,OPTS), applied as W*channelBySampleData.
%   C is a square Hermitian covariance, not a samples-by-channels data array.
%   If C=hnf_robcov(X) for samples-by-variables X, apply X*W'. For native
%   channel-by-sample D use C=D*D'/denominator (conjugate transpose).
%   Options: method='zca' (n-by-n) or 'pca' (rank-by-n);
%   mode='truncate' (default) or 'floor'; relativeTolerance=n*eps;
%   ridge=0 (relative to componentScale(C)); psd='error' or 'clip';
%   negativeTolerance=32*n*eps; symmetry='check'/'symmetrize';
%   symmetryTolerance=32*n*eps.
%   Truncation keeps eigenvalues > relativeTolerance*largest eigenvalue.
%   Floor mode explicitly raises every eigenvalue to that positive floor.
%   CUSED is the spectrally modified covariance. PCA yields W*CUSED*W'=I;
%   rank-deficient ZCA yields a projector, NOT a full identity.
if nargin<2, opts=struct; end
hnf_la_check(C); n=size(C,1); ep=eps(class(C));
o=hnf_la_opts(opts,struct('method','zca','mode','truncate', ...
    'relativeTolerance',min(1,max(1,n)*ep),'ridge',0,'psd','error', ...
    'negativeTolerance',min(1,32*max(1,n)*ep),'symmetry','check', ...
    'symmetryTolerance',min(1,32*max(1,n)*ep)));
hnf_la_scalar(o.relativeTolerance,'relativeTolerance',0,1);
hnf_la_scalar(o.negativeTolerance,'negativeTolerance',0,1);
hnf_la_scalar(o.ridge,'ridge',0,realmax(class(C)));
if o.ridge>0 && cast(o.ridge,'like',C)==0, error('hnf:InvalidOption','ridge is below the input precision range.'); end
if ~any(strcmp(o.method,{'zca','pca'})) || ~any(strcmp(o.mode,{'truncate','floor'})) || ~any(strcmp(o.psd,{'error','clip'}))
    error('hnf:InvalidOption','Invalid method, mode, or psd policy.');
end
if strcmp(o.mode,'floor') && o.relativeTolerance==0, error('hnf:InvalidOption','Floor mode requires a positive relativeTolerance.'); end
[H,scale,defect]=hnf_la_hermitian(C,o.symmetry,o.symmetryTolerance);
[V,D]=eig(H); raw=real(diag(D)); dn=raw+cast(o.ridge,'like',raw);
if any(~isfinite(dn)), error('hnf:UnrepresentableOutput','Regularized spectrum exceeds floating-point range.'); end
if isempty(dn), spectralScale=1; else, spectralScale=max(abs(dn)); end
if strcmp(o.psd,'error') && any(dn < -o.negativeTolerance*spectralScale)
    error('hnf:NotPositiveSemidefinite','Whitening requires a positive-semidefinite covariance.');
end
dn=max(dn,0);
if isempty(dn), top=0; else, top=max(dn); end
cutoff=cast(o.relativeTolerance,'like',dn)*top;
if strcmp(o.mode,'floor')
    if top==0, cutoff=cast(o.relativeTolerance,'like',dn); end
    if cutoff==0, error('hnf:InvalidOption','Requested spectral floor is below the input precision range.'); end
    dn=max(dn,cutoff); keep=true(size(dn));
else
    keep=dn>cutoff;
end
P=(V(:,keep)./sqrt(dn(keep))')';
if strcmp(o.method,'pca'), W=P; else, W=V(:,keep)*P; end
W=hnf_la_rescale(W,1,sqrt(scale));
if nargout>1
    used=dn; used(~keep)=0;
    info=struct('method',o.method,'mode',o.mode,'rank',nnz(keep),'scale',scale, ...
        'ridge',o.ridge,'relativeTolerance',o.relativeTolerance, ...
        'normalizedCutoff',cutoff,'eigenvaluesNormalized',dn,'kept',keep, ...
        'symmetryDefect',defect,'modifiedEigenvalues',nnz(used~=raw));
end
if nargout>2
    used=dn; used(~keep)=0;
    Cused=hnf_la_rescale((V.*used')*V',scale,1);
end
end
