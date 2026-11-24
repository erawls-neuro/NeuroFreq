function [logAbsDet, phase, info] = hnf_roblogdet(A, opts)
% HNF_ROBLOGDET  Stable log(abs(det(A))) and signed/complex determinant phase.
%   [LOGABSDET,PHASE,INFO]=hnf_roblogdet(A,OPTS) avoids forming det(A).
%   Options: structure='general' (LU) or 'positive' (Cholesky, no repair);
%   relativeTolerance=0 (exact zero-pivot criterion); ridge=0 (relative to
%   max component of A); symmetry='check'/'symmetrize';
%   symmetryTolerance=32*n*eps. Symmetry options apply to 'positive' only.
%   General LU uses row equilibration; positive Cholesky uses diagonal
%   congruence scaling. Their log scale contributions are added separately.
%   A singular general matrix gives (-Inf,0), intentionally. det(zeros(0))
%   is 1, so the empty result is (0,1). This is NOT a pseudodeterminant.
if nargin<2, opts=struct; end
hnf_la_check(A); n=size(A,1);
if size(A,2)~=n, error('hnf:Shape','Expected a square matrix.'); end
o=hnf_la_opts(opts,struct('structure','general','relativeTolerance',0,'ridge',0, ...
    'symmetry','check','symmetryTolerance',min(1,32*max(1,n)*eps(class(A)))));
hnf_la_scalar(o.relativeTolerance,'relativeTolerance',0,1);
hnf_la_scalar(o.ridge,'ridge',0,realmax(class(A)));
if o.ridge>0 && cast(o.ridge,'like',A)==0, error('hnf:InvalidOption','ridge is below the input precision range.'); end
if ~any(strcmp(o.structure,{'general','positive'})), error('hnf:InvalidOption','structure must be general or positive.'); end
scale=hnf_la_scale(A); phase=cast(1,'like',A); singular=false;
if strcmp(o.structure,'positive')
    [H,scale]=hnf_la_hermitian(A,o.symmetry,o.symmetryTolerance);
else
    H=A./scale;
end
idx=1:n+1:n*n; H(idx)=H(idx)+cast(o.ridge,'like',H);
if any(~isfinite(H(:))), error('hnf:UnrepresentableOutput','Regularized matrix exceeds floating-point range.'); end
% With no requested ridge, retain original-unit entries until per-row or
% congruence equilibration; one global divisor can erase valid tiny pivots.
if o.ridge==0
    H=A; logScale=cast(0,'like',real(A));
else
    logScale=n*log(scale);
end
if n==0
    logAbsDet=cast(0,'like',real(A)); method='empty';
elseif strcmp(o.structure,'positive')
    diagonal=real(diag(H));
    if any(diagonal<=0), error('hnf:NotPositiveDefinite','Positive log determinant requires a positive diagonal.'); end
    sd=sqrt(diagonal);
    H=(H./sd)./sd';
    if any(~isfinite(H(:))), error('hnf:NotPositiveDefinite','Congruence scaling overflowed; matrix is not numerically positive definite.'); end
    H=H/2+H'/2;
    logScale=logScale+2*sum(log(sd));
    [L,p]=chol(H,'lower');
    if p~=0, error('hnf:NotPositiveDefinite','Positive log determinant requires a positive-definite matrix.'); end
    v=real(diag(L)); method='chol';
    singular=any(v<=sqrt(o.relativeTolerance)*max(v));
    logAbsDet=2*sum(log(v))+logScale;
else
    rowScale=max(abs(real(H)),[],2);
    if ~isreal(H), rowScale=max(rowScale,max(abs(imag(H)),[],2)); end
    rowScale(rowScale==0)=1;
    Heq=H./rowScale;
    erased=(real(H)~=0 & real(Heq)==0) | (imag(H)~=0 & imag(Heq)==0);
    if ~any(erased(:))
        H=Heq; logScale=logScale+sum(log(rowScale));
    end
    % Rare extreme within-row ranges use unscaled LU instead of silently
    % zeroing pivots. Fail explicitly if that fallback itself overflows.
    [~,U,p]=lu(H,'vector'); v=diag(U); method='lu';
    if any(~isfinite(U(:))), error('hnf:UnrepresentableOutput','LU factors exceed floating-point range.'); end
    pivotScale=max(abs(real(v)),abs(imag(v))); pivotScale(pivotScale==0)=1;
    vn=v./pivotScale; a=abs(vn); loga=log(pivotScale)+log(a);
    singular=any(a==0);
    if o.relativeTolerance>0, singular=singular || any(loga<=log(o.relativeTolerance)+max(loga)); end
    if ~singular
        visited=false(n,1); cycles=0;
        for k=1:n
            if ~visited(k)
                cycles=cycles+1; j=k;
                while ~visited(j), visited(j)=true; j=p(j); end
            end
        end
        phase=cast((-1)^(n-cycles),'like',A)*prod(vn./a);
        phase=phase/abs(phase);
        logAbsDet=sum(loga)+logScale;
    else
        logAbsDet=cast(-Inf,'like',real(A));
    end
end
if singular, logAbsDet=cast(-Inf,'like',real(A)); phase=cast(0,'like',A); end
if nargout>2
    info=struct('method',method,'singular',singular,'scale',scale, ...
        'ridge',o.ridge,'relativeTolerance',o.relativeTolerance);
end
end
