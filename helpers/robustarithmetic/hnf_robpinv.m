function [P, info] = hnf_robpinv(A, opts)
% HNF_ROBPINV  Scaled Moore-Penrose inverse; no observation-sized identity.
%   P = hnf_robpinv(A) or hnf_robpinv(A,TOL) retains singular values > TOL
%   in original units, matching pinv's absolute-tolerance convention.
%   Options: absoluteTolerance=[]; relativeTolerance=size*eps;
%   rcondThreshold=1e-12 (double), 1e-5 (single); ridge=0.
%   An explicitly supplied absoluteTolerance uses SVD, so no fast path can
%   bypass it. Relative ridge is defined on A/max-component(A); positive
%   ridge produces a regularized inverse, not a Moore-Penrose inverse.
% The original no-options/numeric-tolerance call zeroes nonfinite entries
% and retains its sparse inverse fast path. Struct options request strict input.
if nargin < 2, opts = []; end
if ~isstruct(opts)
    A(~isfinite(A)) = 0;
    if issparse(A) && isempty(opts) && ~isempty(A) && size(A,1)==size(A,2) && ...
            condest(A) < 1/max(1e-8,numel(A)*eps(class(A)))
        P = A \ speye(size(A,1));
        if nargout > 1, info = struct('method','sparse-solve','rank',size(A,1)); end
        return
    end
    if issparse(A), A = full(A); end
end
if isempty(opts), opts=struct; end
if isnumeric(opts), opts=struct('absoluteTolerance',opts); end
hnf_la_check(A);
threshold=1e-12; if isa(A,'single'), threshold=1e-5; end
o=hnf_la_opts(opts,struct('absoluteTolerance',[],'relativeTolerance',min(1,max(size(A))*eps(class(A))), ...
    'rcondThreshold',threshold,'ridge',0));
hnf_la_scalar(o.relativeTolerance,'relativeTolerance',0,1);
hnf_la_scalar(o.rcondThreshold,'rcondThreshold',0,1);
hnf_la_scalar(o.ridge,'ridge',0,realmax(class(A)));
if o.ridge>0 && cast(o.ridge,'like',A)==0, error('hnf:InvalidOption','ridge is below the input precision range.'); end
if ~isempty(o.absoluteTolerance), hnf_la_scalar(o.absoluteTolerance,'absoluteTolerance',0,realmax('double')); end
[m,n]=size(A); scale=hnf_la_scale(A); An=A./scale;
method='svd'; rankA=0; rc=0; actualTolerance=0;
if isempty(A)
    P=zeros(n,m,'like',A);
elseif m==n && ~isempty(o.absoluteTolerance) && o.ridge==0 && nnz(A)==nnz(diag(A))
    % With an absolute cutoff, tiny independent diagonal directions may
    % matter even beside realmax. Global scaling could underflow them to
    % zero before SVD. Invert precisely those original directions retained.
    d=diag(A); keep=abs(d)>o.absoluteTolerance;
    P=zeros(n,m,'like',A); i=find(keep);
    if ~isempty(i)
        reciprocal=hnf_robdiv(ones(size(i),'like',d),d(i),struct('overflow','allow'));
        if any(~isfinite(reciprocal(:)))
            error('hnf:UnrepresentableOutput','Pseudoinverse exceeds the output floating-point range.');
        end
        P(sub2ind([n m],i,i))=reciprocal;
    end
    rankA=nnz(keep); method='diagonal';
    relative=abs(d./scale);
    if any(relative), rc=double(min(relative)/max(relative)); end
    actualTolerance=double(o.absoluteTolerance)/double(scale);
else
    fast=false;
    if m>=n && o.ridge==0 && isempty(o.absoluteTolerance)
        [Q,R]=qr(An,0); rc=rcond(R);
        if rc>max(o.relativeTolerance,o.rcondThreshold)
            [P,~]=linsolve(R,Q',struct('UT',true));
            fast=all(isfinite(P(:))); method='qr'; rankA=n;
        end
        if ~fast, [U,S,V]=svd(R,'econ'); U=Q*U; end
    elseif m<n && o.ridge==0 && isempty(o.absoluteTolerance)
        [Q,R]=qr(An',0); rc=rcond(R);
        if rc>max(o.relativeTolerance,o.rcondThreshold)
            [T,~]=linsolve(R,Q',struct('UT',true)); P=T';
            fast=all(isfinite(P(:))); method='qr-transpose'; rankA=m;
        end
        if ~fast, [U,S,V]=svd(An,'econ'); end
    else
        [U,S,V]=svd(An,'econ');
    end
    if ~fast
        method='svd'; k=min(m,n); s=diag(S); s=s(1:k); U=U(:,1:k); V=V(:,1:k);
        if isempty(o.absoluteTolerance)
            cutoff=o.relativeTolerance*s(1);
        else
            cutoff=double(o.absoluteTolerance)/double(scale);
        end
        keep=s>cutoff; rankA=nnz(keep); actualTolerance=cutoff;
        if s(1)>0, rc=double(s(end)/s(1)); else, rc=0; end
        if o.ridge>0
            g=s./(s.^2+cast(o.ridge,'like',s));
            P=(V.*g')*U'; method='svd-ridge';
        else
            P=(V(:,keep)./s(keep)')*U(:,keep)';
        end
    end
    P=hnf_la_rescale(P,1,scale);
end
if nargout>1
    info=struct('method',method,'rank',rankA,'reciprocalCondition',rc, ...
        'scale',scale,'relativeTolerance',o.relativeTolerance, ...
        'absoluteTolerance',o.absoluteTolerance,'normalizedSVDCutoff',actualTolerance,'ridge',o.ridge);
end
end
