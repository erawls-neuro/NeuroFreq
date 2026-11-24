function [L, info, Aused] = hnf_robchol(A, opts)
% HNF_ROBCHOL  Quiet Cholesky with bounded, reported diagonal-shift repair.
%   [L,INFO,AUSED] = hnf_robchol(A,OPTS) gives L*L' approximately AUSED.
%   hnf_robchol(A,TOL) retains legacy meaning: add absolute TOL to diagonal.
%   Options: shift=0; repair='shift' or 'none'; symmetry='check' (default)
%   or 'symmetrize'; symmetryTolerance=32*n*eps; relativeFloor=16*n*eps;
%   maxRelativeShift=max(1e-6,128*n*eps); maxIterations=4.
%   The default repair handles roundoff-sized indefiniteness, not arbitrary
%   covariance changes. Larger repairs need explicit maxRelativeShift.
%   relativeFloor is applied only after Cholesky fails; a successful direct
%   factorization is returned unchanged apart from Hermitian roundoff averaging.
%   The repair budget bounds EXTRA shift / max(componentScale(A),shift).
%   No nonfinite-entry substitution and no changes to MATLAB warning state.
% The positional interface zeroes nonfinite entries, symmetrizes, and
% retains the original singular-covariance ridge. Struct options opt into
% the bounded repair and diagnostics described above.
if nargin < 2 || isempty(opts), opts=0; end
if ~isstruct(opts)
    [L,info,Aused] = classic(A,opts,nargout);
    return
end
hnf_la_check(A); n=size(A,1); ep=eps(class(A));
o=hnf_la_opts(opts,struct('shift',0,'repair','shift','symmetry','check', ...
    'symmetryTolerance',min(1,32*max(1,n)*ep),'relativeFloor',min(1,16*max(1,n)*ep), ...
    'maxRelativeShift',max(1e-6,128*max(1,n)*ep),'maxIterations',4));
hnf_la_scalar(o.shift,'shift',0,realmax(class(A)));
if o.shift>0 && cast(o.shift,'like',A)==0, error('hnf:InvalidOption','shift is below the input precision range.'); end
hnf_la_scalar(o.relativeFloor,'relativeFloor',ep,1);
hnf_la_scalar(o.maxRelativeShift,'maxRelativeShift',0,realmax(class(A)));
hnf_la_scalar(o.maxIterations,'maxIterations',1,100);
if o.maxIterations~=fix(o.maxIterations), error('hnf:InvalidOption','maxIterations must be an integer.'); end
if ~any(strcmp(o.repair,{'none','shift'})), error('hnf:InvalidOption','repair must be none or shift.'); end
[H,scale,defect]=hnf_la_hermitian(A,o.symmetry,o.symmetryTolerance);
if o.shift>scale, H=H.*(scale/o.shift); scale=cast(o.shift,'like',real(A)); end
shiftN=cast(o.shift,'like',real(A))/scale;
idx=1:n+1:n*n; H(idx)=H(idx)+shiftN;
[Ln,p]=chol(H,'lower'); autoShift=0; attempts=1; minimumEigenvalue=[];
if p~=0
    if strcmp(o.repair,'none'), error('hnf:NotPositiveDefinite','Cholesky failed and repair is disabled.'); end
    eigenvalues=real(eig(H)); minimumEigenvalue=min(eigenvalues);
    floorN=cast(o.relativeFloor,'like',eigenvalues)*max(1,max(abs(eigenvalues)));
    autoShift=max(floorN-minimumEigenvalue,floorN);
    Hbase=H;
    for k=1:o.maxIterations
        if autoShift>o.maxRelativeShift
            error('hnf:ShiftLimit','Required diagonal shift exceeds maxRelativeShift.');
        end
        H=Hbase; H(idx)=H(idx)+autoShift;
        [Ln,p]=chol(H,'lower'); attempts=attempts+1;
        if p==0, break; end
        autoShift=2*autoShift;
    end
    if p~=0, error('hnf:FactorizationFailed','Bounded Cholesky repair did not succeed.'); end
end
L=Ln.*sqrt(scale);
if any(~isfinite(L(:))), error('hnf:UnrepresentableOutput','Cholesky factor exceeds the output range.'); end
if nargout>1
    added=hnf_la_rescale(autoShift,scale,1);
    totalShift=double(o.shift)+double(added);
    if ~isfinite(totalShift), error('hnf:UnrepresentableOutput','Total diagonal shift exceeds floating-point range.'); end
    info=struct('method','chol','repaired',autoShift>0,'shift',totalShift, ...
        'requestedShift',o.shift,'repairShift',added,'relativeRepairShift',autoShift, ...
        'scale',scale,'symmetryDefect',defect,'attempts',attempts, ...
        'minimumEigenvalueNormalized',minimumEigenvalue);
end
if nargout>2, Aused=hnf_la_rescale(H,scale,1); end
end

function [L,info,Aused] = classic(A,shift,nout)
validateattributes(shift,{'numeric'},{'scalar','real','finite','nonnegative'});
if ~ismatrix(A) || size(A,1)~=size(A,2), error('hnf:Shape','Expected a square matrix.'); end
A(~isfinite(A)) = 0;
% Keep a well-scaled small diagonal in its original units even if a
% different diagonal is near realmax. Average only unequal Hermitian pairs.
idx = A~=A';
if any(idx(:))
    partner = A'; a=A(idx); b=partner(idx);
    midpoint = (a+b)./2;
    risk = ~isfinite(midpoint);
    midpoint(risk) = a(risk)./2 + b(risk)./2;
    A(idx) = midpoint;
end
n=size(A,1); diagIdx=1:n+1:n*n;
scale=1;
if any(real(diag(A))>realmax(class(A))-shift)
    scale=max(shift,max(abs(real(diag(A)))));
    A=A./scale;
end
A(diagIdx)=A(diagIdx)+shift./scale;
[L,p]=chol(A,'lower'); ridge=0;
if p
    d=abs(diag(A)); peak=max(d);
    if peak==0
        ridge=eps(class(A))./scale;
    else
        rel=max(1e-8,10*eps(class(A)));
        ridge=max(eps(cast(0,'like',d)),peak .* (rel.*mean(d./peak)));
    end
    if any(real(diag(A))>realmax(class(A))-ridge)
        A=A./peak;
        ridge=ridge./peak;
        scale=scale.*peak;
    end
    A(diagIdx)=A(diagIdx)+ridge;
    L=chol(A,'lower');
end
L=L.*sqrt(scale);
info=[]; Aused=[];
if nout>1
    info=struct('method','chol','repaired',ridge>0,'shift',shift+ridge*scale, ...
        'requestedShift',shift,'repairShift',ridge*scale,'scale',scale,'attempts',1+(ridge>0));
end
if nout>2, Aused=A.*scale; end
end
