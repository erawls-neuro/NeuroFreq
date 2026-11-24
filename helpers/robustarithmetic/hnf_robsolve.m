function [X, info] = hnf_robsolve(A, B, opts)
% HNF_ROBSOLVE  Solve A*X approximately B with quiet QR/Cholesky/SVD routes.
%   [X,INFO]=hnf_robsolve(A,B,OPTS) accepts hnf_robols options plus
%   structure='general' (default) or 'positive' (Hermitian positive definite).
%   The positive hint uses Cholesky if successful and well conditioned;
%   otherwise ordinary QR/SVD supplies a minimum-norm least-squares answer.
%   No diagonal shift occurs automatically. For an intentional repaired
%   covariance solve, first obtain AUSED from hnf_robchol explicitly.
%   Rectangular/rank-deficient A, WLS and robust IRLS use hnf_robols.
if nargin<3 || isempty(opts), opts=struct; end
if ~isstruct(opts) || ~isscalar(opts), error('hnf:InvalidOptions','Options must be a scalar structure.'); end
structure='general';
if isfield(opts,'structure'), structure=opts.structure; opts=rmfield(opts,'structure'); end
if ~any(strcmp(structure,{'general','positive'})), error('hnf:InvalidOption','structure must be general or positive.'); end
% Validate and solve through the common interface when fitting controls
% are present; this prevents the structural hint bypassing any LS options.
canChol=strcmp(structure,'positive') && isempty(fieldnames(opts));
if canChol
    hnf_la_check(A); hnf_la_check(B);
    if ~strcmp(class(A),class(B)), error('hnf:ClassMismatch','A and B must have the same precision.'); end
    if size(A,1)~=size(B,1), error('hnf:Shape','A and B must have matching rows.'); end
    n=size(A,1);
    [H,scale,defect]=hnf_la_hermitian(A,'check',min(1,32*max(1,n)*eps(class(A))));
    [L,p]=chol(H,'lower');
    if p==0 && ~isempty(H)
        rc=rcond(L); threshold=1e-12; if isa(A,'single'), threshold=1e-5; end
        if rc^2>max(threshold,max(size(A))*eps(class(A)))
            scaleB=ones(1,size(B,2),'like',real(B));
            for j=1:size(B,2), scaleB(j)=hnf_la_scale(B(:,j)); end
            [T,~]=linsolve(L,B./scaleB,struct('LT',true));
            [X,~]=linsolve(L',T,struct('UT',true));
            for j=1:size(X,2), X(:,j)=hnf_la_rescale(X(:,j),scaleB(j),scale); end
            if nargout>1
                info=struct('method','chol','rank',n,'factorReciprocalCondition',rc, ...
                    'usedPseudoinverse',false,'ridge',0,'symmetryDefect',defect);
            end
            return;
        end
    end
end
if nargout>1
    [X,~,~,info]=hnf_robols(A,B,opts);
else
    X=hnf_robols(A,B,opts);
end
end
