function [B, reciprocalCondition, usedPseudoinverse, info] = hnf_robols(A, Y, opts)
% HNF_ROBOLS  Quiet, scaled QR/SVD least squares and optional robust IRLS.
%   [B,RC,USEDPINV,INFO] = hnf_robols(A,Y,OPTS)
%   A is observations-by-predictors; Y is observations-by-responses.
%   No intercept is added. Complex ordinary/weighted LS uses conjugation.
%
%   OPTS fields (all optional):
%     relativeTolerance  SVD cutoff / largest singular value (size*eps).
%     rcondThreshold     QR fast-path threshold (1e-12 double, 1e-5 single).
%     ridge              >=0 penalty on normalized design (default 0).
%                        Equivalent original penalty ridge*scaleA^2, where
%                        scaleA=max(abs(real(A)),abs(imag(A))).
%     weights            nonnegative observation weights; normalized to max 1.
%     missing            'error' (default) or 'omitrows' (listwise A and Y).
%     robust             'none' (default), 'huber', or 'bisquare' (real only).
%     tuning             [] selects 1.345 Huber / 4.685 bisquare.
%     maxIterations      IRLS limit (50); convergenceTolerance (sqrt(eps)).
%   Responses use separate IRLS weights. Huber/bisquare control residual
%   outliers, not leverage contamination; no inferential SE/p-values implied.
%   A zero or rank-deficient design has the minimum-norm truncated-SVD fit.
%   Nonfinite/unsupported inputs and unrepresentable results raise hnf:* errors.
%   INFO is built only when requested; its weights correspond to the final
%   weighted solve (including supplied weights), not a post-fit extra update.
% The published two-input interface uses a 1e-9 rank/QR cutoff.
% Explicit options select the scaled, configurable solver described above.
if nargin < 3 || isempty(opts)
    opts = struct('relativeTolerance',1e-9,'rcondThreshold',1e-9);
    if issparse(A) || issparse(Y)
        % Octave/MATLAB use sparse QR for tall designs. Keep that route
        % without constructing a dense observation-by-predictor matrix.
        [Q,R] = qr(A,0);
        [B,reciprocalCondition,usedPseudoinverse,info] = hnf_robols(full(R),full(Q'*Y),opts);
        return
    end
end
hnf_la_check(A,true); hnf_la_check(Y,true);
if ~strcmp(class(A),class(Y)), error('hnf:ClassMismatch','A and Y must have the same precision.'); end
if size(A,1) ~= size(Y,1), error('hnf:Shape','A and Y must have matching rows.'); end
ep = eps(class(A));
threshold = 1e-12; if isa(A,'single'), threshold = 1e-5; end
o = hnf_la_opts(opts,struct('relativeTolerance',min(1,max(size(A))*ep), ...
    'rcondThreshold',threshold,'ridge',0,'weights',[],'missing','error', ...
    'robust','none','tuning',[],'maxIterations',50,'convergenceTolerance',sqrt(ep)));
hnf_la_scalar(o.relativeTolerance,'relativeTolerance',0,1);
hnf_la_scalar(o.rcondThreshold,'rcondThreshold',0,1);
hnf_la_scalar(o.ridge,'ridge',0,realmax(class(A)));
if o.ridge>0 && cast(o.ridge,'like',A)==0, error('hnf:InvalidOption','ridge is below the input precision range.'); end
hnf_la_scalar(o.maxIterations,'maxIterations',1,1e6);
if o.maxIterations ~= fix(o.maxIterations), error('hnf:InvalidOption','maxIterations must be an integer.'); end
hnf_la_scalar(o.convergenceTolerance,'convergenceTolerance',ep,1);
if ~any(strcmp(o.missing,{'error','omitrows'})) || ~any(strcmp(o.robust,{'none','huber','bisquare'}))
    error('hnf:InvalidOption','Invalid missing or robust policy.');
end
isRobust = ~strcmp(o.robust,'none');
if isRobust && (~isreal(A) || ~isreal(Y)), error('hnf:ComplexRobust','IRLS requires real inputs.'); end
if isempty(o.tuning)
    o.tuning = 1.345; if strcmp(o.robust,'bisquare'), o.tuning = 4.685; end
end
hnf_la_scalar(o.tuning,'tuning',ep,realmax(class(A)));
mOriginal = size(A,1);
if isempty(o.weights)
    weights = ones(mOriginal,1,'like',real(A));
else
    hnf_la_check(o.weights);
    if ~isreal(o.weights) || numel(o.weights) ~= mOriginal || ~isvector(o.weights) || any(o.weights(:)<0)
        error('hnf:InvalidWeights','weights must be a nonnegative vector with one value per row.');
    end
    weights = o.weights(:);
end
rows = all(isfinite(A),2) & all(isfinite(Y),2);
if any(~rows) && strcmp(o.missing,'error'), error('hnf:NonfiniteInput','A or Y contains NaN or Inf.'); end
if any(~rows), A = A(rows,:); Y = Y(rows,:); weights = weights(rows); end
if ~isempty(weights) && max(weights)>0, weights=weights./max(weights); end
weights=cast(weights,'like',real(A));
if ~isempty(weights) && ~any(weights>0), error('hnf:InvalidWeights','At least one retained weight must be positive.'); end
if mOriginal>0 && isempty(weights), error('hnf:NoObservations','No finite observations remain.'); end
% Inactive observations must not set the scale of the actual fitting data.
inactive=weights==0;
if any(inactive), A(inactive,:)=0; Y(inactive,:)=0; end
scaleA = hnf_la_scale(A);
scaleY = ones(1,size(Y,2),'like',real(Y));
for j=1:size(Y,2), scaleY(j)=hnf_la_scale(Y(:,j)); end
An = A./scaleA; Yn = Y./scaleY;
sw = sqrt(weights);
[Bn,reciprocalCondition,usedPseudoinverse,rankA,method] = hnf_la_ls(An.*sw,Yn.*sw,o);
iterations = zeros(1,size(Y,2)); converged = true(1,size(Y,2));
if nargout > 3 || isRobust, finalWeights = repmat(weights,1,size(Y,2)); end
if isRobust && ~isempty(Y)
    ranks = repmat(rankA,1,size(Y,2)); rcs = repmat(reciprocalCondition,1,size(Y,2));
    pinvs = repmat(usedPseudoinverse,1,size(Y,2)); methods = cell(1,size(Y,2));
    for j = 1:size(Y,2)
        b = Bn(:,j); converged(j) = false;
        for it = 1:o.maxIterations
            r = Yn(:,j)-An*b;
            active = weights>0;
            if any(~isfinite(r(active))), error('hnf:UnrepresentableOutput','IRLS residual exceeds floating-point range.'); end
            sigma = median(abs(r(active)-median(r(active))))/0.674489750196082;
            floorScale = ep*max(1,max(abs(Yn(active,j))));
            if sigma<=floorScale
                % Repeated residuals may have zero MAD before the intercept
                % is corrected. An uncentered median prevents all-zero
                % bisquare weights in that common exact-fit limit.
                sigma=median(abs(r(active)))/0.674489750196082;
            end
            sigma = max(sigma,floorScale);
            u = abs(r)./sigma./o.tuning;
            if strcmp(o.robust,'huber')
                rw = ones(size(u),'like',u); large = u>1; rw(large) = 1./u(large);
            else
                rw = max(0,1-min(u,1).^2).^2;
            end
            w = weights.*rw;
            if ~any(w>0), error('hnf:DegenerateWeights','IRLS removed every observation.'); end
            swj = sqrt(w);
            [next,rc,up,rk,meth] = hnf_la_ls(An.*swj,Yn(:,j).*swj,o);
            finalWeights(:,j) = w; ranks(j)=rk; rcs(j)=rc; pinvs(j)=up; methods{j}=meth;
            delta = norm(next-b)/max(1,norm(b)); b = next;
            iterations(j)=it;
            if delta <= o.convergenceTolerance, converged(j)=true; break; end
        end
        Bn(:,j)=b;
    end
    rankA=ranks; reciprocalCondition=min(rcs); usedPseudoinverse=any(pinvs); method=methods;
end
B=Bn;
for j=1:size(B,2), B(:,j)=hnf_la_rescale(Bn(:,j),scaleY(j),scaleA); end
if nargout > 3
    rowsUsed=rows; rowsUsed(rows)=weights>0;
    info = struct('method',{method},'rank',rankA,'reciprocalCondition',reciprocalCondition, ...
        'usedPseudoinverse',usedPseudoinverse,'relativeTolerance',o.relativeTolerance, ...
        'ridge',o.ridge,'designScale',scaleA,'responseScale',scaleY, ...
        'rowsUsed',rowsUsed,'rowsRetained',rows,'weights',finalWeights,'robust',o.robust, ...
        'iterations',iterations,'converged',converged);
end
end
