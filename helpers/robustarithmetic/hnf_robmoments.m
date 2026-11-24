function [state, mu, sd] = hnf_robmoments(X, state, opts)
% HNF_ROBMOMENTS  Mergeable translated/scaled running mean and SD.
%   state = hnf_robmoments(X,[]) starts; X is observations-by-variables.
%   [state,mu,sd] = hnf_robmoments(X,state,opts) adds a batch.
%   [state,mu,sd] = hnf_robmoments([],state,opts) reads without adding data.
%   opts.nonfinite='error' | 'omit' (per variable); opts.normalization=0
%   uses sample SD, =1 uses population SD. opts.overflow='error' | 'allow'
%   | 'saturate' applies to requested SD. Counts are exact doubles up to 2^53.
%   A per-variable origin retains small variations around a large DC offset;
%   scaled centered sums avoid sum(X.^2)-n*mean(X).^2 cancellation. If a
%   difference from the origin would overflow, that column uses origin zero.
%   State is O(nVariables). Chunk size may change rounding. This is numerical
%   stabilization, not an outlier-robust estimator or compensated streaming
%   sum. Extreme signed cancellation can lose low-order mean contributions
%   (e.g., the 1 in [1e16;1;-1e16]); central-variance stability is the target.
%   Do not edit state fields.
if nargin<2, state = []; end
if nargin<3, opts = []; end
o = hnf_aux_options(opts,struct('nonfinite','error','normalization',0,'overflow','error'));
if ~any(strcmp(o.nonfinite,{'error','omit'})) || ...
        ~isscalar(o.normalization) || ~isreal(o.normalization) || ...
        ~any(o.normalization==[0 1]) || ~any(strcmp(o.overflow,{'error','allow','saturate'}))
    error('hnf:Options','Invalid running-moment options.');
end
readOnly = isequal(size(X),[0 0]) && ~isempty(state);
if ~readOnly
    hnf_aux_real(X);
    if ~ismatrix(X), error('hnf:Shape','X must be observations-by-variables.'); end
    p = size(X,2);
end
if isempty(state)
    if isempty(X), error('hnf:Empty','Initialization requires a nonempty batch.'); end
    z = zeros(1,p,'like',X);
    state = struct('version',2,'count',zeros(1,p),'origin',z,'scale',z, ...
        'meanScaled',z,'m2Scaled',z);
else
    validate_state(state);
    if ~readOnly && (size(X,2)~=numel(state.count) || ~strcmp(class(X),class(state.scale)))
        error('hnf:State','Batch variable count and precision must match the state.');
    end
end
if ~readOnly && size(X,1)>0
    good = isfinite(X);
    if strcmp(o.nonfinite,'error') && any(~good(:))
        error('hnf:Nonfinite','Running-moment input must be finite.');
    end
    X(~good) = 0; nb = sum(good,1);
    [~,first] = max(good,[],1);
    ob = X(sub2ind(size(X),first,1:p));
    D = X-ob; D(~good) = 0;
    badDifference = any(~isfinite(D),1);
    ob(badDifference) = 0; D(:,badDifference) = X(:,badDifference);
    sb = max(abs(D),[],1); safe = sb; safe(safe==0) = 1;
    U = D./safe; mb = sum(U,1)./max(nb,1);
    D = U-mb; D(~good) = 0; m2b = sum(D.^2,1);
    na = state.count; nt = na+nb;
    if any(nt>flintmax), error('hnf:CountOverflow','Count exceeds exact double integer range.'); end
    origin = state.origin; origin(na==0) = ob(na==0);
    da = state.origin-origin; da(na==0) = 0; db = ob-origin;
    needZeroOrigin = ~isfinite(db);
    origin(needZeroOrigin) = 0;
    da(needZeroOrigin) = state.origin(needZeroOrigin);
    db(needZeroOrigin) = ob(needZeroOrigin);
    s = max(max(state.scale,sb),max(abs(da),abs(db)));
    safe = s; safe(safe==0) = 1;
    ra = state.scale./safe; rb = sb./safe;
    ma = da./safe+state.meanScaled.*ra; ma(na==0) = 0;
    mb = db./safe+mb.*rb;
    delta = mb-ma; fraction = nb./max(nt,1);
    mergedMean = ma+delta.*fraction;
    mergedM2 = state.m2Scaled.*ra.^2+m2b.*rb.^2+delta.^2.*(na.*fraction);
    if any(~isfinite(mergedMean)) || any(~isfinite(mergedM2))
        error('hnf:StateOverflow','Running centered moments exceed the state precision range.');
    end
    active = nb>0;
    state.count(active) = nt(active); state.origin(active) = origin(active);
    state.scale(active) = s(active); state.meanScaled(active) = mergedMean(active);
    state.m2Scaled(active) = max(mergedM2(active),0);
end
if nargout>1
    mu = state.origin+state.scale.*state.meanScaled; mu(state.count==0) = NaN;
    retry = isinf(mu) & state.scale>0;
    if any(retry)
        mu(retry) = state.scale(retry).*(state.meanScaled(retry)+ ...
            state.origin(retry)./state.scale(retry));
    end
    if any(isinf(mu)), error('hnf:Overflow','Requested mean is not representable.'); end
end
if nargout>2
    denom = state.count-(o.normalization==0);
    sd = state.scale.*sqrt(state.m2Scaled./max(denom,1)); sd(denom<=0) = NaN;
    overflow = isinf(sd);
    if any(overflow)
        if strcmp(o.overflow,'error'), error('hnf:Overflow','Requested SD is not representable.'); end
        if strcmp(o.overflow,'saturate'), sd(overflow)=realmax(class(state.scale)); end
    end
end
end

function validate_state(s)
fields = {'version','count','origin','scale','meanScaled','m2Scaled'};
if ~isstruct(s) || ~isscalar(s) || ~all(isfield(s,fields)) || ~isequal(s.version,2)
    error('hnf:State','Invalid running-moment state.');
end
if ~isa(s.count,'double') || ~isrow(s.count) || any(~isfinite(s.count)) || ...
        any(s.count<0 | s.count~=fix(s.count) | s.count>flintmax)
    error('hnf:State','Invalid running counts.');
end
for k=3:numel(fields)
    A = s.(fields{k}); hnf_aux_real(A);
    if ~isequal(size(A),size(s.count)) || any(~isfinite(A)) || ...
            ~strcmp(class(A),class(s.scale))
        error('hnf:State','Invalid running-moment fields.');
    end
end
if any(s.scale<0) || any(s.m2Scaled<0)
    error('hnf:State','Invalid scale or centered sum of squares.');
end
end
