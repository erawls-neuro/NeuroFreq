function [mu, info] = hnf_robmean(x, dim, opts)
% HNF_ROBMEAN  Stable classical or opt-in robust location along a dimension.
%   MU = HNF_ROBMEAN(X,DIM,OPTS). DIM=[] selects first nonsingleton dimension.
%   X: dense CPU single/double; complex supported for method='mean' only.
%   opts.method: 'mean' (default), 'median', 'trim', 'winsor', 'huber'.
%   opts.nonfinite: 'error' (default), 'omit' (NaN AND Inf), 'propagate'.
%   opts.trim=.1: fraction in EACH tail; trim drops floor(n*trim), whereas
%     winsor clips to interpolated type-7 quantiles (not a trimmed mean).
%   opts.tuning=1.345 (eps<=tuning<=1e6), maxiter=50,
%     tol=sqrt(eps(class(X))): Huber location
%     IRLS with fixed 1.4826022185*MAD scale initialized at the median.
%     A zero initial MAD returns the median, flagged in INFO.zeroScale.
%   Empty/all-omitted slices return NaN. No printing or global warning changes.
%   INFO reports finite counts, counts used, iteration/convergence status.
if nargin < 2, dim = []; end
if nargin < 3, opts = struct; end
native = ischar(opts);
if native
    nanflag=opts;
    if ~any(strcmp(nanflag,{'omitnan','includenan'})), error('hnf:Options','Invalid NaN flag.'); end
    opts=struct('nonfinite','omit');
end
opts = stat_options(opts,struct('method','mean','nonfinite','error', ...
    'trim',.1,'tuning',1.345,'maxiter',50,'tol',sqrt(eps(class(x)))));
opts.method = stat_choice(opts.method,{'mean','median','trim','winsor','huber'},'method');
stat_number(opts.trim,0,.5-eps,'trim'); stat_number(opts.tuning,eps,1e6,'tuning');
stat_number(opts.maxiter,1,1e6,'maxiter',true); stat_number(opts.tol,0,1,'tol');
[a,v,n,s,o,shape,order,bad,raw] = stat_prepare(x,dim,opts.nonfinite,strcmp(opts.method,'mean'));
it = zeros(size(n)); conv = true(size(n)); zs = false(size(n)); used = n;
switch opts.method
    case 'median'
        mu = stat_quantile(raw,v,.5);
    case 'huber'
        [mu,it,conv,zs] = stat_huber(raw,v,opts);
    otherwise
        if any(strcmp(opts.method,{'trim','winsor'}))
            [raw,~,used] = stat_transform(raw,v,opts.method,opts.trim);
            [a,v,~,s,o,~,~,~,raw] = stat_prepare(raw,1,'omit');
        end
        inner = opts; inner.method = 'mean';
        mu = stat_location(a,v,inner);
        mu = stat_location_restore(mu,s,o);
        count = sum(v,1);
        threshold = sqrt(eps(class(raw)))*s;
        repair = count > 0 & abs(real(mu)) <= threshold;
        if ~isreal(raw), repair = repair | (count > 0 & abs(imag(mu)) <= threshold); end
        if any(repair)
            values = raw(:,repair);
            recovered = componentMean(real(values),count(repair));
            if ~isreal(values)
                recovered = complex(recovered,componentMean(imag(values),count(repair)));
            end
            mu(repair) = recovered;
        end
end
if strcmp(opts.nonfinite,'propagate'), mu(bad) = NaN; end
mu = stat_restore(mu,shape,order);
if native
    if isempty(dim), dim=find(size(x)~=1,1); if isempty(dim), dim=1; end, end
    if isreal(x)
        plusInf=any(x==Inf,dim); minusInf=any(x==-Inf,dim);
        mu(plusInf)=Inf; mu(minusInf)=-Inf; mu(plusInf & minusInf)=NaN;
    elseif any(isinf(x(:)) & ~isnan(x(:)))
        % Preserve the independent real/imaginary components of an IEEE mean.
        % If either component is NaN, 'omitnan' omits the entire observation.
        z=x; z(isnan(z))=complex(NaN,NaN);
        mu=complex(hnf_robmean(real(z),dim,nanflag),hnf_robmean(imag(z),dim,nanflag));
    end
    if strcmp(nanflag,'includenan'), mu(any(isnan(x),dim))=NaN; end
end
if nargout > 1
    info = struct('method',opts.method,'nonfinite',opts.nonfinite, ...
        'count',stat_restore(n,shape,order),'usedCount',stat_restore(used,shape,order), ...
        'iterations',stat_restore(it,shape,order),'converged',stat_restore(conv,shape,order), ...
        'zeroScale',stat_restore(zs,shape,order),'invalid',stat_restore(bad,shape,order));
end
end

function y = componentMean(values,count)
% Repair a cancellation-dominated sum with Neumaier compensation. Scale by
% a power of two bounded by the sample count only where overflow threatens.
% Scaling by max(abs(x)) would erase 1e-200 in [realmax -realmax 1e-200].
values = double(values);
limit = realmax ./ max(count,1);
shift = zeros(size(count));
danger = max(abs(values),[],1) > limit;
shift(danger) = ceil(log2(count(danger)));
factor = pow2(1,-shift);
total = zeros(size(count));
correction = zeros(size(count));
for i = 1:size(values,1)
    x = values(i,:).*factor;
    next = total + x;
    dominant = abs(total)>=abs(x);
    residual = (x-next)+total;
    residual(dominant) = (total(dominant)-next(dominant))+x(dominant);
    correction = correction + residual;
    total = next;
end
y = (total+correction)./(count.*factor);
end
