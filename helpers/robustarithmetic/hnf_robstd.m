function [s, info] = hnf_robstd(x, dim, opts, nanflag)
% HNF_ROBSTD  Stable standard deviation, MAD scale, or winsorized deviation.
%   S = HNF_ROBSTD(X,DIM,OPTS). Dense CPU single/double; DIM=[] automatic.
%   Complex X supported for method='std': sqrt(sum(abs(X-MU).^2)/denom).
%   opts.method='std'|'mad'|'winsor' (default 'std').
%   opts.normalization=0 uses max(n-1,1); 1 uses n. Applies to std/winsor.
%   opts.consistency=1.482602218505602 multiplies raw median absolute
%     deviation for method='mad'; this is an asymptotic normal-model factor.
%   opts.trim=.1 clips each tail at type-7 quantiles for method='winsor'.
%   opts.nonfinite='error'|'omit'|'propagate'; omit removes NaN AND Inf.
%   opts.overflow='error'|'allow'|'saturate'; default error. An actually
%     unrepresentable SD cannot be made finite without changing its value.
%   Constants/singletons have SD zero; empty/all-omitted slices give NaN.
if nargin < 2, dim = []; end
if nargin < 3, opts = struct; end
% MATLAB-style hnf_robstd(X,weight,DIM[,nanflag]) accepts scalar weights
% 0/1 for sample or population SD. A struct in position three retains
% the full robust-estimator options of hnf_robstd(X,DIM,OPTS).
native = isnumeric(opts) && nargin>=3;
if native
    w=dim; dim=opts;
    if nargin<4, nanflag='includenan'; end
    if ~any(strcmp(nanflag,{'omitnan','includenan'})), error('hnf:Options','Invalid NaN flag.'); end
    opts=struct('nonfinite','omit','normalization',w,'overflow','allow');
end
opts = stat_options(opts,struct('method','std','normalization',0, ...
    'consistency',1.482602218505602,'trim',.1,'nonfinite','error','overflow','error'));
opts.method = stat_choice(opts.method,{'std','mad','winsor'},'method');
stat_number(opts.normalization,0,1,'normalization',true);
stat_number(opts.consistency,0,Inf,'consistency'); stat_number(opts.trim,0,.5-eps,'trim');
[a,v,n,scale,~,shape,order,bad,raw] = stat_prepare(x,dim,opts.nonfinite,strcmp(opts.method,'std'));
if strcmp(opts.method,'mad')
    s = stat_rawmad(raw,v)*opts.consistency;
else
    if strcmp(opts.method,'winsor')
        raw = stat_transform(raw,v,'winsor',opts.trim);
        [a,v,~,scale] = stat_prepare(raw,1,'omit');
    end
    inner = opts; inner.method = 'std';
    if ~isreal(raw) && strcmp(opts.method,'std')
        % A huge constant real component must not erase a tiny varying
        % imaginary component (and conversely). Center/scale each part.
        re=real(raw); im=imag(raw); re(~v)=NaN; im(~v)=NaN;
        inner.nonfinite='omit'; inner.overflow='allow';
        s=hypot(hnf_robstd(re,1,inner),hnf_robstd(im,1,inner));
    else
        s = stat_residual_scale(a,v,inner).*scale;
    end
end
if strcmp(opts.nonfinite,'propagate'), s(bad) = NaN; end
[s,overflow] = stat_overflow(s,opts.overflow);
s = stat_restore(s,shape,order);
if native
    if isempty(dim), dim=find(size(x)~=1,1); if isempty(dim), dim=1; end, end
    rejected=any(isinf(x) & ~isnan(x),dim);
    if strcmp(nanflag,'includenan'), rejected=rejected | any(isnan(x),dim); end
    s(rejected)=NaN;
end
if nargout > 1
    info = struct('method',opts.method,'nonfinite',opts.nonfinite, ...
        'count',stat_restore(n,shape,order),'normalization',opts.normalization, ...
        'constant',s == 0,'overflow',stat_restore(overflow,shape,order), ...
        'invalid',stat_restore(bad,shape,order));
end
end
