function [z, location, scale, info] = hnf_robzscore(x, dim, opts)
% HNF_ROBZSCORE  Stable standardization without subtracting huge raw values.
%   [Z,LOCATION,SCALE,INFO] = HNF_ROBZSCORE(X,DIM,OPTS), CPU single/double.
%   Complex X supported in standard mode, using real magnitude-based SD.
%   opts.method='standard' (mean/SD) or 'mad' (median/normal-consistent MAD).
%   opts.normalization=0 (sample) or 1 (population); standard mode only.
%   opts.constant='zero' (default), 'nan', 'error' controls zero-scale slices.
%     This includes nonconstant slices with zero MAD in method='mad'.
%   opts.nonfinite='error'|'omit'|'propagate'; omitted entries in Z stay NaN.
%   opts.overflow='error'|'allow'|'saturate' controls unrepresentable Z/SCALE.
%   SCALE may overflow even when Z is representable: scale overflow is only
%   evaluated when that third output is requested. INFO includes finite counts.
if nargin < 2, dim = []; end
if nargin < 3, opts = struct; end
opts = stat_options(opts,struct('method','standard','normalization',0, ...
    'constant','zero','nonfinite','error','overflow','error'));
opts.method = stat_choice(opts.method,{'standard','mad'},'method');
opts.constant = stat_choice(opts.constant,{'zero','nan','error'},'constant');
stat_number(opts.normalization,0,1,'normalization',true);
[a,v,n,s,o,shape,order,bad,raw] = stat_prepare(x,dim,opts.nonfinite,strcmp(opts.method,'standard'));
sopts = struct('method','std','normalization',opts.normalization, ...
    'consistency',1.482602218505602,'trim',.1);
isMad = strcmp(opts.method,'mad');
scaleMultiplier = 1;
if isMad
    [sd,mu] = stat_rawmad(raw,v);
    a = raw; s = ones(size(sd),'like',sd);
    scaleMultiplier = 1.482602218505602;
else
    [sd,mu] = stat_residual_scale(a,v,sopts);
    centeredRaw = raw; centeredRaw(~v) = NaN;
    physicalLocation = hnf_robmean(centeredRaw,1,struct('nonfinite','omit'));
    a = bsxfun(@minus,raw,physicalLocation);
    overflowDifference = isinf(a);
    a = bsxfun(@rdivide,a,s);
    if any(overflowDifference(:))
        fallback = bsxfun(@minus,bsxfun(@rdivide,raw,s),physicalLocation./s);
        a(overflowDifference) = fallback(overflowDifference);
    end
    mu = zeros(size(mu),'like',mu);
end
constant = sd == 0 & n > 0;
if strcmp(opts.constant,'error') && any(constant)
    error('hnf:statistics:ZeroScale','At least one slice has zero scale.');
end
denominator = sd; denominator(constant) = 1;
z = bsxfun(@minus,a,mu);
badDifference = isinf(z);
z = bsxfun(@rdivide,z,denominator);
if any(badDifference(:))
    fallback = bsxfun(@minus,bsxfun(@rdivide,a,denominator),mu./denominator);
    z(badDifference) = fallback(badDifference);
end
z = z/scaleMultiplier;
if strcmp(opts.constant,'zero'), z(:,constant) = 0;
else, z(:,constant) = NaN; end
z(~v) = NaN;
if strcmp(opts.nonfinite,'propagate'), z(:,bad) = NaN; end
[z,zoverflow] = stat_overflow(z,opts.overflow);
z = stat_restore(z,shape,order);
if nargout > 1
    if isMad, location = mu;
    else, location = physicalLocation; end
    if strcmp(opts.nonfinite,'propagate'), location(bad) = NaN; end
    location = stat_restore(location,shape,order);
end
if nargout > 2
    scale = (sd.*s)*scaleMultiplier;
    if strcmp(opts.nonfinite,'propagate'), scale(bad) = NaN; end
    [scale,soverflow] = stat_overflow(scale,opts.overflow);
    scale = stat_restore(scale,shape,order);
end
if nargout > 3
    info = struct('method',opts.method,'nonfinite',opts.nonfinite, ...
        'count',stat_restore(n,shape,order),'constant',stat_restore(constant,shape,order), ...
        'zOverflow',stat_restore(zoverflow,shape,order),'scaleOverflow',stat_restore(soverflow,shape,order));
end
end
