function [s, mu] = stat_residual_scale(a, valid, opts)
count = sum(valid,1);
if strcmp(opts.method,'mad')
    mu = stat_quantile(a,valid,.5);
    s = stat_quantile(abs(bsxfun(@minus,a,mu)),valid,.5)*opts.consistency;
else
    if strcmp(opts.method,'winsor')
        cuts = stat_quantile(a,valid,[opts.trim,1-opts.trim]);
        a = bsxfun(@max,bsxfun(@min,a,cuts(2,:)),cuts(1,:));
        a(~valid) = 0;
    end
    mu = sum(a,1)./max(count,1);
    r = bsxfun(@minus,a,mu); r(~valid) = 0;
    if isempty(r)
        s = nan(1,size(a,2),'like',real(a)); return
    end
    peak = max(abs(r),[],1); peak(peak == 0) = 1;
    r = bsxfun(@rdivide,r,peak);
    if opts.normalization == 0, denominator = max(count-1,1);
    else, denominator = max(count,1); end
    s = peak.*sqrt(sum(abs(r).^2,1)./denominator);
end
s(count == 0) = NaN; mu(count == 0) = NaN;
end
