function [mu, iterations, converged, zeroScale, usedCount] = stat_location(a, valid, opts)
count = sum(valid,1); usedCount = count;
iterations = zeros(1,size(a,2));
converged = true(1,size(a,2)); zeroScale = false(1,size(a,2));
switch opts.method
    case 'mean'
        mu = sum(a,1)./max(count,1);
    case 'median'
        mu = stat_quantile(a,valid,.5);
    case {'trim','winsor'}
        if strcmp(opts.method,'winsor')
            cut = stat_quantile(a,valid,[opts.trim,1-opts.trim]);
            a = bsxfun(@max,bsxfun(@min,a,cut(2,:)),cut(1,:));
            a(~valid) = 0;
            mu = sum(a,1)./max(count,1);
        else
            a(~valid) = Inf; a = sort(a,1);
            discard = floor(count*opts.trim);
            keep = bsxfun(@gt,(1:size(a,1))',discard) & ...
                bsxfun(@le,(1:size(a,1))',count-discard);
            usedCount = count-2*discard;
            a(~keep) = 0;
            mu = sum(a,1)./max(usedCount,1);
        end
    case 'huber'
        mu = stat_quantile(a,valid,.5);
        residual = abs(bsxfun(@minus,a,mu));
        s = stat_quantile(residual,valid,.5)*cast(1.482602218505602,'like',a);
        zeroScale = s == 0 & count > 0;
        active = count > 0 & ~zeroScale;
        converged(active) = false;
        for k = 1:opts.maxiter
            if ~any(active), break; end
            residual = abs(bsxfun(@minus,a,mu));
            threshold = opts.tuning*s;
            denominator = residual; denominator(denominator == 0) = 1;
            w = min(1,bsxfun(@rdivide,threshold,denominator));
            w(residual == 0) = 1;
            w(~valid) = 0;
            candidate = sum(w.*a,1)./max(sum(w,1),realmin(class(a)));
            delta = abs(candidate-mu);
            mu(active) = candidate(active);
            iterations(active) = k;
            tolerance = opts.tol.*s + eps(class(a)).*abs(mu);
            done = active & delta <= tolerance;
            converged(done) = true; active(done) = false;
        end
end
mu(count == 0) = NaN;
end
