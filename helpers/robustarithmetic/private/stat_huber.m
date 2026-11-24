function [mu, iterations, converged, zeroScale] = stat_huber(raw, valid, opts)
% Fixed-MAD Huber IRLS via bounded scores; no underflowing weight-times-outlier.
[s,mu] = stat_rawmad(raw,valid);
count = sum(valid,1);
zeroScale = s == 0 & count > 0;
active = count > 0 & ~zeroScale;
converged = ~active; iterations = zeros(size(count));
cutoff = opts.tuning*1.482602218505602;
denScale = s; denScale(~active) = 1;
if isempty(raw), return; end
inputScale = max(abs(raw),[],1); inputScale(inputScale == 0) = 1;
for k = 1:opts.maxiter
    if ~any(active), break; end
    d = bsxfun(@minus,raw,mu);
    bad = isinf(d);
    d = bsxfun(@rdivide,d,denScale);
    if any(bad(:))
        fallback = bsxfun(@minus,bsxfun(@rdivide,raw,denScale),mu./denScale);
        d(bad) = fallback(bad);
    end
    score = max(-cutoff,min(cutoff,d)); score(~valid) = 0;
    denominator = abs(d); denominator(denominator == 0) = 1;
    w = min(1,cutoff./denominator); w(d == 0) = 1; w(~valid) = 0;
    step = sum(score,1)./max(sum(w,1),realmin(class(raw)));
    candidate = mu + denScale.*step;
    overflow = isinf(candidate);
    if any(overflow)
        fallback = inputScale.*(mu./inputScale + (denScale./inputScale).*step);
        candidate(overflow) = fallback(overflow);
    end
    candidate = max(-inputScale,min(inputScale,candidate));
    delta = abs(candidate-mu);
    tolerance = (opts.tol*1.482602218505602)*s + eps(class(raw))*abs(mu);
    done = active & delta <= tolerance;
    mu(active) = candidate(active); iterations(active) = k;
    converged(done) = true; active(done) = false;
end
mu(count == 0) = NaN;
end
