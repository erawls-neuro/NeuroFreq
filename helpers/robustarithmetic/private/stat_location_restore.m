function y = stat_location_restore(mu, scale, origin)
% Convex-location estimates must remain inside [-scale,+scale].
y = bsxfun(@plus,bsxfun(@times,mu,scale),origin);
bad = ~isfinite(y);
if any(bad(:))
    fallback = bsxfun(@times,bsxfun(@plus,mu,origin./scale),scale);
    y(bad) = fallback(bad);
end
undefined = isnan(mu);
if isreal(y)
    y = bsxfun(@max,bsxfun(@min,y,scale),-scale);
else
    re = bsxfun(@max,bsxfun(@min,real(y),scale),-scale);
    im = bsxfun(@max,bsxfun(@min,imag(y),scale),-scale);
    y = complex(re,im);
end
y(undefined) = NaN;
end
