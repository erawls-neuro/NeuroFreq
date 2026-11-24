function q = stat_quantile(a, valid, p)
% Hyndman-Fan type 7: h = 1 + (n-1)*p. No Statistics Toolbox required.
[n,m] = size(a);
q = nan(numel(p),m,'like',a);
if n == 0 || m == 0, return; end
a(~valid) = Inf;
a = sort(a,1);
count = sum(valid,1);
for k = 1:numel(p)
    h = 1 + max(count-1,0)*p(k);
    lo = floor(h); hi = ceil(h); t = cast(h-lo,'like',a);
    left = a(lo + (0:m-1)*n);
    right = a(hi + (0:m-1)*n);
    v = hnf_roblerp(left,right,t);
    same = left == right;
    v(same) = left(same);
    v = max(left,min(right,v));
    v(count == 0) = NaN;
    q(k,:) = v;
end
end
