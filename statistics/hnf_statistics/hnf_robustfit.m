function b = hnf_robustfit(X,y)
% Coefficient-only default robustfit used by nf_stregress. Octave releases
% without robustfit use a Tukey bisquare M-estimate with leverage adjustment.
if ~exist('OCTAVE_VERSION','builtin') || exist('robustfit','file') == 2
    b = robustfit(X,y);
    return
end
y = y(:);
X = [ones(size(X,1),1), X];
valid = all(isfinite(X),2) & isfinite(y);
X = X(valid,:);
y = y(valid);
[n,p] = size(X);
if n <= p || rank(X) < p
    error('NeuroFreq:RobustDesign','Robust fit needs a full-rank design and residual degrees of freedom.');
end
[Q,~] = qr(X,0);
adjustment = sqrt(1 - min(.9999,sum(Q.^2,2)));
b = X \ y;
for iter = 1:50
    resid = (y - X*b) ./ adjustment;
    ordered = sort(abs(resid));
    % Removing the p-1 smallest residuals accounts for the fitted columns.
    scale = max(median(ordered(p:end))/.6745, 1e-6 * std(y));
    if scale == 0, return; end
    u = resid ./ (4.685*scale);
    % Square root of the bisquare weight is sufficient for weighted LS.
    weight = max(0,1-u.^2);
    design = X .* weight;
    if rank(design) < p
        error('NeuroFreq:RobustDesign','Bisquare weights leave the estimated design rank deficient.');
    end
    next = design \ (weight .* y);
    converged = all(abs(next-b) <= sqrt(eps) .* max(abs(next),abs(b)));
    b = next;
    if converged, break; end
end
end
