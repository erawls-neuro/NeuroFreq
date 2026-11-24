function b = hnf_stats_restore_beta(b,yScale,xScale)
% Restore beta in measurement/predictor units. The quotient of scale factors
% alone can overflow even when the coefficient is representable.
[fraction,be] = log2(b);
[~,ye] = log2(yScale);
[~,xe] = log2(xScale(:));
b = pow2(fraction,be + ye - xe);
end
