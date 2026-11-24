function [x,scale] = hnf_stats_scale(x,scale)
% Equilibrate measurement columns with exact powers of two. Common positive
% scaling preserves correlations, t statistics, and the sign/permutation null.
if nargin < 2
    magnitude = max(abs(x),[],1);
    [~,exponent] = log2(magnitude);
    scale = pow2(exponent - 1);
    scale(magnitude == 0) = 1;
end
x = x ./ scale;
end
