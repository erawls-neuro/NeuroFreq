function [z,constant] = hnf_stats_zscore(x)
% Standardize finite observation columns without squaring large physical
% units. Center the deviations from an observed value first, so a large
% common offset does not absorb the mean of small representable changes.
x = hnf_stats_scale(double(x));
x = x - x(1,:);
x = x - mean(x,1);
s = std(x,0,1);
constant = s == 0;
z = x ./ s;
end
