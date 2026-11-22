function p = hnf_tfce_pvalues(observed, maxima)
% Two-sided spatial FWER for ONE tested effect. Count ties in the tail.
% Include the observed arrangement in the Monte Carlo reference set.
if isempty(maxima) || any(~isfinite(maxima(:)))
    error('NeuroFreq:InvalidNull', 'Permutation maxima must be nonempty and finite.');
end
% Sorting the two vectors once avoids a permutations-by-voxels array.
[levels, order] = sort(abs(observed(:)));
null = sort(maxima(:));
n = numel(null);
p = nan(size(observed));
first = 1;
for j = 1:numel(levels)
    if ~isfinite(levels(j)), continue; end
    % Equivalent label partitions can sum in different orders. Count tiny
    % relative roundoff differences as ties (as in scipy.stats.permutation_test).
    lower = levels(j) * (1-100*eps);
    while first <= n && null(first) < lower
        first = first + 1;
    end
    p(order(j)) = (1 + n - first + 1) / (n + 1);
end
end
