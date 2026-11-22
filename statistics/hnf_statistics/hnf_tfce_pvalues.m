function p = hnf_tfce_pvalues(observed, maxima)
% Two-sided spatial FWER for ONE tested effect. Count ties in the tail.
% Include the observed arrangement in the Monte Carlo reference set.
if isempty(maxima) || any(~isfinite(maxima(:)))
    error('NeuroFreq:InvalidNull', 'Permutation maxima must be nonempty and finite.');
end
% For larger maps, search one sorted null vector in native code. All
% other input types and installations retain the MATLAB implementation.
if numel(observed) >= 128 && isa(observed,'double') && isreal(observed) && ~issparse(observed) && ...
        isa(maxima,'double') && isreal(maxima) && ~issparse(maxima) && ...
        exist('hnf_mex_available','file') == 2 && hnf_mex_available('hnf_tfce_tail_mex')
    try
        p = hnf_tfce_tail_mex(observed,maxima);
        return
    catch failure
        if ~hnf_mex_load_error(failure), rethrow(failure); end
        warning('NeuroFreq:NativeTailLoad','%s Using MATLAB for TFCE tail counts.',failure.message);
    end
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
