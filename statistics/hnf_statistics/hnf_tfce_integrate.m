function enhanced = hnf_tfce_integrate(stat, neighbours, exponents)
% HNF_TFCE_INTEGRATE  Signed EPT-style TFCE on chan x sample or chan x freq x time.
% A point connects to the explicitly listed channel neighbours at the same
% or either adjacent bin of each other axis. As in ept_ChN2, zero ends a
% neighbour row and the row may contain its own channel. Each point counts
% only once. Keep both signs independent across the entire integration.
%
% Use a sparse connected-component decomposition for symmetric montages.
% For large graphs and directed neighbour lists, use a bounded-memory queue.
% Preserve the interpreted implementation for platforms without a native
% build and for direct numerical comparisons under the scoped MEX policy.
if exist('hnf_mex_available','file') == 2 && hnf_mex_available('hnf_tfce_clusters_mex')
    try
        enhanced = hnf_tfce_clusters_mex(full(double(stat)),full(double(neighbours)),full(double(exponents)));
        return
    catch err
        if ~hnf_mex_load_error(err), rethrow(err); end
    end
end
peak = max(abs(stat(:)));
enhanced = zeros(size(stat));
if isempty(peak) || peak == 0, return; end
% Ordinary EEG maps use direct midpoint arithmetic. If an extent/height
% factor or their product can leave normal double range, accumulate log
% magnitudes and apply the sign at the end. Small per-step increments can
% then combine to give a representable cluster integral.
logPeak = log(peak);
extentBound = exponents(1)*log(numel(stat));
heightBound = exponents(2)*max(abs(logPeak),abs(logPeak-log(100)));
useLog = extentBound + heightBound + max(abs(logPeak),abs(logPeak-log(50))) > -log(realmin);
if useLog
    stat = stat / peak;
    enhanced(:) = -Inf;
    peak = 1;
end
dims = [size(stat,1) size(stat,2) size(stat,3)];
channels = cell(dims(1),1);
for c = 1:dims(1)
    row = neighbours(c,:);
    stop = find(row == 0,1);
    if ~isempty(stop), row = row(1:stop-1); end
    channels{c} = unique(row(:));
end
from = repelem((1:dims(1))',cellfun(@numel,channels));
chanGraph = sparse(from,vertcat(channels{:}),true,dims(1),dims(1));

% Limit the explicit graph to about two million nonzeros. The queue path
% has O(number of voxels + channel-neighbour entries) memory use.
nLinks = nnz(chanGraph) * (3*dims(2)-2) * (3*dims(3)-2);
useGraph = issymmetric(chanGraph) && nLinks <= 2e6;
if useGraph
    graph = kron(sample_links(dims(3)),kron(sample_links(dims(2)),chanGraph));
    % Structural diagonal ensures dmperm returns connected components.
    graph = spones(graph + speye(numel(stat)));
else
    queue = zeros(numel(stat),1);
end
delta = peak/50;
for k = 1:50
    threshold = delta/2 + (k-1)*delta;
    for signValue = [1 -1]
        active = signValue*stat >= threshold;
        indices = find(active);
        if isempty(indices), continue; end
        if useLog
            logWeight = exponents(2)*(logPeak+log(threshold)) + (logPeak-log(50));
        else
            heightWeight = signValue * threshold^exponents(2) * delta;
        end
        if useGraph
            [order,~,bounds] = dmperm(graph(indices,indices));
            sizes = diff(bounds);
            nodes = indices(order);
            if useLog
                weights = repelem(exponents(1)*log(sizes)+logWeight,sizes);
                enhanced(nodes) = log_add(enhanced(nodes),reshape(weights,size(nodes)));
            else
                weights = repelem(sizes.^exponents(1),sizes);
                enhanced(nodes) = enhanced(nodes) + reshape(weights,size(nodes))*heightWeight;
            end
        else
            for seed = reshape(indices,1,[])
                if ~active(seed), continue; end
                head = 1; tail = 1; queue(1) = seed; active(seed) = false;
                while head <= tail
                    [c,f,t] = ind2sub(dims,queue(head));
                    fRange = max(1,f-1):min(dims(2),f+1);
                    tRange = max(1,t-1):min(dims(3),t+1);
                    offsets = dims(1)*((fRange(:)-1) + dims(2)*(tRange-1));
                    near = channels{c} + reshape(offsets,1,[]);
                    new = near(active(near));
                    queue(tail+(1:numel(new))) = new;
                    active(new) = false;
                    tail = tail+numel(new); head = head+1;
                end
                nodes = queue(1:tail);
                if useLog
                    enhanced(nodes) = log_add(enhanced(nodes),exponents(1)*log(tail)+logWeight);
                else
                    enhanced(nodes) = enhanced(nodes) + tail^exponents(1)*heightWeight;
                end
            end
        end
    end
end
if useLog, enhanced = sign(stat).*exp(enhanced); end
end

function total = log_add(existing,increment)
larger = max(existing,increment);
correction = exp(min(existing,increment)-larger);
correction(~isfinite(larger)) = 0;
total = larger + log1p(correction);
end

function links = sample_links(n)
if n == 1, links = speye(1); else, links = spdiags(ones(n,3),[-1 0 1],n,n); end
end
