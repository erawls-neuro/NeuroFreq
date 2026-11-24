function z = hnf_pac_permutation_z(phase, weights, observed, nperm, metric, nbin, useMex, cross)
% HNF_PAC_PERMUTATION_Z Score bounded batches using the caller's RNG.
% First axis = trials for ERPAC or time for MI/MVL. For ERPAC, time columns
% of phase repeat for every high frequency and amplitude channel in weights.
% Keep the established per-draw calls and order. A single large randi call
% can consume a different number of random values on some MATLAB runtimes.
samples = size(weights,1);
weights = reshape(double(weights),samples,[]);
phase = reshape(double(phase),samples,[]);
series = size(weights,2);
% MVL scales linearly with its weights. The z score is dimensionless.
% Rescale only extreme columns before the permutation variance squares them.
if strcmp(metric,'mvl')
    scales = max(abs(weights),[],1);
    extreme = isfinite(scales) & (scales > 1e150 | (scales > 0 & scales < 1e-150));
    if any(extreme)
        weights(:,extreme) = weights(:,extreme)./scales(extreme);
        observed(extreme) = observed(extreme)./reshape(scales(extreme),size(observed(extreme)));
    end
end
erpac = strcmp(metric,'erpac');
if erpac
    nativeName = 'hnf_erpac_perm_mex';
    indexRows = samples;
else
    nativeName = 'hnf_pac_shift_scores_mex';
    indexRows = series;
end
native = useMex && hnf_mex_available(nativeName);
% Limit transient indices and scores to roughly 4 MiB each, apart from a
% single very wide draw. Within-channel reductions retain MATLAB mean/std.
% Cross-channel reductions retain the original online update order.
block = max(1,min(nperm,floor(524288/max(series,indexRows))));
if cross
    mu = zeros(1,series);
    m2 = zeros(1,series);
else
    null = zeros(nperm,series);
end
for first = 1:block:nperm
    last = min(first+block-1,nperm);
    count = last-first+1;
    indices = zeros(indexRows,count);
    for p = 1:count
        if erpac
            indices(:,p) = randperm(samples).';
        elseif samples > 1
            indices(:,p) = randi([1,samples-1],1,series).';
        end
    end
    if native
        try
            if erpac
                scores = hnf_erpac_perm_mex(phase,weights,indices);
            else
                scores = hnf_pac_shift_scores_mex(phase,weights,indices,double(nbin));
            end
        catch err
            if ~hnf_mex_load_error(err), rethrow(err); end
            warning('NeuroFreq:PACBatchFallback','Using the fallback for %s: %s',nativeName,err.message);
            native = false; % Reuse exactly the sampled indices below.
        end
    end
    if ~native
        scores = explicitScores(phase,weights,indices,metric,nbin,useMex);
    end
    if cross
        for p = 1:count
            delta = scores(p,:) - mu;
            mu = mu + delta ./ (first+p-1);
            m2 = m2 + delta .* (scores(p,:)-mu);
        end
    else
        null(first:last,:) = scores;
    end
end
if cross
    sd = sqrt(m2 ./ max(1,nperm-1));
else
    mu = mean(null,1);
    sd = std(null,0,1);
end
sd(sd <= 8*eps(abs(mu))) = NaN;
z = reshape((observed(:).'-mu)./sd,size(observed));
end

function scores = explicitScores(phase,weights,indices,metric,nbin,useMex)
% A checkout can have only the earlier, single-permutation MEX kernels.
% Keep using those when installed; otherwise use the MATLAB reference.
[n,series] = size(weights);
count = size(indices,2);
scores = zeros(count,series);
if strcmp(metric,'erpac')
    if ~useMex || ~hnf_mex_available('hnf_circ_corrcl_mex')
        scores = hnf_circ_corrcl_vec(phase,weights,indices);
        return
    end
    expanded = repmat(phase,1,series/size(phase,2));
    for p = 1:count
        scores(p,:) = hnf_clcorr_kernel(expanded,weights(indices(:,p),:),useMex);
    end
    return
end
% Gather columns directly instead of dispatching circshift for each trial.
% Chunk very wide arrays so indexing needs only about 4 MiB of workspace.
width = max(1,floor(524288/n));
sample = (0:n-1).';
for first = 1:width:series
    cols = first:min(first+width-1,series);
    offset = n*(cols-1)+1;
    for p = 1:count
        gather = bsxfun(@plus,mod(bsxfun(@minus,sample,indices(cols,p).'),n),offset);
        shifted = reshape(weights(gather),n,numel(cols));
        if strcmp(metric,'mi')
            scores(p,cols) = hnf_mi_kernel(phase(:,cols),shifted,nbin,useMex);
        else
            scores(p,cols) = hnf_mvl_kernel(phase(:,cols),shifted,useMex);
        end
    end
end
end
