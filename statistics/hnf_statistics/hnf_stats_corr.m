function r = hnf_stats_corr(Y, X, method)
% Cross-correlations between measurement and predictor columns.
% nf_stregress passes real finite data; tied observations share average ranks.
if strcmpi(method,'pearson')
    % Both inputs have already been checked for finite real observations.
    r = hnf_stats_zscore(Y).' * hnf_stats_zscore(X) / (size(X,1)-1);
    % Leave undefined constant-column correlations as NaN.
    r(r > 1) = 1;
    r(r < -1) = -1;
    return
end
if ~exist('OCTAVE_VERSION','builtin')
    r = corr(Y,X,'Type',method);
    return
end
switch lower(method)
    case 'spearman'
        r = corr(ranks(Y,1),ranks(X,1));
    case 'kendall'
        % Tau-b counts each observation pair once. Process blocks of voxels
        % to avoid the N-by-N-by-voxel allocation in Octave's kendall().
        r = zeros(size(Y,2),size(X,2));
        for first = 1:128:size(Y,2)
            cols = first:min(first+127,size(Y,2));
            numer = zeros(numel(cols),size(X,2));
            pairsY = zeros(numel(cols),1);
            pairsX = zeros(1,size(X,2));
            for i = 1:size(Y,1)-1
                dY = sign(Y(i+1:end,cols) - Y(i,cols));
                dX = sign(X(i+1:end,:) - X(i,:));
                numer = numer + dY.' * dX;
                pairsY = pairsY + sum(dY ~= 0,1).';
                pairsX = pairsX + sum(dX ~= 0,1);
            end
            r(cols,:) = numer ./ sqrt(pairsY * pairsX);
        end
    otherwise
        error('NeuroFreq:Correlation','Unknown correlation method %s.',method);
end
end
