function Results = hnf_ept_TFCE_corr(data1, data2, chanlocs, varargin)
%
% ATTENTION: THIS FUNCTION IS NOT THE ORIGINAL ept_TFCE.m
% IT HAS BEEN MODIFIED BY ERIC RAWLS 2023.
% MANY ELEMENTS OF THE CODE WERE DELETED FOR THE NEUROFREQ USE CASE.
%
% ORIGINAL DOCUMENTATION FOLLOWS:
%
% ========================================================================
%
% EEG Permutation Test (ept) using Threshold-Free Cluster Enhancement
% Copyright(C) 2012  Armand Mensen (14.12.2010)
%
% [Description]
% This tests initially computes a r value for each channel and sample
% These values are then enhanced or reduced depending on the size of
% the T-value and the size of the cluster the channel belongs to at
% various thresholds...
%
% TFCE is essentially the sum of the clusters to the power of 0.5 multiplied
% by the height of the current threshold squared
%
% [Input]
% Participant data should be organised into two factors
%     e.g. Dataset 1 is controls while Dataset 2 is patients
% Electrode locations file created using eeglab is required to
% calculate channel neighbours
%
%
% data1 should be a "Participants x Channels x Samples" or a
%  "Participants x Channels x Frequencies X Samples" variable
% - Channels must be in the same order as the corresponding electrodes file
% - Samples should be in chronological order
% data2 must be a matrix with one dimension of size "Participants".
% ElecFile is not required. If not specified, must specify flag_ft = 1
%
%

%defaults
[eh, count, grid, neighbours, parallel, corrMode] = hnf_tfce_options(data1, chanlocs, varargin{:});
n = size(data1,1);
if n < 3 || size(data2,1) ~= n || any(~isfinite(data1(:))) || any(~isfinite(data2(:)))
    error('NeuroFreq:CorrelationData', 'Supply at least three finite, paired rows in data and predictors.');
end
dims = size(data1);
spatial = dims(2:end);
[Z,constantZ] = hnf_stats_zscore(data1(:,:));
[X,constantX] = hnf_stats_zscore(data2(:,:));
if any(constantZ) || any(constantX)
    error('NeuroFreq:UndefinedCorrelation', 'Pearson correlation is undefined for constant inputs.');
end
pointwise = strcmpi(corrMode,'pointwise') || ...
    (strcmpi(corrMode,'auto') && isequal(size(data1),size(data2)) && size(X,2)>1);
if pointwise
    if ~isequal(size(data1),size(data2))
        error('NeuroFreq:CorrelationMaps','Pointwise correlation requires two maps with identical axes.');
    end
    r = sum(X .* Z,1)/(n-1);
    corrMode = 'pointwise';
else
    if ~ismatrix(data2), error('NeuroFreq:CorrelationDesign','Predictor design must be a matrix.'); end
    r = X' * Z / (n-1);
    corrMode = 'predictors';
end
r = max(-1,min(1,r)); % Clamp floating-point excursions around exact +/-1.
observed = zeros(size(r));
null = zeros(count, size(r,1));
for k = 1:size(r,1)
    if pointwise, design = X; else, design = X(:,k); end
    map = squeeze(reshape(r(k,:), [1 spatial]));
    temp = hnf_tfce_transformation(map, neighbours, eh, grid);
    observed(k,:) = temp(:)';
    if parallel
        parfor b = 1:count
            null(b,k) = corr_draw(design, Z, spatial, neighbours, eh, grid);
        end
    else
        for b = 1:count
            null(b,k) = corr_draw(design, Z, spatial, neighbours, eh, grid);
        end
    end
end
p = zeros(size(observed));
for k = 1:size(r,1), p(k,:) = hnf_tfce_pvalues(observed(k,:), null(:,k)); end
shape = [size(r,1), spatial];
Results = struct('Obs',reshape(r,shape),'TFCE_Obs',reshape(observed,shape), ...
    'maxTFCE',sort(null,1),'P_Values',reshape(p,shape),'predictorMode',corrMode);
if size(r,1) == 1 % Preserve the original single-predictor shape.
    Results.Obs = squeeze(Results.Obs);
    Results.TFCE_Obs = squeeze(Results.TFCE_Obs);
    Results.P_Values = squeeze(Results.P_Values);
end
end

function maximum = corr_draw(x, Z, spatial, neighbours, eh, grid)
n = size(x,1);
permuted = x(randperm(n),:);
if size(x,2) == 1, coef = permuted' * Z / (n-1);
else, coef = sum(permuted .* Z,1) / (n-1); end
coef = max(-1,min(1,coef));
map = squeeze(reshape(coef,[1 spatial]));
enhanced = hnf_tfce_transformation(map, neighbours, eh, grid);
maximum = max(abs(enhanced(:)));
end
