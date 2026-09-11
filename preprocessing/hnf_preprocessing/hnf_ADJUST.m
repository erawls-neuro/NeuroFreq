function [art, horiz, vert, blink, disc]=hnf_ADJUST(EEG,developmental)
% HNF_ADJUST  NF Adapted Automatic EEG artifact Detector with Joint Use of Spatial and Temporal features
%
% Inputs:
%   EEG            - current dataset structure or structure array (has to be epoched)
%   developmental  - is the dataset developmental (adjusted-ADJUST?)? 1 if
%                        yes
%
% Outputs:
%   art        - List of artifact ICs
%   horiz      - List of HEM ICs
%   vert       - List of VEM ICs
%   blink      - List of EB ICs
%   disc       - List of GD ICs
%
% Internal modification, 2026-09-10:
%   Skip only detectors whose required spatial features are unavailable.
%   Fit only thresholds used by eligible detectors; compute kurtosis only
%   for the classic blink detector. Preserve the existing decision rules.
%
if length(size(EEG.data))==3
    num_epoch=size(EEG.data,3);
else
    num_epoch=0;
end


% Eric Rawls modified 9/10/2026.
% Separate component count from the number of ICA training channels.
W = EEG.icaweights * EEG.icasphere;
nIC = size(W,1);
% Preserve the channel order used to train ICA.
icaChannels = EEG.icachansind(:).';
if isempty(icaChannels)
    icaChannels = 1:size(EEG.data,1);
end
EEG.icachansind = icaChannels;
% The product has COMPONENT rows, even when ICA is rectangular.
if isempty(EEG.icaact)
    icaData = reshape(EEG.data(icaChannels,:,:), ...
        numel(icaChannels), []);
    EEG.icaact = reshape(W * icaData, ...
        nIC, size(EEG.data,2), size(EEG.data,3));
    clear icaData
end


topografie=EEG.icawinv'; %computes IC topographies
for i=1:size(EEG.icawinv,2) % number of ICs
    ScalingFactor=norm(topografie(i,:));
    topografie(i,:)=topografie(i,:)/ScalingFactor;
    if length(size(EEG.data))==3
        EEG.icaact(i,:,:)=ScalingFactor*EEG.icaact(i,:,:);
    else
        EEG.icaact(i,:)=ScalingFactor*EEG.icaact(i,:);
    end
end
% Eric Rawls modified 9/10/2026.
% Mixing-map columns correspond to ICA training channels, in this order.
icaLocs = reshape(EEG.chanlocs(icaChannels), 1, []);
% Identify channels with usable coordinates.
% The original AND test excluded only wholly empty location records.
hasPosition = true(1,numel(icaLocs));
for el = 1:numel(icaLocs)
    loc = icaLocs(el);
    coordinates = {loc.X, loc.Y, loc.Z, loc.theta, loc.radius};
    hasPosition(el) = all(cellfun(@(v) ...
        isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v), ...
        coordinates));
end
if any(~hasPosition)
    warning('hnf_ADJUST:MissingPositions', ...
        ['Channels %s have incomplete location information and will ' ...
         'be excluded from spatial feature calculations.'], ...
        num2str(icaChannels(~hasPosition)));
end
% Apply EXACTLY the same subset to map columns and channel locations.
% Keep all component rows and retain the complete-map normalization.
spatialTopog = topografie(:,hasPosition);
spatialLocs = icaLocs(hasPosition);
% GDSF - Generic Discontinuity Spatial Feature
GDSF = [];
if numel(spatialLocs) >= 2
    GDSF = hnf_compute_GD(spatialTopog, spatialLocs, nIC);
else
    warning('hnf_ADJUST:MissingGDChannels', ...
        'GD detection skipped: fewer than two usable channel locations.');
end
% SED - Spatial Eye Difference
[SED,medie_left,medie_right] = ...
    hnf_compute_SED(spatialTopog, spatialLocs, nIC);
% SAD - Spatial Average Difference
[SAD,var_front,var_back,~,~] = ...
    hnf_compute_SAD(spatialTopog, spatialLocs, nIC);
%SVD - Spatial Variance Difference between front zone and back zone
diff_var=var_front-var_back;
% Empty SAD/SED outputs mean unavailable regions, not zero-valued features.
% VEM and classic blinks also require the signed eye means returned by SED.
doHEM = developmental~=1 && ~isempty(SED);
doVEM = ~isempty(SAD) && ~isempty(SED);
doBlink = developmental~=1 && doVEM;
doGD = ~isempty(GDSF);
needV = doHEM || doVEM || doGD;

% Classic blinks imply doVEM, so needV also covers their temporal work.
if needV
%epoch dynamic range, variance and kurtosis
if doBlink
    K=zeros(num_epoch,size(EEG.icawinv,2)); %kurtosis
end
Vmax=zeros(num_epoch,size(EEG.icawinv,2)); %variance
for i=1:size(EEG.icawinv,2) % number of ICs
    for j=1:num_epoch
        Vmax(j,i)=var(EEG.icaact(i,:,j));
        if doBlink
            K(j,i)=kurt(EEG.icaact(i,:,j));
        end
    end
end
%TK - Temporal Kurtosis
if doBlink
meanK=zeros(1,size(EEG.icawinv,2));
for i=1:size(EEG.icawinv,2)
    if num_epoch>=100
        meanK(1,i)=hnf_trim_and_mean(K(:,i));
    else
        meanK(1,i)=mean(K(:,i));
    end
end
end % doBlink
%MEV - Maximum Epoch Variance
maxvar=zeros(1,size(EEG.icawinv,2));
meanvar=zeros(1,size(EEG.icawinv,2));
for i=1:size(EEG.icawinv,2)
    if num_epoch>=100
        maxvar(1,i)=hnf_trim_and_max(Vmax(:,i)');
        meanvar(1,i)=hnf_trim_and_mean(Vmax(:,i)');
    else
        maxvar(1,i)=max(Vmax(:,i));
        meanvar(1,i)=mean(Vmax(:,i));
    end
end
% MEV in reviewed formulation:
nuovaV=maxvar./meanvar;
end % needV
% Thresholds computation
% Do not pass empty features to EM or fit unused developmental thresholds.
if doBlink
    soglia_K=hnf_EM(meanK);
end
if doHEM
    soglia_SED=hnf_EM(SED);
end
if doVEM
    soglia_SAD=hnf_EM(SAD);
end
if doGD
    soglia_GDSF=hnf_EM(GDSF);
end
if needV
    soglia_V=hnf_EM(nuovaV);
end
% vectors creation
horiz=[]; vert=[]; blink=[]; disc=[];
% Horizontal eye movements (HEM)
if developmental==1
    cwb=hnf_MARA_TF(EEG);
    horiz=hnf_beall_horizontal(cwb,spatialTopog,spatialLocs,nIC);
elseif doHEM
    horiz=intersect(intersect(find(SED>=soglia_SED),find(medie_left.*medie_right<0)),...
        (find(nuovaV>=soglia_V)));
end
% Vertical eye movements (VEM)
if doVEM
vert=intersect(intersect(find(SAD>=soglia_SAD),find(medie_left.*medie_right>0)),...
    intersect(find(diff_var>0),find(nuovaV>=soglia_V)));
end
% Eye Blink (EB)
if developmental==1
    blink=hnf_beall_blink(cwb,spatialTopog,spatialLocs,nIC);
elseif doBlink
    blink=intersect ( intersect( find(SAD>=soglia_SAD),find(medie_left.*medie_right>0) ) ,...
        intersect ( find(meanK>=soglia_K),find(diff_var>0) ));
end
% Generic Discontinuities (GD)
if doGD
disc=intersect(find(GDSF>=soglia_GDSF),find(nuovaV>=soglia_V));
end
%compute output variable
art = nonzeros( union (union(blink,horiz) , union(vert,disc)) )'; %artifact ICs
% now optionally do adjusted-ADJUST's rescue step
if developmental==1
    % loop through artifact ICs and remove any ICs that have bumps
    remove_indices=[];
    for ii = 1:length(art)
        for jj = 1:length(cwb) %components_with_bumps
            if art(ii) == cwb(jj) %components_with_bumps(jj)
                remove_indices = [remove_indices ii]; %#ok<AGROW>
            end
        end
    end
    art(remove_indices)=[];
    vert = intersect(vert,art);
    disc = intersect(disc,art);
end
