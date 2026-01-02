function TFcpm = nf_cpm( TF, trlvec )
% NF_CPM    Calculates concurrent phaser method from TF structures
%
% GENERAL
% -------
% Calculates concurrent phaser method from TF structures. Tf structures
% must be calculated from nf_tftransform or from any of the nf tf 
% functions. Power is averaged in order to separately extract activity
% phase-locked and non-phase-locked to stimulus presentation. The CPM
% is calculated according to the reference:
%
% Singhal S, Ghosh P, Kumar N, Banerjee A. Parametric separation of 
%   phase-locked and non-phase-locked activity. J Neurophysiol. 2023 
%   Jan 1;129(1):199-210. doi: 10.1152/jn.00467.2022. Epub 2022 Dec 21. 
%   PMID: 36541609.
%
% Note: this function only operates on TF structures containing a phase
% field because we must calculate an analytic signal prior to CPM method.
%
%
% OUTPUT
% ------
% TFcpm - structure output by tftransform.m and tf_fun functions. Contains 
%   a new field .cpm, containing separate phase-locked and non-phase-locked 
%   power estimates.
%
% INPUT
% -----
% 1) TF - structure output by nf_tftransform.m and tf_fun functions
% 2) trlvec - vector describing which trials to average together. Each 
%     unique value in the vector is taken to be a different trial type.
%
%


%argument handling
if nargin<2 || isempty(trlvec)
    disp('No trial vector supplied - averaging all trials together');
    trlvec = ones(1, size(TF.power,ndims(TF.power)));
end
if nargin<1 || isempty(TF)
    error('at least a TF structure is required input');
end
%dimension figuring
if TF.ntrls==1
    error('multiple trials required for CPM.');
end
if TF.nsensor==1
    flagsens=1;
else
    flagsens=0;
end

%first average a version of the TF - this way we have total power
TFcpm = nf_avebase( TF, 'none', [], trlvec, 0 );

%get trial info
conds = unique(trlvec);

%preallocate cpm
if flagsens~=1
    TFcpm.cpm.PLmean = zeros(size(TFcpm.power,1),size(TFcpm.power,2),size(TFcpm.power,3),numel(conds));
    TFcpm.cpm.PLvar = zeros(size(TFcpm.power,1),size(TFcpm.power,2),size(TFcpm.power,3),numel(conds));
    TFcpm.cpm.NPLmean = zeros(size(TFcpm.power,1),size(TFcpm.power,2),size(TFcpm.power,3),numel(conds));
    TFcpm.cpm.NPLvar = zeros(size(TFcpm.power,1),size(TFcpm.power,2),size(TFcpm.power,3),numel(conds));
    TFcpm.cpm.alpha = zeros(size(TFcpm.power,1),size(TFcpm.power,2),size(TFcpm.power,3),numel(conds));
else
    TFcpm.cpm.PLmean = zeros(size(TFcpm.power,1),size(TFcpm.power,2),numel(conds));
    TFcpm.cpm.PLvar = zeros(size(TFcpm.power,1),size(TFcpm.power,2),numel(conds));
    TFcpm.cpm.NPLmean = zeros(size(TFcpm.power,1),size(TFcpm.power,2),numel(conds));
    TFcpm.cpm.NPLvar = zeros(size(TFcpm.power,1),size(TFcpm.power,2),numel(conds));
    TFcpm.cpm.alpha = zeros(size(TFcpm.power,1),size(TFcpm.power,2),numel(conds));
end

for n=1:numel(conds)

    %get current (condition-specific) power/phase
    if flagsens~=1
        currPow = TF.power(:,:,:,find(trlvec==conds(n)));
        currPhase = TF.phase(:,:,:,find(trlvec==conds(n)));
    else
        currPow = TF.power(:,:,find(trlvec==conds(n)));
        currPhase = TF.phase(:,:,find(trlvec==conds(n)));
    end
    % Convert power to amplitude (square root of power)
    amplitude = sqrt(currPow);
    % Combine amplitude and phase to create the analytic signal
    analyticSignal = amplitude .* exp(1i * currPhase);
    %get real/imaginary components
    R_array = 2 * real( analyticSignal );
    I_array = -2 * imag( analyticSignal );
    
    if flagsens~=1
        % estimate alpha
        est_alpha = atan2(mean(R_array,4), mean(I_array,4));
        %calculate CPM
        z = (R_array .* repmat(cos(est_alpha), [1 1 1 size(R_array,4)])) ...
            - (I_array .* repmat(sin(est_alpha), [1 1 1 size(R_array,4)]));
        zz = z.*(z > 0);
        zz = sum(zz,4)./sum(zz~=0,4);
        %get phase-locked
        Ar_mean_est = squeeze(mean(I_array,4)) ./ cos(est_alpha);
        Ar_var_est = (mean(I_array.^2,4) - mean(R_array.^2,4)) ./ cos(2 * est_alpha) - Ar_mean_est.^2;
        %get non-phase-locked
        Br_mean_est = pi * 0.5 * zz;
        Br_var_est = 2 * (mean(R_array.^2,4) .* cos(est_alpha).^2 - mean(I_array.^2,4) .* sin(est_alpha).^2) ./ cos(2 * est_alpha) - Br_mean_est.^2;
        %add to preallocated output
        TFcpm.cpm.PLmean(:,:,:,n) = Ar_mean_est;
        TFcpm.cpm.PLvar(:,:,:,n) = Ar_var_est;
        TFcpm.cpm.NPLmean(:,:,:,n) = Br_mean_est;
        TFcpm.cpm.NPLvar(:,:,:,n) = Br_var_est;
        TFcpm.cpm.alpha(:,:,:,n) = est_alpha;
    else
        % estimate alpha
        est_alpha = atan2(mean(R_array,3), mean(I_array,3));
        %calculate CPM
        z = (R_array .* repmat(cos(est_alpha), [1 1 size(R_array,3)])) ...
            - (I_array .* repmat(sin(est_alpha), [1 1 size(R_array,3)]));
        zz = z.*(z > 0);
        zz = sum(zz,3)./sum(zz~=0,3);
        %get phase-locked
        Ar_mean_est = squeeze(mean(I_array,3)) ./ cos(est_alpha);
        Ar_var_est = (mean(I_array.^2,3) - mean(R_array.^2,3)) ./ cos(2 * est_alpha) - Ar_mean_est.^2;
        %get non-phase-locked
        Br_mean_est = pi * 0.5 * zz;
        Br_var_est = 2 * (mean(R_array.^2,3) .* cos(est_alpha).^2 - mean(I_array.^2,3) .* sin(est_alpha).^2) ./ cos(2 * est_alpha) - Br_mean_est.^2;
        %add to preallocated output
        TFcpm.cpm.PLmean(:,:,n) = Ar_mean_est;
        TFcpm.cpm.PLvar(:,:,n) = Ar_var_est;
        TFcpm.cpm.NPLmean(:,:,n) = Br_mean_est;
        TFcpm.cpm.NPLvar(:,:,n) = Br_var_est;
        TFcpm.cpm.alpha(:,:,n) = est_alpha;
    end
    
end

end



