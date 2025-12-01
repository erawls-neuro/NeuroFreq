function TF = nf_rmerp( TF, trlvec, scale )
% NF_RMERP    Removes the trial-averaged ERP from TF sets
%
% GENERAL
% -------
% Removes the trial-averaged ERP from TF sets. Reconstitutes analytic
% signal by combining power and phase, calculates trial average, and
% removes it, then converts the analytic signal back to power/phase.
%
% OUTPUT
% ------
% TF - structure output by tfUtility.m and tf_fun functions
%
% INPUT
% -----
% 1) TF - structure output by tfUtility.m and tf_fun functions
% 2) trlvec - vector describing which trials to average together. Each
%     unique value in the vector is taken to be a different trial type.
% 3) scale - scale single-trial ERP amplitudes according to the similarity
%     between single-trial EEG and the ERP? default 0 (no)
%
%


if nargin<3 || isempty(scale)
    scale=0;
end
if nargin<2 || isempty(trlvec)
    disp('No trial vector supplied - averaging all trials together');
    trlvec = ones(1, size(TF.power,ndims(TF.power)));
end
if nargin<1 || isempty(TF)
    error('at least a TF structure is required input');
end

%check if already averaged
if TF.ntrls==1
    error('TF is already averaged (ntrls = 1)');
end
if TF.nsensor==1
    flagsens=1;
else
    flagsens=0;
end

%get trial info
conds = unique(trlvec);

%loop thru conditions
for n=1:numel(conds)
    
    if flagsens==0
        %get only condition-specific activity (square root of power)
        condP = sqrt(TF.power(:,:,:,find(trlvec==conds(n))));
        %get phase
        condPh = TF.phase(:,:,:,find(trlvec==conds(n)));
        % Combine amplitude and phase to create the analytic signal
        analyticSignal = condP .* exp(1i * condPh);
        % average analytic over trials to get ERP, repmat to full size
        erp = repmat(squeeze(mean(analyticSignal,4)),[1 1 1 numel(find(trlvec==conds(n)))]);
        %optionally scale single-trial ERPs according to their similarity
        if scale==1
            %get dot product of each trial with ERP
            dtp = squeeze(dot(analyticSignal,erp,3));
            %normalize dot products so they are positive and average over
            %trials to 1 (that is, they recapitulate the average ERP)
            dtp = dtp-min(dtp(:))+eps;
            dtp_m = mean(dtp,3);
            dtp = dtp./repmat(dtp_m,[1 1 size(dtp,3)]);
            dtp_n = zeros(size(dtp,1),size(dtp,2),1,size(dtp,3));
            dtp_n(:,:,1,:) = dtp;
            dtp = repmat(dtp_n,[1 1 size(erp,3) 1]);
            %multiply scaled dot products by ERP
            erp = erp .* dtp;
        end
        % subtract "ERP" (actually filtered but mathematically equal)
        analyticSignal = analyticSignal - erp;
        %add them to power/phase in structure
        TF.erprem.erprempow(:,:,:,find(trlvec==conds(n))) = abs(analyticSignal).^2;
        TF.erprem.erpremphase(:,:,:,find(trlvec==conds(n))) = angle(analyticSignal);
        TF.erprem.erppow(:,:,:,find(trlvec==conds(n))) = abs(erp).^2;
    else
        %get only condition-specific power
        condP = sqrt(TF.power(:,:,find(trlvec==conds(n))));
        %get phase
        condPh = TF.phase(:,:,find(trlvec==conds(n)));
        % Combine amplitude and phase to create the analytic signal
        analyticSignal = condP .* exp(1i * condPh);
        % average analytic over trials to get ERP, repmat to full size
        erp = repmat(squeeze(mean(analyticSignal,3)),[1 1 numel(find(trlvec==conds(n)))]);
        %optionally scale single-trial ERPs according to their similarity
        if scale==1
            %get dot product of each trial with ERP
            dtp = squeeze(dot(analyticSignal,erp,2));
            %normalize dot products so they are positive and average over
            %trials to 1 (that is, they will recapitulate the average ERP)
            dtp = dtp-min(dtp(:)+eps);
            dtp_m = mean(dtp,2);
            dtp = dtp./repmat(dtp_m,[1 size(dtp,2)]);
            dtp_n = zeros(size(dtp,1),1,size(dtp,2));
            dtp_n(:,:,1,:) = dtp;
            dtp = repmat(dtp_n,[1 size(erp,2) 1]);
            %multiply scaled dot products by ERP
            erp = erp .* dtp;
        end
        % subtract "ERP" (actually filtered but mathematically equal)
        analyticSignal = analyticSignal - erp;
        %add them to power/phase in structure
        TF.erprem.erprempow(:,:,find(trlvec==conds(n))) = abs(analyticSignal).^2;
        TF.erprem.erpremphase(:,:,find(trlvec==conds(n))) = angle(analyticSignal);
        TF.erprem.erppow(:,:,find(trlvec==conds(n))) = abs(erp).^2;
    end
    
end


end











