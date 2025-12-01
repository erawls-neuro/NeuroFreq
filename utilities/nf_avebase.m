function TF = nf_avebase( TF, method, blTimes, trlvec, avmode )
% NF_AVEBASE    Averages/baselines TF structures
%
% GENERAL
% -------
% Averages/baselines TF structures from tftransform or from tf_fun. 
% Power is averaged and optionally baseline corrected according to either 
% dB, %-change, or z-score. If the signal contains phase, phase is 
% processed into ITC.
%
% OUTPUT
% ------
% TF - structure output by nf_tftransform.m and tf functions (averaged)
%
% INPUT
% -----
% 1) TF - structure output by tfUtility.m and tf_fun functions
% 2) method - 'none', 'zscore', 'db', 'percent'
% 3) blTimes - [min max] times for baseline, in ms
% 4) trlvec - vector describing which trials to average together. Each 
%     unique value in the vector is taken to be a different trial type.
% 5) avmode - 'mean' or 'median'
%
%


if nargin<5 || isempty(avmode)
    disp('No averaging mode supplied - using mean');
    avmode = 'mean';
end
if nargin<4 || isempty(trlvec)
    disp('No trial vector supplied - averaging all trials together');
    trlvec = ones(1, size(TF.power,ndims(TF.power)));
end
if nargin<3 || isempty(blTimes)
    disp('no times provided, using all times before 0');
    blTimes = find( TF.times<=0 );
else
    blTimes = find( (TF.times>=blTimes(1)) + (TF.times<=blTimes(end)) == 2);
end
if nargin<2 || isempty(method)
    method='none';
else
    method=lower(method);
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

%preallocate power/phase
sz = size(TF.power);
sz=sz(1:end-1);
newpower = zeros([sz,numel(conds)]);
%if we have phase
if isfield(TF,'phase')
    newphase = zeros([sz,numel(conds)]);
end
%if tf is parameterized
if isfield(TF,'SPRiNT') %parameterized
    newosc = zeros([sz,numel(conds)]);
    newap = zeros([sz,numel(conds)]);
end
%if erp was removed
if isfield(TF,'erprem') %erp-removed
    newerprempow = zeros([sz,numel(conds)]);
    newerppow = zeros([sz,numel(conds)]);
    newerpremphase = zeros([sz,numel(conds)]);
end
%loop thru conditions
if isfield(TF,'behavior')
    newBeh = TF.behavior;
    newBeh(numel(conds)+1:end)=[];
end
for n=1:numel(conds)
    if flagsens==0
        %get only condition-specific power
        condP = TF.power(:,:,:,find(trlvec==conds(n)));
        %average power over trials
        newpower(:,:,:,n) = corrP(condP,blTimes,method,flagsens,avmode);
        %get parameterized
        if isfield(TF,'SPRiNT') %parameterized
            condap = TF.SPRiNT.ap_power(:,:,:,find(trlvec==conds(n)));
            condosc = TF.SPRiNT.osc_power(:,:,:,find(trlvec==conds(n)));
            newap(:,:,:,n) = corrP(condap,blTimes,method,flagsens,avmode);
            newosc(:,:,:,n) = corrP(condosc,blTimes,method,flagsens,avmode);
        end
        %get ITC
        if isfield(TF,'phase')
            condphase = TF.phase(:,:,:,find(trlvec==conds(n)));
            newphase(:,:,:,n) = squeeze(abs(mean(exp(1i*condphase),4)));
        end
        %get erp removed
        if isfield(TF,'erprem') %erp-removed
            conderprempow = TF.erprem.erprempow(:,:,:,find(trlvec==conds(n)));
            conderpremphase = TF.erprem.erpremphase(:,:,:,find(trlvec==conds(n)));
            newerprempow(:,:,:,n) = corrP(conderprempow,blTimes,method,flagsens,avmode);
            newerpremphase(:,:,:,n) = squeeze(abs(mean(exp(1i*conderpremphase),4)));
            conderppow = TF.erprem.erppow(:,:,:,find(trlvec==conds(n)));
            newerppow(:,:,:,n) = corrP(conderppow,blTimes,method,flagsens,avmode);
        end
    else
        %get only condition-specific power
        condP = TF.power(:,:,find(trlvec==conds(n)));
        %average power over trials
        newpower(:,:,n) = corrP(condP,blTimes,method,flagsens,avmode);
        %get parameterized
        if isfield(TF,'SPRiNT') %parameterized
            condap = TF.SPRiNT.ap_power(:,:,find(trlvec==conds(n)));
            condosc = TF.SPRiNT.osc_power(:,:,find(trlvec==conds(n)));
            newap(:,:,n) = corrP(condap,blTimes,method,flagsens,avmode);
            newosc(:,:,n) = corrP(condosc,blTimes,method,flagsens,avmode);
        end
        %get ITC
        if isfield(TF,'phase')
            condphase = TF.phase(:,:,find(trlvec==conds(n)));
            newphase(:,:,n) = squeeze(abs(mean(exp(1i*condphase),3)));
        end
    end
    %set trials/erp
    TF.trlerp(n)=numel(find(trlvec==conds(n)));
    %get behavior
    if isfield(TF,'behavior')
        newBeh(n) = fieldfun( @(varargin) nanmean([varargin{:}]), TF.behavior(find(trlvec==conds(n))) );
    end
end
TF.conds = n;
%add to output
TF.ntrls = 1;
TF.power = newpower;
if isfield(TF,'phase')
    TF.phase = newphase;
end
if isfield(TF,'SPRiNT')
    TF.SPRiNT.osc_power = newosc;
    TF.SPRiNT.ap_power = newap;
end
if isfield(TF,'erprem')
    TF.erprem.erprempow = newerprempow;
    TF.erprem.erpremphase = newerpremphase;
    TF.erprem.erppow = newerppow;
end
if isfield(TF,'behavior')
    TF.behavior = newBeh;
end

%get rid of redundant fields
if isfield(TF, 'event')
    TF = rmfield(TF,'event');
end
if isfield(TF, 'epoch')
    TF = rmfield(TF, 'epoch');
end

end




function power = corrP(power,blTimes,method,flagsens,avmode) %average power over trials
tfPow = squeeze(mean( power,ndims(power) ));
if flagsens==0
    if strcmp(avmode,'mean')
        blPow = repmat(squeeze(mean( tfPow(:,:,blTimes), 3)), [1,1,size(power,3)]);
        blPowSTD = repmat(squeeze(std(permute(tfPow(:,:,blTimes),[3,1,2]))), [1,1,size(power,3)]);
    elseif strcmp(avmode,'median')
        blPow = repmat(squeeze(median( tfPow(:,:,blTimes), 3)), [1,1,size(power,3)]);
        blPowSTD = repmat(squeeze(mad(permute(tfPow(:,:,blTimes),[3,1,2]))), [1,1,size(power,3)]);
    end
else
    if strcmp(avmode,'mean')
        blPow = repmat(squeeze(mean( tfPow(:,blTimes), 2)), [1,size(power,2)]);
        blPowSTD = repmat(squeeze(std(permute(tfPow(:,blTimes),[2,1]))), [1,size(power,2)]);
    elseif strcmp(avmode,'median')
        blPow = repmat(squeeze(median( tfPow(:,blTimes), 2)), [1,size(power,2)]);
        blPowSTD = repmat(squeeze(mad(permute(tfPow(:,blTimes),[2,1]))), [1,size(power,2)]);
    end
end
%correct power
switch method
    case 'zscore'
        disp('baselining with z-score');
        tfPow = (tfPow-blPow)./blPowSTD;
    case 'db'
        disp('baselining with decibel');
        tfPow = 10*log10( tfPow./blPow );
    case 'percent'
        disp('baselining using %-change');
        tfPow = 100*(tfPow-blPow)./blPow;
    case 'none'
        disp('NOT baselining power estimates');
end
power = tfPow;

end

