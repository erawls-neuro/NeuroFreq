function reg_out = nf_stregress(EEG,regs,varargin)
% NF_STREGRESS    Performs single-trial analysis of TF transformed EEG
%
% GENERAL
% -------
% performs single-trial analysis of time-domain EEG or time-frequency
% transformed scalp-recorded EEG. Operates on the output of
% nf_tftransform.m or the tf_fun functions. Uses tfce correlations,
% robust regression, or Pearson, Spearman or Kendall correlations.
%
% USAGE
% -----
% Outputs:
% 1) reg - a structure with information about the regression
%     analysis, including weights and pvalues
%
% Inputs:
% 1) TF/EEG - required. an EEGLAB set or TF structure from nf_tftransform
% 2) regs - required. Must be a vector or matrix, with one dimension equal
%     in size to EEG trials.
%
% various name-value pairs:
%
%     'toi': a pair of times [xmin xmax] to be returned. Will analyze all
%       times between these points.
%
%     'foi': pair of frequencies [fmin fmax] to be returned. Will analyze
%       all freqs between these points.
%
%     'method': a string detailing the method for analysis. Supports
%       'tfce', 'robust', 'pearson', 'spearman', 'kendall'
%


% parse inputs
p = inputParser;
validNumMatrix = @(x) isnumeric(x) && isvector(x);
defaultTois = [EEG.times(1) EEG.times(end)];
defaultFois = [];
expectedMethods = {'tfce' 'robust' 'kendall' 'pearson' 'spearman'};
defaultmethod = 'robust';
addRequired(p,'EEG');
addRequired(p,'regs');
addParameter(p,'toi',defaultTois,validNumMatrix);
addParameter(p,'foi',defaultFois,validNumMatrix);
addParameter(p,'method',defaultmethod,@(x) any(validatestring(x,expectedMethods)));
parse(p,EEG,regs,varargin{:});

%first, what type of input?
tIndices = find( (EEG.times>=p.Results.toi(1)) + (EEG.times<=p.Results.toi(end)) == 2);
if isfield(EEG,'freqs') %its TF
    tf=1;
    trials=size(EEG.power,4);
    %trim freqs
    fIndices = find( (EEG.freqs>=p.Results.foi(1)) + (EEG.freqs<=p.Results.foi(end)) == 2);
    power=EEG.power(:,fIndices,tIndices,:);
    freqs=EEG.freqs(fIndices);
else %its EEG
    tf=0;
    trials=EEG.trials;
    power=EEG.data(:,tIndices,:);
end
times=EEG.times(tIndices);

% Find regressor size
[n,m]=size(regs);

%check size match
if n==trials
    reg=regs;
elseif m==trials
    reg=regs';
else
    error('mismatch in regressor size?');
end

%get power and reshape for efficiency
EEGData=power;
OEEGD=EEGData;
if tf==1
    EEGData=permute(EEGData,[4 1 2 3]);
    datashape=size(squeeze(EEGData(1,:,:,:)));
elseif tf==0
    EEGData=permute(EEGData,[3 1 2]);
    datashape=size(squeeze(EEGData(1,:,:)));
end

%get method
method=p.Results.method;

%matricize data for faster results, reshape it later
if ~strcmp(method, 'tfce') %tfce does its own thing
    s1 = size(EEGData);
    EEGData = EEGData(:,:);
    pts = 1:size(EEGData,2);
end

%create output structure
reg_out.ntrls = trials; % # of stimulus bins
if tf==1
    reg_out.times = times; % freq points of interest in the analysis
    reg_out.freqs = freqs; % freq points of interest in the analysis
    reg_out.Fs = EEG.Fs;
elseif tf==0
    reg_out.times = times/1000;
    freqs=[];
    reg_out.Fs = EEG.srate;
end
if isfield(EEG,'chanlocs')
    reg_out.chanlocs = EEG.chanlocs; % electrodes included in the analysis
    reg_out.nsensor=length(reg_out.chanlocs);
end
reg_out.ntrls=1;
reg_out.conds=min(size(regs));
reg_out.scale='stat';
    
%switch method
switch method
    
    %TFCE case
    case 'tfce'
        if size(reg,2)==1
            tfce_out = nf_corr_TFCE(EEGData,reg,reg_out.chanlocs,freqs,times,1000,1);
        else
            error('no glm defined for tfce testing');
        end
        obs=tfce_out.corrcoef;
        tVals=tfce_out.tstat;
        pVals=tfce_out.p_vals;
        type=tfce_out.type;
        ncond=size(reg,2);
        
    %robust regression case
    case 'robust'
        obs=zeros(size(EEGData,2),size(reg,2));
        tVals=obs;
        pVals=obs;
        if isempty(gcp)
            warning off;
        else
            parfevalOnAll(@warning,0,'off','all');
        end
        disp('running single-trial regression...');
        parfor pt = pts
            EEGDataNow = squeeze(EEGData(:,pt));
            [bweights,stat]=robustfit(zscore(reg),zscore(EEGDataNow));
            obs(pt,:)=bweights(2:end);
            tVals(pt,:)=stat.t(2:end);
            pVals(pt,:)=stat.p(2:end);
        end
        if isempty(gcp)
            warning on;
        else
            parfevalOnAll(@warning,0,'on','all');
        end
        %resize things
        obs=squeeze(reshape(obs,[datashape size(reg,2)]));
        tVals=squeeze(reshape(tVals,[datashape size(reg,2)]));
        pVals=squeeze(reshape(pVals,[datashape size(reg,2)]));
        ncond=size(reg,2);
        
        %figure type
        if ~isempty(times)
            e3 = 1;
        else
            e3=0;
        end
        if ~isempty(freqs)
            e2 = 1;
        else
            e2=0;
        end
        if ~isempty(EEG.chanlocs)
            e1 = 1;
        else
            e1=0;
        end
        %figure out type for plotting
        if e1==1 && e2==1 && e3==1
            if size(reg,2)==1
                type = 'corr_chXfsXts';
            else
                type = 'glm_chXfsXts';
            end
        elseif e1==1 && e2==1
            if size(reg,2)==1
                type = 'corr_chXts';
            else
                type = 'glm_chXts';
            end
        elseif e1==1 && e3==1
            if size(reg,2)==1
                type = 'corr_chXts';
            else
                type = 'glm_chXts';
            end
        end
        EEGData=reshape(EEGData,s1);

    %correlation case
    case {'spearman' 'pearson' 'kendall'}
        disp('running single-trial correlation...');
        [obs,pVals] = corr(EEGData,reg,'type',method);
        t = obs .* sqrt((length(reg)-2)./(1-obs.^2));
        obs=reshape(obs,datashape);
        tVals=reshape(t,datashape);
        pVals=reshape(pVals,datashape);
        ncond=size(reg,2);
        
        %figure type
        if ~isempty(times)
            e3 = 1;
        else
            e3=0;
        end
        if ~isempty(freqs)
            e2 = 1;
        else
            e2=0;
        end
        if ~isempty(EEG.chanlocs)
            e1 = 1;
        else
            e1=0;
        end
        %figure out type for plotting
        if e1==1 && e2==1 && e3==1
            type = 'corr_chXfsXts';
        elseif e1==1 && e2==1
            type = 'corr_chXts';
        elseif e1==1 && e3==1
            type = 'corr_chXts';
        end
        EEGData=reshape(EEGData,s1);
end

reg_out.power=OEEGD;
reg_out.corrvar=reg;
reg_out.df=size(EEGData,1)-2;
reg_out.corrcoef=obs;
reg_out.tstat=tVals;
reg_out.p_vals=pVals;
reg_out.method=method;
reg_out.type=type;
reg_out.conds=ncond;


end









