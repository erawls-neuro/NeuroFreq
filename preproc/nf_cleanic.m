%
% cleans epoched adult EEG data with ICA and ICLabel
%
function EEG = nf_cleanic( EEG,method,aggressive )
if nargin<2 || isempty(method)
    method = 'iclabel';
end
if nargin<3 || isempty(aggressive)
    aggressive=0;
end
% if continuous - copy and regepoch
EEG_copy = EEG;
if EEG.trials==1
    % regepochs
    EEG_copy = eeg_regepochs(EEG_copy,...
        'recurrence', 1, 'limits',[0 1], 'rmbase', NaN, 'eventtype', '999'); % insert temporary marker 1 second apart and create epochs
end
EEG_copy = pop_rmbase (EEG_copy, [] );
% threshhold
EEG_copy = pop_eegthresh(EEG_copy, 1, 1:EEG_copy.nbchan, -500,  ...
    500, EEG_copy.xmin, EEG_copy.xmax, 0, 1);
EEG_copy = pop_rejspec( EEG_copy, 1, 'threshold', [-100 30], ...
    'freqlimits', [20 40], 'eegplotreject', 1 ); % remove the worst segments that cannot be repaired by ICA
% filter more aggressively at 1 Hz - helps ICA
EEG_copy = pop_eegfiltnew( EEG_copy, 'locutoff', 1);
% run ICA and reject
EEG = icrun(EEG,EEG_copy,method,aggressive);
if strcmp(method,'iclabel')
    EEG = icrun(EEG,EEG_copy,method,aggressive);
end
end

function EEG = icrun( EEG,EEG_copy,method,aggressive )
%clear existing ICA info if any
EEG.icaact=[];
EEG.icachansind=[];
EEG.icasphere=[];
EEG.icaweights=[];
EEG.icawinv=[];
EEG_copy.icaact=[];
EEG_copy.icachansind=[];
EEG_copy.icasphere=[];
EEG_copy.icaweights=[];
EEG_copy.icawinv=[];
%run ICA
EEG_copy = pop_runica(EEG_copy, 'pca', getrank(EEG_copy), 'extended', 1);
%label ICs
if strcmpi(method,'iclabel')
    EEG_copy = pop_iclabel(EEG_copy, 'default');
    %spare all brain and other
    [~,a]=max(EEG_copy.etc.ic_classification.ICLabel.classifications'); %#ok
    if ~aggressive
        a(a==7)=1;
    end
    rej = find(a>1);
elseif strcmpi(method,'made')
    rej = adjusted_ADJUST(EEG_copy);
end
% copy components
EEG.icachansind=EEG_copy.icachansind;
EEG.icasphere=EEG_copy.icasphere;
EEG.icaweights=EEG_copy.icaweights;
EEG.icawinv=EEG_copy.icawinv;
%delete components
EEG = pop_subcomp(EEG,rej);
end

function tmprank2 = getrank( EEG )
tmpdata = double(EEG.data(:,:));
%Here: alternate computation of the rank by Sven Hoffman
covarianceMatrix = cov(tmpdata', 1);
[~, D] = eig (covarianceMatrix);
rankTolerance = 1e-7;
tmprank2=sum (diag (D) > rankTolerance);
end




