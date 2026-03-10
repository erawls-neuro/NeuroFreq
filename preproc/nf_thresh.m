%
%
% simple little script that thresholds and artifact-rejects epoched EEGs
% using specified input methods
%
%

function [EEG,rej] = nf_thresh(EEG,voltthresh,powthresh,freqrange,times,interp,maxbadchan,frontalchan)

if nargin<8 || isempty(frontalchan)
    [~,frontalchan]=maxk([EEG.chanlocs.X],floor(numel(EEG.chanlocs)/10));
    frontalchan = {EEG.chanlocs(frontalchan).labels};
end
if nargin<7 || isempty(maxbadchan)
    maxbadchan = round(numel(EEG.chanlocs)/10);
end
if nargin<6 || isempty(interp)
    interp = 1;
end
if nargin<5 || isempty(times)
    times=[EEG.xmin EEG.xmax];
end
if nargin<4 || isempty(freqrange)
    freqrange = [20 40];
end
if nargin<3 || isempty(powthresh)
    powthresh=[-100 30];
end
if nargin<2 || isempty(voltthresh)
    voltthresh = 125;
end
if EEG.trials==1
    EEG = eeg_regepochs(EEG,'recurrence', 1, 'limits',[0 1], 'rmbase', NaN, 'eventtype', '999'); % insert temporary marker 1 second apart and create epochs;
end
%thresholding
EEG = pop_rmbase(EEG,[]);
EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -1*voltthresh,  ...
    voltthresh, times(1), times(2), 0, 0);
EEG = pop_rejspec( EEG, 1, 'threshold', powthresh, ...
    'freqlimits', freqrange, 'eegplotreject', 0);
if interp==0
    rej = (EEG.reject.rejthresh+EEG.reject.rejfreq)>0;
elseif interp==1
    chans = EEG.chanlocs;
    rej=zeros(1,EEG.trials);
    for k=1:EEG.trials
        % put trial number in behavior
        EEG.etc.behavior(k).trial = k;
    end
    for k=1:EEG.trials
        % thresholding
        badc = find((EEG.reject.rejthreshE(:,k)+EEG.reject.rejfreqE(:,k))>0);
        % if any frontal channels are bad, the epoch is dead
        if ~isempty(intersect(badc,find(ismember({EEG.chanlocs.labels},frontalchan))))
            rej(k)=1;
            %epoch also dies if there are > max bad channels
        elseif numel(badc)>maxbadchan
            rej(k)=1;
            %if neither, it is a good epoch - we can interpolate
        elseif numel(badc)==numel(EEG.chanlocs)
            rej(k)=1;
        else
            rej(k)=0;
        end
        if rej(k)==0
            % we can continue
            EEG.etc.behavior(k).nLocalBadChan = numel(badc);
            % cut to single epoch
            data = nf_epochbynumber(EEG,k);
            % remove the bad channels
            data = pop_select( data, 'nochannel', badc );
            % interpolate missing channels
            data = pop_interp(data, chans, 'sphericalKang');
            % average reference again following interpolation
            data = pop_reref(data,[]);
            % copy into EEG
            EEG.data(:,:,k) = data.data;
        else
            % we cannot continue
            EEG.etc.behavior(k).nLocalBadChan = 999;
        end
    end
end

%reject
EEG = nf_epochbynumber( EEG, ~rej);

%interpolate missing channels
EEG = pop_interp(EEG, EEG.etc.ogchan, 'sphericalKang');

%average reference again following interpolation
EEG = pop_reref(EEG,[]);

end


