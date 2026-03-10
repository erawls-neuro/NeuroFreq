%Function that simply calls pop_epoch eeglab function but selects only
%those epochs, that are indicated in 'epochs' input.
%AGF, 2013
%
% Modified by ELR, 2014-2022
%

function [EEG, indices, com] = nf_epochbynumber( EEG, epochs, xmin, xmax) 

if nargin<4
    xmax = EEG.xmax;
end
if nargin<3
    xmin = EEG.xmin;
end
if nargin<2
    epochs = 1:EEG.trials;
end

Ev_array=[];

%first make sure the event codes in behavior work...
if isfield(EEG.etc,'behavior') && numel(EEG.etc.behavior)~=EEG.trials
    error('EEG and behavioral codes do not align, please fix!');
end

%Get events at latency 0 in epoch structure
for c = 1 : length(EEG.epoch)
    if iscell(EEG.epoch(c).eventlatency(1)) % Check if data is cell or integer array and read that out
        Ev0=find([EEG.epoch(c).eventlatency{:}]==0);
        Ev0 = Ev0(1);
    else
        Ev0=find([EEG.epoch(c).eventlatency(:)]==0);
        Ev0 = Ev0(1);
    end
    Ev_array=[Ev_array EEG.epoch(c).event(Ev0)];
end

[EEG, indices, com] = pop_epoch( EEG, {},[xmin xmax+(1/EEG.srate)], 'newname', ...
    EEG.setname, 'epochinfo', 'yes', 'eventindices', Ev_array(epochs));

%fix the custom event codes
if isfield(EEG.etc,'behavior')
   EEG.etc.behavior=EEG.etc.behavior(epochs);
end

%finally catch if there is a trial mismatch
if isfield(EEG.etc,'behavior') && numel(EEG.etc.behavior)~=EEG.trials
    error('EEG and behavioral codes do not align, please fix!');
end

end