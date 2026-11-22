function TF = hnf_event_grid(TF,oldTimes,newTimes,keepTrials,selectedSamples)
% Translate flat EEGLAB sample latencies and the matching per-epoch event
% lists when cropping, resampling, or removing trials. Times use the same
% units on both grids. Urevent references always refer to the original EEG.
if isfield(TF,'epoch') && ~isempty(TF.epoch)
    TF.epoch = TF.epoch(keepTrials);
end
if ~isfield(TF,'event') || isempty(TF.event) || ~isfield(TF.event,'latency')
    return
end
oldTimes = double(oldTimes(:));
newTimes = double(newTimes(:));
event = TF.event;
eventMap = zeros(size(event));
position = nan(size(event));
kept = false(size(event));
for i = 1:numel(event)
    latency = event(i).latency;
    if ~isnumeric(latency) || ~isscalar(latency) || ~isfinite(latency), continue; end
    if isfield(event,'epoch') && ~isempty(event(i).epoch)
        oldTrial = event(i).epoch;
    else
        oldTrial = floor((latency-1)/numel(oldTimes))+1;
    end
    newTrial = find(keepTrials == oldTrial,1);
    if isempty(newTrial), continue; end
    sample = latency - (oldTrial-1)*numel(oldTimes);
    if sample<1 || sample>numel(oldTimes), continue; end
    if ~isempty(selectedSamples) && ~ismember(round(sample),selectedSamples), continue; end
    if numel(oldTimes)==1
        t = oldTimes(1);
    else
        t = interp1(1:numel(oldTimes),oldTimes,sample);
    end
    if numel(newTimes)==1
        if ~isequal(t,newTimes(1)), continue; end
        newSample = 1;
    else
        newSample = interp1(newTimes,1:numel(newTimes),t,'linear',NaN);
        if ~isfinite(newSample), continue; end
    end
    kept(i) = true;
    position(i) = newSample;
    event(i).latency = newSample+(newTrial-1)*numel(newTimes);
    if isfield(event,'epoch'), event(i).epoch = newTrial; end
    if isfield(event,'duration') && isnumeric(event(i).duration) && ...
            numel(newTimes)>1 && numel(oldTimes)>1
        event(i).duration = event(i).duration*median(diff(oldTimes))/median(diff(newTimes));
    end
end
eventMap(kept) = 1:nnz(kept);
TF.event = event(kept);
if ~isfield(TF,'epoch') || ~isfield(TF.epoch,'event'), return; end
for e = 1:numel(TF.epoch)
    refs = TF.epoch(e).event;
    if isempty(refs), continue; end
    if iscell(refs), refs = [refs{:}]; end
    subset = kept(refs);
    fields = fieldnames(TF.epoch);
    for j = 1:numel(fields)
        name = fields{j};
        value = TF.epoch(e).(name);
        if strcmp(name,'event')
            TF.epoch(e).event = eventMap(refs(subset));
        elseif strcmp(name,'eventposition') && ~isfield(event,'position')
            if iscell(value), TF.epoch(e).(name) = num2cell(position(refs(subset)));
            else, TF.epoch(e).(name) = position(refs(subset)); end
        elseif strncmp(name,'event',5) && numel(value)==numel(refs) && ~ischar(value)
            TF.epoch(e).(name) = value(subset);
        end
    end
end
end
