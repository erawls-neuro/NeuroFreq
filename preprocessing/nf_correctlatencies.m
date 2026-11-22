function [EEG,sLats] = nf_correctlatencies(EEG, targetEvent, shiftEvent, expectedCount)
% NF_CORRECTLATENCIES  Align TTL events to their measured physical onsets.
%   [EEG, delayMs] = nf_correctlatencies(EEG,'stim','photodiode',400);
%   [EEG, delayMs] = nf_correctlatencies(EEG,'stim',4,400);
% A numeric shift is a fixed number of samples. Returned delays are in ms.
% Otherwise pair each target with the first matching onset in the following
% 100 ms, stopping at the next target, EEGLAB boundary or epoch boundary.
% Retain exact zero delays. Missing measurements use the median observed
% delay. EEGLAB permits fractional-sample latencies after resampling.

types = {EEG.event.type};
sample = [EEG.event.latency];
targets = find(eventMatches(types,targetEvent));
if numel(targets) ~= expectedCount || isempty(targets)
    error('nf_correctlatencies:TargetCount', ...
        'Expected %d target events and found %d.', expectedCount, numel(targets));
end
if isnumeric(shiftEvent) && isscalar(shiftEvent)
    sLats = repmat(double(shiftEvent) * 1000 / EEG.srate, 1, numel(targets));
else
    markers = find(eventMatches(types,shiftEvent));
    bounds = sample(strcmp(types,'boundary'));
    sLats = nan(1,numel(targets));
    [~,order] = sort(sample(targets));
    for n = 1:numel(order)
        k = order(n);
        target = targets(k);
        t = sample(target);
        stop = t + EEG.srate / 10;
        if n < numel(order), stop = min(stop,sample(targets(order(n+1)))); end
        nextBoundary = bounds(bounds >= t);
        if ~isempty(nextBoundary), stop = min(stop,min(nextBoundary)); end
        candidates = markers(sample(markers) >= t & sample(markers) < stop);
        if isfield(EEG.event,'epoch') && ~isempty(EEG.event(target).epoch)
            sameEpoch = arrayfun(@(i) isequal(EEG.event(i).epoch,EEG.event(target).epoch),candidates);
            candidates = candidates(sameEpoch);
        end
        if ~isempty(candidates)
            sLats(k) = (min(sample(candidates))-t) * 1000 / EEG.srate;
        end
    end
    missed = isnan(sLats);
    if all(missed)
        error('nf_correctlatencies:NoMeasurements','No target had a matching physical onset within 100 ms.');
    end
    if any(missed)
        fprintf('[nf_correctlatencies] Substituting the observed median for %d missing onsets.\n',nnz(missed));
        sLats(missed) = median(sLats(~missed));
    end
end
newSample = sample(targets) + sLats*EEG.srate/1000;
if any(newSample < 1 | newSample > EEG.pnts*EEG.trials)
    error('nf_correctlatencies:OutsideRecording', ...
        'A corrected onset falls outside the recording. Check the fixed or imputed delay.');
end
for k = 1:numel(targets)
    EEG.event(targets(k)).latency = newSample(k);
end
EEG = eeg_checkset(EEG,'eventconsistency');
fprintf('[nf_correctlatencies] Corrected %d events; median shift %.3g ms.\n',numel(targets),median(sLats));
end

function mask = eventMatches(types, wanted)
% EEGLAB accepts text labels, numeric codes, and cell arrays of labels.
if ~iscell(wanted)
    if isstring(wanted), wanted = cellstr(wanted); else, wanted = {wanted}; end
end
mask = false(size(types));
for k = 1:numel(wanted)
    label = wanted{k};
    if ischar(label) || isstring(label)
        mask = mask | strcmp(types,label);
    else
        mask = mask | cellfun(@(v) isequal(v,label),types);
    end
end
end
