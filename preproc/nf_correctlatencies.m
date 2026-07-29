function [EEG, latencies] = nf_correctlatencies( EEG, targetEvent, shiftEvent, expectedCount )
% NF_CORRECTLATENCIES  Correct event-code latencies using a second event set.
%
% [EEG, LATENCIES] = NF_correctlatencies(EEG, 'stim', 'photodiode', 400);
%
% Corrects EEG TTL trigger latency using another set of event codes
% The other event codes are expected to come from an event marking device,
% such as Cortech Solutions "Triggy", Cedrus "StimTracker", Black Box
% Toolkit M2/M3, or some other comparable device.
%
% Name/value inputs:
%   targetEvent           a target event code to be corrected
%   shiftEvent            the event code marking the empirical onset
%                             (photodiode or other)
%   expectedCount         expected number of events
%

sEEG = pop_epoch(EEG, targetEvent, [0,.1]);
sLats = zeros(1,expectedCount);
for k = 1:length(sEEG.epoch)
    lt = cell2mat(sEEG.epoch(k).eventlatency(find(strcmp(sEEG.epoch(k).eventtype,shiftEvent))));
    if ~isempty(lt)
        sLats(k) = lt;
    end
end
% now, in case we missed an event during recording, we need to get
% median
sLat_tmp = sLats; % make a copy in preparation for calculating the median
disp([num2str(sum(sLat_tmp(sLat_tmp == 0))) ' external event ' ...
    'signals were missed during recording - substituting the median.']);
sLat_tmp(sLat_tmp == 0) = []; % get rid of missed photodiodes before calculating median
sLats(sLats==0) = median(sLat_tmp); % fill in median if we missed a photodiode somewhere (unlikely)

% now we loop again - but this time we rewrite the latency of the
% events to correct for the photodiode
sDone = 0;
for k = 1:length(EEG.event)
    if strcmp(EEG.event(k).type,targetEvent)
        sDone = sDone + 1;
        disp(['shifting event ' num2str(k) ' by ' ...
            num2str(round((sLats(sDone)/(1000/EEG.srate)))) ' samples.']);
        EEG.event(k).latency = EEG.event(k).latency + ... % original event latency
            round((sLats(sDone)/(1000/EEG.srate))); % shift by a rounded number of samples according to shift
    end
end
EEG = eeg_checkset(EEG, 'event');

% fin!!!
disp([num2str(sDone) ' events successfully corrected for an ' ...
    'average latency shift of ' num2str(median(sLats)) ' ms.']);

end
