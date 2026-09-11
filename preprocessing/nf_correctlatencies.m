function [EEG,sLats] = nf_correctlatencies( EEG, targetEvent, shiftEvent, expectedCount )
% NF_CORRECTLATENCIES  Correct event-code latencies using a second event set.
%
% [EEG, LATENCIES] = nf_correctlatencies(EEG, 'stim', 'photodiode', 400);
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
%                             OR a scalar indicating a standard sample shift 
%                             (only used when shiftEvent is missing)
%   expectedCount         expected number of events
%

% check for what type of event correction - static or empirical
if isnumeric(shiftEvent) & isscalar(shiftEvent)
    disp(['[nf_correctlatencies]: correcting events using a standard delay; assuming ' ...
        'triggy/cedrus/bbtk/etc was disconnected/not available.']);
    % preallocate latencies vector
    sLats = zeros(1,expectedCount) + shiftEvent; 
else % almost all of this is only needed for empirical events
    photo_count = sum(strcmp({EEG.event.type},shiftEvent));
    if photo_count==0
        error('no shift events detected in dataset; assuming triggy was disconnected.');
    end
    % epoch EEG around stimuli closely (0-100 ms)
    sEEG = pop_epoch(EEG, targetEvent, [0,.1]);
    % preallocate latencies vector
    sLats = zeros(1,expectedCount);
    % loop thru every epoch
    for k = 1:length(sEEG.epoch)
        % take the difference in latencies between the target and shift event,
        % in samples, for every event
        lt = cell2mat(sEEG.epoch(k).eventlatency(find(strcmp(sEEG.epoch(k).eventtype,shiftEvent))));
        % if this value exists - add it to the vector
        if ~isempty(lt)
            sLats(k) = lt;
        end
    end
    % now, in case we missed an event during recording, we need to get
    % median
    sLat_tmp = sLats; % make a copy in preparation for calculating the median
    missed = sum(sLat_tmp(sLat_tmp == 0));
    % provide info on missed codes - if any
    if missed>0
        disp(['[nf_correctlatencies]: ' num2str(missed) ' external event ' ...
            'signals were missed during recording - substituting the median.']);
    else
        disp('[nf_correctlatencies]: no missed event signals!');
    end
    % get rid of missed photodiodes before calculating median
    sLat_tmp(sLat_tmp == 0) = [];
    % fill in median if we missed a photodiode somewhere (unlikely)
    sLats(sLats==0) = median(sLat_tmp);
end

% now we loop again - but this time we rewrite the latency of the
% events to correct for the photodiode
sDone = 0;
for k = 1:length(EEG.event)
    if strcmp(EEG.event(k).type,targetEvent)
        sDone = sDone + 1;
        EEG.event(k).latency = EEG.event(k).latency + ... % original event latency
            round((sLats(sDone)/(1000/EEG.srate))); % shift by a rounded number of samples according to shift
    end
end
EEG = eeg_checkset(EEG, 'event');

% fin!!!
disp(['[nf_correctlatencies]: ' num2str(sDone) ' events successfully corrected for an ' ...
    'average latency shift of ' num2str(median(sLats)) ' ms.']);

end
