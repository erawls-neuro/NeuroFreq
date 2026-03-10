function TFout = nf_epoch(TF, types, tlimits, varargin)
% NF_EPOCH  Epoch a continuous NF TF set using EEGLAB-style TF.event latencies.
%
% EEGLAB-LIKE OUTPUT BEHAVIOR (KEY):
% - TFout.event is a flat list of events (only events that occur inside kept epochs).
% - TFout.epoch is a struct array with per-epoch event listings:
%     .event         : indices into TFout.event
%     .eventtype     : cell array of event types
%     .eventlatency  : cell array of latencies (ms) relative to epoch (0 at timelock)
%     .eventposition : cell array of positions (points) within epoch
%     .eventduration : cell array (if duration exists)
%     .eventurevent  : cell array (if urevent exists)
%
% NOTE
% - This function does NOT create TFout.event(i).epoch (per your request).
%
% USAGE
%   TFep = nf_epoch(TF, {'stim'}, [-0.2 0.8])
%   TFep = nf_epoch(TF, 11, [-200 800], 'tunit', 'ms')
%
% REQUIRED
%   TF      : continuous TF struct (time is continuous, no trials)
%   types   : event types to epoch around (cellstr/string/numeric/mixed)
%   tlimits : [tmin tmax] around event
%
% OPTIONS (name-value)
%   'tunit'         : 'seconds' (default) or 'ms'
%   'eventfield'    : 'type' (default)
%   'latencyfield'  : 'latency' (default; assumed in POINTS)
%   'verbose'       : true/false (default true)

% -------------------------
% Parse inputs
% -------------------------
p = inputParser;

addRequired(p, 'TF', @isstruct);
addRequired(p, 'types');
addRequired(p, 'tlimits', @(x) isnumeric(x) && isvector(x) && numel(x) == 2);

addParameter(p, 'tunit', 'seconds', @(x) ischar(x) || isstring(x));
addParameter(p, 'eventfield', 'type', @(x) ischar(x) || isstring(x));
addParameter(p, 'latencyfield', 'latency', @(x) ischar(x) || isstring(x));
addParameter(p, 'verbose', true, @(x) islogical(x) && isscalar(x));

parse(p, TF, types, tlimits, varargin{:});

TFin = p.Results.TF;

% -------------------------
% Validate TF
% -------------------------
[~, rep] = nf_checkset(TFin, 'mode', 'error');

if strcmpi(rep.setType, 'stat')
    error('nf_epoch:BadType', 'nf_epoch does not operate on stat sets.');
end

if ~isfield(TFin, 'event')
    error('nf_epoch:MissingEvent', 'Continuous TF set must contain TF.event.');
end

if ~isstruct(TFin.event)
    error('nf_epoch:BadEvent', 'TF.event must be a struct array.');
end

eventField = char(p.Results.eventfield);
latField = char(p.Results.latencyfield);

if ~isfield(TFin.event, eventField)
    error('nf_epoch:BadEvent', 'TF.event missing eventfield: %s', eventField);
end

if ~isfield(TFin.event, latField)
    error('nf_epoch:BadEvent', 'TF.event missing latencyfield: %s', latField);
end

Fs = double(TFin.Fs);

if Fs <= 0
    error('nf_epoch:BadFs', 'TF.Fs must be > 0.');
end

% -------------------------
% Ensure continuous (no trials) for measured fields
% -------------------------
[powerFields, phaseFields, wpliFields] = local_discover_measured_fields(TFin);

if isempty(powerFields) && isempty(phaseFields) && isempty(wpliFields)
    error('nf_epoch:MissingData', 'No power/phase/wpli numeric fields found to epoch.');
end

local_error_if_epoched(TFin, powerFields, phaseFields, wpliFields);

% -------------------------
% Convert tlimits to seconds (pop_epoch convention)
% -------------------------
tmin = double(tlimits(1));
tmax = double(tlimits(2));

if strcmpi(char(p.Results.tunit), 'ms')
    tmin = tmin / 1000;
    tmax = tmax / 1000;
end

if tmax <= tmin
    error('nf_epoch:BadLimits', 'tlimits must satisfy tmax > tmin.');
end

smin = round(tmin * Fs);
smax = round(tmax * Fs);

L = smax - smin + 1;

if L < 2
    error('nf_epoch:BadEpochLength', 'Epoch length too short after rounding.');
end

nTime = numel(TFin.times);

% -------------------------
% Select matching timelocking events
% -------------------------
ev = TFin.event;

matchMask = false(1, numel(ev));

for iE = 1:numel(ev)
    et = ev(iE).(eventField);
    tf = local_match_type(et, types);
    matchMask(iE) = tf;
end

selIdx = find(matchMask);

if isempty(selIdx)
    error('nf_epoch:NoEvents', 'No events matched requested types.');
end

% -------------------------
% Precompute all event latencies (assumed POINTS)
% -------------------------
allLat = nan(1, numel(ev));

for iE = 1:numel(ev)
    latNow = ev(iE).(latField);

    if isempty(latNow)
        allLat(iE) = NaN;
    else
        allLat(iE) = double(latNow);
    end
end

if any(isfinite(allLat) & (allLat < 0))
    error('nf_epoch:BadLatency', 'Some event latencies are < 0. Expected points (>=1).');
end

if any(isfinite(allLat) & (allLat > (nTime + 1e-6)))
    warning('nf_epoch:LatencyRange', 'Some event latencies exceed length(TF.times). Are latencies really in points?');
end

% -------------------------
% Determine valid epochs (exclude boundary-truncated epochs)
% -------------------------
keepAnchorIdx = [];
epochWindows = [];
anchorLatRound = [];

for iK = 1:numel(selIdx)

    iAnchor = selIdx(iK);

    lat = allLat(iAnchor);

    if ~isfinite(lat)
        continue
    end

    latR = round(lat);

    if latR < 1
        continue
    end

    if latR > nTime
        continue
    end

    w1 = latR + smin;
    w2 = latR + smax;

    if w1 < 1
        continue
    end

    if w2 > nTime
        continue
    end

    keepAnchorIdx = [keepAnchorIdx iAnchor]; %#ok<AGROW>
    epochWindows = [epochWindows; [w1 w2]]; %#ok<AGROW>
    anchorLatRound = [anchorLatRound latR]; %#ok<AGROW>

end

if isempty(keepAnchorIdx)
    error('nf_epoch:NoValidEpochs', 'All matched events were too close to edges for the epoch window.');
end

nEp = size(epochWindows, 1);

if p.Results.verbose == true
    disp(['nf_epoch: creating ' num2str(nEp) ' epochs.']);
end

% -------------------------
% Build output struct core
% -------------------------
TFout = TFin;

TFout.times = (double(smin):double(smax)) ./ Fs;

if local_times_are_ms(TFin.times, Fs) == true
    TFout.times = TFout.times .* 1000;
end

TFout.ntrls = nEp;

% -------------------------
% Epoch all measured numeric fields (ch x f x t) -> (ch x f x t x ep)
% -------------------------
TFout = local_epoch_fields(TFout, powerFields, epochWindows, L, nEp);
TFout = local_epoch_fields(TFout, phaseFields, epochWindows, L, nEp);
TFout = local_epoch_fields(TFout, wpliFields, epochWindows, L, nEp);

% -------------------------
% Build EEGLAB-like TFout.event + TFout.epoch
% -------------------------
TFout.event = struct([]);
TFout.epoch = struct([]);

eventCursor = 0;

for iEp = 1:nEp

    w1 = epochWindows(iEp, 1);
    w2 = epochWindows(iEp, 2);

    inWin = false(1, numel(ev));

    for iE = 1:numel(ev)

        lat = allLat(iE);

        if ~isfinite(lat)
            inWin(iE) = false;
        else
            latR = round(lat);
            if latR >= w1 && latR <= w2
                inWin(iE) = true;
            else
                inWin(iE) = false;
            end
        end

    end

    evIdx = find(inWin);

    % Sort events within epoch by latency
    if ~isempty(evIdx)
        latVec = zeros(1, numel(evIdx));
        for k = 1:numel(evIdx)
            latVec(k) = round(allLat(evIdx(k)));
        end
        [~, ord] = sort(latVec, 'ascend');
        evIdx = evIdx(ord);
    end

    epEventIdxOut = zeros(1, numel(evIdx));
    epEventType = cell(1, numel(evIdx));
    epEventLatencyMs = cell(1, numel(evIdx));
    epEventPosition = cell(1, numel(evIdx));

    epEventDuration = [];
    hasDuration = false;

    epEventUrevent = [];
    hasUrevent = false;

    if ~isempty(evIdx)
        if isfield(ev, 'duration')
            hasDuration = true;
            epEventDuration = cell(1, numel(evIdx));
        end
        if isfield(ev, 'urevent')
            hasUrevent = true;
            epEventUrevent = cell(1, numel(evIdx));
        end
    end

    for k = 1:numel(evIdx)

        iE = evIdx(k);

        latOrig = allLat(iE);
        latR = round(latOrig);

        pos = latR - w1 + 1;

        if pos < 1
            error('nf_epoch:Internal', 'Computed event position < 1.');
        end

        if pos > L
            error('nf_epoch:Internal', 'Computed event position > epoch length.');
        end

        eventCursor = eventCursor + 1;

        eNew = ev(iE);

        if isfield(eNew, latField)
            eNew.latency_orig = eNew.(latField);
        else
            eNew.latency_orig = latOrig;
        end

        % EEGLAB-style concatenated event latency
        eNew.(latField) = (iEp - 1) * L + pos;

        TFout.event(eventCursor) = eNew;

        epEventIdxOut(k) = eventCursor;
        epEventType{k} = eNew.(eventField);

        % Event latency in ms relative to epoch timeline (0 at timelock)
        % timelock event is at ~t = 0 because start is tmin seconds.
        latSecRel = (double(pos - 1) / Fs) + tmin;
        epEventLatencyMs{k} = latSecRel * 1000;

        epEventPosition{k} = pos;

        if hasDuration == true
            epEventDuration{k} = eNew.duration;
        end

        if hasUrevent == true
            epEventUrevent{k} = eNew.urevent;
        end

    end

    ep = struct();

    ep.event = epEventIdxOut;
    ep.eventtype = epEventType;
    ep.eventlatency = epEventLatencyMs;
    ep.eventposition = epEventPosition;

    if hasDuration == true
        ep.eventduration = epEventDuration;
    end

    if hasUrevent == true
        ep.eventurevent = epEventUrevent;
    end

    % Keep a reference to the timelocking event index from original TF.event
    % (EEGLAB often has this implicitly via eventtype/latency; this is harmless metadata)
    ep.timelock_urevent_index = keepAnchorIdx(iEp);
    ep.timelock_latency_points = anchorLatRound(iEp);

    TFout.epoch(iEp) = ep;

end

% -------------------------
% Final consistency check
% -------------------------
[~, ~] = nf_checkset(TFout, 'mode', 'warn');

end

% ======================================================================
% Local helpers
% ======================================================================

function [powerFields, phaseFields, wpliFields] = local_discover_measured_fields(TF)

fn = fieldnames(TF);

powerFields = {};
phaseFields = {};
wpliFields = {};

for i = 1:numel(fn)

    name = fn{i};
    val = TF.(name);

    if ~isnumeric(val)
        continue
    end

    lname = lower(name);

    if contains(lname, 'power')
        powerFields{end + 1, 1} = name; %#ok<AGROW>
    elseif contains(lname, 'phase')
        phaseFields{end + 1, 1} = name; %#ok<AGROW>
    elseif contains(lname, 'wpli')
        wpliFields{end + 1, 1} = name; %#ok<AGROW>
    end

end

powerFields = sort(powerFields);
phaseFields = sort(phaseFields);
wpliFields = sort(wpliFields);

end

function local_error_if_epoched(TF, powerFields, phaseFields, wpliFields)

fields = [powerFields; phaseFields; wpliFields];

for i = 1:numel(fields)

    name = fields{i};
    X = TF.(name);

    if ndims(X) >= 4
        if size(X, 4) > 1
            error('nf_epoch:AlreadyEpoched', ...
                'Field %s indicates trials already exist (dim4>1). nf_epoch expects continuous TF.', name);
        end
    end

    if ndims(X) >= 5
        error('nf_epoch:Unsupported', ...
            'Field %s has 5+ dims; nf_epoch expects continuous TF with <=4 dims.', name);
    end

end

end

function tf = local_match_type(et, types)

tf = false;

if iscell(types)
    list = types;
else
    list = {types};
end

for i = 1:numel(list)

    t = list{i};

    if isequal(et, t)
        tf = true;
        return
    end

    if isnumeric(et) && isnumeric(t)
        if isequal(double(et), double(t))
            tf = true;
            return
        end
    end

    try
        a = char(string(et));
        b = char(string(t));
        if strcmpi(a, b)
            tf = true;
            return
        end
    catch
    end

end

end

function TFout = local_epoch_fields(TFout, fields, epochWindows, L, nEp)

if isempty(fields)
    return
end

nChan = numel(TFout.chanlocs);
nFreq = numel(TFout.freqs);

for iF = 1:numel(fields)

    name = fields{iF};
    X = TFout.(name);

    if ndims(X) == 3
        if size(X, 1) ~= nChan
            error('nf_epoch:BadShape', 'Field %s dim1 must match chanlocs.', name);
        end
        if size(X, 2) ~= nFreq
            error('nf_epoch:BadShape', 'Field %s dim2 must match freqs.', name);
        end
    elseif ndims(X) == 4
        if size(X, 4) ~= 1
            error('nf_epoch:BadShape', 'Field %s dim4 must be singleton for continuous TF.', name);
        end
        X = squeeze(X(:, :, :, 1));
    else
        error('nf_epoch:BadShape', 'Field %s must be 3D (chxfxt) or 4D with singleton dim4.', name);
    end

    if size(X, 3) < max(epochWindows(:, 2))
        error('nf_epoch:BadShape', 'Field %s time dimension shorter than epoch window.', name);
    end

    Y = zeros(nChan, nFreq, L, nEp);

    for iEp = 1:nEp

        w1 = epochWindows(iEp, 1);
        w2 = epochWindows(iEp, 2);

        idx = w1:w2;

        Y(:, :, :, iEp) = X(:, :, idx);

    end

    TFout.(name) = Y;

end

end

function tf = local_times_are_ms(times, Fs)

tf = true;

if numel(times) < 2
    tf = true;
    return
end

dt = median(diff(double(times(:))));

dt_ms = 1000 / Fs;
dt_s = 1 / Fs;

err_ms = abs(dt - dt_ms);
err_s = abs(dt - dt_s);

if err_s < err_ms
    tf = false;
else
    tf = true;
end

end








