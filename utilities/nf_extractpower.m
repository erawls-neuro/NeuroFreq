function [behTab, TF] = nf_extractpower(TF, rois, varargin)
% NF_EXTRACTPOWER  Extract ROI-averaged power from ALL top-level NF power* fields
%                 and append to TF.behavior, returning behavior as a table.
%
% USAGE
% -----
%   rois(1).name  = 'mfTheta';
%   rois(1).chans = {'FCz','Fz','Cz'};
%   rois(1).foi   = [4 8];
%   rois(1).toi   = [200 400];   % ms (auto-detected vs TF.times)
%
%   rois(2).name  = 'postAlpha';
%   rois(2).chans = {'POz','Oz','Pz'};
%   rois(2).foi   = [8 12];
%   rois(2).toi   = [300 800];   % ms
%
%   [behTab, TF] = nf_extractpower(TF, rois);
%
% INPUTS
% ------
% TF:
%   Must contain TF.freqs, TF.times.
%   If TF.chanlocs exists, channel labels can be resolved from it.
%   If TF.behavior exists (struct array or table), new columns are appended.
%   If TF.behavior is missing, an empty behavior struct is created.
%
% rois:
%   struct array with fields:
%     .name  (string/char)     -> ROI name used in output column naming
%     .chans (indices/logical/cellstr/string) -> channels to average
%     .foi   ([fmin fmax] or vector of freqs) -> frequency selection
%     .toi   ([tmin tmax] or vector of times) -> time selection
%
% NAME-VALUE
% ----------
% 'time_units'  : 'auto' (default) | 'ms' | 's'
% 'overwrite'   : 0/1 (default 0). If 0, collisions get suffixed _2, _3, ...
% 'verbose'     : 0/1 (default 0)
% 'powerFields' : cellstr to restrict which top-level power fields are used (default all)
%
% OUTPUTS
% -------
% behTab : table version of TF.behavior with added ROI columns
% TF     : TF with updated TF.behavior (struct array)

% ----------------------------
% Parse options
% ----------------------------
p = inputParser;

addRequired(p, 'TF', @(x) isstruct(x) == 1);
addRequired(p, 'rois');

addParameter(p, 'time_units', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'overwrite', 0, @(x) isscalar(x) && ismember(x, [0 1]));
addParameter(p, 'verbose', 0, @(x) isscalar(x) && ismember(x, [0 1]));
addParameter(p, 'powerFields', {}, @(x) iscell(x) || isstring(x));

parse(p, TF, rois, varargin{:});

opt = p.Results;

if isstring(opt.time_units) == 1
    opt.time_units = char(opt.time_units);
end

% ----------------------------
% Validate TF axes
% ----------------------------
if isfield(TF, 'freqs') == 0
    error('nf_extractpower:BadTF', 'TF.freqs is required.');
end

if isfield(TF, 'times') == 0
    error('nf_extractpower:BadTF', 'TF.times is required.');
end

freqs = TF.freqs(:);
times = double(TF.times(:));

if isempty(freqs) == 1 || isnumeric(freqs) == 0
    error('nf_extractpower:BadTF', 'TF.freqs must be numeric and non-empty.');
end

if isempty(times) == 1 || isnumeric(times) == 0
    error('nf_extractpower:BadTF', 'TF.times must be numeric and non-empty.');
end

% ----------------------------
% Normalize ROI input
% ----------------------------
if iscell(rois) == 1
    rois = [rois{:}];
end

if isstruct(rois) == 0
    error('nf_extractpower:BadROIs', 'rois must be a struct array (or cell array of structs).');
end

if isempty(rois) == 1
    error('nf_extractpower:BadROIs', 'rois is empty.');
end

reqRoiFields = {'name', 'chans', 'foi', 'toi'};

for iR = 1:numel(rois)
    for j = 1:numel(reqRoiFields)
        if isfield(rois(iR), reqRoiFields{j}) == 0
            error('nf_extractpower:BadROIs', 'ROI %d missing field "%s".', iR, reqRoiFields{j});
        end
    end
end

% ----------------------------
% Ensure behavior exists as struct array
% ----------------------------
if isfield(TF, 'behavior') == 1

    if istable(TF.behavior) == 1
        TF.behavior = table2struct(TF.behavior, 'ToScalar', false);
    end

    if isstruct(TF.behavior) == 0
        error('nf_extractpower:BadBehavior', 'TF.behavior must be a struct array or table.');
    end

else
    TF.behavior = struct([]);
end

% ----------------------------
% Find top-level power fields (NF convention)
% ----------------------------
allPowFields = local_find_top_power_fields(TF);

if isempty(opt.powerFields) == 0

    pf = opt.powerFields;

    if isstring(pf) == 1
        pf = cellstr(pf(:));
    end

    if iscell(pf) == 0
        error('nf_extractpower:BadArgs', 'powerFields must be a cellstr or string array.');
    end

    keep = false(numel(allPowFields), 1);

    for i = 1:numel(allPowFields)
        keep(i) = any(strcmpi(allPowFields{i}, pf));
    end

    allPowFields = allPowFields(keep);

end

if isempty(allPowFields) == 1
    error('nf_extractpower:NoPower', 'No top-level power* fields found.');
end

% ----------------------------
% Determine nObs from behavior or data fields
% ----------------------------
nObs = 0;

if isempty(TF.behavior) == 0
    nObs = numel(TF.behavior);
end

if nObs < 1
    % Try common index fields
    if isfield(TF, 'trial_subjectindex') == 1
        nObs = numel(TF.trial_subjectindex);
    end
end

if nObs < 1
    % Fall back to first power field observation dim guess
    X0 = TF.(allPowFields{1});
    nObs = local_guess_nobs_from_field(X0, numel(freqs), numel(times));
end

if nObs < 1
    error('nf_extractpower:BadObs', 'Could not determine number of observations.');
end

% If behavior was missing, create empty behavior rows now
if isempty(TF.behavior) == 1
    TF.behavior = repmat(struct(), 1, nObs);

    if isfield(TF, 'trial_subjectindex') == 1
        if numel(TF.trial_subjectindex) == nObs
            for i = 1:nObs
                TF.behavior(i).trial_subjectindex = TF.trial_subjectindex(i);
            end
        end
    end

    if isfield(TF, 'trial_condindex') == 1
        if numel(TF.trial_condindex) == nObs
            for i = 1:nObs
                TF.behavior(i).trial_condindex = TF.trial_condindex(i);
            end
        end
    end

    if isfield(TF, 'trial_subid') == 1
        if numel(TF.trial_subid) == nObs
            for i = 1:nObs
                TF.behavior(i).trial_subid = TF.trial_subid{i};
            end
        end
    end

end

% ----------------------------
% Pre-resolve ROI indices
% ----------------------------
roiResolved = rois;

for iR = 1:numel(rois)

    roiName = rois(iR).name;

    if isstring(roiName) == 1
        roiName = char(roiName);
    end

    if ischar(roiName) == 0 || isempty(roiName) == 1
        error('nf_extractpower:BadROIs', 'ROI %d name must be a non-empty string.', iR);
    end

    roiResolved(iR).name = roiName;

    roiResolved(iR).chan_idx = local_resolve_chan_idx(TF, rois(iR).chans);
    roiResolved(iR).freq_idx = local_resolve_axis_idx(freqs, rois(iR).foi, 'freq');
    roiResolved(iR).time_idx = local_resolve_time_idx(times, rois(iR).toi, opt.time_units);

    if isempty(roiResolved(iR).chan_idx) == 1
        error('nf_extractpower:BadROIs', 'ROI "%s" resolved to 0 channels.', roiName);
    end

    if isempty(roiResolved(iR).freq_idx) == 1
        error('nf_extractpower:BadROIs', 'ROI "%s" resolved to 0 freqs.', roiName);
    end

    if isempty(roiResolved(iR).time_idx) == 1
        error('nf_extractpower:BadROIs', 'ROI "%s" resolved to 0 times.', roiName);
    end

end

% ----------------------------
% Main extraction loop
% ----------------------------
if opt.verbose == 1
    fprintf('nf_extractpower: nObs=%d, powerFields=%d, rois=%d\n', nObs, numel(allPowFields), numel(roiResolved));
end

for iF = 1:numel(allPowFields)

    powName = allPowFields{iF};
    Xraw = TF.(powName);

    if isnumeric(Xraw) == 0
        continue
    end

    if isempty(Xraw) == 1
        continue
    end

    [chanDim, freqDim, timeDim, obsDim] = local_detect_dims(Xraw, TF, freqs, times, nObs, powName);

    for iR = 1:numel(roiResolved)

        roiName = roiResolved(iR).name;

        chanIdx = roiResolved(iR).chan_idx;
        freqIdx = roiResolved(iR).freq_idx;
        timeIdx = roiResolved(iR).time_idx;

        vals = local_extract_roi_mean(Xraw, chanDim, freqDim, timeDim, obsDim, chanIdx, freqIdx, timeIdx, nObs);

        colBase = [powName '_' roiName];
        colName = matlab.lang.makeValidName(colBase);

        if opt.overwrite == 0
            colName = local_unique_fieldname(TF.behavior, colName);
        end

        for iObs = 1:nObs
            TF.behavior(iObs).(colName) = vals(iObs);
        end

    end

end

% ----------------------------
% Return behavior table
% ----------------------------
behTab = struct2table(TF.behavior);

end

% ====================================================================== %
% Local helpers
% ====================================================================== %

function powFields = local_find_top_power_fields(TF)

fn = fieldnames(TF);
powFields = {};

for i = 1:numel(fn)

    name = fn{i};

    if strcmpi(name, 'freqs') == 1
        continue
    end

    if strcmpi(name, 'times') == 1
        continue
    end

    if strcmpi(name, 'chanlocs') == 1
        continue
    end

    if strcmpi(name, 'behavior') == 1
        continue
    end

    if strcmpi(name, 'scale') == 1
        continue
    end

    if strcmpi(name, 'aggtype') == 1
        continue
    end

    v = TF.(name);

    if isnumeric(v) == 0
        continue
    end

    isPow = false;

    if strncmpi(name, 'power', 5) == 1
        isPow = true;
    end

    if strncmpi(name, 'pow_', 4) == 1
        isPow = true;
    end

    if isPow == 0
        continue
    end

    powFields{end + 1, 1} = name; %#ok<AGROW>

end

% Prefer 'power' first if present
idx = find(strcmpi(powFields, 'power'), 1);

if isempty(idx) == 0
    powFields(idx) = [];
    powFields = [{'power'}; powFields];
end

end

function nObs = local_guess_nobs_from_field(X, nFreq, nTime)

nObs = 0;

if isnumeric(X) == 0
    return
end

sz = size(X);

if numel(sz) < 4
    return
end

% Heuristic: if looks like C x F x T x N
if sz(2) == nFreq && sz(3) == nTime
    nObs = sz(4);
    return
end

% If looks like N x C x F x T
if sz(3) == nFreq && sz(4) == nTime
    nObs = sz(1);
    return
end

end

function chanIdx = local_resolve_chan_idx(TF, chansIn)

if isempty(chansIn) == 1
    error('nf_extractpower:BadChans', 'ROI chans cannot be empty.');
end

if isnumeric(chansIn) == 1

    chanIdx = unique(chansIn(:)');

elseif islogical(chansIn) == 1

    chanIdx = find(chansIn(:) ~= 0);
    chanIdx = chanIdx(:)';

elseif isstring(chansIn) == 1 || iscell(chansIn) == 1

    if isstring(chansIn) == 1
        chansIn = cellstr(chansIn(:));
    end

    if isfield(TF, 'chanlocs') == 0
        error('nf_extractpower:BadChans', 'Channel labels requested but TF.chanlocs is missing.');
    end

    labels = local_chan_labels(TF.chanlocs);

    want = chansIn(:);
    wantIdx = [];

    for i = 1:numel(want)

        nm = want{i};

        if isstring(nm) == 1
            nm = char(nm);
        end

        if ischar(nm) == 0
            continue
        end

        hit = find(strcmpi(labels, nm), 1);

        if isempty(hit) == 0
            wantIdx = [wantIdx hit]; %#ok<AGROW>
        end

    end

    chanIdx = unique(wantIdx);

else

    error('nf_extractpower:BadChans', 'ROI chans must be numeric/logical indices or channel labels.');

end

chanIdx = chanIdx(chanIdx >= 1);

end

function labels = local_chan_labels(chanlocs)

nChan = numel(chanlocs);
labels = cell(nChan, 1);

for i = 1:nChan

    if isstruct(chanlocs(i)) == 1

        if isfield(chanlocs(i), 'labels') == 1
            labels{i} = chanlocs(i).labels;
        elseif isfield(chanlocs(i), 'label') == 1
            labels{i} = chanlocs(i).label;
        else
            labels{i} = sprintf('Ch%d', i);
        end

    else
        labels{i} = sprintf('Ch%d', i);
    end

end

end

function idx = local_resolve_axis_idx(axisVals, spec, which)

if isempty(spec) == 1
    error('nf_extractpower:BadROI', '%s spec is empty.', which);
end

axisVals = double(axisVals(:));

if isnumeric(spec) == 0
    error('nf_extractpower:BadROI', '%s spec must be numeric.', which);
end

spec = double(spec(:)');

if numel(spec) == 2

    lo = min(spec(1), spec(2));
    hi = max(spec(1), spec(2));

    idx = find(axisVals >= lo & axisVals <= hi);

else

    % treat as a list of exact values (within tolerance)
    tol = 1e-8;
    idx = [];

    for i = 1:numel(spec)
        hit = find(abs(axisVals - spec(i)) <= tol);
        idx = [idx; hit(:)]; %#ok<AGROW>
    end

    idx = unique(idx);

end

end

function idx = local_resolve_time_idx(times, toi, time_units)

if isempty(toi) == 1
    error('nf_extractpower:BadROI', 'toi is empty.');
end

if isnumeric(toi) == 0
    error('nf_extractpower:BadROI', 'toi must be numeric.');
end

toi = double(toi(:)');

tAxis = double(times(:));

if strcmpi(time_units, 'ms') == 1
    toi_use = toi ./ 1000;
elseif strcmpi(time_units, 's') == 1
    toi_use = toi;
elseif strcmpi(time_units, 'auto') == 1

    if max(abs(tAxis)) < 20 && max(abs(toi)) > 20
        toi_use = toi ./ 1000;
    else
        toi_use = toi;
    end

else
    error('nf_extractpower:BadArgs', 'time_units must be auto|ms|s.');
end

idx = local_resolve_axis_idx(tAxis, toi_use, 'time');

end

function [chanDim, freqDim, timeDim, obsDim] = local_detect_dims(X, TF, freqs, times, nObs, fname)

sz = size(X);
nd = ndims(X);

if nd < 3
    error('nf_extractpower:BadField', 'Field %s must be at least 3D.', fname);
end

nChan = 0;

if isfield(TF, 'chanlocs') == 1
    nChan = numel(TF.chanlocs);
end

nFreq = numel(freqs);
nTime = numel(times);

% Find dims by size matching (prefer NF convention if possible)
dimsChan = find(sz == nChan);
dimsFreq = find(sz == nFreq);
dimsTime = find(sz == nTime);
dimsObs  = find(sz == nObs);

if isempty(dimsObs) == 1
    error('nf_extractpower:BadField', 'Field %s has no dimension matching nObs=%d.', fname, nObs);
end

% obsDim: prefer last, then first, then the first match
if any(dimsObs == nd)
    obsDim = nd;
elseif any(dimsObs == 1)
    obsDim = 1;
else
    obsDim = dimsObs(1);
end

% freqDim: prefer 2 if matches, else first match not equal obsDim
if any(dimsFreq == 2) && obsDim ~= 2
    freqDim = 2;
else
    freqDim = local_first_not(dimsFreq, obsDim);
end

% timeDim: prefer 3 if matches, else first match not equal obsDim/freqDim
if any(dimsTime == 3) && obsDim ~= 3 && freqDim ~= 3
    timeDim = 3;
else
    timeDim = local_first_not(dimsTime, [obsDim freqDim]);
end

% chanDim: prefer 1 if matches, else first match not equal obsDim/freqDim/timeDim
if nChan > 0 && any(dimsChan == 1) && obsDim ~= 1 && freqDim ~= 1 && timeDim ~= 1
    chanDim = 1;
else
    chanDim = local_first_not(dimsChan, [obsDim freqDim timeDim]);
end

if isempty(freqDim) == 1 || isempty(timeDim) == 1
    error('nf_extractpower:BadField', 'Could not resolve freq/time dims for field %s.', fname);
end

if isempty(chanDim) == 1
    % If chanlocs missing, assume chanDim is the remaining non-obs/freq/time dim
    allDims = 1:nd;
    used = unique([obsDim freqDim timeDim]);
    rem = setdiff(allDims, used);

    if isempty(rem) == 0
        chanDim = rem(1);
    else
        chanDim = 1;
    end
end

end

function d = local_first_not(cands, banned)

d = [];

if isempty(cands) == 1
    return
end

if isempty(banned) == 1
    d = cands(1);
    return
end

for i = 1:numel(cands)
    if ~any(cands(i) == banned)
        d = cands(i);
        return
    end
end

end

function vals = local_extract_roi_mean(X, chanDim, freqDim, timeDim, obsDim, chanIdx, freqIdx, timeIdx, nObs)

nd = ndims(X);

idx = repmat({':'}, 1, nd);

idx{chanDim} = chanIdx;
idx{freqDim} = freqIdx;
idx{timeDim} = timeIdx;
idx{obsDim} = 1:size(X, obsDim);

Xroi = X(idx{:});

perm = 1:nd;

% bring obsDim to front
perm = [obsDim perm(perm ~= obsDim)];

Xp = permute(Xroi, perm);

szp = size(Xp);

if szp(1) ~= nObs
    % nObs mismatch; try to coerce if singleton behavior
    if szp(1) == 1 && nObs > 1
        Xp = repmat(Xp, [nObs ones(1, nd - 1)]);
        szp = size(Xp);
    else
        error('nf_extractpower:BadObs', 'ROI extraction produced obs dim %d but expected %d.', szp(1), nObs);
    end
end

Xp2 = reshape(double(Xp), [nObs prod(szp(2:end))]);

vals = local_nanmean_rows(Xp2);

end

function m = local_nanmean_rows(X)

n = size(X, 1);

m = nan(n, 1);

for i = 1:n

    row = X(i, :);

    ok = isfinite(row);

    if ~any(ok)
        m(i) = NaN;
    else
        m(i) = sum(row(ok)) ./ sum(ok);
    end

end

end

function nameOut = local_unique_fieldname(beh, nameIn)

nameOut = nameIn;

if isempty(beh) == 1
    return
end

if isfield(beh, nameOut) == 0
    return
end

k = 2;

while true

    cand = [nameIn '_' num2str(k)];

    if isfield(beh, cand) == 0
        nameOut = cand;
        return
    end

    k = k + 1;

end

end





