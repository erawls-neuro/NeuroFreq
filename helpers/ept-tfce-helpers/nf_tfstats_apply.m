function out = nf_tfstats_apply(tfIn, runner, varargin)

groups = {'power', 'phase'};
slices = [];
labelMode = 'auto';
metric = 'corrcoeff';

if nargin < 2
    error('nf_tfstats_apply:BadArgs', 'Requires tfIn and runner.');
end

if isstruct(tfIn) == 0
    error('nf_tfstats_apply:BadArgs', 'tfIn must be a struct.');
end

if isa(runner, 'function_handle') == 0
    error('nf_tfstats_apply:BadArgs', 'runner must be a function handle.');
end

if mod(numel(varargin), 2) ~= 0
    error('nf_tfstats_apply:BadArgs', 'Arguments must be name/value pairs.');
end

for i = 1:2:numel(varargin)
    key = varargin{i};
    val = varargin{i + 1};

    if isstring(key) == 1
        key = char(key);
    end

    if ischar(key) == 0
        error('nf_tfstats_apply:BadArgs', 'Option names must be strings.');
    end

    key = lower(key);

    if strcmp(key, 'groups') == 1
        groups = val;
    elseif strcmp(key, 'slices') == 1
        slices = val;
    elseif strcmp(key, 'labelmode') == 1
        labelMode = val;
    elseif strcmp(key, 'metric') == 1
        if isstring(val) == 1
            val = char(val);
        end
        metric = val;
    else
        error('nf_tfstats_apply:BadArgs', 'Unknown option: %s', key);
    end
end

if isstring(groups) == 1
    groups = cellstr(groups(:));
end

if iscell(groups) == 0
    error('nf_tfstats_apply:BadArgs', 'groups must be a cell array of strings.');
end

for iG = 1:numel(groups)
    if isstring(groups{iG}) == 1
        groups{iG} = char(groups{iG});
    end

    if ischar(groups{iG}) == 0
        error('nf_tfstats_apply:BadArgs', 'groups must contain strings.');
    end
end

if isempty(slices) == 0
    if isnumeric(slices) == 0
        error('nf_tfstats_apply:BadArgs', 'slices must be numeric.');
    end
end

if isstring(labelMode) == 1
    labelMode = char(labelMode);
end

if ischar(labelMode) == 0
    error('nf_tfstats_apply:BadArgs', 'labelMode must be a string.');
end

if isstring(metric) == 1
    metric = char(metric);
end

if ischar(metric) == 0
    error('nf_tfstats_apply:BadArgs', 'metric must be a string.');
end

scale = local_detect_scale(tfIn);

out = struct();
out.scale_in = scale;
out.results = struct();

for iG = 1:numel(groups)

    gtoken = groups{iG};
    gkey = matlab.lang.makeValidName(gtoken);

    candidates = local_select_candidates(tfIn, gtoken);

    if isempty(candidates) == 1
        continue
    end

    if isfield(out.results, gkey) == 0
        out.results.(gkey) = struct();
    end

    for iC = 1:numel(candidates)

        fn = candidates(iC).name;
        Xraw = candidates(iC).value;

        X = [];
        metricFound = '';
        outKey = '';
        outPath = '';

        if isnumeric(Xraw) == 1

            if isempty(Xraw) == 1
                continue
            end

            X = Xraw;

            outPath = fn;
            outKey = matlab.lang.makeValidName(fn);

        elseif isstruct(Xraw) == 1

            [ok, realMetricName] = local_find_field_ci(Xraw, metric);

            if ok == 0
                continue
            end

            Xm = Xraw.(realMetricName);

            if isnumeric(Xm) == 0
                continue
            end

            if isempty(Xm) == 1
                continue
            end

            X = Xm;
            metricFound = realMetricName;

            outPath = [fn '.' realMetricName];
            outKey = matlab.lang.makeValidName([fn '__' realMetricName]);

        else

            continue

        end

        [sliceCount, sliceLabels] = local_slice_info(tfIn, X, scale, labelMode);

        useSlices = 1:sliceCount;

        if isempty(slices) == 0
            useSlices = slices(:)';

            if any(useSlices < 1) == 1 || any(useSlices > sliceCount) == 1
                error('nf_tfstats_apply:BadSlice', 'Slice index out of bounds for field %s.', fn);
            end
        end

        res = struct();
        res.group = gtoken;
        res.field = fn;
        res.metric = metricFound;
        res.path = outPath;
        res.scale_in = scale;

        res.slice_labels = sliceLabels;
        res.slice_index = useSlices(:);
        res.tfce = cell(numel(useSlices), 1);

        for iS = 1:numel(useSlices)

            sIdx = useSlices(iS);

            Xslice = local_get_slice(X, scale, sIdx, sliceCount);
            Xobs = local_obs_first(Xslice, scale);

            res.tfce{iS} = runner(Xobs);

        end

        outKey = local_make_unique_key(out.results.(gkey), outKey);
        out.results.(gkey).(outKey) = res;

    end
end

end

% ---------------------------
% Local helpers
% ---------------------------

function scale = local_detect_scale(tf)

scale = '';

if isfield(tf, 'scale') == 1
    if ischar(tf.scale) == 1
        scale = tf.scale;
    end
end

if isempty(scale) == 1
    if isfield(tf, 'aggtype') == 1
        if ischar(tf.aggtype) == 1
            scale = tf.aggtype;
        end
    end
end

if isempty(scale) == 1
    if isfield(tf, 'trial_subjectindex') == 1
        scale = 'singletrial';
    elseif isfield(tf, 'conds') == 1
        scale = 'avg';
    elseif isfield(tf, 'stat') == 1
        scale = 'stat';
    else
        scale = 'avg';
    end
end

end

function candidates = local_select_candidates(tf, gtoken)

candidates = struct('name', {}, 'value', {});
fns = fieldnames(tf);

if isstring(gtoken) == 1
    gtoken = char(gtoken);
end

if isempty(gtoken) == 1
    return
end

wantPower = false;
wantPhase = false;

if strcmpi(gtoken, 'power') == 1
    wantPower = true;
end

if strcmpi(gtoken, 'phase') == 1
    wantPhase = true;
end

k = 0;

% NOTE:
% For gtoken=='power' or 'phase', we DO NOT allow an exact-match early return,
% because NF convention is top-level prefix fields (power*, pow_, phase*, ph_).
% We include an exact-match field if present, but we still scan for siblings.

if isfield(tf, gtoken) == 1
    if wantPower == false && wantPhase == false
        k = 1;
        candidates(k).name = gtoken;
        candidates(k).value = tf.(gtoken);
        return
    else
        k = k + 1;
        candidates(k).name = gtoken;
        candidates(k).value = tf.(gtoken);
    end
end

for i = 1:numel(fns)

    fn = fns{i};

    if strcmpi(fn, gtoken) == 1
        continue
    end

    % Skip obvious header/meta fields
    if strcmpi(fn, 'freqs') == 1
        continue
    end

    if strcmpi(fn, 'times') == 1
        continue
    end

    if strcmpi(fn, 'chanlocs') == 1
        continue
    end

    if strcmpi(fn, 'conds') == 1
        continue
    end

    if strcmpi(fn, 'scale') == 1
        continue
    end

    if strcmpi(fn, 'aggtype') == 1
        continue
    end

    if strcmpi(fn, 'stat') == 1
        continue
    end

    isMatch = false;

    if wantPower == true
        if strncmpi(fn, 'power', 5) == 1
            isMatch = true;
        end

        if strncmpi(fn, 'pow_', 4) == 1
            isMatch = true;
        end
    elseif wantPhase == true
        if strncmpi(fn, 'phase', 5) == 1
            isMatch = true;
        end

        if strncmpi(fn, 'ph_', 3) == 1
            isMatch = true;
        end
    else
        if strncmpi(fn, gtoken, numel(gtoken)) == 1
            isMatch = true;
        end
    end

    if isMatch == false
        continue
    end

    v = tf.(fn);

    if isnumeric(v) == 1 || isstruct(v) == 1
        k = k + 1;
        candidates(k).name = fn;
        candidates(k).value = v;
    end

end

end

function [ok, realName] = local_find_field_ci(S, want)

ok = 0;
realName = '';

if isstruct(S) == 0
    return
end

if ischar(want) == 0
    return
end

wantList = {want};

if strcmpi(want, 'corrcoeff') == 1
    wantList{end + 1} = 'corrcoef';
end

fns = fieldnames(S);

for iW = 1:numel(wantList)

    w = wantList{iW};

    for i = 1:numel(fns)

        nm = fns{i};

        if strcmpi(nm, w) == 1
            ok = 1;
            realName = nm;
            return
        end

    end

end

end

function [sliceCount, sliceLabels] = local_slice_info(tf, X, scale, labelMode)

sliceCount = 1;
sliceLabels = {''};

if isnumeric(X) == 0
    return
end

nd = ndims(X);

if strcmpi(scale, 'singletrial') == 1

    if nd == 5
        sliceCount = size(X, 5);
        sliceLabels = local_default_slice_labels(sliceCount);
    end

    return

end

if strcmpi(scale, 'avg') == 1

    if nd == 5
        sliceCount = size(X, 5);
        sliceLabels = local_cond_labels(tf, sliceCount, labelMode);
        return
    end

    if nd == 4
        sliceCount = size(X, 4);
        sliceLabels = local_cond_labels(tf, sliceCount, labelMode);
        return
    end

    return

end

if strcmpi(scale, 'stat') == 1

    if nd == 5
        sliceCount = size(X, 5);
        sliceLabels = local_stat_labels(tf, sliceCount, labelMode);
        return
    end

    % Support non-aggregated stat maps: C x F x T x K
    if nd == 4
        if size(X, 4) > 1
            sliceCount = size(X, 4);
            sliceLabels = local_stat_labels(tf, sliceCount, labelMode);
            return
        end
    end

    return

end

end

function labels = local_cond_labels(tf, sliceCount, labelMode)

labels = local_default_slice_labels(sliceCount);

if strcmpi(labelMode, 'auto') == 1
    if isfield(tf, 'conds') == 1
        if iscell(tf.conds) == 1
            if numel(tf.conds) == sliceCount
                labels = tf.conds(:);
                return
            end
        end
    end
end

end

function labels = local_stat_labels(tf, sliceCount, labelMode)

labels = local_default_slice_labels(sliceCount);

if strcmpi(labelMode, 'auto') == 1
    if isfield(tf, 'prednames') == 1
        if iscell(tf.prednames) == 1
            if numel(tf.prednames) == sliceCount
                labels = tf.prednames(:);
                return
            end
        end
    end

    if isfield(tf, 'conds') == 1
        if iscell(tf.conds) == 1
            if numel(tf.conds) == sliceCount
                labels = tf.conds(:);
                return
            end
        end
    end
end

end

function labels = local_default_slice_labels(sliceCount)

labels = cell(sliceCount, 1);

for i = 1:sliceCount
    labels{i} = ['slice' num2str(i)];
end

end

function Xslice = local_get_slice(X, scale, sIdx, sliceCount)

Xslice = X;

if isnumeric(X) == 0
    return
end

nd = ndims(X);

if strcmpi(scale, 'singletrial') == 1
    if nd == 5
        if sIdx < 1 || sIdx > sliceCount
            error('nf_tfstats_apply:BadSlice', 'Slice index out of bounds.');
        end

        Xslice = X(:, :, :, :, sIdx);
    end

    return
end

if strcmpi(scale, 'avg') == 1

    if nd == 5
        if sIdx < 1 || sIdx > sliceCount
            error('nf_tfstats_apply:BadSlice', 'Slice index out of bounds.');
        end

        Xslice = X(:, :, :, :, sIdx);
        return
    end

    if nd == 4
        if sIdx < 1 || sIdx > sliceCount
            error('nf_tfstats_apply:BadSlice', 'Slice index out of bounds.');
        end

        Xslice = X(:, :, :, sIdx);
        return
    end

    return

end

if strcmpi(scale, 'stat') == 1

    if nd == 5
        if sIdx < 1 || sIdx > sliceCount
            error('nf_tfstats_apply:BadSlice', 'Slice index out of bounds.');
        end

        Xslice = X(:, :, :, :, sIdx);
        return
    end

    if nd == 4
        if sliceCount > 1
            if sIdx < 1 || sIdx > sliceCount
                error('nf_tfstats_apply:BadSlice', 'Slice index out of bounds.');
            end

            Xslice = X(:, :, :, sIdx);
            return
        end
    end

    return

end

end

function Xobs = local_obs_first(X, scale)

Xobs = X;

if isnumeric(X) == 0
    return
end

nd = ndims(X);

% CRITICAL:
% Only singletrial needs re-ordering (trial -> obs-first).
% avg/stat are assumed to already be in the correct orientation for group stats.
% Do NOT insert singleton obs dimensions (that breaks TFCE channel checks).

if strcmpi(scale, 'singletrial') == 1

    if nd == 4
        Xobs = permute(X, [4 1 2 3]);
        return
    end

    if nd == 5
        Xobs = permute(X, [4 1 2 3 5]);
        return
    end

    error('nf_tfstats_apply:BadShape', 'singletrial data must be 4D or 5D.');

end

end

function keyOut = local_make_unique_key(S, keyIn)

keyOut = keyIn;

if isfield(S, keyOut) == 0
    return
end

idx = 2;

while true

    cand = [keyIn '_' num2str(idx)];

    if isfield(S, cand) == 0
        keyOut = cand;
        return
    end

    idx = idx + 1;

end

end






