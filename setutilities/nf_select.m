function TFout = nf_select(TF, varargin)
% NF_SELECT  Select channels/frequencies/times/epochs from an NF TF struct.
%
% USAGE
%   TFout = nf_select(TF, 'channels', [1 2 3], 'freqs', [4 8 12], 'times', [-200 0 200])
%   TFout = nf_select(TF, 'channels', {'Fz','Cz'}, 'nofreqs', [60], 'notimes', [0])
%   TFout = nf_select(TF, 'epochs', 1:50)
%   TFout = nf_select(TF, 'noepochs', [1 2 3])
%
% INPUTS (name-value)
%   'channels'   : numeric indices OR cellstr labels to KEEP
%   'nochannels' : numeric indices OR cellstr labels to REMOVE
%   'freqs'      : numeric values OR integer indices to KEEP
%   'nofreqs'    : numeric values OR integer indices to REMOVE
%   'times'      : numeric values OR integer indices to KEEP
%   'notimes'    : numeric values OR integer indices to REMOVE
%   'epochs'     : integer indices to KEEP (4th dim, if present)
%   'noepochs'   : integer indices to REMOVE (4th dim, if present)
%
% NOTES
% - Applies selection to ALL top-level numeric fields containing:
%     'power', 'phase', or 'wpli'
% - Also supports TF.scale == 'stat' by subsetting stat maps in
%   TF.power / TF.phase / TF.wpli containers when present.
%

% -------------------------
% Parse inputs
% -------------------------
p = inputParser;

addRequired(p, 'TF', @isstruct);

addParameter(p, 'channels', [], @(x) local_is_chan_selector(x));
addParameter(p, 'nochannels', [], @(x) local_is_chan_selector(x));

addParameter(p, 'freqs', [], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'nofreqs', [], @(x) isnumeric(x) && isvector(x));

addParameter(p, 'times', [], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'notimes', [], @(x) isnumeric(x) && isvector(x));

addParameter(p, 'epochs', [], @(x) local_is_epoch_selector(x));
addParameter(p, 'noepochs', [], @(x) local_is_epoch_selector(x));

parse(p, TF, varargin{:});

TFout = TF;

% -------------------------
% Light validation
% -------------------------
try
    [~, ~] = nf_checkset(TFout, 'mode', 'warn'); 
catch
end

if isfield(TFout, 'scale')
    if ischar(TFout.scale) || isstring(TFout.scale)
        if strcmpi(char(TFout.scale), 'stat')
            TFout = local_select_stat_set(TFout, p.Results);
            return
        end
    end
end

if ~isfield(TFout, 'freqs') || ~isfield(TFout, 'times') || ~isfield(TFout, 'Fs')
    error('nf_select:MissingCore', 'TF must contain freqs, times, and Fs.');
end

if ~isfield(TFout, 'chanlocs')
    error('nf_select:MissingChanlocs', 'TF must contain chanlocs for channel selection.');
end

nChan = numel(TFout.chanlocs);
nFreq = numel(TFout.freqs);
nTime = numel(TFout.times);

% -------------------------
% Resolve channel indices
% -------------------------
chanKeep = (1:nChan);

if ~isempty(p.Results.channels)
    chanKeep = local_resolve_channels(TFout.chanlocs, p.Results.channels);
end

if ~isempty(p.Results.nochannels)
    chanDrop = local_resolve_channels(TFout.chanlocs, p.Results.nochannels);
    chanKeep = setdiff(chanKeep, chanDrop);
end

if isempty(chanKeep)
    error('nf_select:EmptySelection', 'Channel selection resulted in empty set.');
end

% -------------------------
% Resolve frequency indices
% -------------------------
freqKeep = (1:nFreq);

if ~isempty(p.Results.freqs)
    freqKeep = local_resolve_axis(TFout.freqs, p.Results.freqs, 'freqs');
end

if ~isempty(p.Results.nofreqs)
    freqDrop = local_resolve_axis(TFout.freqs, p.Results.nofreqs, 'freqs');
    freqKeep = setdiff(freqKeep, freqDrop);
end

if isempty(freqKeep)
    error('nf_select:EmptySelection', 'Frequency selection resulted in empty set.');
end

% -------------------------
% Resolve time indices
% -------------------------
timeKeep = (1:nTime);

if ~isempty(p.Results.times)
    timeKeep = local_resolve_axis(TFout.times, p.Results.times, 'times');
end

if ~isempty(p.Results.notimes)
    timeDrop = local_resolve_axis(TFout.times, p.Results.notimes, 'times');
    timeKeep = setdiff(timeKeep, timeDrop);
end

if isempty(timeKeep)
    error('nf_select:EmptySelection', 'Time selection resulted in empty set.');
end

% -------------------------
% Resolve epoch indices (4th dim)
% -------------------------
epochKeep = [];
epochSelectionRequested = false;

if ~isempty(p.Results.epochs)
    epochSelectionRequested = true;
end

if ~isempty(p.Results.noepochs)
    epochSelectionRequested = true;
end

nEpoch = local_infer_nEpoch(TFout,nChan,nFreq,nTime);
epochKeep = 1:nEpoch;
if epochSelectionRequested
    if isempty(nEpoch)
        error('nf_select:EpochsRequestedButMissing', ...
            'Requested epochs/noepochs but no eligible TF numeric field had a 4th (epoch) dimension.');
    end

    epochKeep = (1:nEpoch);

    if ~isempty(p.Results.epochs)
        epochKeep = local_resolve_epoch_axis(nEpoch, p.Results.epochs, 'epochs');
    end

    if ~isempty(p.Results.noepochs)
        epochDrop = local_resolve_epoch_axis(nEpoch, p.Results.noepochs, 'epochs');
        epochKeep = setdiff(epochKeep, epochDrop);
    end

    if isempty(epochKeep)
        error('nf_select:EmptySelection', 'Epoch selection resulted in empty set.');
    end
end

% -------------------------
% Apply to core axes
% -------------------------
if epochSelectionRequested || numel(timeKeep) ~= nTime
    TFout = hnf_event_grid(TFout,TF.times,TF.times(timeKeep),epochKeep,timeKeep);
end
if isfield(TFout,'nsensor'), TFout.nsensor = numel(chanKeep); end
if isfield(TFout,'phaseReference') && strcmp(TFout.phaseReference,'epoch-start-carrier')
    TFout.phaseReferenceSeconds = TFout.phaseReferenceSeconds(timeKeep);
end
if numel(timeKeep)>1 && all(diff(timeKeep)==timeKeep(2)-timeKeep(1))
    TFout.Fs = TFout.Fs/(timeKeep(2)-timeKeep(1));
end
TFout.freqs = TFout.freqs(freqKeep);
TFout.times = TFout.times(timeKeep);
TFout.chanlocs = TFout.chanlocs(chanKeep);

% If user keeps an epochs vector in the struct, subset it too (common patterns).
if epochSelectionRequested == true

    aligned = {'behavior','trlvec','trial_subjectindex','trial_condindex','trial_subid','trial_fileindex','trial_condlabel'};
    for i = 1:numel(aligned)
        name = aligned{i};
        if isfield(TFout,name) && numel(TFout.(name))==nEpoch
            TFout.(name) = TFout.(name)(epochKeep);
        end
    end

    if isfield(TFout, 'trialinfo')
        if isa(TFout.trialinfo, 'table')
            if height(TFout.trialinfo) == nEpoch
                TFout.trialinfo = TFout.trialinfo(epochKeep, :);
            end
        end
    end

    if isfield(TFout,'conds') && isscalar(TFout.conds) && TFout.conds==nEpoch && isfield(TFout,'ntrls') && TFout.ntrls==1
        TFout.conds = numel(epochKeep);
        if isfield(TFout,'trlerp'), TFout.trlerp = TFout.trlerp(epochKeep); end
    elseif isfield(TFout,'ntrls')
        TFout.ntrls = numel(epochKeep);
    end
end

% -------------------------
% Apply to all relevant numeric top-level fields
% -------------------------
fn = fieldnames(TFout);

for iF = 1:numel(fn)

    name = fn{iF};
    val = TFout.(name);

    if ~isnumeric(val)
        continue
    end

    lname = lower(name);

    if isempty(regexp(lname,'^(power|phase|itpc|wpli|complex)($|_)','once')), continue; end
    [X,compact] = hnf_tf_grid(val,nChan,nFreq,nTime,nEpoch,name);
    X = X(chanKeep,freqKeep,timeKeep,epochKeep);
    if compact, X = squeeze(X); end
    TFout.(name) = X;

end

disp('[nf_select]: Time-Frequency data selection complete.');

end

% ======================================================================
% Local helpers
% ======================================================================

function tf = local_is_chan_selector(x)

tf = false;

if isempty(x)
    tf = true;
    return
end

if isnumeric(x) && isvector(x)
    tf = true;
    return
end

if iscell(x)
    if all(cellfun(@(c) ischar(c) || isstring(c), x))
        tf = true;
        return
    end
end

end

function tf = local_is_epoch_selector(x)

tf = false;

if isempty(x)
    tf = true;
    return
end

if isnumeric(x) && isvector(x)
    if all(x == round(x))
        tf = true;
        return
    end
end

end

function idx = local_resolve_channels(chanlocs, sel)

labels = {chanlocs.labels};

if isnumeric(sel)

    idx = unique(sel(:)');

    if any(idx < 1) || any(idx > numel(labels))
        error('nf_select:BadChannels', 'Channel indices out of range.');
    end

elseif iscell(sel)

    sel = cellfun(@char, sel, 'UniformOutput', false);

    idx = [];

    for i = 1:numel(sel)

        m = find(strcmpi(labels, sel{i}));

        if isempty(m)
            error('nf_select:BadChannels', 'Channel label not found: %s', sel{i});
        end

        idx = [idx m]; %#ok<AGROW>

    end

    idx = unique(idx);

else

    error('nf_select:BadChannels', 'channels selector must be numeric or cellstr.');

end

end

function idx = local_resolve_axis(axisVec, req, label)

axisVec = double(axisVec(:));
req = double(req(:));

n = numel(axisVec);

% Decide: indices vs values
isIndexCandidate = false;

if all(req == round(req))
    if all(req >= 1) && all(req <= n)
        isIndexCandidate = true;
    end
end

if isIndexCandidate == true

    % Disambiguate: if requested integers look like actual axis VALUES, treat as values
    % Otherwise treat as indices.
    tol = local_axis_tol(axisVec);

    matchCount = 0;

    for i = 1:numel(req)
        d = min(abs(axisVec - req(i)));
        if d <= tol
            matchCount = matchCount + 1;
        end
    end

    if matchCount < numel(req)
        idx = unique(req(:)') ;
        return
    end

end

% Value-based matching (nearest within tolerance)
tol = local_axis_tol(axisVec);

idx = [];

for i = 1:numel(req)

    [d, k] = min(abs(axisVec - req(i)));

    if d > tol
        error('nf_select:BadAxis', 'Requested %s value not found (within tol): %.6f', label, req(i));
    end

    idx = [idx k]; %#ok<AGROW>

end

idx = unique(idx);

end

function tol = local_axis_tol(axisVec)

if numel(axisVec) < 2
    tol = 1e-12;
    return
end

d = abs(diff(axisVec));
d = d(~isnan(d));
d = d(d > 0);

if isempty(d)
    tol = 1e-12;
    return
end

tol = min(d) / 2;

end

function idx = local_resolve_epoch_axis(nEpoch, req, label)

req = double(req(:));

idx = unique(req(:)');

if any(idx < 1) || any(idx > nEpoch)
    error('nf_select:BadEpochs', 'Requested %s indices out of range (1..%d).', label, nEpoch);
end

end

function nEpoch = local_infer_nEpoch(TF, nChan, nFreq, nTime)

nEpoch = [];
names = fieldnames(TF);
for i = 1:numel(names)
    if isnumeric(TF.(names{i})) && ~isempty(regexpi(names{i},'^(power|phase|itpc|wpli|complex)($|_)','once'))
        nEpoch = numel(TF.(names{i}))/(nChan*nFreq*nTime);
        if nEpoch~=fix(nEpoch), error('nf_select:BadShape','Invalid TF data dimensions.'); end
        return
    end
end
end

function TFout = local_select_stat_set(TFout, opts)

% Minimal stat subsetting:
% - requires freqs/times/chanlocs to exist if you want those selections.
% - supports epochs/noepochs when a 4th dim exists in stat maps.

if ~isfield(TFout, 'chanlocs')
    error('nf_select:StatMissingChanlocs', 'Stat set missing chanlocs.');
end

nChan = numel(TFout.chanlocs);

chanKeep = (1:nChan);

if ~isempty(opts.channels)
    chanKeep = local_resolve_channels(TFout.chanlocs, opts.channels);
end

if ~isempty(opts.nochannels)
    chanDrop = local_resolve_channels(TFout.chanlocs, opts.nochannels);
    chanKeep = setdiff(chanKeep, chanDrop);
end

if isempty(chanKeep)
    error('nf_select:EmptySelection', 'Channel selection resulted in empty set.');
end

freqKeep = [];
timeKeep = [];

if isfield(TFout, 'freqs') && ~isempty(TFout.freqs)
    nFreq = numel(TFout.freqs);
    freqKeep = (1:nFreq);

    if ~isempty(opts.freqs)
        freqKeep = local_resolve_axis(TFout.freqs, opts.freqs, 'freqs');
    end

    if ~isempty(opts.nofreqs)
        freqDrop = local_resolve_axis(TFout.freqs, opts.nofreqs, 'freqs');
        freqKeep = setdiff(freqKeep, freqDrop);
    end

    TFout.freqs = TFout.freqs(freqKeep);
end

if isfield(TFout, 'times') && ~isempty(TFout.times)
    nTime = numel(TFout.times);
    timeKeep = (1:nTime);

    if ~isempty(opts.times)
        timeKeep = local_resolve_axis(TFout.times, opts.times, 'times');
    end

    if ~isempty(opts.notimes)
        timeDrop = local_resolve_axis(TFout.times, opts.notimes, 'times');
        timeKeep = setdiff(timeKeep, timeDrop);
    end

    TFout.times = TFout.times(timeKeep);
end

TFout.chanlocs = TFout.chanlocs(chanKeep);

if isempty(freqKeep) || isempty(timeKeep)
    error('nf_select:EmptySelection','Stat selection must contain frequencies and times.');
end
% nf_stregress writes metric arrays directly under each measure and primary
% aliases at the root. Older sets have one more level of containers.
group = isfield(TFout,'aggtype') && strcmpi(TFout.aggtype,'stat');
[TFout,axis4] = local_stat_maps(TFout,[nChan nFreq nTime],chanKeep,freqKeep,timeKeep,opts,group);
if ~isempty(opts.epochs) || ~isempty(opts.noepochs)
    if isempty(axis4), error('nf_select:EpochDimMissing','No selectable fourth map axis.'); end
    if isfield(TFout,'prednames') && ~isempty(TFout.prednames)
        TFout.prednames = TFout.prednames(axis4);
        if isfield(TFout,'regressors'), TFout.regressors = TFout.regressors(:,axis4); end
    else
        if isfield(TFout,'epochs'), TFout.epochs = TFout.epochs(axis4); end
        if isfield(TFout,'trialinfo'), TFout.trialinfo = TFout.trialinfo(axis4,:); end
        if isfield(TFout,'ntrls'), TFout.ntrls = numel(axis4); end
    end
end
end

function [S,selected] = local_stat_maps(S,dims,c,f,t,opts,group)
selected = [];
names = fieldnames(S);
metrics = {'corrcoef','tstat','p_vals','tfce_tstat','beta'};
for i = 1:numel(names)
    name = names{i};
    X = S.(name);
    if isstruct(X) && isscalar(X) && ~strcmp(name,'chanlocs')
        [S.(name),child] = local_stat_maps(X,dims,c,f,t,opts,group);
        if ~isempty(child), selected = child; end
    elseif isnumeric(X) && any(strcmp(name,metrics)) && ~isempty(X)
        offset = double(group);
        shape = [size(X,1+offset),size(X,2+offset),size(X,3+offset)];
        if ~isequal(shape,dims), error('nf_select:BadShape','Stat map %s disagrees with its axes.',name); end
        n = size(X,4+offset);
        kept = 1:n;
        if ~isempty(opts.epochs), kept = local_resolve_epoch_axis(n,opts.epochs,'epochs'); end
        if ~isempty(opts.noepochs), kept = setdiff(kept,local_resolve_epoch_axis(n,opts.noepochs,'noepochs')); end
        if isempty(kept), error('nf_select:EmptySelection','No observations or predictors remain.'); end
        selected = kept;
        if group, S.(name) = X(:,c,f,t,kept); else, S.(name) = X(c,f,t,kept); end
    end
end
end
