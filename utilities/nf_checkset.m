function [TF, report] = nf_checkset(TF, varargin)
% NF_CHECKSET  Validate an NF TF struct for internal consistency.
%
% USAGE
%   [TF, report] = nf_checkset(TF)
%   [TF, report] = nf_checkset(TF, 'mode', 'warn')
%
% OPTIONS
%   'mode'   : 'error' (default), 'warn', or 'silent'
%   'fix'    : false (default). If true, applies safe reshapes for 1ch legacy
%
% OUTPUT
%   TF      : returned TF (possibly minimally fixed if 'fix' true)
%   report  : struct with fields:
%             .errors, .warnings, .setType, .powerFields, .phaseFields, .wpliFields
%
% Checks (high-level)
% - Required core fields (for non-stat): freqs, times, Fs, chanlocs
% - Discovers ALL top-level numeric power/phase/wpli fields and ensures
%   consistent dimensionality with freqs/times/chanlocs.
% - Checks behavior length requirements when trials/conds exist.
% - Checks event struct validity when present.
% - Detects set type: stat, avg, singletrial, continuous
%

% -------------------------
% Parse options
% -------------------------
p = inputParser;

addRequired(p, 'TF', @isstruct);

addParameter(p, 'mode', 'error', @(x) ischar(x) || isstring(x));
addParameter(p, 'fix', false, @(x) islogical(x) && isscalar(x));

parse(p, TF, varargin{:});

mode = lower(char(p.Results.mode));
doFix = p.Results.fix;

report = struct();
report.errors = {};
report.warnings = {};
report.setType = '';
report.powerFields = {};
report.phaseFields = {};
report.wpliFields = {};

% -------------------------
% Helper lambdas
% -------------------------
addErr = @(msg) local_add_msg('error', msg);
addWarn = @(msg) local_add_msg('warn', msg);

% -------------------------
% Detect stat sets
% -------------------------
if isfield(TF, 'scale')
    if ischar(TF.scale) || isstring(TF.scale)
        if strcmpi(char(TF.scale), 'stat')
            report.setType = 'stat';
            local_check_stat();
            local_finalize();
            return
        end
    end
end

% -------------------------
% Non-stat core checks
% -------------------------
local_require_field('freqs');
local_require_field('times');
local_require_field('Fs');
local_require_field('chanlocs');

if isempty(report.errors)
    if ~isnumeric(TF.freqs) || ~isvector(TF.freqs)
        addErr('TF.freqs must be a numeric vector.');
    end
    if ~isnumeric(TF.times) || ~isvector(TF.times)
        addErr('TF.times must be a numeric vector.');
    end
    if ~isnumeric(TF.Fs) || ~isscalar(TF.Fs) || TF.Fs <= 0
        addErr('TF.Fs must be a positive scalar.');
    end
    if ~isstruct(TF.chanlocs)
        addErr('TF.chanlocs must be a struct array.');
    end
end

if isempty(report.errors)
    if any(diff(double(TF.freqs(:))) <= 0)
        addWarn('TF.freqs is not strictly increasing.');
    end
    if any(diff(double(TF.times(:))) == 0)
        addWarn('TF.times contains duplicates.');
    end
end

% -------------------------
% Discover measured fields
% -------------------------
[powerFields, phaseFields, wpliFields] = local_discover_fields(TF);

report.powerFields = powerFields;
report.phaseFields = phaseFields;
report.wpliFields = wpliFields;

if isempty(powerFields) && isempty(phaseFields) && isempty(wpliFields)
    addErr('No numeric top-level power/phase/wpli fields were found.');
end

if ~isempty(report.errors)
    local_finalize();
    return
end

nChan = numel(TF.chanlocs);
nFreq = numel(TF.freqs);
nTime = numel(TF.times);

% chanlocs labels
if ~isfield(TF.chanlocs, 'labels')
    addErr('TF.chanlocs missing labels field.');
else
    labs = {TF.chanlocs.labels};
    if numel(labs) ~= nChan
        addErr('TF.chanlocs labels count mismatch.');
    end
    if numel(unique(labs)) ~= numel(labs)
        addWarn('TF.chanlocs labels are not unique.');
    end
end

% -------------------------
% Validate all measured fields
% -------------------------
fields = [powerFields; phaseFields; wpliFields];

baseShape = [];

for iF = 1:numel(fields)

    name = fields{iF};
    X = TF.(name);

    if ~isnumeric(X)
        addErr(['Field ' name ' is not numeric.']);
        continue
    end

    if ndims(X) < 3
        addErr(['Field ' name ' has <3 dims; expected ch x f x t x ...']);
        continue
    end

    % Allow fix: if user stored 1ch as f x t x trials
    if doFix == true
        if ndims(X) == 3
            if size(X, 1) == nFreq && size(X, 2) == nTime
                X = reshape(X, 1, size(X, 1), size(X, 2), size(X, 3));
                TF.(name) = X;
                addWarn(['Fixed 1ch legacy shape for field ' name '.']);
            end
        end
    end

    if size(X, 1) ~= nChan
        addErr(['Field ' name ' dim1 must match chanlocs (' num2str(nChan) ').']);
        continue
    end

    if size(X, 2) ~= nFreq
        addErr(['Field ' name ' dim2 must match freqs (' num2str(nFreq) ').']);
        continue
    end

    if size(X, 3) ~= nTime
        addErr(['Field ' name ' dim3 must match times (' num2str(nTime) ').']);
        continue
    end

    if isempty(baseShape)
        baseShape = size(X);
    else
        if size(X, 1) ~= baseShape(1)
            addErr(['Field ' name ' channel dim mismatch vs other fields.']);
        end
        if size(X, 2) ~= baseShape(2)
            addErr(['Field ' name ' freq dim mismatch vs other fields.']);
        end
        if size(X, 3) ~= baseShape(3)
            addErr(['Field ' name ' time dim mismatch vs other fields.']);
        end
    end

    if any(~isfinite(X(:)))
        addWarn(['Field ' name ' contains NaN/Inf values.']);
    end

end

% -------------------------
% Determine setType (avg / singletrial / continuous)
% -------------------------
report.setType = local_detect_type_from_fields(TF, powerFields, phaseFields, wpliFields);

% Behavior checks
if strcmpi(report.setType, 'singletrial')
    nTr = local_infer_trials(TF, powerFields, phaseFields, wpliFields);
    if ~isfield(TF, 'behavior')
        addErr('Single-trial TF set requires TF.behavior with one entry per trial.');
    else
        if ~isstruct(TF.behavior)
            addErr('TF.behavior must be a struct array.');
        else
            if numel(TF.behavior) ~= nTr
                addErr(['TF.behavior length must equal trials (' num2str(nTr) ').']);
            end
        end
    end
elseif strcmpi(report.setType, 'avg')
    if isfield(TF, 'conds')
        if isfield(TF, 'behavior')
            if isstruct(TF.behavior)
                if numel(TF.behavior) ~= TF.conds
                    addErr('Avg TF: TF.behavior length must match TF.conds.');
                end
            else
                addErr('Avg TF: TF.behavior must be struct array if present.');
            end
        end
    end
end

% Event checks (optional)
if isfield(TF, 'event')
    if ~isstruct(TF.event)
        addErr('TF.event must be a struct array.');
    else
        if numel(TF.event) > 0
            if ~isfield(TF.event, 'type')
                addWarn('TF.event missing field "type".');
            end
            if ~isfield(TF.event, 'latency')
                addWarn('TF.event missing field "latency".');
            else
                lat = [TF.event.latency];
                lat = double(lat);
                if any(~isfinite(lat))
                    addWarn('TF.event latency contains NaN/Inf.');
                else
                    if any(lat < 1) || any(lat > nTime)
                        addWarn('Some TF.event.latency values fall outside [1, nTime].');
                    end
                end
            end
        end
    end
end

local_finalize();
return

% ======================================================================
% Nested helper: stat set checks
% ======================================================================
    function local_check_stat()

        containers = {'power','phase','wpli'};
        found = false;

        for iC = 1:numel(containers)

            cname = containers{iC};

            if ~isfield(TF, cname)
                continue
            end

            if ~isstruct(TF.(cname))
                continue
            end

            meas = fieldnames(TF.(cname));

            for iM = 1:numel(meas)

                m = meas{iM};
                S = TF.(cname).(m);

                if ~isstruct(S)
                    continue
                end

                if isfield(S, 'corrcoef') && isfield(S, 'tstat') && isfield(S, 'p_vals')
                    found = true;
                end

            end

        end

        if found == false
            addErr('Stat TF set: could not find any measures with corrcoef/tstat/p_vals in power/phase/wpli containers.');
        end

    end

% ======================================================================
% Nested helper: finalize behavior
% ======================================================================
    function local_finalize()

        if strcmpi(mode, 'silent') ~= 1

            if ~isempty(report.errors)
                disp('nf_checkset: ERRORS');
                for i = 1:numel(report.errors)
                    disp(['  - ' report.errors{i}]);
                end
            end

            if ~isempty(report.warnings)
                disp('nf_checkset: WARNINGS');
                for i = 1:numel(report.warnings)
                    disp(['  - ' report.warnings{i}]);
                end
            end

            disp(['nf_checkset: setType = ' report.setType]);

        end

        if strcmpi(mode, 'error') == 1
            if ~isempty(report.errors)
                error('nf_checkset:InvalidSet', report.errors{1});
            end
        end

    end

% ======================================================================
% Nested helper: add message
% ======================================================================
    function local_add_msg(kind, msg)

        if strcmpi(kind, 'error')
            report.errors{end + 1, 1} = msg;
        else
            report.warnings{end + 1, 1} = msg;
        end

    end

% ======================================================================
% Nested helper: require field
% ======================================================================
    function local_require_field(fname)

        if ~isfield(TF, fname)
            addErr(['Missing required field: ' fname]);
        end

    end

end

% ======================================================================
% Local helpers (non-nested)
% ======================================================================

function [powerFields, phaseFields, wpliFields] = local_discover_fields(TF)

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

function setType = local_detect_type_from_fields(TF, powerFields, phaseFields, wpliFields)

setType = 'continuous';

fields = [powerFields; phaseFields; wpliFields];

hasTrials = false;

for i = 1:numel(fields)

    X = TF.(fields{i});

    if ndims(X) >= 4
        if size(X, 4) > 1
            hasTrials = true;
        end
    end

end

if hasTrials == true
    setType = 'singletrial';
end

if isfield(TF, 'ntrls') && isfield(TF, 'conds')
    if isnumeric(TF.ntrls) && TF.ntrls == 1
        % check if looks like averaged: last dim equals conds for first power field
        if ~isempty(powerFields)
            X = TF.(powerFields{1});
            if ndims(X) == 4
                if size(X, 4) == TF.conds
                    setType = 'avg';
                end
            end
            if ndims(X) == 3
                if size(X, 3) == TF.conds
                    setType = 'avg';
                end
            end
        end
    end
end

end

function nTr = local_infer_trials(TF, powerFields, phaseFields, wpliFields)

nTr = [];

fields = [powerFields; phaseFields; wpliFields];

for i = 1:numel(fields)

    X = TF.(fields{i});

    if ndims(X) >= 4
        nTrNow = size(X, 4);
        if isempty(nTr)
            nTr = nTrNow;
        else
            if nTrNow ~= nTr
                error('nf_checkset:TrialMismatch', 'Trial counts differ across measured fields.');
            end
        end
    end

end

if isempty(nTr)
    nTr = 1;
end

end










