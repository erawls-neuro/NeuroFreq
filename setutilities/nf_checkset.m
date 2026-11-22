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
    if ~isnumeric(TF.freqs) || ~isvector(TF.freqs) || any(~isfinite(TF.freqs(:)))
        addErr('TF.freqs must be a numeric vector.');
    end
    if ~isnumeric(TF.times) || ~isvector(TF.times) || any(~isfinite(TF.times(:)))
        addErr('TF.times must be a numeric vector.');
    end
    if ~isnumeric(TF.Fs) || ~isscalar(TF.Fs) || ~isreal(TF.Fs) || ~isfinite(TF.Fs) || TF.Fs <= 0
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

    nObs = numel(X)/(nChan*nFreq*nTime);
    try
        [X,compact] = hnf_tf_grid(X,nChan,nFreq,nTime,nObs,name);
    catch err
        addErr(err.message);
        continue
    end
    if doFix && compact, TF.(name) = X; end
    if isempty(baseShape)
        baseShape = [nChan nFreq nTime nObs];
    elseif nObs~=baseShape(4)
        addErr(['Observation count differs in ' name '.']);
    end

    if any(~isfinite(X(:)))
        addWarn(['Field ' name ' contains NaN/Inf values.']);
    end

end

% -------------------------
% Determine setType (avg / singletrial / continuous)
% -------------------------
if isempty(baseShape), local_finalize(); return; end
nTr = baseShape(4);
report.setType = 'continuous';
if isfield(TF,'conds') && isfield(TF,'ntrls') && isequal(TF.ntrls,1) && isequal(TF.conds,nTr)
    report.setType = 'avg';
elseif nTr>1 || (isfield(TF,'epoch') && ~isempty(TF.epoch))
    report.setType = 'singletrial';
end
if strcmpi(report.setType,'singletrial')
    if isfield(TF,'behavior') && (~isstruct(TF.behavior) || numel(TF.behavior)~=nTr)
        addErr('TF.behavior must contain one structure per trial when present.');
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
                    if any(lat < 1) || any(lat > nTime*nTr)
                        addWarn('Some TF.event.latency values fall outside the concatenated sample range.');
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

        containers = fieldnames(TF);
        found = false;

        for iC = 1:numel(containers)

            cname = containers{iC};

            if ~isfield(TF, cname)
                continue
            end

            if ~isstruct(TF.(cname)) || ~isscalar(TF.(cname))
                continue
            end

            if all(isfield(TF.(cname),{'corrcoef','tstat','p_vals'}))
                found = true;
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

    if ~isempty(regexp(lname,'^(power|pow)($|_)','once'))
        powerFields{end + 1, 1} = name; %#ok<AGROW>
    elseif ~isempty(regexp(lname,'^(phase|ph|itpc)($|_)','once'))
        phaseFields{end + 1, 1} = name; %#ok<AGROW>
    elseif ~isempty(regexp(lname,'^wpli($|_)','once'))
        wpliFields{end + 1, 1} = name; %#ok<AGROW>
    end

end

powerFields = sort(powerFields);
phaseFields = sort(phaseFields);
wpliFields = sort(wpliFields);

end
