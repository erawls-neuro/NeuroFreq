function behav = nf_makebehavior(tab, lags)
% NF_MAKEBEHAVIOR    Make a NeuroFreq-style 'behavior' file from a table of task stimuli and behaviors.
%
% GENERAL
% -------
% Make a NeuroFreq-style 'behavior' file from a table of task conditions, stimuli, and
% behaviors. Table must have column headers and a number of rows equal to
% the number of trials in the input EEG.
%
% Extension:
%   Optional 'lags' vector (e.g., [-2 -1 1 2]) creates lagged versions of EVERY
%   variable in the table, NaN-padded where unavailable.
%   Example: accuracy -> accuracy_1back (lag=-1), accuracy_2forward (lag=+2)
%
% OUTPUT
% ------
% behav - behavior in NF format (struct array, 1 per trial).
%
% INPUT
% -----
% 1) tab  - table of task conditions, stimuli, and behaviors OR a filename to a readable table.
% 2) lags - optional vector of integer lags (e.g., [-2 -1 1 2]). Default: [] (no lag features).
%

if nargin < 1 || isempty(tab)
    error('tab is required.');
end

if nargin < 2 || isempty(lags)
    lags = [];
end

if isempty(lags) == 0
    if ~isnumeric(lags) || ~isvector(lags)
        error('lags must be a numeric vector of integer offsets (e.g., [-2 -1 1 2]).');
    end

    if any(isnan(lags(:))) || any(isinf(lags(:)))
        error('lags must not contain NaN or Inf.');
    end

    if any(abs(lags(:) - round(lags(:))) > 0)
        error('lags must contain integer values only.');
    end

    lags = unique(round(lags(:).'), 'stable');

    if any(lags == 0)
        error('lags must not contain 0 (0-lag would duplicate the original variable).');
    end
end

if istable(tab) == 1
    behavior = tab;
elseif ischar(tab) == 1 || isstring(tab) == 1
    filepath = char(tab);

    if exist(filepath, 'file') == 2
        behavior = readtable(filepath);
    else
        error('File not found: %s', filepath);
    end
else
    error('tab must be a table or a filename.');
end

nTrials = size(behavior, 1);

if isempty(lags) == 0
    vars = behavior.Properties.VariableNames;

    for j = 1:numel(vars)
        vname = vars{j};

        for k = 1:numel(lags)
            lag = lags(k);

            newName = local_make_lagged_name(vname, lag);

            if any(strcmp(behavior.Properties.VariableNames, newName))
                error('Lagged variable name collision: %s already exists in table.', newName);
            end

            behavior.(newName) = local_make_lagged_column(behavior.(vname), lag, nTrials);
        end
    end
end

vars = behavior.Properties.VariableNames;

behav = [];

for h = 1:nTrials
    for j = 1:numel(vars)
        behav(h).(vars{j}) = table2array(behavior(h, j)); %#ok<AGROW>
    end
end

end

% =====================================================================
% Local: lagged column naming
% =====================================================================
function newName = local_make_lagged_name(baseName, lag)

if lag < 0
    steps = abs(lag);
    suffix = sprintf('_%dback', steps);
else
    steps = lag;
    suffix = sprintf('_%dforward', steps);
end

newName = [baseName suffix];

end

% =====================================================================
% Local: build NaN-padded lagged column
% =====================================================================
function out = local_make_lagged_column(x, lag, nTrials)

% Ensure column vector for indexing operations
if isrow(x)
    x = x.';
end

% Handle table columns that are cell arrays, strings, categoricals, etc.
% Requirement: pad with NaN. We'll produce a numeric column when possible,
% otherwise a cell column padded with {NaN}.
if isnumeric(x) == 1 || islogical(x) == 1
    out = nan(nTrials, 1);

    if lag < 0
        k = abs(lag);
        if nTrials > k
            out((k + 1):end) = double(x(1:(end - k)));
        end
    else
        k = lag;
        if nTrials > k
            out(1:(end - k)) = double(x((k + 1):end));
        end
    end

elseif iscell(x) == 1
    out = cell(nTrials, 1);

    for i = 1:nTrials
        out{i, 1} = NaN;
    end

    if lag < 0
        k = abs(lag);
        if nTrials > k
            out((k + 1):end) = x(1:(end - k));
        end
    else
        k = lag;
        if nTrials > k
            out(1:(end - k)) = x((k + 1):end);
        end
    end

elseif isstring(x) == 1

    % Strings cannot be NaN-padded naturally; if user wants lags for string columns,
    % we convert to cellstr-like via cell and pad with NaN (cell).
    outCell = cell(nTrials, 1);

    for i = 1:nTrials
        outCell{i, 1} = NaN;
    end

    if lag < 0
        k = abs(lag);
        if nTrials > k
            outCell((k + 1):end) = cellstr(x(1:(end - k)));
        end
    else
        k = lag;
        if nTrials > k
            outCell(1:(end - k)) = cellstr(x((k + 1):end));
        end
    end

    out = outCell;

elseif iscategorical(x) == 1
    out = cell(nTrials, 1);

    for i = 1:nTrials
        out{i, 1} = NaN;
    end

    if lag < 0
        k = abs(lag);
        if nTrials > k
            out((k + 1):end) = cellstr(x(1:(end - k)));
        end
    else
        k = lag;
        if nTrials > k
            out(1:(end - k)) = cellstr(x((k + 1):end));
        end
    end

else
    % Fallback: store as cell with NaN padding
    out = cell(nTrials, 1);

    for i = 1:nTrials
        out{i, 1} = NaN;
    end

    if lag < 0
        k = abs(lag);
        if nTrials > k
            for i = (k + 1):nTrials
                out{i, 1} = x(i - k);
            end
        end
    else
        k = lag;
        if nTrials > k
            for i = 1:(nTrials - k)
                out{i, 1} = x(i + k);
            end
        end
    end
end

end













