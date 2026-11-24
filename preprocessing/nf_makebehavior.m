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
    error('nf_makebehavior:TabRequired','tab is required.');
end

if nargin < 2 || isempty(lags)
    lags = [];
end

if isempty(lags) == 0
    if ~isnumeric(lags) || ~isvector(lags)
        error('nf_makebehavior:LagsNotNumeric',...
            'lags must be a numeric vector of integer offsets (e.g., [-2 -1 1 2]).');
    end
    if any(isnan(lags(:))) || any(isinf(lags(:)))
        error('nf_makebehavior:LagsContainsNaNOrInf',...
            'lags must not contain NaN or Inf.');
    end
    if any(abs(lags(:) - round(lags(:))) > 0)
        error('nf_makebehavior:LagsContainsNonIntegers',...
            'lags must contain integer values only.');
    end

    lags = unique(round(lags(:).'), 'stable');

    if any(lags == 0)
        error('nf_makebehavior:LagsCannotContain0',...
            'lags must not contain 0 (0-lag would duplicate the original variable).');
    end
end

if isa(tab, 'table') == 1
    behavior = tab;
elseif ischar(tab) == 1 || isstring(tab) == 1
    filepath = char(tab);

    if exist(filepath, 'file') == 2
        [~,~,extension] = fileparts(filepath);
        if exist('OCTAVE_VERSION','builtin')
            behavior = hnf_octave_readtable(filepath);
        elseif strcmpi(extension,'.tsv')
            behavior = readtable(filepath,'FileType','text','Delimiter',sprintf('\t'));
        else
            behavior = readtable(filepath);
        end
    else
        error('File not found: %s', filepath);
    end
else
    error('nf_makebehavior:TabMustBeTableorFile',...
        'tab must be a table or a filename.');
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
                error('nf_makebehavior:LaggedVariableNameCollision',...
                    'Lagged variable name collision: %s already exists in table.', newName);
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

disp('[nf_makebehavior]: Behavioral file complete.');

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

% Shift the trial dimension. A table variable may store several values
% per trial, all of which must travel together. Padding uses numeric NaN
% or cells containing NaN for nonnumeric types.
dims = size(x);
if isnumeric(x) || islogical(x)
    out = nan(dims);
    values = double(x);
elseif ischar(x)
    values = cellstr(x);
    out = repmat({NaN},nTrials,1);
else
    if iscell(x)
        values = x;
    elseif isstring(x) || isa(x, 'categorical')
        values = reshape(cellstr(x(:)),dims);
    else
        values = num2cell(x);
    end
    out = repmat({NaN},dims);
end
target = max(1,1-lag):min(nTrials,nTrials-lag);
subs = repmat({':'},1,ndims(values));
source = subs;
subs{1} = target;
source{1} = target+lag;
out(subs{:}) = values(source{:});
end
