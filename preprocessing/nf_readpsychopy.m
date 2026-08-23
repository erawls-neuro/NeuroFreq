function [trials, info] = nf_readpsychopy(inputData, varargin)
% NF_READPSYCHOPY  Read the trial rows from a PsychoPy-wide data file.
%
% TRIALS = NF_READPSYCHOPY(INPUT) accepts a PsychoPy CSV, TSV, TXT, XLS,
% or XLSX file, or an already imported MATLAB table. It finds the dominant
% PsychoPy trial-loop handler and returns only rows declared by that loop.
%
% [TRIALS, INFO] = NF_READPSYCHOPY(...) also returns the selected handler,
% source-row indices, candidate handlers, recovered run metadata, likely
% setup rows, and index diagnostics.
%
% Name/value inputs:
%   'TrialVariable'      Exact loop-index variable or loop prefix. For
%                        example, 'trials.thisN' or 'trials'. Default:
%                        infer the dominant non-practice loop.
%   'TrialRows'          Explicit raw-table row indices or a logical mask.
%                        This bypasses loop inference. Default: [].
%   'ExcludeRows'        Raw-table rows to omit after detection. Default: [].
%   'MinimumTrials'      Minimum credible loop size. Default: 2.
%   'AmbiguityRatio'     Competing-loop score ratio that triggers an error.
%                        Default: 0.85.
%   'PropagateMetadata'  Copy constant pre-task run parameters found only
%                        before the trial rows into the table. Default: true.
%   'OmitFlaggedRows'    Omit rows flagged as possible setup/seed rows.
%                        Default: false. Inspect INFO.flagged first.
%   'DropEmptyVariables' Remove variables empty across all retained trials.
%                        Default: false.
%   'Verbose'            Print a compact import summary. Default: true.
%
% A response, RT, correctness value, or component stop time is never
% required for trial membership. Declared no-response and partial trials
% are retained. Suspicious rows are reported rather than silently removed.
%
% Example:
%   [trials, info] = nf_readpsychopy('sub-3005_task-its.csv');
%   analysisTrials = nf_readpsychopy('sub-3005_task-its.csv', ...
%       'OmitFlaggedRows', true);

parser = inputParser;
parser.FunctionName = mfilename;
addRequired(parser, 'inputData', @nf_valid_input);
addParameter(parser, 'TrialVariable', '', @nf_valid_text_scalar);
addParameter(parser, 'TrialRows', [], @nf_valid_rows_option);
addParameter(parser, 'ExcludeRows', [], @nf_valid_rows_option);
addParameter(parser, 'MinimumTrials', 2, @nf_valid_positive_integer);
addParameter(parser, 'AmbiguityRatio', 0.85, @nf_valid_ratio);
addParameter(parser, 'PropagateMetadata', true, @nf_valid_logical_scalar);
addParameter(parser, 'OmitFlaggedRows', false, @nf_valid_logical_scalar);
addParameter(parser, 'DropEmptyVariables', false, @nf_valid_logical_scalar);
addParameter(parser, 'Verbose', true, @nf_valid_logical_scalar);
parse(parser, inputData, varargin{:});
options = parser.Results;
if isstring(options.TrialVariable) && ...
        ~ismissing(options.TrialVariable) && ...
        strlength(options.TrialVariable) == 0
    options.TrialVariable = '';
end

[raw, source] = nf_read_input(inputData);
if height(raw) == 0
    error('nf_readpsychopy:EmptyInput', ...
        'The input contains variables but no data rows.');
end
if width(raw) == 0
    error('nf_readpsychopy:EmptyInput', ...
        'The input contains data rows but no variables.');
end

presence = nf_presence_matrix(raw);
if ~isempty(options.TrialRows)
    declaredMask = nf_rows_to_mask(options.TrialRows, height(raw), ...
        'TrialRows');
    detection = nf_explicit_rows_detection(declaredMask);
else
    [declaredMask, detection] = nf_detect_trial_rows(raw, presence, ...
        options);
end

declaredRows = find(declaredMask);
if numel(declaredRows) < options.MinimumTrials
    error('nf_readpsychopy:TooFewTrials', ...
        ['Only %d trial rows were identified, fewer than the requested ' ...
        'MinimumTrials value of %d.'], ...
        numel(declaredRows), options.MinimumTrials);
end

flagged = nf_flag_suspicious_rows(raw, presence, declaredRows);
retainedMask = declaredMask;
explicitExcludeMask = nf_rows_to_mask(options.ExcludeRows, height(raw), ...
    'ExcludeRows');
retainedMask(explicitExcludeMask) = false;
if options.OmitFlaggedRows
    retainedMask(flagged.sourceRows) = false;
end

retainedRows = find(retainedMask);
if isempty(retainedRows)
    error('nf_readpsychopy:NoTrialsRemain', ...
        'No trial rows remain after the requested exclusions.');
end
flagged.retainedSourceRows = flagged.sourceRows( ...
    ismember(flagged.sourceRows, retainedRows));
flagged.excludedSourceRows = flagged.sourceRows( ...
    ~ismember(flagged.sourceRows, retainedRows));

trials = raw(retainedRows, :);
[metadata, metadataNames, metadataSourceRows] = ...
    nf_recover_run_metadata(raw, presence, declaredRows);
if options.PropagateMetadata
    trials = nf_propagate_run_metadata(trials, raw, metadataNames, ...
        metadataSourceRows);
end

trialPresence = nf_presence_matrix(trials);
emptyVariableMask = ~any(trialPresence, 1);
emptyVariables = trials.Properties.VariableNames(emptyVariableMask);
if options.DropEmptyVariables
    trials(:, emptyVariableMask) = [];
end

indexDiagnostics = nf_index_diagnostics(raw, detection.selectedVariable, ...
    declaredRows);
warnings = nf_build_warnings(detection, flagged, indexDiagnostics);

info = struct();
info.schemaVersion = '1.0.0';
info.source = source;
info.method = detection.method;
info.confidence = detection.confidence;
info.selectedVariable = detection.selectedVariable;
info.selectedLoop = detection.selectedLoop;
info.candidates = detection.candidates;
info.inputRows = height(raw);
info.inputVariables = width(raw);
info.declaredTrialRows = declaredRows;
info.sourceRows = retainedRows;
info.excludedSourceRows = find(declaredMask & ~retainedMask);
info.firstSourceRow = retainedRows(1);
info.lastSourceRow = retainedRows(end);
info.sourceSegments = nf_source_segments(retainedRows);
info.outputRows = height(trials);
info.outputVariables = width(trials);
info.flagged = flagged;
info.recoveredMetadata = metadata;
info.metadataVariables = metadataNames;
info.metadataSourceRows = metadataSourceRows;
info.metadataPropagated = logical(options.PropagateMetadata);
info.emptyVariables = emptyVariables;
info.emptyVariablesDropped = logical(options.DropEmptyVariables);
info.index = indexDiagnostics;
info.warnings = warnings;

if options.Verbose
    nf_print_summary(info, options);
end

end

function valid = nf_valid_input(value)
validCharacterInput = ischar(value) && ...
    (isrow(value) || isempty(value));
validStringInput = isstring(value) && isscalar(value) && ...
    ~ismissing(value);
valid = istable(value) || validCharacterInput || validStringInput;
end

function valid = nf_valid_text_scalar(value)
validCharacterInput = ischar(value) && ...
    (isrow(value) || isempty(value));
validStringInput = isstring(value) && isscalar(value) && ...
    ~ismissing(value);
valid = validCharacterInput || validStringInput;
end

function valid = nf_valid_rows_option(value)
valid = isempty(value) || ((isnumeric(value) || islogical(value)) && ...
    isreal(value) && isvector(value));
end

function valid = nf_valid_positive_integer(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value >= 1 && value == round(value);
end

function valid = nf_valid_ratio(value)
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value > 0 && value <= 1;
end

function valid = nf_valid_logical_scalar(value)
if islogical(value) && isscalar(value)
    valid = true;
    return
end
valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && (value == 0 || value == 1);
end

function [raw, source] = nf_read_input(inputData)
if istable(inputData)
    raw = inputData;
    source = '<MATLAB table>';
    return
end

fileName = char(inputData);
if ~isfile(fileName)
    error('nf_readpsychopy:FileNotFound', ...
        'The input file does not exist: %s', fileName);
end

[~, ~, extension] = fileparts(fileName);
extension = lower(extension);
supported = {'.csv', '.tsv', '.txt', '.xls', '.xlsx'};
if ~any(strcmp(extension, supported))
    error('nf_readpsychopy:UnsupportedFile', ...
        ['Unsupported file type "%s". Export PsychoPy data as CSV, TSV, ' ...
        'TXT, XLS, or XLSX, or pass an imported table. PsychoPy PSYDAT ' ...
        'and LOG files are not rectangular trial tables.'], extension);
end

try
    importOptions = detectImportOptions(fileName, ...
        'VariableNamingRule', 'preserve');
catch firstError
    try
        importOptions = detectImportOptions(fileName);
    catch
        rethrow(firstError)
    end
    if isprop(importOptions, 'VariableNamingRule')
        importOptions.VariableNamingRule = 'preserve';
    end
end

try
    raw = readtable(fileName, importOptions);
catch readError
    error('nf_readpsychopy:ReadFailed', ...
        'MATLAB could not import "%s": %s', fileName, readError.message);
end
source = fileName;
end

function presence = nf_presence_matrix(inputTable)
rowCount = height(inputTable);
variableCount = width(inputTable);
presence = false(rowCount, variableCount);
names = inputTable.Properties.VariableNames;
for variableIndex = 1:variableCount
    presence(:, variableIndex) = ...
        nf_value_present(inputTable.(names{variableIndex}));
end
end

function present = nf_value_present(value)
rowCount = size(value, 1);
if isnumeric(value) || islogical(value)
    missing = ismissing(value);
    present = any(~missing, 2);
    return
end
if isstring(value)
    presentValues = ~ismissing(value) & strlength(strtrim(value)) > 0;
    present = any(presentValues, 2);
    return
end
if ischar(value)
    present = any(~isspace(value), 2);
    return
end
if iscell(value)
    present = false(rowCount, 1);
    for rowIndex = 1:rowCount
        rowValue = value(rowIndex, :);
        for columnIndex = 1:numel(rowValue)
            if nf_cell_value_present(rowValue{columnIndex})
                present(rowIndex) = true;
                break
            end
        end
    end
    return
end

try
    missing = ismissing(value);
    present = any(~missing, 2);
catch
    present = true(rowCount, 1);
end
end

function present = nf_cell_value_present(value)
if isempty(value)
    present = false;
    return
end
if ischar(value)
    present = ~isempty(strtrim(value));
    return
end
if isstring(value)
    present = any(~ismissing(value) & strlength(strtrim(value)) > 0);
    return
end
if isnumeric(value) || islogical(value)
    missing = ismissing(value);
    present = any(~missing(:));
    return
end

try
    missing = ismissing(value);
    present = any(~missing(:));
catch
    present = true;
end
end

function mask = nf_rows_to_mask(rows, rowCount, optionName)
mask = false(rowCount, 1);
if isempty(rows)
    return
end
if islogical(rows)
    if numel(rows) ~= rowCount
        error('nf_readpsychopy:InvalidRows', ...
            '%s must contain one logical value per raw-table row.', ...
            optionName);
    end
    mask = rows(:);
    return
end
if any(~isfinite(rows)) || any(rows ~= round(rows)) || ...
        any(rows < 1) || any(rows > rowCount)
    error('nf_readpsychopy:InvalidRows', ...
        '%s must contain finite integer raw-table rows from 1 through %d.', ...
        optionName, rowCount);
end
mask(unique(rows(:))) = true;
end

function detection = nf_explicit_rows_detection(mask)
detection = struct();
detection.method = 'explicit rows';
detection.confidence = 'specified';
detection.selectedVariable = '';
detection.selectedLoop = '';
detection.candidates = struct([]);
detection.fallback = false;
detection.scoreRatio = 0;
detection.rowCount = sum(mask);
end

function [mask, detection] = nf_detect_trial_rows(raw, presence, options)
names = raw.Properties.VariableNames;
normalizedNames = cellfun(@nf_normalize_name, names, ...
    'UniformOutput', false);

if ~isempty(options.TrialVariable)
    variableIndex = nf_resolve_trial_variable(names, normalizedNames, ...
        char(options.TrialVariable));
    mask = presence(:, variableIndex);
    if sum(mask) < options.MinimumTrials
        error('nf_readpsychopy:TooFewTrials', ...
            'TrialVariable "%s" identifies only %d rows.', ...
            names{variableIndex}, sum(mask));
    end
    detection = nf_explicit_variable_detection(raw, presence, ...
        variableIndex);
    return
end

candidates = nf_loop_candidates(raw, presence, normalizedNames, ...
    options.MinimumTrials);
if isempty(candidates)
    [mask, detection] = nf_fallback_detection(raw, presence, ...
        normalizedNames, options);
    return
end

groups = nf_group_loop_candidates(candidates, presence);
eligible = true(1, numel(groups));
hasQualified = any([groups.hasQualified]);
if hasQualified
    eligible = [groups.hasQualified];
end
hasNonPractice = any(eligible & ~[groups.isPractice]);
if hasNonPractice
    eligible = eligible & ~[groups.isPractice];
end
eligibleIndices = find(eligible);
eligibleScores = [groups(eligibleIndices).score];
[~, order] = sort(eligibleScores, 'descend');
orderedIndices = eligibleIndices(order);
bestGroupIndex = orderedIndices(1);

scoreRatio = 0;
if numel(orderedIndices) > 1
    secondGroupIndex = orderedIndices(2);
    scoreRatio = groups(secondGroupIndex).score / ...
        groups(bestGroupIndex).score;
    if scoreRatio >= options.AmbiguityRatio
        nf_throw_loop_ambiguity(groups(bestGroupIndex), ...
            groups(secondGroupIndex), candidates);
    end
end

bestGroup = groups(bestGroupIndex);
mask = bestGroup.mask;
detection = struct();
detection.method = 'PsychoPy loop handler';
detection.confidence = nf_confidence_label(scoreRatio, ...
    numel(orderedIndices));
detection.selectedVariable = ...
    candidates(bestGroup.representative).variable;
detection.selectedLoop = candidates(bestGroup.representative).loop;
detection.candidates = nf_candidate_report(candidates, groups);
detection.fallback = false;
detection.scoreRatio = scoreRatio;
detection.rowCount = sum(mask);
end

function variableIndex = nf_resolve_trial_variable(names, normalizedNames, query)
exactMatch = find(strcmp(names, query));
if isscalar(exactMatch)
    variableIndex = exactMatch;
    return
end

caseInsensitiveMatch = find(strcmpi(names, query));
if isscalar(caseInsensitiveMatch)
    variableIndex = caseInsensitiveMatch;
    return
end

normalizedQuery = nf_normalize_name(query);
normalizedMatch = find(strcmp(normalizedNames, normalizedQuery));
if isscalar(normalizedMatch)
    variableIndex = normalizedMatch;
    return
end
if numel(normalizedMatch) > 1
    error('nf_readpsychopy:AmbiguousTrialVariable', ...
        ['TrialVariable "%s" matches multiple variables after header ' ...
        'normalization. Use one exact preserved variable name.'], query);
end

loopVariable = [normalizedQuery '.thisn'];
prefixMatch = find(strcmp(normalizedNames, loopVariable));
if isscalar(prefixMatch)
    variableIndex = prefixMatch;
    return
end
if numel(prefixMatch) > 1
    error('nf_readpsychopy:AmbiguousTrialVariable', ...
        'TrialVariable loop prefix "%s" matches multiple variables.', query);
end

if isempty(exactMatch) && isempty(caseInsensitiveMatch) && ...
        isempty(prefixMatch)
    error('nf_readpsychopy:UnknownTrialVariable', ...
        ['TrialVariable "%s" does not match a variable or a loop prefix. ' ...
        'Use a preserved table variable name such as "trials.thisN".'], ...
        query);
end
error('nf_readpsychopy:AmbiguousTrialVariable', ...
    'TrialVariable "%s" does not resolve to exactly one variable.', query);
end

function detection = nf_explicit_variable_detection(raw, presence, ...
        variableIndex)
names = raw.Properties.VariableNames;
normalizedName = nf_normalize_name(names{variableIndex});
loopName = nf_loop_prefix(normalizedName);
detection = struct();
detection.method = 'specified PsychoPy loop handler';
detection.confidence = 'specified';
detection.selectedVariable = names{variableIndex};
detection.selectedLoop = loopName;
detection.candidates = struct([]);
detection.fallback = false;
detection.scoreRatio = 0;
detection.rowCount = sum(presence(:, variableIndex));
end

function candidates = nf_loop_candidates(raw, presence, normalizedNames, ...
        minimumTrials)
names = raw.Properties.VariableNames;
candidates = struct('variableIndex', {}, 'variable', {}, 'loop', {}, ...
    'mask', {}, 'rowCount', {}, 'numericFraction', {}, ...
    'integerFraction', {}, 'stepFraction', {}, 'monotoneFraction', {}, ...
    'uniqueFraction', {}, 'companions', {}, 'contentSupport', {}, ...
    'isPractice', {}, 'quality', {});

for variableIndex = 1:numel(names)
    normalizedName = normalizedNames{variableIndex};
    if isempty(regexp(normalizedName, '(^|[.])thisn$', 'once'))
        continue
    end

    mask = presence(:, variableIndex);
    rowCount = sum(mask);
    if rowCount < minimumTrials
        continue
    end

    loopName = nf_loop_prefix(normalizedName);
    values = nf_numeric_vector(raw.(names{variableIndex}));
    selectedValues = values(mask);
    finiteValues = isfinite(selectedValues);
    numericFraction = mean(finiteValues);
    if any(finiteValues)
        integerFraction = mean(abs(selectedValues(finiteValues) - ...
            round(selectedValues(finiteValues))) < 1e-10);
        uniqueFraction = numel(unique(selectedValues(finiteValues))) / ...
            sum(finiteValues);
    else
        integerFraction = 0;
        uniqueFraction = 0;
    end
    if numericFraction < 0.80 || integerFraction < 0.80
        continue
    end

    finiteSequence = selectedValues;
    finiteSequence = finiteSequence(isfinite(finiteSequence));
    if numel(finiteSequence) > 1
        differences = diff(finiteSequence);
        stepFraction = mean(abs(differences - 1) < 1e-10);
        monotoneFraction = mean(differences >= 0);
    else
        stepFraction = 0;
        monotoneFraction = 0;
    end

    companions = nf_count_companions(loopName, mask, presence, ...
        normalizedNames);
    contentSupport = nf_count_matching_masks(mask, presence);
    isPractice = nf_is_practice_loop(loopName);
    namedBonus = ~isempty(loopName);
    quality = 5 * numericFraction + 2 * integerFraction + ...
        2 * stepFraction + monotoneFraction + uniqueFraction + ...
        0.5 * companions + 0.5 * namedBonus - 2 * isPractice;

    candidate = struct();
    candidate.variableIndex = variableIndex;
    candidate.variable = names{variableIndex};
    candidate.loop = loopName;
    candidate.mask = mask;
    candidate.rowCount = rowCount;
    candidate.numericFraction = numericFraction;
    candidate.integerFraction = integerFraction;
    candidate.stepFraction = stepFraction;
    candidate.monotoneFraction = monotoneFraction;
    candidate.uniqueFraction = uniqueFraction;
    candidate.companions = companions;
    candidate.contentSupport = contentSupport;
    candidate.isPractice = isPractice;
    candidate.quality = quality;
    candidates(end + 1) = candidate; %#ok<AGROW>
end
end

function loopName = nf_loop_prefix(normalizedName)
loopName = regexprep(normalizedName, '(^|[.])thisn$', '');
if endsWith(loopName, '.')
    loopName = loopName(1:end - 1);
end
end

function companionCount = nf_count_companions(loopName, mask, presence, ...
        normalizedNames)
suffixes = {'thistrialn', 'thisrepn', 'thisindex'};
companionCount = 0;
for suffixIndex = 1:numel(suffixes)
    if isempty(loopName)
        companionName = suffixes{suffixIndex};
    else
        companionName = [loopName '.' suffixes{suffixIndex}];
    end
    variableIndex = find(strcmp(normalizedNames, companionName), 1);
    if ~isempty(variableIndex) && isequal(presence(:, variableIndex), mask)
        companionCount = companionCount + 1;
    end
end
end

function count = nf_count_matching_masks(mask, presence)
count = 0;
for variableIndex = 1:size(presence, 2)
    if isequal(presence(:, variableIndex), mask)
        count = count + 1;
    end
end
end

function practice = nf_is_practice_loop(loopName)
pattern = ['(^|[.])(practice|practise|prac|training|train|' ...
    'tutorial|demo|example)(trials?|blocks?)?([.]|$)'];
practice = ~isempty(regexp(loopName, pattern, 'once'));
end

function groups = nf_group_loop_candidates(candidates, presence)
groups = struct('mask', {}, 'candidateIndices', {}, 'rowCount', {}, ...
    'representative', {}, 'isPractice', {}, 'hasQualified', {}, ...
    'score', {});
for candidateIndex = 1:numel(candidates)
    groupIndex = 0;
    for existingIndex = 1:numel(groups)
        if isequal(groups(existingIndex).mask, ...
                candidates(candidateIndex).mask)
            groupIndex = existingIndex;
            break
        end
    end
    if groupIndex == 0
        group = struct();
        group.mask = candidates(candidateIndex).mask;
        group.candidateIndices = candidateIndex;
        group.rowCount = candidates(candidateIndex).rowCount;
        group.representative = candidateIndex;
        group.isPractice = candidates(candidateIndex).isPractice;
        group.hasQualified = ~isempty(candidates(candidateIndex).loop);
        group.score = 0;
        groups(end + 1) = group; %#ok<AGROW>
    else
        groups(groupIndex).candidateIndices(end + 1) = candidateIndex;
        groups(groupIndex).hasQualified = ...
            groups(groupIndex).hasQualified || ...
            ~isempty(candidates(candidateIndex).loop);
        currentRepresentative = groups(groupIndex).representative;
        if candidates(candidateIndex).quality > ...
                candidates(currentRepresentative).quality
            groups(groupIndex).representative = candidateIndex;
        end
    end
end

for groupIndex = 1:numel(groups)
    candidateIndices = groups(groupIndex).candidateIndices;
    namedMask = false(1, numel(candidateIndices));
    for index = 1:numel(candidateIndices)
        namedMask(index) = ...
            ~isempty(candidates(candidateIndices(index)).loop);
    end
    namedIndices = candidateIndices(namedMask);
    if isempty(namedIndices)
        groups(groupIndex).isPractice = false;
    else
        groups(groupIndex).isPractice = ...
            all([candidates(namedIndices).isPractice]);
    end
    representative = groups(groupIndex).representative;
    support = numel(groups(groupIndex).candidateIndices);
    contentSupport = nf_count_matching_masks(groups(groupIndex).mask, ...
        presence);
    groups(groupIndex).score = groups(groupIndex).rowCount + ...
        5 * candidates(representative).quality + ...
        2 * support + min(contentSupport, 25);
end
end

function nf_throw_loop_ambiguity(firstGroup, secondGroup, candidates)
firstNames = nf_group_variable_names(firstGroup, candidates);
secondNames = nf_group_variable_names(secondGroup, candidates);
error('nf_readpsychopy:AmbiguousTrialLoop', ...
    ['Two different PsychoPy loop masks are similarly plausible: [%s] ' ...
    '(%d rows) and [%s] (%d rows). Specify TrialVariable so NeuroFreq ' ...
    'does not silently choose between task, practice, or nested loops.'], ...
    firstNames, firstGroup.rowCount, secondNames, secondGroup.rowCount);
end

function namesText = nf_group_variable_names(group, candidates)
names = cell(1, numel(group.candidateIndices));
for index = 1:numel(group.candidateIndices)
    names{index} = ...
        candidates(group.candidateIndices(index)).variable;
end
namesText = strjoin(names, ', ');
end

function report = nf_candidate_report(candidates, groups)
report = struct('variable', {}, 'loop', {}, 'rowCount', {}, ...
    'quality', {}, 'groupScore', {}, 'practice', {});
for candidateIndex = 1:numel(candidates)
    groupScore = NaN;
    for groupIndex = 1:numel(groups)
        if any(groups(groupIndex).candidateIndices == candidateIndex)
            groupScore = groups(groupIndex).score;
            break
        end
    end
    item = struct();
    item.variable = candidates(candidateIndex).variable;
    item.loop = candidates(candidateIndex).loop;
    item.rowCount = candidates(candidateIndex).rowCount;
    item.quality = candidates(candidateIndex).quality;
    item.groupScore = groupScore;
    item.practice = candidates(candidateIndex).isPractice;
    report(end + 1) = item; %#ok<AGROW>
end
end

function [mask, detection] = nf_fallback_detection(raw, presence, ...
        normalizedNames, options)
names = raw.Properties.VariableNames;
rawRowCount = height(raw);
groups = struct('mask', {}, 'names', {}, 'rowCount', {}, ...
    'support', {}, 'substantiveSupport', {}, 'score', {});

for variableIndex = 1:width(raw)
    mask = presence(:, variableIndex);
    count = sum(mask);
    if count < options.MinimumTrials
        continue
    end

    groupIndex = 0;
    for existingIndex = 1:numel(groups)
        if isequal(groups(existingIndex).mask, mask)
            groupIndex = existingIndex;
            break
        end
    end
    if groupIndex == 0
        group = struct();
        group.mask = mask;
        group.names = {names{variableIndex}}; %#ok
        group.rowCount = count;
        group.support = 1;
        group.substantiveSupport = ...
            nf_is_substantive_name(normalizedNames{variableIndex});
        group.score = 0;
        groups(end + 1) = group; %#ok<AGROW>
    else
        groups(groupIndex).names{end + 1} = names{variableIndex};
        groups(groupIndex).support = groups(groupIndex).support + 1;
        groups(groupIndex).substantiveSupport = ...
            groups(groupIndex).substantiveSupport + ...
            nf_is_substantive_name(normalizedNames{variableIndex});
    end
end

credible = false(1, numel(groups));
for groupIndex = 1:numel(groups)
    groups(groupIndex).score = groups(groupIndex).rowCount * ...
        (1 + log1p(groups(groupIndex).substantiveSupport)) + ...
        groups(groupIndex).support;
    credible(groupIndex) = groups(groupIndex).substantiveSupport >= 2;
end
partial = [groups.rowCount] < rawRowCount;
if any(credible & partial)
    credible = credible & partial;
end

credibleIndices = find(credible);
if isempty(credibleIndices)
    error('nf_readpsychopy:NoTrialLoop', ...
        ['No credible PsychoPy trial-loop index or dominant trial-row ' ...
        'pattern was found. Specify TrialVariable or TrialRows.']);
end

scores = [groups(credibleIndices).score];
[~, order] = sort(scores, 'descend');
orderedIndices = credibleIndices(order);
bestIndex = orderedIndices(1);
scoreRatio = 0;
if numel(orderedIndices) > 1
    secondIndex = orderedIndices(2);
    scoreRatio = groups(secondIndex).score / groups(bestIndex).score;
    if scoreRatio >= options.AmbiguityRatio
        error('nf_readpsychopy:AmbiguousRowPattern', ...
            ['No loop-index variable was found, and two row patterns are ' ...
            'similarly plausible (%d and %d rows). Specify TrialRows.'], ...
            groups(bestIndex).rowCount, groups(secondIndex).rowCount);
    end
end

mask = groups(bestIndex).mask;
detection = struct();
detection.method = 'dominant trial-row pattern';
detection.confidence = nf_confidence_label(scoreRatio, ...
    numel(orderedIndices));
detection.selectedVariable = '';
detection.selectedLoop = '';
detection.candidates = nf_fallback_candidate_report(groups);
detection.fallback = true;
detection.scoreRatio = scoreRatio;
detection.rowCount = sum(mask);
end

function substantive = nf_is_substantive_name(normalizedName)
nonSubstantive = {'participant', 'session', 'date', 'expname', ...
    'expversion', 'psychopyversion', 'framerate', 'expstart', 'notes'};
if any(strcmp(normalizedName, nonSubstantive))
    substantive = false;
    return
end
excludedSuffix = ...
    '([.]started|[.]stopped|[.]duration|[.]rt|[.]keys|[.]corr)$';
substantive = isempty(regexp(normalizedName, excludedSuffix, 'once'));
end

function report = nf_fallback_candidate_report(groups)
report = struct('variables', {}, 'rowCount', {}, 'support', {}, ...
    'substantiveSupport', {}, 'score', {});
for groupIndex = 1:numel(groups)
    item = struct();
    item.variables = groups(groupIndex).names;
    item.rowCount = groups(groupIndex).rowCount;
    item.support = groups(groupIndex).support;
    item.substantiveSupport = groups(groupIndex).substantiveSupport;
    item.score = groups(groupIndex).score;
    report(end + 1) = item; %#ok<AGROW>
end
end

function label = nf_confidence_label(scoreRatio, candidateCount)
if candidateCount <= 1 || scoreRatio < 0.50
    label = 'high';
elseif scoreRatio < 0.70
    label = 'moderate';
else
    label = 'guarded';
end
end

function numeric = nf_numeric_vector(value)
rowCount = size(value, 1);
numeric = NaN(rowCount, 1);
if isnumeric(value) || islogical(value)
    if size(value, 2) == 1
        numeric = double(value);
    end
    return
end
if isstring(value)
    if size(value, 2) == 1
        numeric = str2double(value);
    end
    return
end
if ischar(value)
    numeric = str2double(cellstr(value));
    return
end
if iscategorical(value)
    if size(value, 2) == 1
        numeric = str2double(string(value));
    end
    return
end
if iscell(value) && size(value, 2) == 1
    for rowIndex = 1:rowCount
        item = value{rowIndex};
        if isnumeric(item) && isscalar(item)
            numeric(rowIndex) = double(item);
        elseif islogical(item) && isscalar(item)
            numeric(rowIndex) = double(item);
        elseif ischar(item) || (isstring(item) && isscalar(item))
            numeric(rowIndex) = str2double(string(item));
        end
    end
end
end

function normalized = nf_normalize_name(name)
normalized = lower(char(name));
normalized = regexprep(normalized, '[^a-z0-9]+', '.');
normalized = regexprep(normalized, '^[.]|[.]$', '');
end

function flagged = nf_flag_suspicious_rows(raw, presence, sourceRows)
localRowCount = numel(sourceRows);
variableCount = width(raw);
evidence = false(localRowCount, variableCount);
names = raw.Properties.VariableNames;

for variableIndex = 1:variableCount
    normalizedName = nf_normalize_name(names{variableIndex});
    if ~nf_flag_eligible_name(normalizedName)
        continue
    end
    localPresence = presence(sourceRows, variableIndex);
    presenceFraction = mean(localPresence);
    rareMissing = false(localRowCount, 1);
    if presenceFraction >= 0.90 && presenceFraction < 1
        rareMissing = ~localPresence;
    end
    sentinel = nf_sentinel_rows(raw.(names{variableIndex}), sourceRows);
    evidence(:, variableIndex) = rareMissing | sentinel;
end

evidenceCount = sum(evidence, 2);
localFlagMask = evidenceCount >= 2;
flaggedLocalRows = find(localFlagMask);
flagged = struct();
flagged.sourceRows = sourceRows(localFlagMask);
flagged.declaredTrialRows = flaggedLocalRows;
flagged.evidenceCount = evidenceCount(localFlagMask);
flagged.variables = cell(numel(flaggedLocalRows), 1);
for flaggedIndex = 1:numel(flaggedLocalRows)
    variableIndices = evidence(flaggedLocalRows(flaggedIndex), :);
    flagged.variables{flaggedIndex} = names(variableIndices);
end
flagged.omittedByDefault = false;
end

function eligible = nf_flag_eligible_name(normalizedName)
excludedNames = {'participant', 'session', 'date', 'expname', ...
    'expversion', 'psychopyversion', 'framerate', 'expstart', 'notes', ...
    'thisn', 'thistrialn', 'thisrepn', 'thisindex', 'thisrow.t'};
if any(strcmp(normalizedName, excludedNames))
    eligible = false;
    return
end
if ~isempty(regexp(normalizedName, ...
        '(^|[.])this(n|trialn|repn|index)$', 'once'))
    eligible = false;
    return
end
excludedSuffix = ...
    '([.]started|[.]stopped|[.]duration|[.]rt|[.]keys|[.]corr)$';
eligible = isempty(regexp(normalizedName, excludedSuffix, 'once'));
end

function sentinel = nf_sentinel_rows(value, sourceRows)
sentinel = false(numel(sourceRows), 1);
sentinels = {'na', 'n/a', 'not applicable', 'dummy', 'practice', ...
    'training', 'demo', 'seed', 'setup'};
for localIndex = 1:numel(sourceRows)
    rowIndex = sourceRows(localIndex);
    rowTexts = nf_row_text_values(value, rowIndex);
    for textIndex = 1:numel(rowTexts)
        normalizedText = strtrim(rowTexts{textIndex});
        if any(strcmpi(normalizedText, sentinels))
            sentinel(localIndex) = true;
            break
        end
    end
end
end

function texts = nf_row_text_values(value, rowIndex)
texts = {};
if isstring(value)
    rowValue = value(rowIndex, :);
    for itemIndex = 1:numel(rowValue)
        if ~ismissing(rowValue(itemIndex))
            texts{end + 1} = char(rowValue(itemIndex)); %#ok<AGROW>
        end
    end
    return
end
if ischar(value)
    texts = {value(rowIndex, :)};
    return
end
if iscategorical(value)
    rowValue = value(rowIndex, :);
    for itemIndex = 1:numel(rowValue)
        if ~ismissing(rowValue(itemIndex))
            texts{end + 1} = char(string(rowValue(itemIndex))); %#ok<AGROW>
        end
    end
    return
end
if iscell(value)
    rowValue = value(rowIndex, :);
    for itemIndex = 1:numel(rowValue)
        item = rowValue{itemIndex};
        if ischar(item)
            texts{end + 1} = item; %#ok<AGROW>
        elseif isstring(item) && isscalar(item) && ~ismissing(item) 
            texts{end + 1} = char(item); %#ok<AGROW>
        end
    end
end
end

function [metadata, metadataNames, sourceRows] = ...
        nf_recover_run_metadata(raw, presence, trialRows)
metadata = table();
metadataNames = {};
sourceRows = [];
names = raw.Properties.VariableNames;

for variableIndex = 1:width(raw)
    normalizedName = nf_normalize_name(names{variableIndex});
    if ~nf_metadata_eligible_name(normalizedName)
        continue
    end
    if any(presence(trialRows, variableIndex))
        continue
    end
    presentRows = find(presence(:, variableIndex));
    if isempty(presentRows)
        continue
    end
    if any(presentRows >= trialRows(1))
        continue
    end
    value = raw.(names{variableIndex});
    if size(value, 2) ~= 1 || ...
            ~nf_rows_have_one_value(value, presentRows)
        continue
    end

    name = names{variableIndex};
    metadata.(name) = raw.(name)(presentRows(1), :);
    metadataNames{end + 1} = names{variableIndex}; %#ok<AGROW>
    sourceRows(end + 1) = presentRows(1); %#ok<AGROW>
end
end

function eligible = nf_metadata_eligible_name(normalizedName)
excludedSuffix = ...
    '([.]started|[.]stopped|[.]duration|[.]rt|[.]keys|[.]corr)$';
if ~isempty(regexp(normalizedName, excludedSuffix, 'once'))
    eligible = false;
    return
end
if ~isempty(regexp(normalizedName, ...
        '(^|[.])this(n|trialn|repn|index)$', 'once')) || ...
        strcmp(normalizedName, 'thisrow.t')
    eligible = false;
    return
end
eligible = true;
end

function oneValue = nf_rows_have_one_value(value, rows)
reference = value(rows(1), :);
oneValue = true;
for rowIndex = 2:numel(rows)
    comparison = value(rows(rowIndex), :);
    if ~isequaln(reference, comparison)
        oneValue = false;
        return
    end
end
end

function trials = nf_propagate_run_metadata(trials, raw, names, sourceRows)
for metadataIndex = 1:numel(names)
    name = names{metadataIndex};
    value = raw.(name);
    trials.(name) = repmat(value(sourceRows(metadataIndex), :), ...
        height(trials), 1);
end
end

function diagnostics = nf_index_diagnostics(raw, selectedVariable, ...
        sourceRows)
diagnostics = struct();
diagnostics.numeric = false;
diagnostics.startValue = NaN;
diagnostics.endValue = NaN;
diagnostics.resetSourceRows = [];
diagnostics.duplicateSourceRows = [];
diagnostics.gapSourceRows = [];
if isempty(selectedVariable)
    return
end

values = nf_numeric_vector(raw.(selectedVariable));
selectedValues = values(sourceRows);
if any(~isfinite(selectedValues))
    return
end
diagnostics.numeric = true;
diagnostics.startValue = selectedValues(1);
diagnostics.endValue = selectedValues(end);
if numel(selectedValues) < 2
    return
end

differences = diff(selectedValues);
diagnostics.resetSourceRows = sourceRows(find(differences < 0) + 1);
diagnostics.duplicateSourceRows = sourceRows(find(differences == 0) + 1);
diagnostics.gapSourceRows = sourceRows(find(differences > 1) + 1);
end

function segments = nf_source_segments(sourceRows)
if isempty(sourceRows)
    segments = zeros(0, 2);
    return
end
sourceRows = sourceRows(:);
segmentStarts = [1; find(diff(sourceRows) > 1) + 1];
segmentEnds = [segmentStarts(2:end) - 1; numel(sourceRows)];
segments = [sourceRows(segmentStarts) sourceRows(segmentEnds)];
end

function warnings = nf_build_warnings(detection, flagged, indexDiagnostics)
warnings = {};
if detection.fallback
    warnings{end + 1} = ...
        ['No PsychoPy thisN loop variable was available; rows were ' ...
        'selected from the dominant shared nonmissing-data pattern.'];
end
if ~isempty(flagged.sourceRows)
    if ~isempty(flagged.excludedSourceRows)
        warnings{end + 1} = sprintf( ...
            ['%d possible setup/seed rows were omitted because ' ...
            'of the requested row exclusions.'], ...
            numel(flagged.excludedSourceRows));
    end
    if ~isempty(flagged.retainedSourceRows)
        warnings{end + 1} = sprintf( ...
            ['%d selected rows have at least two setup/seed indicators. ' ...
            'They were retained; inspect INFO.flagged before exclusion.'], ...
            numel(flagged.retainedSourceRows));
    end
end
if ~isempty(indexDiagnostics.resetSourceRows)
    warnings{end + 1} = sprintf( ...
        ['The selected trial index resets at %d source rows. All declared ' ...
        'rows were retained as possible appended runs.'], ...
        numel(indexDiagnostics.resetSourceRows));
end
if ~isempty(indexDiagnostics.duplicateSourceRows)
    warnings{end + 1} = sprintf( ...
        ['The selected trial index repeats at %d source rows. Confirm that ' ...
        'the handler emits one row per intended trial.'], ...
        numel(indexDiagnostics.duplicateSourceRows));
end
if ~isempty(indexDiagnostics.gapSourceRows)
    warnings{end + 1} = sprintf( ...
        ['The selected trial index has %d forward gaps. This can indicate ' ...
        'truncation or deliberately omitted handler entries.'], ...
        numel(indexDiagnostics.gapSourceRows));
end
end

function nf_print_summary(info, options)
fprintf(['[nf_readpsychopy]: retained %d of %d raw rows using ' ...
    '%s'], info.outputRows, info.inputRows, info.method);
if ~isempty(info.selectedVariable)
    fprintf(' (%s)', info.selectedVariable);
end
fprintf('.\n');
if ~isempty(info.flagged.retainedSourceRows)
    fprintf(['[nf_readpsychopy]: %d possible setup/seed rows were flagged and ' ...
        'retained; inspect INFO.flagged.\n'], ...
        numel(info.flagged.retainedSourceRows));
end
if ~isempty(info.flagged.excludedSourceRows)
    fprintf(['[nf_readpsychopy]: %d possible setup/seed rows were flagged and ' ...
        'omitted by the requested exclusions.\n'], ...
        numel(info.flagged.excludedSourceRows));
end
if ~isempty(info.metadataVariables) && options.PropagateMetadata
    fprintf('[nf_readpsychopy]: Propagated run metadata: %s.\n', ...
        strjoin(info.metadataVariables, ', '));
end
end
