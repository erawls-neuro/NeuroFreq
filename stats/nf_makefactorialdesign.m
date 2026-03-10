function [X, designTab, info] = nf_makefactorialdesign(TF, varargin)
% NF_MAKEFACTORIALDESIGN  Build a factorial design matrix/table from TF.behavior.
%
% SUMMARY
% -------
% Converts a NeuroFreq-style behavior struct (often TF.behavior) into a
% fixed-effect design matrix for GLM and a design table compatible with LME.
%
% RULES
% -----
% 1) Base predictors:
%    - numeric scalar fields become predictors
%    - z-score all non-binary numeric predictors
%    - binary predictors are NOT z-scored
% 2) Interactions:
%    - builds all multiplicative interactions up to maxOrder
%    - subject ID never participates in interactions
% 3) Subject ID:
%    - optionally detected or specified via subjectField
%    - subject column can be included in designTab
%    - subject is excluded from X by default (for LME random intercept use)
%
% USAGE
% -----
%   [X, designTab, info] = nf_makefactorialdesign(TF)
%   [X, designTab, info] = nf_makefactorialdesign(TF, 'maxOrder', 2)
%   [X, designTab, info] = nf_makefactorialdesign(TF, 'subjectField', 'subID', 'requireSubjectID', true)
%
% INPUT
% -----
%   TF  - TF struct containing TF.behavior, OR a behavior struct array, OR a table.
%
% NAME-VALUE PAIRS
% ----------------
%   'standardize'               (default: true)
%       Z-score all numerical regressors?
%
%   'includeIntercept'          (default: false)
%       Adds a leading intercept column of ones to X and designTab.
%
%   'includeBinaryInteractions' (default: true)
%       If true, binary predictors may participate in interactions.
%       If false, only z-scored (continuous) predictors participate.
%
%   'maxOrder'                  (default: [])
%       Maximum interaction order. [] means "full" (up to eligible predictors).
%
%   'maxTerms'                  (default: 5000)
%       Hard cap on number of columns in X (including interactions & intercept).
%
%   'binaryTol'                 (default: 1e-12)
%       Numeric tolerance for detecting binary values.
%
%   'subjectField'              (default: '')
%       Field name storing subject ID, e.g. 'subID'.
%       If empty, function will attempt auto-detection from common names.
%
%   'requireSubjectID'          (default: false)
%       If true, error if no subject field is found. Recommended for LME.
%
%   'subjectAsCategorical'      (default: true)
%       If true, stores subject ID in designTab as categorical (recommended for LME).
%
%   'keepSubjectInDesignTab'    (default: true)
%       If true, subject column is included in designTab (even though excluded from X).
%
%   'keepSubjectInX'            (default: false)
%       If true, subject column is included in X (not recommended; usually wrong).
%
%   'excludeFields'             (default: {})
%       Cellstr of field names to exclude from predictor building entirely.
%
% OUTPUT
% ------
%   X         - numeric matrix (nObs x nFixedTerms) for GLM / TFCE-GLM
%   designTab - table including subject ID (optional) + fixed terms, for LME
%   info      - struct metadata:
%               .nObs
%               .subjectField
%               .hasSubject
%               .subID
%               .fixedNames
%               .baseNames
%               .interactionNames
%               .interactionOrders
%               .isBinary
%               .mu
%               .sigma
%               .eligibleInteractionIdx
%

parser = inputParser();
parser.FunctionName = mfilename();

addParameter(parser, 'standardize', false);
addParameter(parser, 'includeIntercept', false);
addParameter(parser, 'includeBinaryInteractions', true);
addParameter(parser, 'maxOrder', []);
addParameter(parser, 'maxTerms', 5000);
addParameter(parser, 'binaryTol', 1e-12);

addParameter(parser, 'subjectField', '');
addParameter(parser, 'requireSubjectID', false);
addParameter(parser, 'subjectAsCategorical', true);
addParameter(parser, 'keepSubjectInDesignTab', true);
addParameter(parser, 'keepSubjectInX', false);

addParameter(parser, 'excludeFields', {});

parse(parser, varargin{:});

standardize = parser.Results.standardize;
includeIntercept = parser.Results.includeIntercept;
includeBinaryInteractions = parser.Results.includeBinaryInteractions;
maxOrder = parser.Results.maxOrder;
maxTerms = parser.Results.maxTerms;
binaryTol = parser.Results.binaryTol;

subjectField = parser.Results.subjectField;
requireSubjectID = parser.Results.requireSubjectID;
subjectAsCategorical = parser.Results.subjectAsCategorical;
keepSubjectInDesignTab = parser.Results.keepSubjectInDesignTab;
keepSubjectInX = parser.Results.keepSubjectInX;

excludeFields = parser.Results.excludeFields;

if ~iscell(excludeFields)
    error('nf_makefactorialdesign:BadExcludeFields', 'excludeFields must be a cell array of strings.');
end

behaviorTable = local_get_behavior_table(TF);

nObs = size(behaviorTable, 1);
if nObs < 1
    error('nf_makefactorialdesign:EmptyBehavior', 'Behavior has zero rows.');
end

% -----------------------------
% Detect / extract subject ID
% -----------------------------
[hasSubject, subjectFieldResolved] = local_resolve_subject_field(behaviorTable, subjectField);

if requireSubjectID == true && hasSubject == false
    error('nf_makefactorialdesign:MissingSubjectID', ...
        'No subject ID field was found. Provide ''subjectField'' (e.g., ''subID'') or add it to TF.behavior.');
end

subID = [];
if hasSubject == true
    subID = local_extract_subject_vector(behaviorTable, subjectFieldResolved);
end

% -----------------------------
% Exclude fields from modeling
% -----------------------------
predictorTable = behaviorTable;

excluded = {};
if hasSubject == true
    predictorTable(:, subjectFieldResolved) = [];
    excluded = [excluded, {subjectFieldResolved}];
end

if ~isempty(excludeFields)
    for iE = 1:numel(excludeFields)
        nm = excludeFields{iE};
        if isstring(nm)
            nm = char(nm);
        end
        if ~ischar(nm)
            error('nf_makefactorialdesign:BadExcludeName', 'excludeFields elements must be strings.');
        end
        if any(strcmp(predictorTable.Properties.VariableNames, nm))
            predictorTable(:, nm) = [];
            excluded = [excluded, {nm}]; %#ok<AGROW>
        end
    end
end

% -----------------------------
% Sanitize predictor names
% -----------------------------
rawNames = predictorTable.Properties.VariableNames;
nVars = numel(rawNames);

if nVars < 1
    error('nf_makefactorialdesign:NoPredictors', ...
        'After excluding fields, there are no predictors left to model.');
end

baseNames = cell(1, nVars);
for iVar = 1:nVars
    baseNames{iVar} = local_sanitize_varname(rawNames{iVar});
end

% -----------------------------
% Build base predictors
% -----------------------------
[baseX, isBinary, mu, sigma] = local_build_base_predictors(predictorTable, baseNames, binaryTol, standardize);

% -----------------------------
% Determine interaction-eligible predictors
% -----------------------------
eligibleIdx = local_get_interaction_eligible_indices(isBinary, includeBinaryInteractions);
nEligible = numel(eligibleIdx);

if isempty(maxOrder)
    maxOrder = nEligible;
end

if ~isscalar(maxOrder) || ~isnumeric(maxOrder) || maxOrder < 1
    error('nf_makefactorialdesign:BadMaxOrder', 'maxOrder must be a positive scalar integer or [].');
end

maxOrder = floor(maxOrder);

if maxOrder > nEligible
    maxOrder = nEligible;
end

% -----------------------------
% Build interactions
% -----------------------------
interactionX = [];
interactionNames = {};
interactionOrders = [];

if maxOrder >= 2
    [interactionX, interactionNames, interactionOrders] = local_build_interactions( ...
        baseX(:, eligibleIdx), ...
        baseNames(eligibleIdx), ...
        maxOrder ...
    );
end

% -----------------------------
% Assemble X (fixed effects only)
% -----------------------------
X = baseX;
fixedNames = baseNames;

if ~isempty(interactionX)
    X = [X, interactionX];
    fixedNames = [fixedNames, interactionNames];
end

if includeIntercept == true
    X = [ones(nObs, 1), X];
    fixedNames = ['Intercept', fixedNames];
end

% Optionally include subject in X (not recommended)
if keepSubjectInX == true
    if hasSubject == false
        error('nf_makefactorialdesign:NoSubjectForX', ...
            'keepSubjectInX=true, but no subject field exists.');
    end
    subjNum = double(grp2idx(categorical(subID)));
    X = [subjNum(:), X];
    fixedNames = ['subID', fixedNames];
end

if size(X, 2) > maxTerms
    error('nf_makefactorialdesign:TooManyTerms', ...
        ['Design has %d columns (maxTerms=%d). ' ...
         'Reduce predictors, set ''maxOrder'' smaller, or raise ''maxTerms''.'], ...
        size(X, 2), maxTerms);
end

% -----------------------------
% Assemble designTab for LME
% -----------------------------
designTab = array2table(X, 'VariableNames', fixedNames);

if keepSubjectInDesignTab == true && hasSubject == true
    if subjectAsCategorical == true
        designTab = addvars(designTab, categorical(subID), 'Before', 1, 'NewVariableNames', subjectFieldResolved);
    else
        designTab = addvars(designTab, subID, 'Before', 1, 'NewVariableNames', subjectFieldResolved);
    end
end

% -----------------------------
% Info output
% -----------------------------
info = struct();
info.nObs = nObs;
info.subjectField = subjectFieldResolved;
info.hasSubject = hasSubject;
info.subID = subID;

info.excludedFields = excluded;

info.fixedNames = fixedNames;
info.baseNames = baseNames;

info.isBinary = isBinary;
info.mu = mu;
info.sigma = sigma;

info.interactionNames = interactionNames;
info.interactionOrders = interactionOrders;
info.eligibleInteractionIdx = eligibleIdx;

if hasSubject == true
    info.suggestedRandomIntercept = ['(1|' subjectFieldResolved ')'];
else
    info.suggestedRandomIntercept = '';
end

end

% ====================================================================== %
% Helpers
% ====================================================================== %

function behaviorTable = local_get_behavior_table(TF)

if istable(TF)
    behaviorTable = TF;
    return
end

if isstruct(TF)
    if isfield(TF, 'behavior')
        behavior = TF.behavior;
        if ~isstruct(behavior)
            error('nf_makefactorialdesign:BadBehavior', 'TF.behavior must be a struct array.');
        end
        behaviorTable = struct2table(behavior);
        return
    end

    if numel(TF) > 1
        behaviorTable = struct2table(TF);
        return
    end
end

error('nf_makefactorialdesign:BadInput', ...
    'Input must be a TF struct with TF.behavior, a behavior struct array, or a table.');

end

function [hasSubject, subjectFieldResolved] = local_resolve_subject_field(behaviorTable, subjectField)

hasSubject = false;
subjectFieldResolved = '';

vars = behaviorTable.Properties.VariableNames;

if ~isempty(subjectField)
    if isstring(subjectField)
        subjectField = char(subjectField);
    end
    if ~ischar(subjectField)
        error('nf_makefactorialdesign:BadSubjectField', 'subjectField must be a string.');
    end

    if any(strcmp(vars, subjectField))
        hasSubject = true;
        subjectFieldResolved = subjectField;
        return
    end

    % case-insensitive match
    for iV = 1:numel(vars)
        if strcmpi(vars{iV}, subjectField)
            hasSubject = true;
            subjectFieldResolved = vars{iV};
            return
        end
    end

    hasSubject = false;
    subjectFieldResolved = '';
    return
end

% auto-detect
candidates = { ...
    'subID', ...
    'subid', ...
    'subject', ...
    'Subject', ...
    'participant', ...
    'Participant', ...
    'ID', ...
    'id' ...
    };

for iC = 1:numel(candidates)
    nm = candidates{iC};
    if any(strcmp(vars, nm))
        hasSubject = true;
        subjectFieldResolved = nm;
        return
    end
end

% case-insensitive auto-detect
for iC = 1:numel(candidates)
    nm = candidates{iC};
    for iV = 1:numel(vars)
        if strcmpi(vars{iV}, nm)
            hasSubject = true;
            subjectFieldResolved = vars{iV};
            return
        end
    end
end

end

function subID = local_extract_subject_vector(behaviorTable, subjectFieldResolved)

col = behaviorTable.(subjectFieldResolved);

if iscell(col)
    if numel(col) ~= size(behaviorTable, 1)
        error('nf_makefactorialdesign:BadSubjectColumn', 'Subject field has unexpected size.');
    end
    subID = col;
    return
end

if isstring(col)
    subID = cellstr(col);
    return
end

if ischar(col)
    subID = cellstr(col);
    return
end

if iscategorical(col)
    subID = cellstr(string(col));
    return
end

if isnumeric(col) || islogical(col)
    if size(col, 2) ~= 1
        error('nf_makefactorialdesign:BadSubjectColumn', 'Subject field must be Nx1.');
    end
    subID = col(:);
    return
end

error('nf_makefactorialdesign:BadSubjectColumn', ...
    'Subject field "%s" must be numeric/logical/string/cellstr/categorical.', subjectFieldResolved);

end

function [baseX, isBinary, mu, sigma] = local_build_base_predictors(behaviorTable, baseNames, binaryTol, standardize)

nObs = size(behaviorTable, 1);
nVars = size(behaviorTable, 2);

baseX = nan(nObs, nVars);
isBinary = false(1, nVars);
mu = nan(1, nVars);
sigma = nan(1, nVars);

for iVar = 1:nVars

    col = behaviorTable{:, iVar};

    if iscell(col)
        error('nf_makefactorialdesign:NonNumericField', ...
            'Behavior field "%s" is a cell array. Convert to numeric scalars before design build.', ...
            baseNames{iVar});
    end

    if isstring(col) || ischar(col) || iscategorical(col)
        error('nf_makefactorialdesign:NonNumericField', ...
            'Behavior field "%s" is non-numeric (string/char/categorical). Convert to numeric scalars first.', ...
            baseNames{iVar});
    end

    if islogical(col)
        isBinary(iVar) = true;
        baseX(:, iVar) = double(col);
        mu(iVar) = nanmean(baseX(:, iVar));
        sigma(iVar) = nanstd(baseX(:, iVar), 0);
        continue
    end

    if ~isnumeric(col)
        error('nf_makefactorialdesign:NonNumericField', ...
            'Behavior field "%s" must be numeric scalar values per observation.', baseNames{iVar});
    end

    if size(col, 2) ~= 1
        error('nf_makefactorialdesign:NonScalarField', ...
            'Behavior field "%s" must be Nx1 numeric values per observation.', baseNames{iVar});
    end

    x = double(col);

    isBinary(iVar) = local_is_binary_numeric(x, binaryTol);

    if isBinary(iVar) == true
        baseX(:, iVar) = x;
        mu(iVar) = nanmean(x);
        sigma(iVar) = nanstd(x, 0);
    else
        mu(iVar) = nanmean(x);
        sigma(iVar) = nanstd(x, 0);

        if isnan(sigma(iVar)) || sigma(iVar) == 0
            baseX(:, iVar) = zeros(size(x));
        elseif standardize == true
            baseX(:, iVar) = (x - mu(iVar)) ./ sigma(iVar);
        else
            baseX(:, iVar) = x;
        end
    end

end

end

function tf = local_is_binary_numeric(x, binaryTol)

xv = x(~isnan(x));

if isempty(xv)
    tf = false;
    return
end

u = unique(xv);

if numel(u) ~= 2
    tf = false;
    return
end

u = sort(u(:));

is01 = abs(u(1) - 0) <= binaryTol && abs(u(2) - 1) <= binaryTol;
is12 = abs(u(1) - 1) <= binaryTol && abs(u(2) - 2) <= binaryTol;

if is01 || is12
    tf = true;
    return
end

tf = false;

end

function eligibleIdx = local_get_interaction_eligible_indices(isBinary, includeBinaryInteractions)

if includeBinaryInteractions == true
    eligibleIdx = 1:numel(isBinary);
else
    eligibleIdx = find(isBinary == false);
end

eligibleIdx = eligibleIdx(:)';

end

function [interactionX, interactionNames, interactionOrders] = local_build_interactions(baseXEligible, eligibleNames, maxOrder)

nObs = size(baseXEligible, 1);
p = size(baseXEligible, 2);

interactionX = [];
interactionNames = {};
interactionOrders = [];

if p < 2
    return
end

for order = 2:maxOrder

    combos = nchoosek(1:p, order);

    for iC = 1:size(combos, 1)

        idx = combos(iC, :);

        prodVec = ones(nObs, 1);

        for k = 1:numel(idx)
            prodVec = prodVec .* baseXEligible(:, idx(k));
        end

        name = eligibleNames{idx(1)};

        for k = 2:numel(idx)
            name = [name '_X_' eligibleNames{idx(k)}]; %#ok<AGROW>
        end

        interactionX = [interactionX, prodVec]; %#ok<AGROW>
        interactionNames = [interactionNames, {name}]; %#ok<AGROW>
        interactionOrders = [interactionOrders, order]; %#ok<AGROW>

    end

end

end

function cleanName = local_sanitize_varname(nameIn)

if isstring(nameIn)
    nameIn = char(nameIn);
end

if ~ischar(nameIn)
    error('nf_makefactorialdesign:BadVarName', 'Variable names must be char or string.');
end

cleanName = nameIn;

cleanName = regexprep(cleanName, '\s+', '_');
cleanName = regexprep(cleanName, '[^A-Za-z0-9_]', '_');

if isempty(cleanName)
    cleanName = 'Var';
end

if ~isletter(cleanName(1))
    cleanName = ['V_' cleanName];
end

cleanName = matlab.lang.makeValidName(cleanName);

end
