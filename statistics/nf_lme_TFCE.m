function tfce = nf_lme_TFCE(dataX, designTab, varargin)
% NF_LME_TFCE  TFCE-corrected mixed-effects regression for NeuroFreq TF data (TABLE-ONLY design).
%
% DESIGN (TABLE ONLY)
% -------------------
% - designTab MUST be a table.
% - designTab MUST include a subject ID column (auto-detected or via 'subjectField').
% - All other fixed effects MUST be numeric/logical Nx1 columns.
% - Applies the same LME+TFCE model to EVERY TOP-LEVEL power field: power*, pow_*
%
% INPUTS
% ------
% tfce = nf_lme_TFCE(dataX, designTab, chanlocs, faxis, taxis, nperm, par, ...)
%
% dataX may be:
%   1) numeric array with an observation dimension matching height(designTab)
%   2) NeuroFreq TF struct: LME is applied to EVERY top-level power field (power*/pow_*)
%
% Name-value options
% ------------------
%   'subjectField'     : subject ID column name in designTab (default: auto)
%   'requireSubjectID' : true/false (default true)
%   'randomSlopes'     : numeric indices OR predictor names (cellstr) into fixed effects
%   'fitMethod'        : 'REML' (default) or 'ML'
%   'uncorrelated'     : true (default)
%   'flag_ft'          : 0/1 passed to ept_TFCE_lme when chanlocs empty (default 1)
%   'nperm'            : override permutations
%   'par'              : override parallel flag 0/1
%
% NOTES ON AVG vs SINGLETRIAL
% ---------------------------
% - Single-trial aggregated: observations = trials (height(designTab) == nTrialsAll)
% - Avg aggregated: observations = subjects (height(designTab) == nSubjects)
%   If data has an extra slice dim (e.g., condition in dim 5), this function runs
%   the same model separately for each slice.
%

% -----------------------------
% Defaults
% -----------------------------
randomSlopes = [];
fitMethod = 'REML';
uncorrelated = true;

subjectField = '';
requireSubjectID = true;

chanlocs = [];
faxis = [];
taxis = [];
nperm = 1000;
par = 1;

flag_ft = 1;

% -----------------------------
% Basic checks
% -----------------------------
if exist('fitlme', 'file') ~= 2
    error('nf_lme_TFCE:NoFitLME', 'fitlme not found. Requires Statistics and Machine Learning Toolbox.');
end

if istable(designTab) == 0
    error('nf_lme_TFCE:BadDesign', 'designTab must be a table (table-only design is enforced).');
end

if isempty(dataX) == 1
    error('nf_lme_TFCE:BadData', 'dataX must be non-empty.');
end

% -----------------------------
% Parse positional inputs: chanlocs, faxis, taxis, nperm, par
% -----------------------------
nvStart = 1;

if numel(varargin) >= 1
    if ischar(varargin{1}) == 0 && isstring(varargin{1}) == 0
        chanlocs = varargin{1};
        nvStart = 2;
    end
end

if numel(varargin) >= 2
    if ischar(varargin{2}) == 0 && isstring(varargin{2}) == 0
        faxis = varargin{2};
        nvStart = 3;
    end
end

if numel(varargin) >= 3
    if ischar(varargin{3}) == 0 && isstring(varargin{3}) == 0
        taxis = varargin{3};
        nvStart = 4;
    end
end

if numel(varargin) >= 4
    if ischar(varargin{4}) == 0 && isstring(varargin{4}) == 0
        nperm = varargin{4};
        nvStart = 5;
    end
end

if numel(varargin) >= 5
    if ischar(varargin{5}) == 0 && isstring(varargin{5}) == 0
        par = varargin{5};
        nvStart = 6;
    end
end

% -----------------------------
% Parse name-value options
% -----------------------------
if numel(varargin) >= nvStart
    nv = varargin(nvStart:end);

    if mod(numel(nv), 2) ~= 0
        error('nf_lme_TFCE:BadArgs', 'Name-value inputs must come in pairs.');
    end

    for k = 1:2:numel(nv)

        key = nv{k};
        val = nv{k + 1};

        if isstring(key) == 1
            key = char(key);
        end

        if ischar(key) == 0
            error('nf_lme_TFCE:BadArgs', 'Name-value key must be a string.');
        end

        key = lower(key);

        if strcmp(key, 'randomslopes') == 1
            randomSlopes = val;
        elseif strcmp(key, 'fitmethod') == 1
            fitMethod = val;
        elseif strcmp(key, 'uncorrelated') == 1
            uncorrelated = val;
        elseif strcmp(key, 'subjectfield') == 1
            subjectField = val;
        elseif strcmp(key, 'requiresubjectid') == 1
            requireSubjectID = val;
        elseif strcmp(key, 'flag_ft') == 1
            flag_ft = val;
        elseif strcmp(key, 'par') == 1
            par = val;
        elseif strcmp(key, 'nperm') == 1
            nperm = val;
        else
            error('nf_lme_TFCE:BadArgs', 'Unknown option: %s', key);
        end

    end
end

% -----------------------------
% Design extraction (TABLE ONLY)
% -----------------------------
[subID, subjectFieldResolved] = local_extract_subject_from_design(designTab, subjectField, requireSubjectID);

[dataY, fixedNames] = local_extract_fixed_from_design(designTab, subjectFieldResolved);

nObs = height(designTab);

if numel(subID) ~= nObs
    error('nf_lme_TFCE:DesignMismatch', 'designTab height must match length of subject ID column.');
end

if size(dataY, 1) ~= nObs
    error('nf_lme_TFCE:DesignMismatch', 'Fixed-effect rows must match height(designTab).');
end

if any(~isfinite(dataY(:)))
    error('nf_lme_TFCE:BadDesign', 'Fixed-effect design contains NaN/Inf. Clean or impute before LME.');
end

if numel(unique(subID)) < 2
    error('nf_lme_TFCE:TooFewSubjects', 'Need at least 2 unique subjects for mixed-effects modeling.');
end

% Drop constant predictors (including an intercept-like all-ones column)
[dataY, fixedNames, dropped] = local_drop_constant_predictors(dataY, fixedNames);

if isempty(fixedNames) == 1
    error('nf_lme_TFCE:NoFixed', 'No usable fixed-effect predictors remain after dropping constant columns.');
end

randomSlopesIdx = local_random_slopes_to_indices(randomSlopes, fixedNames, size(dataY, 2));

% =========================================================================
% STRUCT MODE: apply to every top-level power field
% =========================================================================
if isstruct(dataX) == 1

    tfIn = dataX;

    if isempty(chanlocs) == 1
        if isfield(tfIn, 'chanlocs') == 1
            chanlocs = tfIn.chanlocs;
        end
    end

    if isempty(faxis) == 1
        if isfield(tfIn, 'freqs') == 1
            faxis = tfIn.freqs;
        end
    end

    if isempty(taxis) == 1
        if isfield(tfIn, 'times') == 1
            taxis = tfIn.times;
        end
    end

    powerFields = local_find_top_level_power_fields(tfIn);

    if isempty(powerFields) == 1
        error('nf_lme_TFCE:NoPowerFields', 'No top-level power fields found (power*/pow_*).');
    end

    tfce = struct();

    tfce.scale = 'stat';
    tfce.method = 'lme_tfce';

    tfce.designTab = designTab;
    tfce.subjectField = subjectFieldResolved;
    tfce.subID = subID;

    tfce.fixedNames = fixedNames;
    tfce.prednames = fixedNames;

    tfce.dropped_predictors = dropped;

    tfce.randomSlopes = randomSlopesIdx;
    tfce.fitMethod = fitMethod;
    tfce.uncorrelated = uncorrelated;
    tfce.flag_ft = flag_ft;

    tfce.nperm = nperm;
    tfce.par = par;

    tfce.power_fields = powerFields;

    if isempty(faxis) == 0
        tfce.freqs = faxis;
    end

    if isempty(taxis) == 0
        tfce.times = taxis;
    end

    if isempty(chanlocs) == 0
        tfce.chanlocs = chanlocs;
    end

    for iF = 1:numel(powerFields)

        fname = powerFields{iF};
        Xraw = tfIn.(fname);

        if isnumeric(Xraw) == 0
            continue
        end

        if isempty(Xraw) == 1
            continue
        end

        Xobs = local_force_obs_dim_first(Xraw, nObs);
        Xobs = local_force_obs_first_to_4d_or_5d(Xobs);

        if ndims(Xobs) == 4

            r = local_run_lme_once(Xobs, dataY, subID, chanlocs, faxis, taxis, nperm, par, flag_ft, fixedNames, randomSlopesIdx, fitMethod, uncorrelated);
            tfce.(matlab.lang.makeValidName(fname)) = r;

        elseif ndims(Xobs) == 5

            nSlices = size(Xobs, 5);

            outS = struct();
            outS.slice_index = (1:nSlices).';
            outS.slice_labels = local_default_slice_labels(nSlices);
            outS.slices = cell(nSlices, 1);

            for s = 1:nSlices

                Xs = Xobs(:, :, :, :, s);

                r = local_run_lme_once(Xs, dataY, subID, chanlocs, faxis, taxis, nperm, par, flag_ft, fixedNames, randomSlopesIdx, fitMethod, uncorrelated);
                outS.slices{s} = r;

            end

            tfce.(matlab.lang.makeValidName(fname)) = outS;

        else

            error('nf_lme_TFCE:BadShape', 'Unexpected dimensionality after obs-first coercion for field %s.', fname);

        end

    end

    return

end

% =========================================================================
% NUMERIC MODE: run once on numeric dataX (still TABLE design)
% =========================================================================
if isnumeric(dataX) == 0
    error('nf_lme_TFCE:BadData', 'dataX must be numeric or a TF struct.');
end

Xobs = local_force_obs_dim_first(dataX, nObs);
Xobs = local_force_obs_first_to_4d_or_5d(Xobs);

if ndims(Xobs) ~= 4
    error('nf_lme_TFCE:BadShape', 'Numeric dataX must resolve to 4D (obs x chan x freq x time).');
end

tfce = local_run_lme_once(Xobs, dataY, subID, chanlocs, faxis, taxis, nperm, par, flag_ft, fixedNames, randomSlopesIdx, fitMethod, uncorrelated);

tfce.designTab = designTab;
tfce.fixedNames = fixedNames;
tfce.prednames = fixedNames;

tfce.subjectField = subjectFieldResolved;
tfce.subID = subID;

tfce.dropped_predictors = dropped;

disp('[nf_lme_TFCE]: tfce-corrected linear mixed effect model computation complete');

end

% ====================================================================== %
% Local helpers
% ====================================================================== %

function r = local_run_lme_once(Xobs4, dataY, subID, chanlocs, faxis, taxis, nperm, par, flag_ft, fixedNames, randomSlopesIdx, fitMethod, uncorrelated)

if isempty(Xobs4) == 1
    error('nf_lme_TFCE:BadData', 'Empty data passed to LME core.');
end

if ndims(Xobs4) ~= 4
    error('nf_lme_TFCE:BadData', 'LME core expects 4D: obs x chan x freq x time.');
end

if size(Xobs4, 1) ~= numel(subID)
    error('nf_lme_TFCE:ObsMismatch', 'dataX obs count does not match subID length.');
end

nChan = size(Xobs4, 2);
nFreq = size(Xobs4, 3);
nTime = size(Xobs4, 4);

if isempty(chanlocs) == 1
    if nChan > 1
        error('nf_lme_TFCE:ChanlocRequired', 'Multi-channel data requires chanlocs (or pre-average channels to 1).');
    end
else
    if numel(chanlocs) ~= nChan
        error('nf_lme_TFCE:ChanlocMismatch', 'chanlocs count (%d) does not match data channels (%d).', numel(chanlocs), nChan);
    end
end

if isempty(faxis) == 0
    if numel(faxis) ~= nFreq
        error('nf_lme_TFCE:FreqMismatch', 'faxis length (%d) does not match data freqs (%d).', numel(faxis), nFreq);
    end
end

if isempty(taxis) == 0
    if numel(taxis) ~= nTime
        error('nf_lme_TFCE:TimeMismatch', 'taxis length (%d) does not match data times (%d).', numel(taxis), nTime);
    end
end

tfc = ept_TFCE_lme( ...
    Xobs4, ...
    dataY, ...
    subID, ...
    chanlocs, ...
    'nperm', nperm, ...
    'par', par, ...
    'flag_ft', flag_ft, ...
    'randomslopes', randomSlopesIdx, ...
    'fitmethod', fitMethod, ...
    'uncorrelated', uncorrelated, ...
    'fixednames', fixedNames );

r = struct();

r.beta = tfc.Beta;
r.tstat = tfc.Obs;
r.p_vals = tfc.P_Values;

r.fixedNames = tfc.fixedNames;

r.nperm = nperm;
r.par = par;

r.fitMethod = fitMethod;
r.uncorrelated = uncorrelated;
r.randomSlopes = randomSlopesIdx;
r.flag_ft = flag_ft;

if isempty(faxis) == 0
    r.freqs = faxis;
end

if isempty(taxis) == 0
    r.times = taxis;
end

if isempty(chanlocs) == 0
    r.chanlocs = chanlocs;
end

if isfield(tfc, 'type') == 1
    r.type = tfc.type;
else
    r.type = 'lme';
end

end

function powerFields = local_find_top_level_power_fields(S)

fns = fieldnames(S);
powerFields = {};

for i = 1:numel(fns)

    nm = fns{i};

    if strcmpi(nm, 'freqs') == 1
        continue
    end

    if strcmpi(nm, 'times') == 1
        continue
    end

    if strcmpi(nm, 'chanlocs') == 1
        continue
    end

    if strcmpi(nm, 'behavior') == 1
        continue
    end

    if strcmpi(nm, 'scale') == 1
        continue
    end

    if strcmpi(nm, 'aggtype') == 1
        continue
    end

    if strcmpi(nm, 'stat') == 1
        continue
    end

    v = S.(nm);

    if isnumeric(v) == 0
        continue
    end

    isPow = false;

    if strncmpi(nm, 'power', 5) == 1
        isPow = true;
    end

    if strncmpi(nm, 'pow_', 4) == 1
        isPow = true;
    end

    if isPow == 0
        continue
    end

    powerFields{end + 1, 1} = nm; %#ok<AGROW>

end

idx = find(strcmpi(powerFields, 'power'), 1);

if isempty(idx) == 0
    powerFields(idx) = [];
    powerFields = [{'power'}; powerFields];
end

end

function X = local_force_obs_dim_first(X, nObs)

if size(X, 1) == nObs
    return
end

nd = ndims(X);

matchDims = [];

for d = 1:nd
    if size(X, d) == nObs
        matchDims = [matchDims d]; %#ok<AGROW>
    end
end

if isempty(matchDims) == 1
    error('nf_lme_TFCE:ObsMismatch', 'Could not find an observation dimension matching nObs=%d in data.', nObs);
end

if numel(matchDims) > 1
    error('nf_lme_TFCE:ObsAmbiguous', 'Multiple dimensions match nObs=%d; cannot auto-permute safely.', nObs);
end

dObs = matchDims(1);

perm = 1:nd;
perm(dObs) = [];
perm = [dObs perm];

X = permute(X, perm);

end

function X = local_force_obs_first_to_4d_or_5d(X)

nd = ndims(X);

if nd == 2
    X = reshape(X, [size(X, 1) 1 1 size(X, 2)]);
    return
end

if nd == 3
    X = reshape(X, [size(X, 1) 1 size(X, 2) size(X, 3)]);
    return
end

if nd == 4
    return
end

if nd == 5
    return
end

error('nf_lme_TFCE:BadShape', 'Data has unsupported dimensionality after obs-first coercion: ndims=%d.', nd);

end

function idx = local_random_slopes_to_indices(randomSlopes, fixedNames, nFixed)

idx = [];

if isempty(randomSlopes) == 1
    return
end

if isnumeric(randomSlopes) == 1
    idx = unique(randomSlopes(:)');
else

    if ischar(randomSlopes) == 1 || isstring(randomSlopes) == 1
        randomSlopes = {char(randomSlopes)};
    end

    if iscell(randomSlopes) == 0
        error('nf_lme_TFCE:BadRandomSlopes', 'randomSlopes must be numeric indices OR cellstr predictor names.');
    end

    idx = [];

    for iR = 1:numel(randomSlopes)

        nm = randomSlopes{iR};

        if isstring(nm) == 1
            nm = char(nm);
        end

        hit = find(strcmp(fixedNames, nm), 1);

        if isempty(hit) == 1
            error('nf_lme_TFCE:BadRandomSlopes', 'randomSlopes name "%s" not found in fixed effects.', nm);
        end

        idx = [idx hit]; %#ok<AGROW>

    end

    idx = unique(idx);

end

if any(idx < 1) || any(idx > nFixed)
    error('nf_lme_TFCE:BadRandomSlopes', 'randomSlopes contains indices outside the fixed-effects range.');
end

end

function [subID, subjField] = local_extract_subject_from_design(T, subjectField, requireSubjectID)

vars = T.Properties.VariableNames;
subjField = '';

if isempty(subjectField) == 0

    if isstring(subjectField) == 1
        subjectField = char(subjectField);
    end

    if ischar(subjectField) == 0
        error('nf_lme_TFCE:BadSubjectField', 'subjectField must be a string.');
    end

    if any(strcmp(vars, subjectField))
        subjField = subjectField;
    else
        for iV = 1:numel(vars)
            if strcmpi(vars{iV}, subjectField) == 1
                subjField = vars{iV};
            end
        end
    end

end

if isempty(subjField) == 1

    candidates = {'subID', 'subid', 'subject', 'Subject', 'participant', 'Participant', 'ID', 'id'};

    for iC = 1:numel(candidates)
        nm = candidates{iC};
        if any(strcmp(vars, nm))
            subjField = nm;
        end
    end

    if isempty(subjField) == 1
        for iC = 1:numel(candidates)
            nm = candidates{iC};
            for iV = 1:numel(vars)
                if strcmpi(vars{iV}, nm) == 1
                    subjField = vars{iV};
                end
            end
        end
    end

end

if isempty(subjField) == 1
    if requireSubjectID == 1
        error('nf_lme_TFCE:NoSubjectID', 'No subject ID field found in design table.');
    else
        error('nf_lme_TFCE:NoSubjectID', 'Subject ID field not found, and requireSubjectID=false is not supported here.');
    end
end

col = T.(subjField);

if isnumeric(col) == 1 || islogical(col) == 1
    subID = double(col(:));
    return
end

if iscategorical(col) == 1
    subID = double(grp2idx(col));
    return
end

if isstring(col) == 1
    subID = double(grp2idx(categorical(col)));
    return
end

if iscell(col) == 1
    subID = double(grp2idx(categorical(col)));
    return
end

error('nf_lme_TFCE:BadSubjectIDType', 'Subject field "%s" is not a supported type.', subjField);

end

function [dataY, fixedNames] = local_extract_fixed_from_design(T, subjField)

vars = T.Properties.VariableNames;

keepVars = {};
fixedNames = {};

for iV = 1:numel(vars)

    nm = vars{iV};

    if strcmp(nm, subjField) == 1
        continue
    end

    col = T.(nm);

    if isnumeric(col) == 1 || islogical(col) == 1

        if size(col, 2) ~= 1
            error('nf_lme_TFCE:BadFixed', 'Fixed-effect variable "%s" must be Nx1 numeric/logical.', nm);
        end

        keepVars = [keepVars {nm}]; %#ok<AGROW>
        fixedNames = [fixedNames {nm}]; %#ok<AGROW>

    else

        error('nf_lme_TFCE:BadFixed', 'Non-numeric fixed-effect column "%s" detected. Convert to numeric first.', nm);

    end

end

if isempty(keepVars) == 1
    error('nf_lme_TFCE:NoFixed', 'No numeric fixed-effect predictors found in design table.');
end

dataY = T{:, keepVars};
dataY = double(dataY);

end

function [Y2, names2, dropped] = local_drop_constant_predictors(Y, names)

tol = 1e-12;

keep = true(1, size(Y, 2));
dropped = {};

for j = 1:size(Y, 2)

    col = Y(:, j);

    if all(abs(col - col(1)) < tol)
        keep(j) = false;
        dropped{end + 1, 1} = names{j}; %#ok<AGROW>
    end

end

Y2 = Y(:, keep);
names2 = names(keep);

end

function labels = local_default_slice_labels(nSlices)

labels = cell(nSlices, 1);

for i = 1:nSlices
    labels{i} = ['slice' num2str(i)];
end

end






