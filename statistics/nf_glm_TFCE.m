function tfce = nf_glm_TFCE(dataX, dataY, chanlocs, faxis, taxis, nperm, par)
%
% NF_GLM_TFCE  TFCE-corrected GLM regression on subject-level data.
%
% DROP-IN REPLACEMENT GUARANTEE (IMPORTANT)
% ----------------------------------------
% If dataX is numeric, behavior matches the prior function:
%   - same inputs and defaults
%   - same output fields and meanings
%
% EXTENSION (NEW)
% ---------------
% If dataX is a NeuroFreq TF struct (avg or singletrial),
% this wrapper will:
%   - detect TF mode (avg vs singletrial)
%   - run TFCE-GLM across ALL power/phase fields
%   - in avg mode, run per condition
%   - return results in:
%       tfce.power.<field>.results.<cond>
%       tfce.phase.<field>.results.<cond>
%     or .results.all in singletrial mode
%
% Design input dataY can be:
%   - numeric matrix (observations x predictors)
%   - table with named predictor columns (recommended)
%
% OUTPUT (numeric path)
% ---------------------
% tfce.beta       : beta maps for each regressor
% tfce.tstat      : t-stat maps (TFCE operates on these)
% tfce.p_vals     : TFCE-corrected p-values
% tfce.z_vals     : signed TFCE z-values from corrected p-values
% tfce.fixedNames : regressor labels
% tfce.designTab  : original design table (if provided)
% tfce.droppedVars: nonnumeric variables dropped from table design
%

if nargin < 7 || isempty(par)
    par = 0;
end

if nargin < 6 || isempty(nperm)
    nperm = 1000;
end

if nargin < 5 || isempty(taxis)
    taxis = [];
    disp('[nf_glm_TFCE]: no time axis supplied - assuming data have no time dimension.');
end

if nargin < 4 || isempty(faxis)
    faxis = [];
    disp('[nf_glm_TFCE]: no frequency axis supplied - assuming data have no frequency dimension.');
end

if nargin < 3 || isempty(chanlocs)
    chanlocs = [];
    disp('[nf_glm_TFCE]: no channel locations supplied - assuming single-channel data.');
end

if nargin < 2 || isempty(dataX) || isempty(dataY)
    error('nf_glm_TFCE:DataInputRequired',...
        'at least data must be supplied...see help for dimensions and other inputs.');
end

% -------------------------------------------------------------------------
% NEW: struct-aware wrapper path (numeric path remains original below)
% -------------------------------------------------------------------------
if isstruct(dataX) == 1
    tfce = local_glm_tfset_wrapper(dataX, dataY, chanlocs, faxis, taxis, nperm, par);
    return
end

% -------------------------------------------------------------------------
% ORIGINAL NUMERIC PATH (drop-in behavior)
% -------------------------------------------------------------------------
e1 = 0;
e2 = 0;
e3 = 0;

designTab = [];
droppedVars = {};

if istable(dataY)
    designTab = dataY;
    [designMat, fixedNames, droppedVars] = local_table_to_numeric_matrix(designTab);
else
    designMat = dataY;
    fixedNames = local_default_names(size(designMat, 2));
end

if isempty(designMat)
    error('Design matrix is empty after processing.');
end

if ~isnumeric(designMat)
    error('Design matrix must be numeric after processing.');
end

% SAFE IMPROVEMENT:
% If numeric dataX is not obs-first but design rows define obs count, try to fix.
nObs = size(designMat, 1);

if size(dataX, 1) ~= nObs
    dataX = local_force_obs_dim_first(dataX, nObs);
end

if ~isempty(chanlocs)
    if par == 0
        tfc = ept_TFCE_glm(dataX, designMat, chanlocs, 'nperm', nperm);
    else
        tfc = ept_TFCE_glm_par(dataX, designMat, chanlocs, 'nperm', nperm);
    end
else
    if par == 0
        tfc = ept_TFCE_glm(dataX, designMat, [], 'nperm', nperm, 'flag_ft', 1);
    else
        tfc = ept_TFCE_glm_par(dataX, designMat, [], 'nperm', nperm, 'flag_ft', 1);
    end
end

tfce = struct();

tfce.data = dataX;

tfce.corrvar = dataY;
tfce.designTab = designTab;

tfce.fixedNames = fixedNames;
tfce.droppedVars = droppedVars;

tfce.obs = squeeze(mean(dataX, 1));
tfce.sd = squeeze(std(dataX, 0, 1));

tfce.beta = tfc.Beta;
tfce.tstat = tfc.Obs;
tfce.p_vals = tfc.P_Values;

% TFCE-corrected p-values -> signed z-values
p = tfce.p_vals;

pmin = 1 ./ (nperm + 1);

if pmin < eps %#ok
    pmin = eps;
end

pmax = 1 - eps;

p = max(p, pmin);
p = min(p, pmax);

sgn = sign(tfce.tstat);

tfce.z_vals = sgn .* norminv(1 - p ./ 2);

if ~isempty(taxis)
    e3 = 1;
    tfce.times = taxis;
end

if ~isempty(faxis)
    e2 = 1;
    tfce.freqs = faxis;
end

if ~isempty(chanlocs)
    e1 = 1;
    tfce.chanlocs = chanlocs;
end

if e1 == 1 && e2 == 1 && e3 == 1
    tfce.type = 'glm_chXfsXts';
elseif e1 == 1 && e2 == 1
    tfce.type = 'glm_chXts';
elseif e1 == 1 && e3 == 1
    tfce.type = 'glm_chXts';
elseif e2 == 1 && e3 == 1
    tfce.type = 'glm_tf';
elseif e3 == 1
    tfce.type = 'glm_ts';
elseif e2 == 1
    tfce.type = 'glm_ts';
else
    tfce.type = 'glm_ts';
end

disp('[nf_glm_TFCE]: tfce-corrected glm computation complete'); 

end

% ====================================================================== %
% NEW: TF-STRUCT WRAPPER
% ====================================================================== %
function tfce = local_glm_tfset_wrapper(tfset, dataY, chanlocs, faxis, taxis, nperm, par)

% Fill axes from tfset if caller passed empty
if isempty(chanlocs) == 1
    if isfield(tfset, 'chanlocs') == 1
        chanlocs = tfset.chanlocs;
    end
end

if isempty(faxis) == 1
    if isfield(tfset, 'freqs') == 1
        faxis = tfset.freqs;
    end
end

if isempty(taxis) == 1
    if isfield(tfset, 'times') == 1
        taxis = tfset.times;
    end
end

% Prepare design (once)
designTab = [];
droppedVars = {};

if istable(dataY)
    designTab = dataY;
    [designMat, fixedNames, droppedVars] = local_table_to_numeric_matrix(designTab);
else
    designMat = dataY;
    fixedNames = local_default_names(size(designMat, 2));
end

if isempty(designMat) == 1
    error('Design matrix is empty after processing.');
end

if isnumeric(designMat) == 0
    error('Design matrix must be numeric after processing.');
end

% Detect mode
mode = local_detect_tf_mode(tfset);

if strcmpi(mode, 'stat') == 1
    error('nf_glm_TFCE:BadInput', 'Input TF struct appears to be stat output. GLM requires raw subject/trial-level data.');
end

% Gather fields
[groupList, fieldMap] = local_collect_power_phase_fields(tfset);

if isempty(groupList) == 1
    error('nf_glm_TFCE:NoFields', 'No power/phase fields found in TF struct.');
end

% Top-level output
tfce = struct();

tfce.method = 'tfce';
tfce.type = 'glm_tfset';
tfce.tf_mode = mode;

tfce.corrvar = dataY;
tfce.designTab = designTab;

tfce.fixedNames = fixedNames;
tfce.droppedVars = droppedVars;

tfce.source_scale = '';
tfce.source_aggtype = '';

if isfield(tfset, 'scale') == 1
    if ischar(tfset.scale) == 1
        tfce.source_scale = tfset.scale;
    end
end

if isfield(tfset, 'aggtype') == 1
    if ischar(tfset.aggtype) == 1
        tfce.source_aggtype = tfset.aggtype;
    end
end

if isempty(taxis) == 0
    tfce.times = taxis;
end

if isempty(faxis) == 0
    tfce.freqs = faxis;
end

if isempty(chanlocs) == 0
    tfce.chanlocs = chanlocs;
end

% Condition labels (avg mode)
condsRaw = {};
condsValid = {};

if strcmpi(mode, 'avg') == 1
    if isfield(tfset, 'conds') == 1
        if iscell(tfset.conds) == 1
            condsRaw = tfset.conds(:)';
        end
    end
end

if isempty(condsRaw) == 0
    condsValid = cell(size(condsRaw));
    for iC = 1:numel(condsRaw)
        condsValid{iC} = matlab.lang.makeValidName(condsRaw{iC});
    end
end

% Run across groups/fields
for iG = 1:numel(groupList)

    grp = groupList{iG};
    fns = fieldMap.(grp);

    if strcmpi(grp, 'power') == 1
        tfce.power = struct();
    end

    if strcmpi(grp, 'phase') == 1
        tfce.phase = struct();
    end

    for iF = 1:numel(fns)

        fn = fns{iF};
        X = tfset.(grp).(fn);

        if isnumeric(X) == 0
            continue
        end

        if isempty(X) == 1
            continue
        end

        fieldOut = struct();

        if strcmpi(mode, 'avg') == 1

            [Xsf, nObs, nChan, nFreq, nTime, nConds] = local_standardize_avg_subjectfirst(X, tfset, chanlocs, faxis, taxis);

            if size(designMat, 1) ~= nObs
                error('nf_glm_TFCE:DesignMismatch', ...
                    'Design rows (%d) do not match observations (%d) for %s.%s.', ...
                    size(designMat, 1), nObs, grp, fn);
            end

            if isempty(condsRaw) == 1
                condsRaw = local_default_conds(nConds);
                condsValid = cell(size(condsRaw));
                for iC = 1:numel(condsRaw)
                    condsValid{iC} = matlab.lang.makeValidName(condsRaw{iC});
                end
            end

            fieldOut.conds_raw = condsRaw;
            fieldOut.conds = condsValid;
            fieldOut.results = struct();

            for c = 1:nConds

                Xc = Xsf(:, :, :, :, c);
                Xc = reshape(Xc, nObs, nChan, nFreq, nTime);

                r = local_glm_core_numeric(Xc, dataY, chanlocs, faxis, taxis, nperm, par);

                r = local_prune_heavy_fields_glm(r);

                fieldOut.results.(condsValid{c}) = r;

            end

        elseif strcmpi(mode, 'singletrial') == 1

            [Xobs, nObs, nChan, nFreq, nTime] = local_standardize_singletrial_obsfirst(X, tfset, chanlocs, faxis, taxis);

            if size(designMat, 1) ~= nObs
                error('nf_glm_TFCE:DesignMismatch', ...
                    'Design rows (%d) do not match trial observations (%d) for %s.%s.', ...
                    size(designMat, 1), nObs, grp, fn);
            end

            Xobs = reshape(Xobs, nObs, nChan, nFreq, nTime);

            r = local_glm_core_numeric(Xobs, dataY, chanlocs, faxis, taxis, nperm, par);

            r = local_prune_heavy_fields_glm(r);

            fieldOut.results = struct();
            fieldOut.results.all = r;

        else
            error('nf_glm_TFCE:BadMode', 'Unknown TF mode: %s', mode);
        end

        if strcmpi(grp, 'power') == 1
            tfce.power.(fn) = fieldOut;
        end

        if strcmpi(grp, 'phase') == 1
            tfce.phase.(fn) = fieldOut;
        end

    end

end

disp(['Recommended methods text: GLM with significance testing' ...
      ' with correction for multiple comparisons was conducted using ' ...
      'threshold-free cluster enhancement (TFCE; Mensen & Khatami, 2013) ' ...
      'using functions from the ept-TFCE toolbox (Mensen & Khatami, 2013).']);

end

function r = local_glm_core_numeric(dataX, dataY, chanlocs, faxis, taxis, nperm, par)

e1 = 0;
e2 = 0;
e3 = 0;

designTab = [];
droppedVars = {};

if istable(dataY)
    designTab = dataY;
    [designMat, fixedNames, droppedVars] = local_table_to_numeric_matrix(designTab);
else
    designMat = dataY;
    fixedNames = local_default_names(size(designMat, 2));
end

if isempty(designMat) == 1
    error('Design matrix is empty after processing.');
end

if isnumeric(designMat) == 0
    error('Design matrix must be numeric after processing.');
end

nObs = size(designMat, 1);

if size(dataX, 1) ~= nObs
    dataX = local_force_obs_dim_first(dataX, nObs);
end

if ~isempty(chanlocs)
    if par == 0
        tfc = ept_TFCE_glm(dataX, designMat, chanlocs, 'nperm', nperm);
    else
        tfc = ept_TFCE_glm_par(dataX, designMat, chanlocs, 'nperm', nperm);
    end
else
    if par == 0
        tfc = ept_TFCE_glm(dataX, designMat, [], 'nperm', nperm, 'flag_ft', 1);
    else
        tfc = ept_TFCE_glm_par(dataX, designMat, [], 'nperm', nperm, 'flag_ft', 1);
    end
end

r = struct();

r.data = dataX;

r.corrvar = dataY;
r.designTab = designTab;

r.fixedNames = fixedNames;
r.droppedVars = droppedVars;

r.obs = squeeze(mean(dataX, 1));
r.sd = squeeze(std(dataX, 0, 1));

r.beta = tfc.Beta;
r.tstat = tfc.Obs;
r.p_vals = tfc.P_Values;

p = r.p_vals;

pmin = 1 ./ (nperm + 1);

if pmin < eps %#ok
    pmin = eps;
end

pmax = 1 - eps;

p = max(p, pmin);
p = min(p, pmax);

sgn = sign(r.tstat);

r.z_vals = sgn .* norminv(1 - p ./ 2);

if ~isempty(taxis)
    e3 = 1;
    r.times = taxis;
end

if ~isempty(faxis)
    e2 = 1;
    r.freqs = faxis;
end

if ~isempty(chanlocs)
    e1 = 1;
    r.chanlocs = chanlocs;
end

if e1 == 1 && e2 == 1 && e3 == 1
    r.type = 'glm_chXfsXts';
elseif e1 == 1 && e2 == 1
    r.type = 'glm_chXts';
elseif e1 == 1 && e3 == 1
    r.type = 'glm_chXts';
elseif e2 == 1 && e3 == 1
    r.type = 'glm_tf';
elseif e3 == 1
    r.type = 'glm_ts';
elseif e2 == 1
    r.type = 'glm_ts';
else
    r.type = 'glm_ts';
end

end

function r = local_prune_heavy_fields_glm(r)

if isfield(r, 'data') == 1
    r.data = [];
end

if isfield(r, 'corrvar') == 1
    r.corrvar = [];
end

if isfield(r, 'designTab') == 1
    r.designTab = [];
end

if isfield(r, 'obs') == 1
    r.obs = [];
end

if isfield(r, 'sd') == 1
    r.sd = [];
end

end

function mode = local_detect_tf_mode(tfset)

mode = 'unknown'; %#ok

if isfield(tfset, 'scale') == 1
    if ischar(tfset.scale) == 1
        if strcmpi(tfset.scale, 'stat') == 1
            mode = 'stat';
            return
        end
    end
end

if isfield(tfset, 'aggtype') == 1
    if ischar(tfset.aggtype) == 1
        if strcmpi(tfset.aggtype, 'stat') == 1
            mode = 'stat';
            return
        end
    end
end

if isfield(tfset, 'conds') == 1
    if iscell(tfset.conds) == 1
        if isempty(tfset.conds) == 0
            mode = 'avg';
            return
        end
    end
end

if isfield(tfset, 'ntrial') == 1
    mode = 'singletrial';
    return
end

if isfield(tfset, 'ntrials') == 1
    mode = 'singletrial';
    return
end

if isfield(tfset, 'trial_subjectindex') == 1
    mode = 'singletrial';
    return
end

mode = 'singletrial';

end

function [groupList, fieldMap] = local_collect_power_phase_fields(tfset)

groupList = {};
fieldMap = struct();

if isfield(tfset, 'power') == 1
    if isstruct(tfset.power) == 1
        fns = fieldnames(tfset.power);
        if isempty(fns) == 0
            groupList{end + 1, 1} = 'power'; 
            fieldMap.power = fns;
        end
    end
end

if isfield(tfset, 'phase') == 1
    if isstruct(tfset.phase) == 1
        fns = fieldnames(tfset.phase);
        if isempty(fns) == 0
            groupList{end + 1, 1} = 'phase'; 
            fieldMap.phase = fns;
        end
    end
end

end

function [Xsf, nObs, nChan, nFreq, nTime, nConds] = local_standardize_avg_subjectfirst(X, tfset, chanlocs, faxis, taxis)

if isempty(chanlocs) == 0
    nChan = numel(chanlocs);
else
    nChan = size(X, 1);
end

if isempty(faxis) == 0
    nFreq = numel(faxis);
else
    nFreq = size(X, 2);
end

if isempty(taxis) == 0
    nTime = numel(taxis);
else
    nTime = size(X, 3);
end

nConds = [];

if isfield(tfset, 'conds') == 1
    if iscell(tfset.conds) == 1
        if isempty(tfset.conds) == 0
            nConds = numel(tfset.conds);
        end
    end
end

if ndims(X) == 4

    if isempty(nConds) == 1
        nConds = size(X, 4);
    end

    if size(X, 1) ~= nChan
        error('nf_glm_TFCE:BadShape', 'Avg field channel dimension mismatch.');
    end

    if size(X, 2) ~= nFreq
        error('nf_glm_TFCE:BadShape', 'Avg field frequency dimension mismatch.');
    end

    if size(X, 3) ~= nTime
        error('nf_glm_TFCE:BadShape', 'Avg field time dimension mismatch.');
    end

    if size(X, 4) ~= nConds
        error('nf_glm_TFCE:BadShape', 'Avg field condition dimension mismatch.');
    end

    Xsf = reshape(X, 1, nChan, nFreq, nTime, nConds);

elseif ndims(X) == 5

    if isempty(nConds) == 1
        nConds = size(X, 5);
    end

    sz = size(X);

    isA = false;
    isB = false;

    if sz(2) == nChan && sz(3) == nFreq && sz(4) == nTime && sz(5) == nConds
        isA = true;
    end

    if sz(1) == nChan && sz(2) == nFreq && sz(3) == nTime && sz(4) == nConds
        isB = true;
    end

    if isA == 1
        Xsf = X;
    elseif isB == 1
        Xsf = permute(X, [5 1 2 3 4]);
    else
        error('nf_glm_TFCE:BadShape', ...
            'Unrecognized avg field layout. Supported: sub x chan x freq x time x cond OR chan x freq x time x cond x sub.');
    end

else
    error('nf_glm_TFCE:BadShape', 'Avg fields must be 4D or 5D.');
end

nObs = size(Xsf, 1);

end

function [Xobs, nObs, nChan, nFreq, nTime] = local_standardize_singletrial_obsfirst(X, ~, chanlocs, faxis, taxis)

if isempty(chanlocs) == 0
    nChan = numel(chanlocs);
else
    nChan = size(X, 1);
end

if isempty(faxis) == 0
    nFreq = numel(faxis);
else
    nFreq = size(X, 2);
end

if isempty(taxis) == 0
    nTime = numel(taxis);
else
    nTime = size(X, 3);
end

if ndims(X) ~= 4
    error('nf_glm_TFCE:BadShape', 'Single-trial fields must be 4D (chan x freq x time x trial).');
end

if size(X, 1) ~= nChan
    error('nf_glm_TFCE:BadShape', 'Single-trial field channel dimension mismatch.');
end

if size(X, 2) ~= nFreq
    error('nf_glm_TFCE:BadShape', 'Single-trial field frequency dimension mismatch.');
end

if size(X, 3) ~= nTime
    error('nf_glm_TFCE:BadShape', 'Single-trial field time dimension mismatch.');
end

Xobs = permute(X, [4 1 2 3]);

nObs = size(Xobs, 1);

end

function conds = local_default_conds(nConds)

conds = cell(1, nConds);

for i = 1:nConds
    conds{i} = ['cond' num2str(i)];
end

end

function X = local_force_obs_dim_first(X, nObs)

if isempty(X) == 1
    error('nf_glm_TFCE:BadData', 'dataX is empty.');
end

if isnumeric(X) == 0
    error('nf_glm_TFCE:BadData', 'dataX must be numeric.');
end

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
    error('nf_glm_TFCE:ObsMismatch', ...
        'Could not find an observation dimension matching nObs=%d in dataX.', nObs);
end

if numel(matchDims) > 1
    error('nf_glm_TFCE:ObsAmbiguous', ...
        'Multiple dimensions match nObs=%d in dataX; cannot auto-permute safely.', nObs);
end

dObs = matchDims(1);

perm = 1:nd;
perm(dObs) = [];
perm = [dObs perm];

X = permute(X, perm);

end

% ====================================================================== %
% Local helpers (same as your current)
% ====================================================================== %
function [X, names, droppedVars] = local_table_to_numeric_matrix(T)

vars = T.Properties.VariableNames;

keepIdx = [];
names = {};
droppedVars = {};

for iV = 1:numel(vars)

    col = T.(vars{iV});

    if isnumeric(col) || islogical(col)

        if size(col, 2) ~= 1
            error('Table variable "%s" must be Nx1 numeric.', vars{iV});
        end

        keepIdx = [keepIdx iV]; %#ok<AGROW>
        names = [names vars(iV)]; %#ok<AGROW>

    else

        droppedVars = [droppedVars vars(iV)]; %#ok<AGROW>

    end

end

if isempty(keepIdx)
    error('No numeric predictors found in design table.');
end

X = T{:, keepIdx};
X = double(X);

end

function names = local_default_names(p)

names = cell(1, p);

for iP = 1:p
    names{iP} = ['X' num2str(iP)];
end

end







