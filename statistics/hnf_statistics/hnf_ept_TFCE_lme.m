function Results = hnf_ept_TFCE_lme(data1, data2, subID, chanlocs, varargin)
% EPT_TFCE_LME  TFCE-corrected mass-univariate mixed-effects regression.
%
% Outputs:
%   Results.Obs      = fixed-effect t-stat maps
%   Results.Beta     = fixed-effect beta maps
%   Results.P_Values = TFCE-corrected p-values
%   Results.fixedNames

% Defaults
E_H = [0.66 2];
nPerm = 1000;
flag_ft = 0;
par = 0;
ChN = [];

randomSlopes = [];
fitMethod = 'REML';
uncorrelated = true;

fixedNamesIn = {};

% Required checks
if nargin < 3
    error('hnf_ept_TFCE_lme:BadArgs', 'Must supply data1, data2, and subID.');
end

if isempty(data1) == 1 || isempty(data2) == 1 || isempty(subID) == 1
    error('hnf_ept_TFCE_lme:BadArgs', 'data1, data2, and subID must be non-empty.');
end

if isnumeric(data1) == 0 || isnumeric(data2) == 0 || isnumeric(subID) == 0
    error('hnf_ept_TFCE_lme:BadArgs', 'data1, data2, and subID must be numeric.');
end

if exist('OCTAVE_VERSION','builtin')
    if exist('fitlmematrix','file') == 2
        backend = 'Octave statistics fitlmematrix';
    else
        backend = 'R lme4 lmer';
    end
elseif exist('fitlme', 'file') == 2
    backend = 'MATLAB fitlme';
else
    error('hnf_ept_TFCE_lme:NoFitLME', 'fitlme not found. Requires Statistics and Machine Learning Toolbox.');
end

subID = subID(:);
if ~isreal(data1) || ~isreal(data2) || any(~isfinite(data1(:))) || any(~isfinite(data2(:))) || any(~isfinite(subID))
    error('NeuroFreq:InvalidLME','Supply finite real measurements, predictors and subject identifiers.');
end

if size(data1, 1) ~= numel(subID)
    error('hnf_ept_TFCE_lme:ObsMismatch', 'data1 must be observations-first and match length(subID).');
end

if size(data2, 1) ~= numel(subID)
    error('hnf_ept_TFCE_lme:ObsMismatch', 'data2 rows must match length(subID).');
end

% Parse name-value args
if isempty(varargin) == 0

    if mod(numel(varargin), 2) ~= 0
        error('hnf_ept_TFCE_lme:BadArgs', 'Name-value inputs must come in pairs.');
    end

    for i = 1:2:numel(varargin)

        Param = varargin{i};
        Value = varargin{i + 1};

        if isstring(Param) == 1
            Param = char(Param);
        end

        if ischar(Param) == 0
            error('hnf_ept_TFCE_lme:BadArgs', 'Flag arguments must be strings.');
        end

        Param = lower(char(Param));

        if strcmp(Param, 'e_h') == 1
            E_H = Value;
        elseif strcmp(Param, 'nperm') == 1
            nPerm = Value;
        elseif strcmp(Param, 'flag_ft') == 1
            flag_ft = Value;
        elseif strcmp(Param, 'chn') == 1
            ChN = Value;
        elseif strcmp(Param, 'par') == 1
            par = Value;
        elseif strcmp(Param, 'randomslopes') == 1
            randomSlopes = Value;
        elseif strcmp(Param, 'fitmethod') == 1
            fitMethod = Value;
        elseif strcmp(Param, 'uncorrelated') == 1
            uncorrelated = Value;
        elseif strcmp(Param, 'fixednames') == 1
            fixedNamesIn = Value;
        else
            error('hnf_ept_TFCE_lme:BadArgs', ['Unknown parameter setting: ' Param]);
        end

    end

end

% Chanloc logic
if nargin < 4 || isempty(chanlocs)
    flag_ft = 1;
end

% Calculate channel neighbours if needed
if ~isscalar(nPerm) || ~isfinite(nPerm) || nPerm < 1 || fix(nPerm) ~= nPerm
    error('NeuroFreq:InvalidPermutations','nperm must be a positive integer.');
end
if flag_ft == 0 && isempty(ChN) == 1
    if numel(chanlocs) == 1, ChN = 1; else, ChN = ept_ChN2(chanlocs); end
end

DataX = double(data1);
Regs = double(data2);

nObs = size(DataX, 1);

siD = size(DataX);
siTail = siD(2:end);

Dmat = reshape(DataX, [nObs, prod(siTail)]);

[Regs2, fixedNames, permInfo] = local_prepare_design(Regs, subID, fixedNamesIn);

% Normalize extreme physical units before optimization. Positive column
% scaling preserves both full and diagonal random-effects structures.
[Regs2,xScale] = hnf_stats_scale(Regs2);
[Dmat,yScale] = hnf_stats_scale(Dmat);

% Fit observed beta + t maps
disp('fitting observed mixed-effects models...')
[B_Obs, T_Obs] = local_lme_stats_mass( ...
    Dmat, ...
    Regs2, ...
    subID, ...
    fixedNames, ...
    randomSlopes, ...
    fitMethod, ...
    uncorrelated, ...
    par, backend );

B_Obs = hnf_stats_restore_beta(B_Obs,yScale,xScale);
B_Obs = reshape(B_Obs, [size(B_Obs, 1), siTail]);
T_Obs = reshape(T_Obs, [size(T_Obs, 1), siTail]);

TFCE_Obs = zeros(size(T_Obs));

for p = 1:size(T_Obs, 1)
    TFCE_Obs(p, :, :, :) = hnf_tfce_transformation(squeeze(T_Obs(p, :, :, :)), ChN, E_H, flag_ft);
end

% Permutation loop: store max abs TFCE per regressor
maxTFCE = zeros(nPerm, size(Regs2, 2));

prog = 0;
fprintf(1, 'calculating permutations: %3d%%', prog);

for iPerm = 1:nPerm

    RegsP = local_permute_design_subjectaware(Regs2, subID, permInfo);

    [~, T_Perm] = local_lme_stats_mass( ...
        Dmat, ...
        RegsP, ...
        subID, ...
        fixedNames, ...
        randomSlopes, ...
        fitMethod, ...
        uncorrelated, ...
        par, backend );

    T_Perm = reshape(T_Perm, [size(T_Perm, 1), siTail]);

    for p = 1:size(T_Perm, 1)
        TFCE_Perm = hnf_tfce_transformation(squeeze(T_Perm(p, :, :, :)), ChN, E_H, flag_ft);
        maxTFCE(iPerm, p) = max(abs(TFCE_Perm(:)));
    end

    prog = 100 * (iPerm / nPerm);
    fprintf(1, '\b\b\b\b%3.0f%%', prog);

end

fprintf(1, '\n');

% Corrected p-values (FWER): compare voxelwise TFCE to perm max distribution
P_Values = zeros(size(TFCE_Obs));

for p = 1:size(TFCE_Obs, 1)

    if permInfo.randomizable(p)
        P_Values(p,:) = hnf_tfce_pvalues(TFCE_Obs(p,:), maxTFCE(:,p));
    else
        P_Values(p,:) = NaN;
        warning('NeuroFreq:FixedPermutationColumn', ...
            'No allowed subject reassignment changes %s; its randomization p-value is undefined.',fixedNames{p});
    end

end

P_Values = reshape(P_Values, size(TFCE_Obs));

% Output
Results.fitBackend = backend;
Results.Obs = T_Obs;
Results.Beta = B_Obs;
Results.TFCE_Obs = TFCE_Obs;
Results.maxTFCE = sort(maxTFCE, 1);
Results.P_Values = P_Values;
Results.fixedNames = fixedNames;
Results.nullHypothesis = 'Joint null of all randomized fixed effects; one spatial family per coefficient.';
Results.type = local_infer_type_from_shape(flag_ft, siTail);

end

% ====================================================================== %
% Local helpers
% ====================================================================== %

function [Regs2, fixedNames, permInfo] = local_prepare_design(Regs, subID, fixedNamesIn)

if isempty(fixedNamesIn) == 0
    if iscell(fixedNamesIn) == 0
        error('hnf_ept_TFCE_lme:BadFixedNames', 'fixedNames must be a cell array of strings.');
    end
    if numel(fixedNamesIn) ~= size(Regs, 2)
        error('hnf_ept_TFCE_lme:BadFixedNames', 'fixedNames length must match #cols in Regs.');
    end
end

keepCols = true(1, size(Regs, 2));

for j = 1:size(Regs, 2)

    col = Regs(:, j);

    if all(col == col(1))
        if col(1) == 1
            keepCols(j) = false;
        end
    end

end

Regs2 = Regs(:, keepCols);

if isempty(Regs2) == 1
    error('hnf_ept_TFCE_lme:BadDesign', 'Design matrix contains only an intercept-like column; no predictors remain.');
end

if isempty(fixedNamesIn) == 0
    fixedNames = fixedNamesIn(keepCols);
else
    fixedNames = cell(1, size(Regs2, 2));
    for j = 1:size(Regs2, 2)
        fixedNames{j} = ['X' num2str(j)];
    end
end

uSubs = unique(subID);
nSub = numel(uSubs);

betweenSub = false(1, size(Regs2, 2));

for j = 1:size(Regs2, 2)

    isBetween = true;

    for s = 1:nSub

        sid = uSubs(s);
        idx = subID == sid;

        v = Regs2(idx, j);

        if numel(v) > 1
            if any(v ~= v(1))
                isBetween = false;
                break
            end
        end

    end

    betweenSub(j) = isBetween;

end

permInfo = struct();
permInfo.betweenSub = betweenSub;
permInfo.randomizable = true(size(betweenSub));
% Mixed designs exchange whole trial profiles within equal-size blocks.
% If a between-subject column is constant in each such block, no allowed
% randomization tests that effect (for example, all trial counts unique).
if any(betweenSub) && ~all(betweenSub)
    [~,first,block] = unique(subID);
    counts = accumarray(block,1);
    for j = find(betweenSub)
        movable = false;
        for count = unique(counts(:))'
            values = Regs2(first(counts==count),j);
            movable = movable || any(values ~= values(1));
        end
        permInfo.randomizable(j) = movable;
    end
end

end

function RegsP = local_permute_design_subjectaware(Regs, subID, permInfo)

% Preserve design rows, including computed interactions. All between-subject
% values move together. In a mixed design, move whole predictor profiles
% between subjects with equal trial counts, then shuffle trial rows jointly.
% Interpretation requires exchangeability at each level actually randomized.
ids = unique(subID);
blocks = arrayfun(@(id) find(subID == id), ids, 'UniformOutput', false);
counts = cellfun(@numel,blocks);
between = permInfo.betweenSub;
RegsP = Regs;
if all(between)
    source = randperm(numel(ids));
    for s = 1:numel(ids)
        RegsP(blocks{s},:) = repmat(Regs(blocks{source(s)}(1),:), counts(s),1);
    end
else
    source = 1:numel(ids);
    if any(between)
        for count = unique(counts(:))'
            group = find(counts == count);
            source(group) = group(randperm(numel(group)));
        end
    end
    for s = 1:numel(ids)
        donor = blocks{source(s)};
        RegsP(blocks{s},:) = Regs(donor(randperm(numel(donor))),:);
    end
end
end

function [B, T] = local_lme_stats_mass(Dmat, Regs, subID, fixedNames, randomSlopes, fitMethod, uncorrelated, par, backend)

if ~strcmp(backend,'MATLAB fitlme')
    [B,T] = hnf_lme_octave(Dmat,Regs,subID,randomSlopes,fitMethod,uncorrelated,backend);
    return
end

nVox = size(Dmat, 2);
nFix = size(Regs, 2);

B = nan(nFix, nVox);
T = nan(nFix, nVox);

if par == 1

    parfor v = 1:nVox
        y = Dmat(:, v);
        [bNow, tNow] = local_lme_stats_single(y, Regs, subID, fixedNames, randomSlopes, fitMethod, uncorrelated);
        B(:, v) = bNow;
        T(:, v) = tNow;
    end

else

    for v = 1:nVox
        y = Dmat(:, v);
        [bNow, tNow] = local_lme_stats_single(y, Regs, subID, fixedNames, randomSlopes, fitMethod, uncorrelated);
        B(:, v) = bNow;
        T(:, v) = tNow;
    end

end

end

function [bVec, tVec] = local_lme_stats_single(y, Regs, subID, fixedNames, randomSlopes, fitMethod, uncorrelated)

nFix = size(Regs, 2);

bVec = nan(nFix, 1);
tVec = nan(nFix, 1);

% Internal formula names keep arbitrary table labels (including 'y') from
% replacing the response, the group identifier, or MATLAB formula syntax.
tbl = table(y,categorical(subID),'VariableNames',{'response','groupID'});
terms = arrayfun(@(j) sprintf('predictor%d',j),1:nFix,'UniformOutput',false);
for j = 1:nFix, tbl.(terms{j}) = Regs(:,j); end
random = [{'1'}, terms(randomSlopes)];
formula = ['response ~ 1 + ' strjoin(terms,' + ') ' + (' strjoin(random,' + ') '|groupID)'];
pattern = 'FullCholesky';
if uncorrelated, pattern = 'Diagonal'; end
lme = fitlme(tbl,formula,'FitMethod',fitMethod,'CovariancePattern',pattern);

coefTbl = lme.Coefficients;

for j = 1:nFix

    rowName = terms{j};

    hit = find(strcmp(coefTbl.Name, rowName), 1);

    if isempty(hit) == 1
        bVec(j) = NaN;
        tVec(j) = NaN;
    else
        bVec(j) = coefTbl.Estimate(hit);
        tVec(j) = coefTbl.tStat(hit);
    end

end

bVec = bVec(:);
tVec = tVec(:);

end

function type = local_infer_type_from_shape(flag_ft, siTail)

if flag_ft == 1
    if numel(siTail) == 2
        type = 'glm_tf';
    elseif isscalar(siTail)
        type = 'glm_ts';
    else
        type = 'glm_tf';
    end
else
    if numel(siTail) == 3
        type = 'glm_chXfsXts';
    elseif numel(siTail) == 2
        type = 'glm_chXts';
    else
        type = 'glm_chXts';
    end
end

end








