function Results = ept_TFCE_lme(data1, data2, subID, chanlocs, varargin)
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
    error('ept_TFCE_lme:BadArgs', 'Must supply data1, data2, and subID.');
end

if isempty(data1) == 1 || isempty(data2) == 1 || isempty(subID) == 1
    error('ept_TFCE_lme:BadArgs', 'data1, data2, and subID must be non-empty.');
end

if isnumeric(data1) == 0 || isnumeric(data2) == 0 || isnumeric(subID) == 0
    error('ept_TFCE_lme:BadArgs', 'data1, data2, and subID must be numeric.');
end

if exist('fitlme', 'file') ~= 2
    error('ept_TFCE_lme:NoFitLME', 'fitlme not found. Requires Statistics and Machine Learning Toolbox.');
end

subID = subID(:);

if size(data1, 1) ~= numel(subID)
    error('ept_TFCE_lme:ObsMismatch', 'data1 must be observations-first and match length(subID).');
end

if size(data2, 1) ~= numel(subID)
    error('ept_TFCE_lme:ObsMismatch', 'data2 rows must match length(subID).');
end

% Parse name-value args
if isempty(varargin) == 0

    if mod(numel(varargin), 2) ~= 0
        error('ept_TFCE_lme:BadArgs', 'Name-value inputs must come in pairs.');
    end

    for i = 1:2:numel(varargin)

        Param = varargin{i};
        Value = varargin{i + 1};

        if isstring(Param) == 1
            Param = char(Param);
        end

        if ischar(Param) == 0
            error('ept_TFCE_lme:BadArgs', 'Flag arguments must be strings.');
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
            error('ept_TFCE_lme:BadArgs', ['Unknown parameter setting: ' Param]);
        end

    end

end

% Chanloc logic
if isempty(chanlocs) == 1
    flag_ft = 1;
end

% Calculate channel neighbours if needed
if flag_ft == 0 && isempty(ChN) == 1
    disp('calculating channel neighbours...')
    ChN = ept_ChN2(chanlocs);
end

DataX = double(data1);
Regs = double(data2);

nObs = size(DataX, 1);

siD = size(DataX);
siTail = siD(2:end);

Dmat = reshape(DataX, [nObs, prod(siTail)]);

[Regs2, fixedNames, permInfo] = local_prepare_design(Regs, subID, fixedNamesIn);

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
    par );

B_Obs = reshape(B_Obs, [size(B_Obs, 1), siTail]);
T_Obs = reshape(T_Obs, [size(T_Obs, 1), siTail]);

TFCE_Obs = zeros(size(T_Obs));

for p = 1:size(T_Obs, 1)
    TFCE_Obs(p, :, :, :) = tfce_transformation(squeeze(T_Obs(p, :, :, :)), ChN, E_H, flag_ft);
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
        par );

    T_Perm = reshape(T_Perm, [size(T_Perm, 1), siTail]);

    for p = 1:size(T_Perm, 1)
        TFCE_Perm = tfce_transformation(squeeze(T_Perm(p, :, :, :)), ChN, E_H, flag_ft);
        maxTFCE(iPerm, p) = max(abs(TFCE_Perm(:)));
    end

    prog = 100 * (iPerm / nPerm);
    fprintf(1, '\b\b\b\b%3.0f%%', prog);

end

fprintf(1, '\n');

% Corrected p-values (FWER): compare voxelwise TFCE to perm max distribution
P_Values = zeros(size(TFCE_Obs));

for p = 1:size(TFCE_Obs, 1)

    dist = maxTFCE(:, p);
    dist = dist(:);

    distSorted = sort(dist);

    absObsVec = abs(TFCE_Obs(p, :));
    absObsVec = absObsVec(:);

    edges = [-inf; distSorted];

    [~, bin] = histc(absObsVec, edges);

    nLessEq = bin - 1;

    pVec = (nPerm - nLessEq + 1) ./ (nPerm + 1);

    P_Values(p, :) = reshape(pVec, [1 size(TFCE_Obs, 2) * size(TFCE_Obs, 3) * size(TFCE_Obs, 4)]);

end

P_Values = reshape(P_Values, size(TFCE_Obs));

% Output
Results.Obs = T_Obs;
Results.Beta = B_Obs;
Results.TFCE_Obs = TFCE_Obs;
Results.maxTFCE = sort(maxTFCE, 1);
Results.P_Values = P_Values;
Results.fixedNames = fixedNames;
Results.type = local_infer_type_from_shape(flag_ft, siTail);

end

% ====================================================================== %
% Local helpers
% ====================================================================== %

function [Regs2, fixedNames, permInfo] = local_prepare_design(Regs, subID, fixedNamesIn)

tol = 1e-12;

if isempty(fixedNamesIn) == 0
    if iscell(fixedNamesIn) == 0
        error('ept_TFCE_lme:BadFixedNames', 'fixedNames must be a cell array of strings.');
    end
    if numel(fixedNamesIn) ~= size(Regs, 2)
        error('ept_TFCE_lme:BadFixedNames', 'fixedNames length must match #cols in Regs.');
    end
end

keepCols = true(1, size(Regs, 2));

for j = 1:size(Regs, 2)

    col = Regs(:, j);

    if all(abs(col - col(1)) < tol)
        if abs(col(1) - 1) < tol
            keepCols(j) = false;
        end
    end

end

Regs2 = Regs(:, keepCols);

if isempty(Regs2) == 1
    error('ept_TFCE_lme:BadDesign', 'Design matrix contains only an intercept-like column; no predictors remain.');
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
        idx = find(subID == sid);

        v = Regs2(idx, j);

        if numel(v) > 1
            if std(v) > tol
                isBetween = false;
                break
            end
        end

    end

    betweenSub(j) = isBetween;

end

permInfo = struct();
permInfo.betweenSub = betweenSub;

end

function RegsP = local_permute_design_subjectaware(Regs, subID, permInfo)

RegsP = Regs;

uSubs = unique(subID);
nSub = numel(uSubs);

betweenSub = permInfo.betweenSub;

for j = 1:size(Regs, 2)

    col = Regs(:, j);

    if all(abs(col - col(1)) < 1e-12)
        continue
    end

    if betweenSub(j) == 1

        subjVals = zeros(nSub, 1);

        for s = 1:nSub
            sid = uSubs(s);
            idx = find(subID == sid);
            subjVals(s) = col(idx(1));
        end

        subjVals = subjVals(randperm(nSub));

        for s = 1:nSub
            sid = uSubs(s);
            idx = find(subID == sid);
            RegsP(idx, j) = subjVals(s);
        end

    else

        for s = 1:nSub
            sid = uSubs(s);
            idx = find(subID == sid);
            RegsP(idx, j) = col(idx(randperm(numel(idx))));
        end

    end

end

end

function [B, T] = local_lme_stats_mass(Dmat, Regs, subID, fixedNames, randomSlopes, fitMethod, uncorrelated, par)

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

tbl = table();
tbl.y = y;
tbl.subID = categorical(subID);

for j = 1:nFix
    tbl.(fixedNames{j}) = Regs(:, j);
end

fixedStr = 'y ~ 1';

for j = 1:nFix
    fixedStr = [fixedStr ' + ' fixedNames{j}];
end

randVars = {};

if isempty(randomSlopes) == 0
    for k = 1:numel(randomSlopes)
        idx = randomSlopes(k);
        if idx >= 1 && idx <= nFix
            randVars{end + 1, 1} = fixedNames{idx}; %#ok<AGROW>
        end
    end
end

if isempty(randVars) == 1

    if uncorrelated == 1
        randStr = ' + (1||subID)';
    else
        randStr = ' + (1|subID)';
    end

else

    slopeStr = randVars{1};

    for k = 2:numel(randVars)
        slopeStr = [slopeStr ' + ' randVars{k}];
    end

    if uncorrelated == 1
        randStr = [' + (1 + ' slopeStr '||subID)'];
    else
        randStr = [' + (1 + ' slopeStr '|subID)'];
    end

end

formula = [fixedStr randStr];

try
    lme = fitlme(tbl, formula, 'FitMethod', fitMethod);
catch
    return
end

coefTbl = lme.Coefficients;

for j = 1:nFix

    rowName = fixedNames{j};

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
    elseif numel(siTail) == 1
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








