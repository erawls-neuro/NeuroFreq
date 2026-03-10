function Results = ept_TFCE_glm(data1, data2, chanlocs, varargin)
%
% EPT_TFCE_GLM  TFCE-corrected mass-univariate GLM using OLS + permutations.
%
% OUTPUT
% ------
% Results.Beta      : beta maps for each regressor
% Results.Obs       : t-stat maps for each regressor (what TFCE operates on)
% Results.TFCE_Obs  : TFCE-transformed observed t-stat maps
% Results.maxTFCE   : max abs TFCE under null for each perm/regressor
% Results.P_Values  : TFCE-corrected p-values
% Results.df        : degrees of freedom (n - rank(X))
%
% Notes:
% - Uses inv(X'X) via pinv for robustness
% - Permutes rows of X (regressors) which preserves X'X, enabling fast perms
%

% -----------------------------
% Defaults
% -----------------------------
E_H = [0.66 2];
nPerm = 1000;
flag_ft = 0;
ChN = [];

% -----------------------------
% Check inputs
% -----------------------------
if nargin < 2
    error('at least data are required');
end

% -----------------------------
% Parse varargin
% -----------------------------
if ~isempty(varargin)
    if mod(numel(varargin), 2) ~= 0
        error('Name-value arguments must come in pairs.');
    end
end

for iArg = 1:2:numel(varargin)

    Param = varargin{iArg};
    Value = varargin{iArg + 1};

    if ~ischar(Param)
        error('Flag arguments must be strings');
    end

    Param = lower(Param);

    if strcmp(Param, 'e_h')
        E_H = Value;
    elseif strcmp(Param, 'nperm')
        nPerm = Value;
    elseif strcmp(Param, 'flag_ft')
        flag_ft = Value;
    elseif strcmp(Param, 'chn')
        ChN = Value;
    else
        error(['Unknown parameter setting: ' Param]);
    end

end

% -----------------------------
% Set things
% -----------------------------
DataX = double(data1);
Regs = double(data2);
e_loc = chanlocs;

% -----------------------------
% Info / error checking
% -----------------------------
nCh = size(DataX, 2);

if flag_ft == 0
    if ~isequal(nCh, length(e_loc))
        error('Number of channels in data does not equal that of locations file');
    end
end

tic;

% -----------------------------
% Channel neighbours
% -----------------------------
if flag_ft == 0 && isempty(ChN)
    disp('calculating channel neighbours...')
    ChN = ept_ChN2(e_loc);
end

% -----------------------------
% Preallocate
% -----------------------------
nRegs = size(Regs, 2);
maxTFCE = zeros(nPerm, nRegs);

% -----------------------------
% Compute observed GLM statistics
% -----------------------------
disp('calculating observed GLM statistics...')

D = DataX;
siD = size(D);
siD = siD(2:end);

Dmat = D(:, :);

[Beta_Obs, T_Obs, df, invXX, diagC] = local_glm_beta_t(Dmat, Regs);

Beta_Obs = reshape(Beta_Obs, [size(Beta_Obs, 1), siD]);
T_Obs = reshape(T_Obs, [size(T_Obs, 1), siD]);

TFCE_Obs = zeros(size(T_Obs));

for p = 1:size(T_Obs, 1)
    TFCE_Obs(p, :, :, :) = tfce_transformation( ...
        squeeze(T_Obs(p, :, :, :)), ...
        ChN, ...
        E_H, ...
        flag_ft );
end

% -----------------------------
% Permutations
% -----------------------------
prog = 1;
fprintf(1, 'calculating permutations: %3d%%\n', prog);

for iPerm = 1:nPerm

    regs1 = Regs(randperm(size(Dmat, 1)), :);

    [~, T_Perm] = local_glm_beta_t_fast(Dmat, regs1, df, invXX, diagC);

    T_Perm = reshape(T_Perm, [size(T_Perm, 1), siD]);

    for p = 1:size(T_Perm, 1)
        TFCE_Perm = tfce_transformation( ...
            squeeze(T_Perm(p, :, :, :)), ...
            ChN, ...
            E_H, ...
            flag_ft );

        maxTFCE(iPerm, p) = max(abs(TFCE_Perm(:)));
    end

    prog = 100 * (iPerm / nPerm);
    fprintf(1, '\b\b\b\b%3.0f%%', prog);

end

fprintf(1, '\n');

% -----------------------------
% Corrected p-values
% -----------------------------
edges = [maxTFCE; max(abs(TFCE_Obs(:, :)), [], 2)'];

sz = size(TFCE_Obs);
sz(1) = [];

P_Values = zeros(size(TFCE_Obs));

for r = 1:size(TFCE_Obs, 1)

    edgeVec = sort(edges(:, r));

    [~, bin] = histc(abs(TFCE_Obs(r, :)), edgeVec);

    P_Values(r, :) = 1 - bin ./ (nPerm + 2);

end

P_Values = reshape(P_Values, [r, sz]);

% -----------------------------
% Output
% -----------------------------
Results = struct();

Results.Beta = Beta_Obs;
Results.Obs = T_Obs;
Results.TFCE_Obs = TFCE_Obs;
Results.maxTFCE = sort(maxTFCE);
Results.P_Values = P_Values;
Results.df = df;

toc;

[min_P, idx] = min(Results.P_Values(:));
max_Obs = Results.Obs(idx);

display(['Peak significance: T = ' num2str(max_Obs) ', p = ' num2str(min_P)]);

end

% ====================================================================== %
% Local helpers
% ====================================================================== %

function [B, T, df, invXX, diagC] = local_glm_beta_t(Dmat, X)

n = size(Dmat, 1);

if size(X, 1) ~= n
    error('Regressor matrix rows must match data observations.');
end

rankX = rank(X);
df = n - rankX;

if df <= 0
    error('Invalid degrees of freedom in GLM. Check design rank/size.');
end

XtX = X' * X;
invXX = pinv(XtX);

diagC = diag(invXX);
diagC(diagC < 0) = 0;

B = invXX * (X' * Dmat);

R = Dmat - (X * B);

SSE = sum(R .^ 2, 1);

mse = SSE ./ df;

se = sqrt(bsxfun(@times, diagC, mse));

T = B ./ se;

T(~isfinite(T)) = NaN;

end

function [B, T] = local_glm_beta_t_fast(Dmat, Xperm, df, invXX, diagC)

B = invXX * (Xperm' * Dmat);

R = Dmat - (Xperm * B);

SSE = sum(R .^ 2, 1);

mse = SSE ./ df;

se = sqrt(bsxfun(@times, diagC, mse));

T = B ./ se;

T(~isfinite(T)) = NaN;

end











