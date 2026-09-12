function reg_out = nf_stregress(TF, regs, varargin)
% NF_STREGRESS    Single-trial regression/correlation on NeuroFreq TF sets
%
% DESIGN
% ------
% 1) INPUT MUST BE an NF TF structure (NOT EEG sets).
% 2) Applies selected POWER method to EVERY numeric top-level field named like power*
%    (power or pow_ prefix), i.e., top-level sibling fields.
% 3) Applies Cohen-style within-channel phase modulation to EVERY numeric top-level field
%    named like phase* (phase or ph_ prefix), i.e., top-level sibling fields.
%
% Regressors (regs) may be:
%   - numeric or logical matrix (trials x K) or (K x trials)
%   - table with K numeric/logical columns; row count must match trials
%     (table variable names are stored in reg_out.prednames)

% ----------------------------
% Parse options
% ----------------------------
p = inputParser;

expectedMethods = {'robust' 'pearson' 'spearman' 'kendall' 'tfce'};

addRequired(p, 'TF', @(x) isstruct(x) == 1);
addRequired(p, 'regs');
addParameter(p, 'toi', [], @(x) isempty(x) || (isnumeric(x) && isvector(x) && numel(x) == 2));
addParameter(p, 'foi', [], @(x) isempty(x) || (isnumeric(x) && isvector(x) && numel(x) == 2));
addParameter(p, 'method', 'robust', @(x) any(strcmpi(char(x), expectedMethods)));
addParameter(p, 'permutations', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'perm_seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'tfce_par', 0, @(x) isscalar(x) && ismember(x, [0 1]));
addParameter(p, 'progress', 1, @(x) isscalar(x) && ismember(x, [0 1]));
addParameter(p, 'progress_detail', 1, @(x) isscalar(x) && ismember(x, [0 1 2]));

parse(p, TF, regs, varargin{:});

opt = p.Results;

method = lower(char(opt.method));
nPerm = round(opt.permutations);

if strcmpi(method, 'tfce') == 1
    if nPerm < 1
        error('nf_stregress:BadPerm', 'TFCE requires permutations >= 1.');
    end
end

if isempty(opt.perm_seed) == 0
    rng(opt.perm_seed);
end

% ----------------------------
% Progress callback
% ----------------------------
progress_cb = [];
progress_detail = opt.progress_detail;

if opt.progress == 1
    if exist('nf_progress', 'file') == 2
        progress_cb = @(s) nf_progress(s);
    else
        progress_cb = @(s) local_progress(s);
    end
end

% ----------------------------
% Require minimal NF TF header
% ----------------------------
if isfield(TF, 'times') == 0
    error('nf_stregress:BadTF', 'TF.times is required.');
end

if isfield(TF, 'freqs') == 0
    error('nf_stregress:BadTF', 'TF.freqs is required.');
end

times_all = TF.times(:);
freqs_all = TF.freqs(:);

if isempty(times_all) == 1 || isnumeric(times_all) == 0
    error('nf_stregress:BadTF', 'TF.times must be numeric and non-empty.');
end

if isempty(freqs_all) == 1 || isnumeric(freqs_all) == 0
    error('nf_stregress:BadTF', 'TF.freqs must be numeric and non-empty.');
end

% ----------------------------
% Find top-level power/phase fields (NF convention)
% ----------------------------
powerFields = nf_find_top_numeric_fields_local(TF, 'power');
powerFields = nf_force_first_local(powerFields, 'power');

if isempty(powerFields) == 1
    error('nf_stregress:NoPower', 'No numeric top-level power fields found.');
end

phaseFields = nf_find_top_numeric_fields_local(TF, 'phase');
phaseFields = nf_force_first_local(phaseFields, 'phase');

% ----------------------------
% Determine trials/chans from first power field
% ----------------------------
firstPow = powerFields{1};
Xfirst = TF.(firstPow);

[~, C, F, T, R] = nf_tf_field_to_4d_local(Xfirst);

if R < 3
    error('nf_stregress:TooFewTrials', 'Need at least 3 trials for single-trial analysis.');
end

if F ~= numel(freqs_all)
    error('nf_stregress:BadLayout', 'Power field %s freq dim does not match numel(TF.freqs).', firstPow);
end

if T ~= numel(times_all)
    error('nf_stregress:BadLayout', 'Power field %s time dim does not match numel(TF.times).', firstPow);
end

% ----------------------------
% toi/foi indices
% ----------------------------
if isempty(opt.toi) == 1
    tIdx = 1:numel(times_all);
else
    if opt.toi(1) > opt.toi(2)
        error('nf_stregress:BadTOI', 'toi must be [tmin tmax] with tmin <= tmax');
    end

    tIdx = find((times_all >= opt.toi(1)) & (times_all <= opt.toi(2)));
end

if isempty(tIdx) == 1
    error('nf_stregress:BadTOI', 'No times in requested toi');
end

if isempty(opt.foi) == 1
    fIdx = 1:numel(freqs_all);
else
    if opt.foi(1) > opt.foi(2)
        error('nf_stregress:BadFOI', 'foi must be [fmin fmax] with fmin <= fmax');
    end

    fIdx = find((freqs_all >= opt.foi(1)) & (freqs_all <= opt.foi(2)));
end

if isempty(fIdx) == 1
    error('nf_stregress:BadFOI', 'No freqs in requested foi');
end

times = times_all(tIdx);
freqs = freqs_all(fIdx);

% ----------------------------
% Regressor parsing (numeric/logical or table) + names
% ----------------------------

if istable(regs) == 1

    if height(regs) ~= R
        error('nf_stregress:BadRegs', 'Table regs must have height == number of trials (%d).', R);
    end

    prednames = regs.Properties.VariableNames(:);

    Xtmp = table2array(regs);

    if isnumeric(Xtmp) == 0 && islogical(Xtmp) == 0
        error('nf_stregress:BadRegs', 'Table regs must contain numeric/logical columns only.');
    end

    regMat = double(Xtmp);

elseif (isnumeric(regs) || islogical(regs)) && ismatrix(regs)

    [nr, nc] = size(regs);

    if nr == R
        regMat = double(regs);
    elseif nc == R
        regMat = double(regs');
    else
        error('nf_stregress:BadRegs', 'regs must have one dimension equal to number of trials (R=%d).', R);
    end

    Kguess = size(regMat, 2);

    prednames = cell(Kguess, 1);

    for k = 1:Kguess
        prednames{k} = ['reg' num2str(k)];
    end

else

    error('nf_stregress:BadRegs', 'regs must be a numeric/logical matrix or a table.');

end

if isempty(regMat) == 1
    error('nf_stregress:BadRegs', 'No regressors found.');
end

if size(regMat, 1) ~= R
    error('nf_stregress:BadRegs', 'Internal regressor matrix must be trials x K.');
end

goodReg = all(isfinite(regMat), 2);

if sum(goodReg) < 3
    error('nf_stregress:BadRegs', 'Too few valid trials after removing NaN/Inf in regressors');
end

reg_used = regMat(goodReg, :);
R_used = size(reg_used, 1);

if size(reg_used, 1) ~= R_used
    error('nf_stregress:BadRegs', 'Internal regressor parsing failure.');
end

% ----------------------------
% chanlocs for TFCE if needed
% ----------------------------
chanlocs = [];

if isfield(TF, 'chanlocs') == 1
    chanlocs = TF.chanlocs;
end

if strcmpi(method, 'tfce') == 1
    if isempty(chanlocs) == 1
        if C > 1
            error('nf_stregress:TFCEChanlocs', 'TFCE GLM requires TF.chanlocs when multiple channels are present.');
        end
    end
end

% ----------------------------
% Progress header
% ----------------------------
if isempty(progress_cb) == 0
    msg = sprintf('nf_stregress: method=%s, perms=%d, trials=%d (used=%d), powerFields=%d, phaseFields=%d', ...
        method, nPerm, R, R_used, numel(powerFields), numel(phaseFields));
    progress_cb(msg);
end

% ----------------------------
% Output header
% ----------------------------
reg_out = struct();

reg_out.times = times;
reg_out.freqs = freqs;
reg_out.scale = 'stat';

if isfield(TF, 'Fs') == 1
    reg_out.Fs = TF.Fs;
end

if isempty(chanlocs) == 0
    reg_out.chanlocs = chanlocs;
end

reg_out.method = method;
reg_out.permutations = nPerm;
reg_out.perm_seed = opt.perm_seed;
reg_out.tfce_par = opt.tfce_par;

reg_out.ntrls_in = R;
reg_out.ntrls_used = R_used;
reg_out.good_trials = goodReg;

reg_out.regressors = reg_used;
reg_out.prednames = prednames;

reg_out.power_fields = powerFields;
reg_out.phase_fields = phaseFields;

% ----------------------------
% POWER: run method on every power field (top-level)
% ----------------------------
for i = 1:numel(powerFields)

    fname = powerFields{i};
    X = TF.(fname);

    if isnumeric(X) == 0
        continue
    end

    if isempty(progress_cb) == 0
        msg = sprintf('POWER %d/%d: %s', i, numel(powerFields), fname);
        progress_cb(msg);
    end

    [X4, ~, Fx, Tx, Rx] = nf_tf_field_to_4d_local(X);

    if Rx ~= R
        error('nf_stregress:TrialMismatch', 'Power field %s has %d trials; expected %d', fname, Rx, R);
    end

    if Fx ~= numel(freqs_all)
        error('nf_stregress:FreqMismatch', 'Power field %s does not match TF.freqs', fname);
    end

    if Tx ~= numel(times_all)
        error('nf_stregress:TimeMismatch', 'Power field %s does not match TF.times', fname);
    end

    X4 = X4(:, fIdx, tIdx, :);
    X4 = X4(:, :, :, goodReg);

    if size(X4, 4) ~= R_used
        error('nf_stregress:TrialMismatch', 'After filtering, power field %s has %d trials; expected %d', fname, size(X4, 4), R_used);
    end

    if any(~isfinite(X4(:)))
        error('nf_stregress:NonFinite', 'Non-finite values found in power field %s. Clean data before regression.', fname);
    end

    if strcmpi(method, 'tfce') == 1
        res = nf_regress_tf_map_tfce_local(X4, reg_used, chanlocs, freqs, times, nPerm, opt.tfce_par, progress_cb);
    else
        res = nf_regress_tf_map_local(X4, reg_used, method, nPerm, progress_cb, progress_detail);
    end

    outName = matlab.lang.makeValidName(fname);
    reg_out.(outName) = res;

    if strcmpi(fname, 'power') == 1
        reg_out.corrcoef = res.corrcoef;
        reg_out.tstat = res.tstat;
        reg_out.p_vals = res.p_vals;
        reg_out.df = res.df;
        reg_out.type = res.type;
        reg_out.primary_measure = 'power';
    end

    if isempty(progress_cb) == 0
        msg = sprintf('POWER %d/%d: %s (done)', i, numel(powerFields), fname);
        progress_cb(msg);
    end

end

% ----------------------------
% PHASE: Cohen modulation for every phase field (top-level)
% ----------------------------
for i = 1:numel(phaseFields)

    pfname = phaseFields{i};
    Phi = TF.(pfname);

    if isnumeric(Phi) == 0
        continue
    end

    if isempty(progress_cb) == 0
        msg = sprintf('PHASE %d/%d: %s', i, numel(phaseFields), pfname);
        progress_cb(msg);
    end

    [P4, ~, Fx, Tx, Rx] = nf_tf_field_to_4d_local(Phi);

    if Rx ~= R
        error('nf_stregress:TrialMismatch', 'Phase field %s has %d trials; expected %d', pfname, Rx, R);
    end

    if Fx ~= numel(freqs_all)
        error('nf_stregress:FreqMismatch', 'Phase field %s does not match TF.freqs', pfname);
    end

    if Tx ~= numel(times_all)
        error('nf_stregress:TimeMismatch', 'Phase field %s does not match TF.times', pfname);
    end

    P4 = P4(:, fIdx, tIdx, :);
    P4 = P4(:, :, :, goodReg);

    if size(P4, 4) ~= R_used
        error('nf_stregress:TrialMismatch', 'After filtering, phase field %s has %d trials; expected %d', pfname, size(P4, 4), R_used);
    end

    if any(~isfinite(P4(:)))
        error('nf_stregress:NonFinite', 'Non-finite values found in phase field %s. Clean data before phase analysis.', pfname);
    end

    pres = nf_phase_cohen_tf_local(P4, reg_used, nPerm, progress_cb, progress_detail, pfname);

    outName = matlab.lang.makeValidName(pfname);
    reg_out.(outName) = pres;

    if isempty(progress_cb) == 0
        msg = sprintf('PHASE %d/%d: %s (done)', i, numel(phaseFields), pfname);
        progress_cb(msg);
    end

end

if isempty(progress_cb) == 0
    progress_cb('nf_stregress: done');
end

disp('[nf_stregress]: Single-trial regression EEG analysis complete.');

end

% =========================================================
% Local helpers
% =========================================================

function fields = nf_find_top_numeric_fields_local(S, which)

fn = fieldnames(S);
fields = {};

wantPower = false;
wantPhase = false;

if strcmpi(which, 'power') == 1
    wantPower = true;
end

if strcmpi(which, 'phase') == 1
    wantPhase = true;
end

for i = 1:numel(fn)

    name = fn{i};

    if strcmpi(name, 'times') == 1
        continue
    end

    if strcmpi(name, 'freqs') == 1
        continue
    end

    if strcmpi(name, 'chanlocs') == 1
        continue
    end

    if strcmpi(name, 'conds') == 1
        continue
    end

    if strcmpi(name, 'scale') == 1
        continue
    end

    if strcmpi(name, 'aggtype') == 1
        continue
    end

    v = S.(name);

    if isnumeric(v) == 0
        continue
    end

    isMatch = false;

    if wantPower == true
        if strncmpi(name, 'power', 5) == 1
            isMatch = true;
        end
        if strncmpi(name, 'pow_', 4) == 1
            isMatch = true;
        end
    elseif wantPhase == true
        if strncmpi(name, 'phase', 5) == 1
            isMatch = true;
        end
        if strncmpi(name, 'ph_', 3) == 1
            isMatch = true;
        end
    end

    if isMatch == false
        continue
    end

    fields{end + 1, 1} = name; %#ok<AGROW>

end

end

function fields = nf_force_first_local(fields, firstName)

idx = find(strcmpi(fields, firstName), 1);

if isempty(idx) == 0
    fields(idx) = [];
    fields = [{firstName}; fields];
end

end

function [X4, C, F, T, R] = nf_tf_field_to_4d_local(X)

nd = ndims(X);

if nd == 4

    X4 = X;

    sz = size(X4);
    C = sz(1);
    F = sz(2);
    T = sz(3);
    R = sz(4);

elseif nd == 3

    sz = size(X);

    F = sz(1);
    T = sz(2);
    R = sz(3);

    C = 1;

    X4 = reshape(X, [1 F T R]);

else

    error('nf_stregress:BadShape', 'Unsupported TF field dimensionality: ndims=%d', nd);

end

end

function local_progress(msg)

ts = datestr(now, 'HH:MM:SS'); %#ok
fprintf('[%s] %s\n', ts, msg);

end

% =========================================================
% Regression backends
% =========================================================

function res = nf_regress_tf_map_tfce_local(Y, reg_used, chanlocs, freqs, times, nPerm, tfce_par, progress_cb)

[~, ~, ~, R] = size(Y);

dataX = permute(Y, [4 1 2 3]);  % R x C x F x T

if isempty(progress_cb) == 0
    progress_cb('TFCE: running nf_glm_TFCE...');
end

tfce = nf_glm_TFCE(dataX, reg_used, chanlocs, freqs, times, nPerm, tfce_par);

if isfield(tfce, 'z_vals') == 1
    zmap = tfce.z_vals;
else
    zmap = nf_tfce_p_to_z_local(tfce.p_vals, tfce.tstat, nPerm);
end

df = R - (size(reg_used, 2) + 1);

res = struct();
res.corrcoef = zmap;
res.tstat = zmap;
res.p_vals = tfce.p_vals;
res.df = df;

res.type = 'glm_chXfsXts';
res.method = 'tfce';

res.beta = tfce.beta;
res.tfce_tstat = tfce.tstat;

if isfield(tfce, 'fixedNames') == 1
    res.fixedNames = tfce.fixedNames;
end

end

function z = nf_tfce_p_to_z_local(p, tstat, nPerm)

pmin = 1 ./ (nPerm + 1);

if pmin < eps %#ok
    pmin = eps;
end

pmax = 1 - eps;

p = max(p, pmin);
p = min(p, pmax);

sgn = sign(tstat);

z = sgn .* norminv(1 - p ./ 2);

end

function res = nf_regress_tf_map_local(Y, reg_used, method, nPerm, progress_cb, progress_detail)

[C, F, T, R] = size(Y);
Kreg = size(reg_used, 2);

P = C * F * T;

Y2 = permute(Y, [4 1 2 3]);
Y2 = reshape(Y2, [R P]);

obs = nan(P, Kreg);

if isempty(progress_cb) == 0
    if strcmpi(method, 'robust') == 1
        progress_cb('ROBUST: computing observed stats...');
    else
        progress_cb('CORR: computing observed stats...');
    end
end

if strcmpi(method, 'robust') == 1

    XZ0 = nf_zscore_mat_local(reg_used);

    pool = gcp('nocreate');

    if isempty(pool) == 1
        usePar = 0;
    else
        usePar = 1;
    end

    if usePar == 1

        if isa(pool, 'parallel.ThreadPool') == 0
            try
                spmd
                    warning('off','all');
                    dbclear if error
                    dbclear if warning
                end
            catch
            end
        end

        parfor p = 1:P

            y = Y2(:, p);
            y = y(:);

            yZ = nf_zscore_vec_local(y);

            try
                b = robustfit(XZ0, yZ);
                b = b(:);
                if numel(b) ~= (Kreg + 1)
                    error('BadB');
                end
                obs(p, :) = b(2:(Kreg + 1)).';
            catch
            end

        end

    else

        warning off;

        for p = 1:P

            y = Y2(:, p);
            y = y(:);

            yZ = nf_zscore_vec_local(y);

            try
                b = robustfit(XZ0, yZ);
                b = b(:);
                if numel(b) ~= (Kreg + 1)
                    error('BadB');
                end
                obs(p, :) = b(2:(Kreg + 1)).';
            catch
            end

        end

        warning on;

    end

else

    try
        obs = corr(Y2, reg_used, 'type', method);
    catch
        error('nf_stregress:corrFail', 'corr failed for method %s', method);
    end

end

if nPerm > 0

    zObs = nan(size(obs));

    if strcmpi(method, 'robust') == 1

        if isempty(progress_cb) == 0
            if progress_detail >= 1
                progress_cb('ROBUST perms: computing null via permutations (progress is coarse).');
            end
        end

        pool = gcp('nocreate');

        if isempty(pool) == 1
            usePar = 0;
        else
            usePar = 1;
        end

        if usePar == 1

            if isa(pool, 'parallel.ThreadPool') == 0
                try
                    spmd
                        warning('off','all');
                        dbclear if error
                        dbclear if warning
                    end
                catch
                end
            end

            parfor p = 1:P

                obsNow = obs(p, :);

                if any(isnan(obsNow))
                    continue
                end

                y = Y2(:, p);
                y = y(:);

                yZ = nf_zscore_vec_local(y);

                mu = zeros(1, Kreg);
                M2 = zeros(1, Kreg);
                nOk = 0;

                for pp = 1:nPerm

                    permIdx = randperm(R);
                    Xp = reg_used(permIdx, :); %#ok

                    xNow = nan(1, Kreg);

                    try
                        XZp = nf_zscore_mat_local(Xp);
                        b0 = robustfit(XZp, yZ);
                        b0 = b0(:);
                        if numel(b0) ~= (Kreg + 1)
                            error('BadB');
                        end
                        xNow = b0(2:(Kreg + 1)).';
                    catch
                    end

                    if any(isnan(xNow))
                        continue
                    end

                    nOk = nOk + 1;

                    if nOk == 1
                        mu = xNow;
                        M2 = zeros(1, Kreg);
                    else
                        d = xNow - mu;
                        mu = mu + d ./ nOk;
                        M2 = M2 + d .* (xNow - mu);
                    end

                end

                if nOk > 1
                    sd = sqrt(M2 ./ (nOk - 1));
                else
                    sd = zeros(1, Kreg);
                end

                sd(sd == 0) = eps;

                zObs(p, :) = (obsNow - mu) ./ sd;

            end

        else

            warning('off','all');

            for p = 1:P

                obsNow = obs(p, :);

                if any(isnan(obsNow))
                    continue
                end

                y = Y2(:, p);
                y = y(:);

                yZ = nf_zscore_vec_local(y);

                mu = zeros(1, Kreg);
                M2 = zeros(1, Kreg);
                nOk = 0;

                for pp = 1:nPerm

                    permIdx = randperm(R);
                    Xp = reg_used(permIdx, :);

                    xNow = nan(1, Kreg);

                    try
                        XZp = nf_zscore_mat_local(Xp);
                        b0 = robustfit(XZp, yZ);
                        b0 = b0(:);
                        if numel(b0) ~= (Kreg + 1)
                            error('BadB');
                        end
                        xNow = b0(2:(Kreg + 1)).';
                    catch
                    end

                    if any(isnan(xNow))
                        continue
                    end

                    nOk = nOk + 1;

                    if nOk == 1
                        mu = xNow;
                        M2 = zeros(1, Kreg);
                    else
                        d = xNow - mu;
                        mu = mu + d ./ nOk;
                        M2 = M2 + d .* (xNow - mu);
                    end

                end

                if nOk > 1
                    sd = sqrt(M2 ./ (nOk - 1));
                else
                    sd = zeros(1, Kreg);
                end

                sd(sd == 0) = eps;

                zObs(p, :) = (obsNow - mu) ./ sd;

            end

            warning('on','all');

        end

    else

        blockSize = 50000;

        if blockSize > P
            blockSize = P;
        end

        if isempty(progress_cb) == 0
            if progress_detail >= 1
                nBlocks = ceil(P ./ blockSize);
                msg = sprintf('CORR perms: %d blocks (blockSize=%d)', nBlocks, blockSize);
                progress_cb(msg);
            end
        end

        for b0 = 1:blockSize:P

            b1 = b0 + blockSize - 1;

            if b1 > P
                b1 = P;
            end

            idx = b0:b1;

            if isempty(progress_cb) == 0
                if progress_detail >= 2
                    nBlocks = ceil(P ./ blockSize);
                    bIdx = ceil(b0 ./ blockSize);
                    msg = sprintf('CORR perms: block %d/%d', bIdx, nBlocks);
                    progress_cb(msg);
                end
            end

            obsBlock = obs(idx, :);

            mu = zeros(numel(idx), Kreg);
            M2 = zeros(numel(idx), Kreg);

            for pp = 1:nPerm

                permIdx = randperm(R);
                Xp = reg_used(permIdx, :);

                nullBlock = corr(Y2(:, idx), Xp, 'type', method);

                if pp == 1
                    mu = nullBlock;
                    M2 = zeros(size(nullBlock));
                else
                    d = nullBlock - mu;
                    mu = mu + d ./ pp;
                    M2 = M2 + d .* (nullBlock - mu);
                end

            end

            if nPerm > 1
                sd = sqrt(M2 ./ (nPerm - 1));
            else
                sd = zeros(size(M2));
            end

            sd(sd == 0) = eps;

            zObs(idx, :) = (obsBlock - mu) ./ sd;

        end

    end

    obs = zObs;

    pVals = 2 .* (1 - normcdf(abs(obs)));
    tVals = obs;

else

    dfTmp = R - 2;

    if strcmpi(method, 'pearson') == 1
        denom = 1 - obs .^ 2;
        denom(denom <= 0) = eps;
        tVals = obs .* sqrt(dfTmp ./ denom);
        pVals = 2 .* (1 - tcdf(abs(tVals), dfTmp));
    else
        tVals = nan(size(obs));
        pVals = nan(size(obs));
    end

end

corrcoef = reshape(obs, [C F T Kreg]);
tstat = reshape(tVals, [C F T Kreg]);
p_out = reshape(pVals, [C F T Kreg]);

if Kreg == 1
    type = 'corr_chXfsXts';
else
    type = 'glm_chXfsXts';
end

if strcmpi(method, 'robust') == 1
    df = R - (Kreg + 1);
else
    df = R - 2;
end

res = struct();
res.corrcoef = corrcoef;
res.tstat = tstat;
res.p_vals = p_out;
res.df = df;
res.type = type;
res.method = method;

end

function res = nf_phase_cohen_tf_local(Phi, reg_used, nPerm, progress_cb, progress_detail, pfname)

[C, F, T, R] = size(Phi);
Kreg = size(reg_used, 2);

out = zeros(C, F, T, Kreg);
pVals = nan(C, F, T, Kreg);

V = zeros(size(reg_used));

for k = 1:Kreg
    V(:, k) = tiedrank(reg_used(:, k));
end

if isempty(progress_cb) == 0
    if progress_detail >= 1
        msg = sprintf('PHASE cohen: %s (C=%d, K=%d, perms=%d)', pfname, C, Kreg, nPerm);
        progress_cb(msg);
    end
end

for k = 1:Kreg

    v = V(:, k);

    if any(~isfinite(v))
        error('nf_stregress:BadRegs', 'Non-finite values in rank-transformed regressor');
    end

    if isempty(progress_cb) == 0
        if progress_detail >= 2
            msg = sprintf('PHASE cohen: %s reg %d/%d', pfname, k, Kreg);
            progress_cb(msg);
        end
    end

    for ch = 1:C

        if isempty(progress_cb) == 0
            if progress_detail >= 2
                msg = sprintf('PHASE cohen: %s reg %d/%d ch %d/%d', pfname, k, Kreg, ch, C);
                progress_cb(msg);
            end
        end

        phiCh = squeeze(Phi(ch, :, :, :));      % F x T x R
        phiCh = reshape(phiCh, [F T R]);

        E = exp(1i .* phiCh);                   % F x T x R
        E2 = reshape(E, [F * T, R]);            % (F*T) x R
        E2 = E2.';                              % R x (F*T)

        obsVec = abs((v.' * E2) ./ R);          % 1 x (F*T)
        obsVec = obsVec(:).';

        if nPerm > 0

            mu = zeros(1, F * T);
            M2 = zeros(1, F * T);

            for pp = 1:nPerm

                permIdx = randperm(R);
                vp = v(permIdx);

                xVec = abs((vp.' * E2) ./ R);
                xVec = xVec(:).';

                if pp == 1
                    mu = xVec;
                    M2 = zeros(1, F * T);
                else
                    d = xVec - mu;
                    mu = mu + d ./ pp;
                    M2 = M2 + d .* (xVec - mu);
                end

            end

            if nPerm > 1
                sd = sqrt(M2 ./ (nPerm - 1));
            else
                sd = zeros(1, F * T);
            end

            sd(sd == 0) = eps;

            zVec = (obsVec - mu) ./ sd;
            zVec = zVec(:).';

            outVec = zVec;

            pVec = 2 .* (1 - normcdf(abs(zVec)));
            pVec = pVec(:).';

            pVals(ch, :, :, k) = reshape(pVec, [F T]);

        else

            outVec = obsVec;

        end

        out(ch, :, :, k) = reshape(outVec, [F T]);

    end

end

res = struct();

res.corrcoef = out;
res.tstat = out;

if nPerm > 0
    res.p_vals = pVals;
else
    res.p_vals = nan(size(out));
end

res.df = R - 2;

if Kreg == 1
    res.type = 'corr_chXfsXts';
else
    res.type = 'glm_chXfsXts';
end

res.method = 'phase_cohen';

end

function z = nf_zscore_vec_local(x)

x = x(:);

mu = mean(x);
sd = std(x);

if ~(isfinite(sd) && sd > 0)
    z = zeros(size(x));
else
    z = (x - mu) ./ sd;
end

end

function Z = nf_zscore_mat_local(X)

mu = mean(X, 1);
sd = std(X, 0, 1);

Z = zeros(size(X));

for c = 1:size(X, 2)

    if ~(isfinite(sd(c)) && sd(c) > 0)
        Z(:, c) = 0;
    else
        Z(:, c) = (X(:, c) - mu(c)) ./ sd(c);
    end

end

end



