function TF = nf_tfspecparam(TF, varargin)
% NF_TFSPECPARAM    Computes a time-resolved spectral parameterization
%
% GENERAL
% -------
% Time resolved spectral parameterization (Wilson et al., 2022) rewritten
% for use with TF structures output by nf_tftransform.m. Computes a
% time-resolved spectral parameterization using nf_specparam.m.
%
% If nf_specparam fails (errors) for ANY time slice within a given trial,
% that entire trial is removed from ALL trial-resolved data fields after
% parameterization completes (including TF.behavior if present).
%
% Accepts TF transforms calculated using any method. Automatically detects
% if power values are on a linear scale and converts to log scale.
%
% OUTPUT
% ------
% TF.specparam fields (all log10 scale)
%   Aperiodic:
%     TF.specparam.rsquare     : C×T×R
%     TF.specparam.mae         : C×T×R
%     TF.specparam.ap_offset   : C×T×R
%     TF.specparam.ap_knee     : C×T×R (NaN if 'fixed')
%     TF.specparam.ap_exponent : C×T×R
%   Peaks (ALL peaks per slice, parfor-safe numeric arrays; padded with NaN):
%     TF.specparam.pk_n        : C×T×R    (#peaks detected)
%     TF.specparam.pk_f        : C×K×T×R  (Hz)        K = maxPeaks
%     TF.specparam.pk_aLog10   : C×K×T×R  (log10 power height)
%     TF.specparam.pk_sd       : C×K×T×R  (Hz)
%
% INPUT
% -----
%   TF.power : C×F×T[×R] power
%   TF.scale : 'linear'|'log10'  (if 'linear', converted to log10)
%   TF.freqs : F×1 (Hz)
%
% Name-Value
% ----------
%   'peakWidthLims' : [min max] Hz (default [0.5 12])
%   'maxPeaks'      : max #peaks per slice (default numel(TF.freqs))
%   'minPeakHeight' : minimum peak height (log10 power units if scale->log10) (default 0)
%   'aPeriodicMode' : 'fixed'|'knee' (default 'fixed')
%   'peakThreshold' : noise SD threshold for detection (default 2.0)
%   'peakType'      : 'gaussian'|'cauchy' (default 'gaussian')
%   'threshAfter'   : re-threshold when no optim (default 1)
%   'optim'         : use Optimization Toolbox if available (default 1)
%   'FoD'           : 'fit' returns model fit; 'data' returns modeled data (default 'fit')
%


p = inputParser;

valid2Scalar = @(x) numel(x) == 2 && isnumeric(x) && isvector(x);
expectedPeakType = @(x) any(validatestring(x, {'gaussian', 'cauchy'}));
expectedAPModes = @(x) any(validatestring(x, {'knee', 'fixed'}));
expectedFoD = @(x) any(validatestring(x, {'fit', 'data'}));

addRequired(p, 'TF');

addParameter(p, 'peakWidthLims', [0.5 12], valid2Scalar);
addParameter(p, 'maxPeaks', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && x >= 0));
addParameter(p, 'minPeakHeight', 0, @isscalar);
addParameter(p, 'aPeriodicMode', 'fixed', expectedAPModes);
addParameter(p, 'peakThreshold', 2.0, @isscalar);
addParameter(p, 'peakType', 'gaussian', expectedPeakType);
addParameter(p, 'threshAfter', 1, @(x) isscalar(x) && ismember(x, [0 1]));
addParameter(p, 'optim', 1, @(x) isscalar(x) && ismember(x, [0 1]));
addParameter(p, 'FoD', 'fit', expectedFoD);

parse(p, TF, varargin{:});
opt = p.Results;

if ~isfield(TF, 'freqs')
    error('TF.freqs is required.');
end

if ~isfield(TF, 'scale')
    error('TF.scale is required.');
end

if ~isfield(TF, 'power')
    error('TF.power is required.');
end

if ~isfield(TF, 'nsensor')
    error('TF.nsensor is required.');
end

if ~isfield(TF, 'ntrls')
    error('TF.ntrls is required.');
end

freqs = TF.freqs(:);

if isempty(opt.maxPeaks)
    opt.maxPeaks = numel(freqs);
end

if ~(strcmpi(TF.scale, 'linear') || strcmpi(TF.scale, 'log10'))
    error('TF.scale must be ''linear'' or ''log10''.');
end

powFields = nf_find_power_fields_local(TF);

if strcmpi(TF.scale, 'linear')
    for i = 1:numel(powFields)
        name = powFields{i};
        X = TF.(name);
        X(X <= 0) = eps;
        TF.(name) = log10(X);
    end
    TF.scale = 'log10';
end

powFields = nf_find_power_fields_local(TF);

badTrialsAll = false(1, TF.ntrls);

for i = 1:numel(powFields)

    srcName = powFields{i};

    if nf_is_param_power_field_local(srcName) == 1
        continue
    end

    surf = TF.(srcName);

    [ap_log, osc_log, specparam, badTrialsThis] = nf_tfspecparam_core_local(surf, freqs, opt, strcmpi(srcName, 'power'));

    TF.([srcName '_osc']) = osc_log;
    TF.([srcName '_ap']) = ap_log;

    if isfield(TF,'conds')
        if numel(badTrialsThis) ~= TF.conds
            error('Bad-trial mask length mismatch: got %d, expected %d.', numel(badTrialsThis), TF.conds);
        end
    elseif numel(badTrialsThis) ~= TF.ntrls
        error('Bad-trial mask length mismatch: got %d, expected %d.', numel(badTrialsThis), TF.ntrls);
    end

    badTrialsAll = badTrialsAll | badTrialsThis(:).';

    if strcmpi(srcName, 'power')
        TF.specparam = specparam;
        TF.specparam.method = 'nf_tfspecparam';
        TF.specparam.options = opt;
        TF.specparam.power_osc_field = 'power_osc';
        TF.specparam.power_ap_field = 'power_ap';
    end

end

if any(badTrialsAll)
    keepIdx = ~badTrialsAll;
    removedIdx = find(badTrialsAll);

    TF = nf_remove_trials_local(TF, keepIdx);

    if isfield(TF, 'specparam') && isstruct(TF.specparam)
        TF.specparam.removed_trials = removedIdx(:);
        TF.specparam.n_removed_trials = numel(removedIdx);
        TF.specparam.removal_reason = 'nf_specparam error in at least one time slice';
    end
end

TF = nf_powerfront_local(TF);

end

function [ap_log, osc_log, specparam, badTrials] = nf_tfspecparam_core_local(surf, freqs, opt, keepspecparam)

Fref = numel(freqs);

nd = ndims(surf);

if nd == 2
    error('Power field must be 3D or 4D.');
end

if nd == 3

    if size(surf, 2) == Fref
        surf = reshape(surf, size(surf, 1), size(surf, 2), size(surf, 3), 1);
    elseif size(surf, 1) == Fref
        surf = reshape(surf, 1, size(surf, 1), size(surf, 2), size(surf, 3));
    else
        error('Cannot infer power field layout relative to freqs.');
    end

end

if nd == 4
    if size(surf, 2) ~= Fref
        error('Second dimension of power must match numel(freqs).');
    end
end

[C, F, T, R] = size4d_local(surf);

if F ~= Fref
    error('Second dimension of power (F=%d) differs from numel(freqs) (%d).', F, Fref);
end

K = max(1, round(opt.maxPeaks));

if K > 200
    warning('maxPeaks is %d. Peak arrays will be large: C×K×T×R. Consider reducing maxPeaks.', K);
end

ap_log = zeros(C, F, T, R, 'like', surf);
osc_log = zeros(C, F, T, R, 'like', surf);

rsq = zeros(C, T, R);
mae = zeros(C, T, R);

ap_offset = nan(C, T, R);
ap_knee = nan(C, T, R);
ap_exponent = nan(C, T, R);

pk_n = zeros(C, T, R);
pk_f = nan(C, K, T, R);
pk_aLog10 = nan(C, K, T, R);
pk_sd = nan(C, K, T, R);

badTrials = false(1, R);

freqRow = double(freqs).';

fprintf(1, 'Time-frequency spectral parameterization progress:   0%%');

for ch = 1:C

    ap_ch = zeros(F, T, R, 'like', surf);
    osc_ch = zeros(F, T, R, 'like', surf);

    rsq_ch = zeros(T, R);
    mae_ch = zeros(T, R);

    ap_offset_ch = nan(T, R);
    ap_knee_ch = nan(T, R);
    ap_exp_ch = nan(T, R);

    pk_n_ch = zeros(T, R);
    pk_f_ch = nan(K, T, R);
    pk_a_ch = nan(K, T, R);
    pk_sd_ch = nan(K, T, R);

    bad_ch = false(1, R);

    % Critical: make sure workers are not in "stop on error" debug state.
    spmd
        warning('off','all');
        dbclear if error
        dbclear if warning
    end

    parfor r = 1:R

        trialFailed = 0;

        ap_r = nan(F, T, 'like', surf);
        osc_r = nan(F, T, 'like', surf);

        rsq_r = nan(T, 1);
        mae_r = nan(T, 1);

        ap_offset_r = nan(T, 1);
        ap_knee_r = nan(T, 1);
        ap_exp_r = nan(T, 1);

        pk_n_r = nan(T, 1);
        pk_f_r = nan(K, T, 'like', surf);
        pk_a_r = nan(K, T, 'like', surf);
        pk_sd_r = nan(K, T, 'like', surf);

        spec_log10_TR = squeeze(surf(ch, :, :, r)).';

        for t = 1:T

            dataX = spec_log10_TR(t, :);

            try
                parm = nf_specparam(dataX, freqRow, ...
                    'aPeriodicMode', opt.aPeriodicMode, ...
                    'minPeakHeight', opt.minPeakHeight, ...
                    'peakThreshold', opt.peakThreshold, ...
                    'peakType', opt.peakType, ...
                    'peakWidthLims', opt.peakWidthLims, ...
                    'threshAfter', opt.threshAfter, ...
                    'optim', opt.optim, ...
                    'plt', 0); %#ok
            catch
                trialFailed = 1;
                break
            end

            if strcmpi(opt.FoD, 'fit')
                ap_r(:, t) = parm.aPeriodicFit(:);
                osc_r(:, t) = parm.PeriodicFit(:);
            else
                ap_r(:, t) = parm.aPeriodicData(:);
                osc_r(:, t) = parm.PeriodicData(:);
            end

            rsq_r(t) = parm.RSquare;
            mae_r(t) = parm.MAE;

            apv = parm.aperiodicParms(:).';

            if strcmpi(parm.options.aPeriodicMode, 'fixed')
                ap_offset_r(t) = apv(1);
                ap_knee_r(t) = NaN;
                ap_exp_r(t) = apv(2);
            else
                ap_offset_r(t) = apv(1);
                ap_knee_r(t) = apv(2);
                ap_exp_r(t) = apv(3);
            end

            pk = parm.peakParms;
            n = size(pk, 1);

            pk_n_r(t) = n;

            if n > 0
                nuse = min(n, K);
                pk_f_r(1:nuse, t) = pk(1:nuse, 1);
                pk_a_r(1:nuse, t) = pk(1:nuse, 2);
                pk_sd_r(1:nuse, t) = pk(1:nuse, 3);
            end

        end

        bad_ch(r) = logical(trialFailed);

        ap_ch(:, :, r) = ap_r;
        osc_ch(:, :, r) = osc_r;

        rsq_ch(:, r) = rsq_r;
        mae_ch(:, r) = mae_r;

        ap_offset_ch(:, r) = ap_offset_r;
        ap_knee_ch(:, r) = ap_knee_r;
        ap_exp_ch(:, r) = ap_exp_r;

        pk_n_ch(:, r) = pk_n_r;
        pk_f_ch(:, :, r) = pk_f_r;
        pk_a_ch(:, :, r) = pk_a_r;
        pk_sd_ch(:, :, r) = pk_sd_r;

    end

    badTrials = badTrials | bad_ch;

    ap_log(ch, :, :, :) = reshape(ap_ch, [1 F T R]);
    osc_log(ch, :, :, :) = reshape(osc_ch, [1 F T R]);

    rsq(ch, :, :) = reshape(rsq_ch, [1 T R]);
    mae(ch, :, :) = reshape(mae_ch, [1 T R]);

    ap_offset(ch, :, :) = reshape(ap_offset_ch, [1 T R]);
    ap_knee(ch, :, :) = reshape(ap_knee_ch, [1 T R]);
    ap_exponent(ch, :, :) = reshape(ap_exp_ch, [1 T R]);

    pk_n(ch, :, :) = reshape(pk_n_ch, [1 T R]);
    pk_f(ch, :, :, :) = reshape(pk_f_ch, [1 K T R]);
    pk_aLog10(ch, :, :, :) = reshape(pk_a_ch, [1 K T R]);
    pk_sd(ch, :, :, :) = reshape(pk_sd_ch, [1 K T R]);

    fprintf(1, '\b\b\b\b%3d%%', round(100 * ch / C));

end

fprintf(1, '\n');

if keepspecparam == 1
    specparam = struct();
    specparam.rsquare = rsq;
    specparam.mae = mae;
    specparam.ap_offset = ap_offset;
    specparam.ap_knee = ap_knee;
    specparam.ap_exponent = ap_exponent;
    specparam.pk_n = pk_n;
    specparam.pk_f = pk_f;
    specparam.pk_aLog10 = pk_aLog10;
    specparam.pk_sd = pk_sd;
else
    specparam = struct();
end

end

function TF = nf_remove_trials_local(TF, keepIdx)

if isfield(TF,'conds')
    nOld = TF.conds;
elseif isfield(TF,'ntrls') && TF.ntrls>1
    nOld = TF.ntrls;
else
    error('Either TF.conds (averaged) or TF.trials > 1 required for removal.');
end

keepIdx = keepIdx(:).';
nNew = sum(keepIdx);

fn = fieldnames(TF);

for i = 1:numel(fn)

    name = fn{i};
    val = TF.(name);

    if isnumeric(val)

        nd = ndims(val);
        sz = size(val);

        if nd >= 3
            if sz(nd) == nOld
                idx = repmat({':'}, 1, nd);
                idx{nd} = keepIdx;
                TF.(name) = val(idx{:});
            end
        end

    elseif istable(val)

        if height(val) == nOld
            TF.(name) = val(keepIdx, :);
        end

    elseif isstruct(val)

        if strcmpi(name, 'behavior')
            TF.(name) = TF.(name)(keepIdx);
        end

        if strcmpi(name, 'specparam') || strcmpi(name, 'cpm')
            fn2 = fieldnames(TF.(name));
            for j = 1:numel(fn2)
                name2 = fn2{j};
                val2 = TF.(name).(name2);
                if isnumeric(val2)
                    nd2 = ndims(val2);
                    sz2 = size(val2);
                    if sz2(nd2) == nOld
                        idx2 = repmat({':'}, 1, nd2);
                        idx2{nd2} = keepIdx;
                        TF.(name).(name2) = val2(idx2{:});
                    end
                end
            end
        end
    end

    if isfield(TF,'conds')
        TF.conds = nNew;
    else
        TF.ntrls = nNew;
    end

end
end


function tf = nf_is_param_power_field_local(name)

tf = 0;

if endsWith(lower(name), '_ap')
    tf = 1;
end

if endsWith(lower(name), '_osc')
    tf = 1;
end

end

function powFields = nf_find_power_fields_local(TF)

fn = fieldnames(TF);

powFields = {};

for i = 1:numel(fn)

    name = fn{i};
    val = TF.(name);

    if isnumeric(val)
        if contains(lower(name), 'power')
            powFields{end + 1, 1} = name; %#ok
        end
    end

end

idx = find(strcmpi(powFields, 'power'));

if ~isempty(idx)
    powFields(idx) = [];
    powFields = [{'power'}; powFields];
end

end

function TF = nf_powerfront_local(TF)

fn = fieldnames(TF);

pow = {};
other = {};

for i = 1:numel(fn)

    name = fn{i};
    val = TF.(name);

    if isnumeric(val)
        if contains(lower(name), 'power')
            pow{end + 1, 1} = name; %#ok
        else
            other{end + 1, 1} = name; %#ok
        end
    else
        other{end + 1, 1} = name; %#ok
    end

end

idx = find(strcmpi(pow, 'power'));

if ~isempty(idx)
    pow(idx) = [];
    pow = [{'power'}; pow];
end

order = [pow; other];

TF = orderfields(TF, order);

end

function [C, F, T, R] = size4d_local(X)

sz = size(X);

C = sz(1);
F = sz(2);

T = 1;
R = 1;

if numel(sz) >= 3
    T = sz(3);
end

if numel(sz) >= 4
    R = sz(4);
end

end





