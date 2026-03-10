function TF = nf_avebase( TF, method, blTimes, trlvec, avmode )
% NF_AVEBASE    Averages/baselines TF structures
%
% GENERAL
% -------
% Averages/baselines TF structures from nf_tftransform.m
% Power is averaged and optionally baseline corrected according to either
% dB, %-change, or z-score. If the signal contains phase, phase is
% processed into ITC.
%
% OUTPUT
% ------
% TF - structure output by nf_tftransform.m (averaged)
%
% INPUT
% -----
% 1) TF - structure output by nf_tftransform.m
% 2) method - 'none', 'zscore', 'db', 'percent'
% 3) blTimes - [min max] times for baseline, in ms
% 4) trlvec - vector describing which trials to average together. Each
%     unique value in the vector is taken to be a different trial type.
% 5) avmode - 'mean' or 'median'
%

if nargin < 5 || isempty(avmode)
    disp('No averaging mode supplied - using mean');
    avmode = 'mean';
end

if nargin < 4 || isempty(trlvec)
    disp('No trial vector supplied - averaging all trials together');
    trlvec = ones(1, size(TF.power, ndims(TF.power)));
end

if nargin < 3 || isempty(blTimes)
    disp('no times provided, using all times before 0');
    blIdx = find(TF.times <= 0);
else
    blIdx = find((TF.times >= blTimes(1)) + (TF.times <= blTimes(end)) == 2);
end

if nargin < 2 || isempty(method)
    method = 'none';
else
    method = lower(method);
end

if nargin < 1 || isempty(TF)
    error('at least a TF structure is required input');
end

if ~isfield(TF, 'power')
    error('TF.power is required');
end

if ~isfield(TF, 'times')
    error('TF.times is required');
end

if ~isfield(TF, 'scale')
    error('TF.scale is required');
end

if ~(strcmpi(TF.scale, 'linear') || strcmpi(TF.scale, 'log10'))
    error('TF.scale must be ''linear'' or ''log10''.');
end

if isempty(blIdx)
    error('Baseline index is empty. Check TF.times and blTimes.');
end

trlvec = trlvec(:).';

if ~isfield(TF, 'ntrls')
    error('TF.ntrls is required');
end

if TF.ntrls == 1
    error('TF is already averaged (ntrls = 1)');
end

if ~isfield(TF, 'nsensor')
    error('TF.nsensor is required');
end

if TF.nsensor == 1
    flagsens = 1;
else
    flagsens = 0;
end

nTrials = size(TF.power, ndims(TF.power));

if numel(trlvec) ~= nTrials
    error('trlvec length (%d) does not match number of trials in TF.power (%d).', numel(trlvec), nTrials);
end

conds = unique(trlvec);

powFields = nf_find_power_fields_local(TF);

newPow = struct();

for i = 1:numel(powFields)

    name = powFields{i};
    X = TF.(name);

    if size(X, ndims(X)) ~= nTrials
        error('Power field %s has %d trials; expected %d.', name, size(X, ndims(X)), nTrials);
    end

    sz = size(X);
    sz = sz(1:end-1);
    newPow.(name) = zeros([sz, numel(conds)], 'like', X);

end

if isfield(TF, 'phase')
    szp = size(TF.phase);
    szp = szp(1:end-1);
    newphase = zeros([szp, numel(conds)], 'like', TF.phase);
end

if isfield(TF, 'behavior')
    newBeh = TF.behavior;
    newBeh(numel(conds) + 1:end) = [];
end

for n = 1:numel(conds)

    idx = find(trlvec == conds(n));

    if isscalar(idx)
        warning(['condition ' num2str(n) ' has only one trial...']);
        flagst = 1;
    else
        flagst = 0;
    end

    for i = 1:numel(powFields)

        name = powFields{i};
        X = TF.(name);

        if flagsens == 0
            condP = X(:,:,:,idx);
            newPow.(name)(:,:,:,n) = corrP(condP, blIdx, method, flagsens, avmode, TF.scale, flagst);
        else
            condP = X(:,:,idx);
            newPow.(name)(:,:,n) = corrP(condP, blIdx, method, flagsens, avmode, TF.scale, flagst);
        end

    end

    if isfield(TF, 'phase')
        if flagsens == 0
            condphase = TF.phase(:,:,:,idx);
            newphase(:,:,:,n) = squeeze(abs(mean(exp(1i .* condphase), 4)));
        else
            condphase = TF.phase(:,:,idx);
            newphase(:,:,n) = squeeze(abs(mean(exp(1i .* condphase), 3)));
        end
    end

    TF.trlerp(n) = numel(idx);

    if isfield(TF, 'behavior')
        newBeh(n) = fieldfun(@(varargin) nanmean([varargin{:}]), TF.behavior(idx));
    end

end

TF.ntrls = 1;
TF.conds = numel(conds);

for i = 1:numel(powFields)
    name = powFields{i};
    TF.(name) = newPow.(name);
end

if isfield(TF, 'phase')
    TF.phase = newphase;
end

if isfield(TF, 'behavior')
    TF.behavior = newBeh;
end

if isfield(TF, 'event')
    TF = rmfield(TF, 'event');
end

if isfield(TF, 'epoch')
    TF = rmfield(TF, 'epoch');
end

TF = nf_powerfront_local(TF);

end

function power = corrP(power, blIdx, method, flagsens, avmode, powScale, flagst)

if ~flagst
    tfPow = squeeze(mean(power, ndims(power)));
else
    tfPow = power;
end

if flagsens == 0

    if strcmp(avmode, 'mean')
        blPow = squeeze(mean(tfPow(:,:,blIdx), 3));
        blPow = repmat(blPow, [1 1 size(power, 3)]);
        blPowSTD = std(permute(tfPow(:,:,blIdx), [3 1 2]));
        blPowSTD = squeeze(blPowSTD);
        blPowSTD = repmat(blPowSTD, [1 1 size(power, 3)]);
    else
        blPow = squeeze(median(tfPow(:,:,blIdx), 3));
        blPow = repmat(blPow, [1 1 size(power, 3)]);
        blPowSTD = mad(permute(tfPow(:,:,blIdx), [3 1 2]));
        blPowSTD = squeeze(blPowSTD);
        blPowSTD = repmat(blPowSTD, [1 1 size(power, 3)]);
    end

else

    if strcmp(avmode, 'mean')
        blPow = squeeze(mean(tfPow(:, blIdx), 2));
        blPow = repmat(blPow, [1 size(power, 2)]);
        blPowSTD = std(permute(tfPow(:, blIdx), [2 1]));
        blPowSTD = squeeze(blPowSTD);
        blPowSTD = repmat(blPowSTD, [1 size(power, 2)]);
    else
        blPow = squeeze(median(tfPow(:, blIdx), 2));
        blPow = repmat(blPow, [1 size(power, 2)]);
        blPowSTD = mad(permute(tfPow(:, blIdx), [2 1]));
        blPowSTD = squeeze(blPowSTD);
        blPowSTD = repmat(blPowSTD, [1 size(power, 2)]);
    end

end

blPowSTD(blPowSTD == 0) = eps;

switch method

    case 'zscore'
        disp('baselining with z-score');
        tfPow = (tfPow - blPow) ./ blPowSTD;

    case 'db'
        disp('baselining with decibel');

        if strcmpi(powScale, 'log10')
            tfPow = 10 .* (tfPow - blPow);
        else
            blPow(blPow <= 0) = eps;
            tfPow(tfPow <= 0) = eps;
            tfPow = 10 .* log10(tfPow ./ blPow);
        end

    case 'percent'
        disp('baselining using %-change');

        if strcmpi(powScale, 'log10')
            P = 10 .^ tfPow;
            Pb = 10 .^ blPow;
            Pb(Pb == 0) = eps;
            tfPow = 100 .* (P - Pb) ./ Pb;
        else
            blPow(blPow == 0) = eps;
            tfPow = 100 .* (tfPow - blPow) ./ blPow;
        end

    case 'none'
        disp('NOT baselining power estimates');

    otherwise
        error('Unknown baselining method: %s', method);

end

power = tfPow;

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

function out = fieldfun(fun, S)

if nargin < 2
    error('fieldfun requires a function handle and a struct array.');
end

if ~isa(fun, 'function_handle')
    error('fieldfun first input must be a function handle.');
end

if ~isstruct(S)
    error('fieldfun second input must be a struct array.');
end

f = fieldnames(S);

out = struct();

for i = 1:numel(f)

    fname = f{i};
    vals = {S.(fname)};

    out.(fname) = fun(vals{:});

end

end

function m = nanmean(x)

if nargin < 1
    error('nanmean requires an input.');
end

if isempty(x)
    m = NaN;
    return
end

if isfloat(x) == 0
    x = double(x);
end

if exist('mean', 'builtin') || exist('mean', 'file')
    try
        m = mean(x, 'omitnan');
        return
    catch
    end
end

x = x(:);
x = x(isfinite(x));

if isempty(x)
    m = NaN;
else
    m = mean(x);
end

end















