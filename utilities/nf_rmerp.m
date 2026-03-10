function TF = nf_rmerp(TF, trlvec, mode)
% NF_RMERP    Removes the trial-averaged ERP from TF sets (analytic-domain)
%
% GENERAL
% -------
% Reconstitutes analytic signal (complex) from TF.power + TF.phase, computes
% the condition ERP (trial mean in complex domain), and separates each trial
% into evoked vs induced components by subtracting a per-trial ERP estimate.
%
% MODES (third argument)
% ---------------------
%   'standard'  : per-trial ERP is exactly the condition mean ERP
%   'dot'       : per-trial ERP is condition ERP scaled by dot-product similarity;
%                 scaling is normalized so mean(scale)=1
%   'dotderiv'  : per-trial ERP is fit as a*mu + b*dmu (mu = condition ERP),
%                 with constraints ENFORCED per (chan,freq):
%                    mean(a)=1 and mean(b)=0
%
% INPUT
% -----
%   TF     : TF structure output by nf_tftransform.m with fields:
%              TF.power, TF.phase, TF.scale ('linear'|'log10'), TF.ntrls, TF.nsensor
%   trlvec : vector (1 x nTrials) defining condition labels per trial
%            If empty, all trials treated as one condition.
%   mode   : string: 'standard' | 'dot' | 'dotderiv'
%
% OUTPUT (added fields)
% ---------------------
%   TF.power_evoked   : same size as TF.power
%   TF.power_induced  : same size as TF.power
%   TF.phase_induced  : same size as TF.phase
%   TF.erprem         : metadata struct
%
% NOTES
% -----
% - Guarantees (within floating point) that in 'dot' and 'dotderiv', the mean
%   over trials of the per-trial ERP exactly equals the condition ERP (mu).
% - Does NOT support averaged TF (TF.ntrls==1).
%

% Inputs / defaults
if nargin < 3 || isempty(mode)
    mode = 'standard';
end

if nargin < 2 || isempty(trlvec)
    disp('No trial vector supplied - averaging all trials together');
    trlvec = ones(1, size(TF.power, ndims(TF.power)));
end

if nargin < 1 || isempty(TF)
    error('at least a TF structure is required input');
end

if isstring(mode) == 1
    mode = char(mode);
end

if ischar(mode) == 0
    error('mode must be a string: ''standard'', ''dot'', or ''dotderiv''.');
end

mode = lower(mode);

if ~(strcmp(mode, 'standard') == 1 || strcmp(mode, 'dot') == 1 || strcmp(mode, 'dotderiv') == 1)
    error('Unknown mode: %s. Allowed: ''standard'', ''dot'', ''dotderiv''.', mode);
end

% Validate TF structure
if ~isfield(TF, 'power')
    error('TF.power is required');
end

if ~isfield(TF, 'phase')
    error('TF.phase is required for ERP removal (analytic reconstruction)');
end

if ~isfield(TF, 'scale')
    error('TF.scale is required (''linear'' or ''log10'')');
end

if ~(strcmpi(TF.scale, 'linear') == 1 || strcmpi(TF.scale, 'log10') == 1)
    error('TF.scale must be ''linear'' or ''log10''.');
end

if ~isfield(TF, 'ntrls')
    error('TF.ntrls is required');
end

if TF.ntrls == 1
    error('TF is already averaged (ntrls = 1)');
end

if ~isfield(TF, 'nsensor')
    error('TF.nsensor is required');
end

trlvec = trlvec(:).';

nTrials = size(TF.power, ndims(TF.power));

if numel(trlvec) ~= nTrials
    error('trlvec length (%d) does not match number of trials in TF.power (%d).', numel(trlvec), nTrials);
end

if TF.nsensor == 1
    flagsens = 1;
else
    flagsens = 0;
end

conds = unique(trlvec);

% Allocate outputs
TF.power_evoked = zeros(size(TF.power), 'like', TF.power);
TF.power_induced = zeros(size(TF.power), 'like', TF.power);
TF.phase_induced = zeros(size(TF.phase), 'like', TF.phase);

% Main loop across conditions
for n = 1:numel(conds)

    idx = find(trlvec == conds(n));

    if flagsens == 0

        % multi-sensor: C x F x T x Rcond
        if strcmpi(TF.scale, 'linear') == 1
            condAmp = sqrt(TF.power(:,:,:,idx));
        else
            condAmp = 10 .^ (TF.power(:,:,:,idx) ./ 2);
        end

        condPh = TF.phase(:,:,:,idx);

        analyticSignal = condAmp .* exp(1i .* condPh);

        mu = mean(analyticSignal, 4);

        if strcmp(mode, 'standard') == 1

            erp = repmat(mu, [1 1 1 numel(idx)]);

        elseif strcmp(mode, 'dot') == 1

            erp = local_fit_erp_dot_4d(analyticSignal, mu);

        elseif strcmp(mode, 'dotderiv') == 1

            erp = local_fit_erp_dotderiv_4d(analyticSignal, mu, TF.times);

        else

            error('Unexpected mode branch.');

        end

        analyticInduced = analyticSignal - erp;

        powInd = abs(analyticInduced) .^ 2;
        powEvp = abs(erp) .^ 2;

        if strcmpi(TF.scale, 'log10') == 1
            powInd(powInd <= 0) = eps;
            powEvp(powEvp <= 0) = eps;
            powInd = log10(powInd);
            powEvp = log10(powEvp);
        end

        TF.power_induced(:,:,:,idx) = powInd;
        TF.power_evoked(:,:,:,idx) = powEvp;
        TF.phase_induced(:,:,:,idx) = angle(analyticInduced);

    else

        % single-sensor: F x T x Rcond
        if strcmpi(TF.scale, 'linear') == 1
            condAmp = sqrt(TF.power(:,:,idx));
        else
            condAmp = 10 .^ (TF.power(:,:,idx) ./ 2);
        end

        condPh = TF.phase(:,:,idx);

        analyticSignal = condAmp .* exp(1i .* condPh);

        mu = mean(analyticSignal, 3);

        if strcmp(mode, 'standard') == 1

            erp = repmat(mu, [1 1 numel(idx)]);

        elseif strcmp(mode, 'dot') == 1

            erp = local_fit_erp_dot_3d(analyticSignal, mu);

        elseif strcmp(mode, 'dotderiv') == 1

            erp = local_fit_erp_dotderiv_3d(analyticSignal, mu, TF.times);

        else

            error('Unexpected mode branch.');

        end

        analyticInduced = analyticSignal - erp;

        powInd = abs(analyticInduced) .^ 2;
        powEvp = abs(erp) .^ 2;

        if strcmpi(TF.scale, 'log10') == 1
            powInd(powInd <= 0) = eps;
            powEvp(powEvp <= 0) = eps;
            powInd = log10(powInd);
            powEvp = log10(powEvp);
        end

        TF.power_induced(:,:,idx) = powInd;
        TF.power_evoked(:,:,idx) = powEvp;
        TF.phase_induced(:,:,idx) = angle(analyticInduced);

    end

end

% Metadata
TF.erprem = struct();
TF.erprem.method = 'nf_rmerp';
TF.erprem.mode = mode;
TF.erprem.trlvec = trlvec(:);
TF.erprem.power_evoked_field = 'power_evoked';
TF.erprem.power_induced_field = 'power_induced';
TF.erprem.phase_induced_field = 'phase_induced';

TF = nf_powerfront_local(TF);

end

% Local: DOT scaling (4D: C x F x T x R)
% Enforces mean(scale)=1 per (C,F) => mean(erp)=mu exactly.
function erp = local_fit_erp_dot_4d(X, mu)

[C, F, T, R] = size(X);

erp = zeros(C, F, T, R, 'like', X);

for ch = 1:C
    for f = 1:F

        xcf = squeeze(X(ch, f, :, :));    % T x R
        muf = squeeze(mu(ch, f, :));      % T x 1

        d = muf' * xcf;                   % 1 x R (complex)
        d = abs(d);                       % similarity magnitude

        d = d - min(d) + eps;             % nonnegative, avoid zeros
        m = mean(d);

        if m == 0
            m = eps;
        end

        w = d ./ m;                       % mean(w)=1

        erp(ch, f, :, :) = repmat(muf, [1 R]) .* reshape(w, [1 R]);

    end
end

end

% Local: DOT scaling (3D: F x T x R)
% Enforces mean(scale)=1 per F => mean(erp)=mu exactly.
function erp = local_fit_erp_dot_3d(X, mu)

[F, T, R] = size(X);

erp = zeros(F, T, R, 'like', X);

for f = 1:F

    xf = squeeze(X(f, :, :));      % T x R
    muf = squeeze(mu(f, :)).';     % T x 1

    d = muf' * xf;                 % 1 x R (complex)
    d = abs(d);

    d = d - min(d) + eps;
    m = mean(d);

    if m == 0
        m = eps;
    end

    w = d ./ m;                    % mean(w)=1

    erp(f, :, :) = (repmat(muf, [1 R]) .* reshape(w, [1 R])).';

end

% Above creates R x T; permute back to F x T x R
erp = permute(erp, [1 2 3]);

end

% Local: DOT+DERIV fit (4D: C x F x T x R)
% Fits ERP_r(t) = a_r * mu(t) + b_r * dmu(t).
% Enforces per (C,F): mean(a)=1 and mean(b)=0 EXACTLY => mean(ERP)=mu.
function erp = local_fit_erp_dotderiv_4d(X, mu, times)

[C, F, T, R] = size(X);

erp = zeros(C, F, T, R, 'like', X);

dmu = local_time_derivative(mu, times);   % C x F x T

for ch = 1:C
    for f = 1:F

        xcf = squeeze(X(ch, f, :, :));        % T x R
        muf = squeeze(mu(ch, f, :));          % T x 1
        dmf = squeeze(dmu(ch, f, :));         % T x 1

        [a, b] = local_fit_ab_constrained(xcf, muf, dmf);

        % Build per-trial ERP
        for r = 1:R
            erp(ch, f, :, r) = (a(r) .* muf) + (b(r) .* dmf);
        end

    end
end

end

% Local: DOT+DERIV fit (3D: F x T x R)
% Same constraints: mean(a)=1 and mean(b)=0 => mean(ERP)=mu.
function erp = local_fit_erp_dotderiv_3d(X, mu, times)

[F, T, R] = size(X);

erp = zeros(F, T, R, 'like', X);

dmu = local_time_derivative_2d(mu, times);   % F x T

for f = 1:F

    xf = squeeze(X(f, :, :));        % T x R
    muf = squeeze(mu(f, :)).';       % T x 1
    dmf = squeeze(dmu(f, :)).';      % T x 1

    [a, b] = local_fit_ab_constrained(xf, muf, dmf);

    for r = 1:R
        erp(f, :, r) = (a(r) .* muf + b(r) .* dmf).';
    end

end

end

% Local: fit a,b per trial then enforce mean(a)=1 and mean(b)=0
% x is T x R, mu and dmu are T x 1.
% Uses complex least squares with 2 regressors, then constraint fix.
function [a, b] = local_fit_ab_constrained(x, mu, dmu)

[~, R] = size(x);

Phi = [mu(:) dmu(:)];               % T x 2 (complex allowed)
G = Phi' * Phi;                     % 2 x 2
Gt = Phi';                          % 2 x T

a = zeros(R, 1);
b = zeros(R, 1);

usePinv = false;

if rcond(double(G)) < 1e-12
    usePinv = true;
end

for r = 1:R

    y = x(:, r);

    if usePinv == false
        ab = G \ (Gt * y);
    else
        ab = pinv(G) * (Gt * y);
    end

    a(r) = ab(1);
    b(r) = ab(2);

end

% Enforce exact recapitulation constraints
% 1) Scale so mean(a)=1
ma = mean(a);

if ma == 0
    ma = eps;
end

k = 1 ./ ma;

a = a .* k;
b = b .* k;

% 2) Center b so mean(b)=0
mb = mean(b);
b = b - mb;

% After this:
% mean( a*mu + b*dmu ) = mean(a)*mu + mean(b)*dmu = 1*mu + 0*dmu = mu

end

% Local: derivative for 3D mu: C x F x T
function dmu = local_time_derivative(mu, times)

sz = size(mu);

if numel(sz) ~= 3
    error('mu must be 3D (C x F x T).');
end

C = sz(1);
F = sz(2);
T = sz(3);

dmu = zeros(size(mu), 'like', mu);

if nargin < 2 || isempty(times) == 1
    dt = 1;
else
    times = times(:);
    if numel(times) ~= T
        dt = 1;
    else
        dt = median(diff(times));
        if dt == 0
            dt = 1;
        end
    end
end

for ch = 1:C
    for f = 1:F

        v = squeeze(mu(ch, f, :));

        if T == 1
            dv = zeros(1, 1, 'like', v);
        elseif T == 2
            dv = [(v(2) - v(1)) / dt; (v(2) - v(1)) / dt];
        else
            dv = zeros(T, 1, 'like', v);
            dv(1) = (v(2) - v(1)) / dt;
            dv(T) = (v(T) - v(T - 1)) / dt;

            for t = 2:(T - 1)
                dv(t) = (v(t + 1) - v(t - 1)) / (2 * dt);
            end
        end

        dmu(ch, f, :) = reshape(dv, [1 1 T]);

    end
end

end

% Local: derivative for 2D mu: F x T
function dmu = local_time_derivative_2d(mu, times)

[F, T] = size(mu);

dmu = zeros(size(mu), 'like', mu);

if nargin < 2 || isempty(times) == 1
    dt = 1;
else
    times = times(:);
    if numel(times) ~= T
        dt = 1;
    else
        dt = median(diff(times));
        if dt == 0
            dt = 1;
        end
    end
end

for f = 1:F

    v = squeeze(mu(f, :)).';

    if T == 1
        dv = zeros(1, 1, 'like', v);
    elseif T == 2
        dv = [(v(2) - v(1)) / dt; (v(2) - v(1)) / dt];
    else
        dv = zeros(T, 1, 'like', v);
        dv(1) = (v(2) - v(1)) / dt;
        dv(T) = (v(T) - v(T - 1)) / dt;

        for t = 2:(T - 1)
            dv(t) = (v(t + 1) - v(t - 1)) / (2 * dt);
        end
    end

    dmu(f, :) = dv(:).';

end

end

% Local: push power fields to front (same logic you used)
function TF = nf_powerfront_local(TF)

fn = fieldnames(TF);

pow = {};
other = {};

for i = 1:numel(fn)

    name = fn{i};
    val = TF.(name);

    if isnumeric(val)
        if contains(lower(name), 'power')
            pow{end + 1, 1} = name; %#ok<AGROW>
        else
            other{end + 1, 1} = name; %#ok<AGROW>
        end
    else
        other{end + 1, 1} = name; %#ok<AGROW>
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





