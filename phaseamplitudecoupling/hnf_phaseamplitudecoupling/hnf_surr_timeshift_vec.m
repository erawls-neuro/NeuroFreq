function Xs = hnf_surr_timeshift_vec(X, minShift, use_mex)
% SURR_TIMESHIFT_VEC    Circular time-shift surrogates
%
%   Xs = hnf_surr_timeshift_vec(X)
%   Xs = hnf_surr_timeshift_vec(X, minShift)
%   Xs = hnf_surr_timeshift_vec(X, minShift, use_mex)
%
% Generates circular time-shifted surrogates for real-valued time series.
%
% INPUT
%   X        : real-valued array, arbitrary size.
%              First dimension = samples (time).
%              All trailing dimensions define independent series.
%   minShift : optional, minimum shift in samples (default = 1).
%              For each series a random integer shift k is drawn
%              uniformly from [minShift, nSamples-1]. If nSamples <= 1,
%              X is returned unchanged.
%   use_mex  : optional logical flag (default = true).
%              If true and surr_timeshift_mex is available, the MEX
%              implementation is used. Otherwise, a pure-MATLAB
%              fallback is used.
%
% OUTPUT
%   Xs       : surrogate array, same size as X.
%
% Notes:
%   - Each column after flattening trailing dims is treated as an
%     independent time series.
%   - Temporal autocorrelation and marginal distribution of each
%     series are preserved exactly; only the alignment along the
%     time axis is changed.
%

% basic checks
if ~isreal(X)
    error('surr_timeshift_vec: X must be real-valued');
end

sz = size(X);
if numel(sz) < 1
    error('surr_timeshift_vec: input must have at least one dimension');
end

nSamples = sz(1);

if nargin < 2 || isempty(minShift)
    minShift = 1;
end

if nargin < 3 || isempty(use_mex)
    use_mex = true;
end

if ~isscalar(minShift) || minShift < 0 || minShift ~= round(minShift)
    error('surr_timeshift_vec: minShift must be a non-negative integer scalar');
end

if ~isscalar(use_mex) || (~islogical(use_mex) && ~ismember(use_mex, [0 1]))
    error('surr_timeshift_vec: use_mex must be a logical scalar');
end

if nSamples <= 1
    Xs = X;
    return
end

maxShift = nSamples - 1;
if minShift > maxShift
    warning('surr_timeshift_vec: minShift > nSamples-1; using minShift = 0');
    minShift = 0;
end

% number of independent series (flatten trailing dims)
if isscalar(sz)
    nVars = 1;
else
    nVars = prod(sz(2:end));
end

% reshape to [nSamples x nVars]
Xmat = reshape(double(X), [nSamples, nVars]);

% Match circshift(x, +k, 1) in both backends. Draw shifts only once.
shifts = randi([minShift, maxShift], 1, nVars);
if use_mex && exist('hnf_surr_timeshift_mex', 'file') == 3
    try
        Xs = reshape(hnf_surr_timeshift_mex(Xmat, shifts), sz);
        return
    catch failure
        warning('NeuroFreq:ShiftFallback', 'Using MATLAB shifts: %s', failure.message);
    end
end
Xs_mat = zeros(nSamples, nVars);
for v = 1:nVars
    Xs_mat(:,v) = circshift(Xmat(:,v), shifts(v), 1);
end
Xs = reshape(Xs_mat, sz);
end
