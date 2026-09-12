function Xs = hnf_surr_timeshift_vec(X, minShift, use_mex)
% SURR_TIMESHIFT_VEC    Circular time-shift surrogates
%
%   Xs = surr_timeshift_vec(X)
%   Xs = surr_timeshift_vec(X, minShift)
%   Xs = surr_timeshift_vec(X, minShift, use_mex)
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

if maxShift == 0
    % all shifts must be zero
    Xs_mat = Xmat;
else
    % draw random circular shifts for each series
    shiftRange = maxShift - minShift + 1;
    shifts = randi(shiftRange, [1, nVars]);
    shifts = shifts + (minShift - 1);
    
    % try MEX if requested and available
    use_mex_here = false;
    if use_mex
        if exist('surr_timeshift_mex', 'file') == 3
            use_mex_here = true;
        end
    end
    
    if use_mex_here
        % MEX implementation: column-wise circular shift in C
        try
            Xs_mat = surr_timeshift_mex(Xmat, shifts);
        catch
            warning('surr_timeshift_vec: using mex failed - using matlab-native (slower).');
            % pure MATLAB fallback: per-column circular shift
            Xs_mat = zeros(nSamples, nVars);
            idx_base = (0:(nSamples - 1))';
            for v = 1:nVars
                k = shifts(v);
                if k == 0
                    Xs_mat(:, v) = Xmat(:, v);
                else
                    idx = idx_base + k;
                    idx = 1 + mod(idx, nSamples);
                    Xs_mat(:, v) = Xmat(idx, v);
                end
            end
        end
    else
        % pure MATLAB fallback: per-column circular shift
        Xs_mat = zeros(nSamples, nVars);
        idx_base = (0:(nSamples - 1))';
        for v = 1:nVars
            k = shifts(v);
            if k == 0
                Xs_mat(:, v) = Xmat(:, v);
            else
                idx = idx_base + k;
                idx = 1 + mod(idx, nSamples);
                Xs_mat(:, v) = Xmat(idx, v);
            end
        end
    end
end

Xs = reshape(Xs_mat, sz);

end





