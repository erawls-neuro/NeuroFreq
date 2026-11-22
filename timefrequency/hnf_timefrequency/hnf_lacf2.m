function lacf = hnf_lacf2(x, mlag)
% lacf2 -- Type II local autocorrelation function (positive lags only)
%          Optionally uses lacf2_mex if present.

validateattributes(x, {'numeric'}, {'vector','nonempty','nonsparse'});
x = double(x(:));
N = length(x);

if nargin < 2 || isempty(mlag)
    mlag = N;
end

validateattributes(mlag, {'numeric'}, {'scalar','integer','real','positive','finite'});
if mlag > N
    error('mlag must be <= length(x).');
end

if exist('hnf_lacf2_mex', 'file') == 3
    lacf = hnf_lacf2_mex(x, double(mlag));
    return;
end

% MATLAB fallback: original Jeff O''Neill logic
lacf = zeros(mlag, N);

for t = 1:N
    mtau = min(mlag, N - t + 1);
    lacf(1:mtau, t) = conj(x(t)) .* x(t:t + mtau - 1);
end

end
