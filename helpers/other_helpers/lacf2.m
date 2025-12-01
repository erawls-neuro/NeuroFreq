function lacf = lacf2(x, mlag)
% lacf2 -- Type II local autocorrelation function (positive lags only)
%          Optionally uses lacf2_mex if present.

x = x(:);
N = length(x);

if nargin < 2 || isempty(mlag)
    mlag = N;
end

if mlag > N
    error('mlag must be <= length(x).');
end

% if issingle(x)
%     x=double(x);
% end

if exist('lacf2_mex', 'file') == 3
    lacf = lacf2_mex(double(x), mlag);
    return;
end

% MATLAB fallback: original Jeff O''Neill logic
lacf = zeros(mlag, N);

for t = 1:N
    mtau = min(mlag, N - t + 1);
    lacf(1:mtau, t) = conj(x(t)) .* x(t:t + mtau - 1);
end

end