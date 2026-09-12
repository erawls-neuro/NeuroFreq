function [rho, pval] = hnf_circ_corrcl_vec(alpha, x)
% CIRC_CORRCL_VEC  Vectorized circular-linear correlations.
%
% [rho, pval] = circ_corrcl_vec(alpha, x)
%   Circular-linear correlation coefficient between one circular
%   (alpha, in radians) and one linear (x) variable.
%
%   Inputs:
%     alpha   array of angles in radians, arbitrary size
%     x       array of linear values, same size as alpha
%
%   Behavior:
%     - First dimension is treated as the sample / observation dimension.
%     - All remaining dimensions define independent variable locations.
%     - For each index in dims 2..end, a circ–lin correlation is computed
%       between alpha(:, idx) and x(:, idx).
%
%   Outputs:
%     rho     circular-linear correlation coefficients, same size as alpha
%             but with the first dimension collapsed to 1
%     pval    p-values from Zar's chi-square approximation, same size as rho
%
%   References:
%     Zar, J. H. (1999). Biostatistical Analysis. (eq. 27.47)
%


% Basic checks
if ~isequal(size(alpha), size(x))
    error('circ_corrcl_vec: alpha and x must have the same size');
end
if ~isreal(alpha)
    error('circ_corrcl_vec: alpha must be real-valued (angles in radians)');
end
if ~isreal(x)
    error('circ_corrcl_vec: x must be real-valued');
end
sz = size(alpha);
if numel(sz) < 1
    error('circ_corrcl_vec: alpha and x must have at least one dimension');
end
% Number of samples (observations) = size along first dimension
nSamples = sz(1);
if nSamples < 2
    error('circ_corrcl_vec: not enough samples along first dimension');
end
% checks complete - you passed!

% Number of variables = product of trailing dimensions
if isscalar(sz)
    nVars = 1;
else
    nVars = prod(sz(2:end));
end

% Reshape to [nSamples x nVars]
X = reshape(double(x),     [nSamples, nVars]);
A = reshape(double(alpha), [nSamples, nVars]);

% Circular components
S = sin(A);
C = cos(A);

% Center all variables along sample dimension (dim 1)
mx = mean(X, 1);
ms = mean(S, 1);
mc = mean(C, 1);
Xc = bsxfun(@minus, X, mx);
Sc = bsxfun(@minus, S, ms);
Cc = bsxfun(@minus, C, mc);

% Linear correlations:
% r_xy = sum(xc .* yc) ./ sqrt(sum(xc.^2) .* sum(yc.^2))
% All sums are along dim 1 (samples)
num_xs = sum(Xc .* Sc, 1);
den_xs = sqrt(sum(Xc.^2, 1) .* sum(Sc.^2, 1));
rxs    = num_xs ./ den_xs;
num_xc = sum(Xc .* Cc, 1);
den_xc = sqrt(sum(Xc.^2, 1) .* sum(Cc.^2, 1));
rxc    = num_xc ./ den_xc;
num_cs = sum(Cc .* Sc, 1);
den_cs = sqrt(sum(Cc.^2, 1) .* sum(Sc.^2, 1));
rcs    = num_cs ./ den_cs;

% Protect against zero-variance cases
zero_xs = den_xs == 0;
zero_xc = den_xc == 0;
zero_cs = den_cs == 0;
if any(zero_xs)
    rxs(zero_xs) = NaN;
end
if any(zero_xc)
    rxc(zero_xc) = NaN;
end
if any(zero_cs)
    rcs(zero_cs) = NaN;
end

% Circular-linear correlation (Zar, eq. 27.47)
% rho = sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs) / (1 - rcs^2))
den_rho = 1 - rcs.^2;
invalid_den = den_rho <= 0;
if any(invalid_den)
    den_rho(invalid_den) = NaN;
end

num_rho = (rxc.^2 + rxs.^2 - 2 .* rxc .* rxs .* rcs);
frac    = num_rho ./ den_rho;

% Small numerical negatives can pop up; clamp them
frac(frac < 0) = 0;
rho_vec = sqrt(frac);

% p-value: chi-square with 2 d.o.f.
% For df = 2, SF(z) = exp(-z/2), where z = nSamples * rho^2
z = nSamples .* rho_vec.^2;
pval_vec = exp(-0.5 .* z);

% Reshape back: same trailing dims, first dim collapsed to 1
out_sz    = sz;
out_sz(1) = 1;
rho       = reshape(rho_vec,  out_sz);
pval      = reshape(pval_vec, out_sz);

end







