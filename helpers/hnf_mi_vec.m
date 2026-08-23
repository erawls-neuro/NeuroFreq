function mi = hnf_mi_vec(phas, pow, nbin)
% MI_VEC    Fully vectorized Tort-modulation-index (MI) computation.
%
% Computes, for each variable (all dims except the first):
%
%   1) Bin phases into nbin bins over [-pi, pi).
%   2) For each bin, compute mean amplitude (pow).
%   3) Normalize to probability distribution P over bins.
%   4) Compute KL divergence of P vs uniform Q, and normalize:
%
%      MI = (1 / log(nbin)) * sum_b P(b) * log(P(b) / Q(b))
%         = (1 / log(nbin)) * sum_b P(b) * log(P(b) * nbin)
%
% FIRST dimension is treated as the sample dimension.
%
% INPUTS
%   phas : phase array (radians), arbitrary size
%   pow  : amplitude/power array, same size as phas
%          size(phas,1) = number of samples / time points
%   nbin : number of phase bins (scalar)
%
% OUTPUT
%   mi   : MI values, same size as phas/pow, but with the first
%          dimension collapsed to 1
%
% NOTES
%   - No permutations or z-scoring here: this is the base MI kernel.
%   - Designed to be MEX-friendly: reshape + accumarray only.
%

% Basic checks
if nargin < 3 || isempty(nbin)
    error('mi_vec: nbin must be provided as a positive integer');
end

if ~isequal(size(phas), size(pow))
    error('mi_vec: phas and pow must have the same size');
end
if ~isreal(phas)
    error('mi_vec: phas must be real-valued (angles in radians)');
end
if ~isscalar(nbin) || nbin < 1 || nbin ~= round(nbin)
    error('mi_vec: nbin must be a positive integer scalar');
end
sz = size(phas);
if numel(sz) < 1
    error('mi_vec: inputs must have at least one dimension');
end
nSamples = sz(1);
if nSamples < 2
    warning('mi_vec: fewer than 2 samples; MI will be ill-defined');
end
% checks complete - you passed!

% Number of variables = product of trailing dimensions
if isscalar(sz)
    nVars = 1;
else
    nVars = prod(sz(2:end));
end

% Cast to double for numerical stability / MEX friendliness
phas_d = double(phas);
pow_d  = double(pow);

% Reshape to [nSamples x nVars]
phas_mat = reshape(phas_d, [nSamples, nVars]);
pow_mat  = reshape(pow_d,  [nSamples, nVars]);

% Wrap phases to [-pi, pi) via complex exponential
phas_wrapped = angle(exp(1i .* phas_mat));

% Compute bin indices for each sample and variable
binWidth = 2 * pi / nbin;

% bin index = floor((phi + pi) / binWidth) + 1, then clamp
binIdx = floor((phas_wrapped + pi) ./ binWidth) + 1;
binIdx(binIdx < 1)    = 1;
binIdx(binIdx > nbin) = nbin;

% Accumulate amplitude within each (bin, variable) pair
pow_vec = pow_mat(:);           % [nSamples * nVars x 1]
bin_vec = binIdx(:);            % same size

% Variable index per sample: 1..nVars, each repeated nSamples times
var_mat = repmat(1:nVars, nSamples, 1);
var_vec = var_mat(:);

% Linear index: 1..(nbin * nVars)
key = (bin_vec - 1) .* nVars + var_vec;

% Sum of amplitudes per (bin, var)
sumAmp = accumarray(key, pow_vec, [nbin * nVars, 1], @sum, 0);

% Count of samples per (bin, var)
countAmp = accumarray(key, 1, [nbin * nVars, 1], @sum, 0);

% Reshape back to [nbin x nVars]
sumAmp_mat   = reshape(sumAmp,   [nbin, nVars]);
countAmp_mat = reshape(countAmp, [nbin, nVars]);

% Mean amplitude per bin; empty bins get 0
amplBin = sumAmp_mat ./ countAmp_mat;
amplBin(countAmp_mat == 0) = 0;

% Convert to probability distribution over bins: P(b, var)
colSum = sum(amplBin, 1);           % [1 x nVars]
zeroCols         = (colSum == 0);   % no amplitude at all
colSum(zeroCols) = NaN;
amplP = bsxfun(@rdivide, amplBin, colSum);   % [nbin x nVars]
amplP(isnan(amplP)) = 0;                     % bins with 0/NaN → P = 0

% KL divergence to uniform and normalized MI
% Q(b) = 1 / nbin, so log(P/Q) = log(P * nbin)
term = zeros(size(amplP));
mask = amplP > 0;
term(mask) = amplP(mask) .* log(amplP(mask) .* nbin);
distKL = sum(term, 1);               % [1 x nVars]
mi_vec = distKL ./ log(nbin);        % Tort normalization

% Columns with zero amplitude across all bins → NaN
mi_vec(zeroCols) = NaN;

% Reshape back to original trailing dimensions
out_sz    = sz;
out_sz(1) = 1;
mi        = reshape(mi_vec, out_sz);

end





