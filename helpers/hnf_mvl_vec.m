function mvl = hnf_mvl_vec(phas, pow)
% MVL_VEC    Fully vectorized Mean Vector Length (MVL) computation.
%
% MVL = | mean( pow .* exp(1i * phas), 1 ) |
%
% Computed over the FIRST dimension (samples), for arbitrary trailing dims.
%
% INPUTS
%   phas : phase array (radians), any size
%   pow  : amplitude/power array, same size as phas
%          First dimension is treated as the sample dimension.
%
% OUTPUT
%   mvl  : MVL values, same size as phas/pow but with the first
%          dimension collapsed to 1.
%
% NOTES
%   - No NaN-handling beyond what MATLAB mean() does by default.
%   - Designed to be MEX-friendly: no loops, no I/O, pure arithmetic.
%


% Basic checks
if ~isequal(size(phas), size(pow))
    error('mvl_vec: phas and pow must have the same size');
end
if ~isreal(phas)
    error('mvl_vec: phas must be real-valued (angles in radians)');
end
if ~isreal(pow)
    error('mvl_vec: pow must be real-valued');
end
sz = size(phas);
if numel(sz) < 1
    error('mvl_vec: inputs must have at least one dimension');
end
% Number of samples along first dimension
nSamples = sz(1);
if nSamples < 1
    error('mvl_vec: first dimension must have at least one sample');
end
% check complete - you passed!

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

% Complex weighted vectors
z = pow_mat .* exp(1i .* phas_mat);

% Average across samples (first dimension)
z_mean = mean(z, 1);                  % [1 x nVars]

% MVL is magnitude of the mean vector
mvl_vec = abs(z_mean);                % [1 x nVars]

% Reshape back to original trailing dimensions
out_sz    = sz;
out_sz(1) = 1;
mvl       = reshape(mvl_vec, out_sz);

end








