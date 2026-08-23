function PAC = nf_mvlpac_within(phas, pow, nperm, use_mex)
% NF_MVLPAC_WITHIN    Phase-amplitude coupling via Mean Vector Length (MVL)
%
% GENERAL
% -------
% Averages phase-amplitude coupling using the MVL measure (Canolty 2006),
% with optional z-scoring via permutation as in Cohen 2014.
% MVL is computed over time samples for each channel × low-freq ×
% high-freq × trial combination.
%
% USAGE
% -----
%   PAC = nf_mvlpac(phas, pow)
%   PAC = nf_mvlpac(phas, pow, nperm)
%   PAC = nf_mvlpac(phas, pow, nperm, use_mex)
%
% OUTPUT:
%   PAC  - MVL or zMVL values with size:
%          channels x loFreqs x hiFreqs x trials
%
% INPUTS:
%   phas    - phase matrix, size:
%             channels x loFreqs x times x trials  (radians, real)
%   pow     - power matrix, size:
%             channels x hiFreqs x times x trials (real)
%   nperm   - number of permutations for z-calculation
%             if 0 or omitted, returns raw MVL
%   use_mex - optional logical (default = true)
%             Controls usage of MVL and time-shift MEX kernels.
%

% checks on inputs
if nargin < 2 || isempty(phas) || isempty(pow)
    error('nf_mvlpac_within: at least phase and power matrices are required as input');
end
if nargin < 3 || isempty(nperm)
    nperm = 0;
    disp('nf_mvlpac_within: nperm = 0, returning raw MVL (no permutation z-scoring)');
end
if nargin < 4 || isempty(use_mex)
    use_mex = true;
end

if ~isscalar(nperm) || nperm < 0 || nperm ~= round(nperm)
    error('nf_mvlpac_within: nperm must be a non-negative integer scalar');
end
if ~isscalar(use_mex) || (~islogical(use_mex) && ~ismember(use_mex, [0 1]))
    error('nf_mvlpac_within: use_mex must be a logical scalar');
end
if ~isreal(phas)
    error('nf_mvlpac_within: phas must be real-valued angles in radians');
end
if ~isreal(pow)
    error('nf_mvlpac_within: pow must be real-valued');
end

% Dimension checks
nChan1  = size(phas, 1);
nChan2  = size(pow, 1);
loF     = size(phas, 2);
hiF     = size(pow, 2);
nTimes1 = size(phas, 3);
nTimes2 = size(pow, 3);
nTrial1 = size(phas, 4);
nTrial2 = size(pow, 4);

if nChan1 ~= nChan2 || nTimes1 ~= nTimes2 || nTrial1 ~= nTrial2
    error('nf_mvlpac_within: dimension mismatch - channels/times/trials must match for phase and power');
end

% checks done - you passed!

nChan  = nChan1;
nTimes = nTimes1;
nTrial = nTrial1;

% Preallocate output: chan x loF x hiF x trial
PAC = zeros(nChan, loF, hiF, nTrial);

% Progress bar
prog = 0;
fprintf(1, 'MVL progress: %3d%%', prog);

% Main computation
for f1 = 1:loF
    
    % Phase for this low frequency: [chan x 1 x time x trial]
    phase_lo = phas(:, f1, :, :);
    
    % Reshape to [time x chan x trial] for MVL kernel
    phase_tmp     = permute(phase_lo, [3, 1, 4, 2]);
    phase_samples = reshape(phase_tmp, [nTimes, nChan, nTrial]);
    
    for f2 = 1:hiF
        
        % Power for this high frequency: [chan x 1 x time x trial]
        power_hi = pow(:, f2, :, :);
        
        % [time x chan x trial]
        power_tmp     = permute(power_hi, [3, 1, 4, 2]);
        power_samples = reshape(power_tmp, [nTimes, nChan, nTrial]);
        
        % Observed MVL using MEX-aware kernel
        mvl_3d = mvl_kernel(phase_samples, power_samples, use_mex);   % [1 x nChan x nTrial]
        obsMVL = reshape(mvl_3d, [nChan, nTrial]);                    % [nChan x nTrial]
        
        if nperm > 0
            
            permutedPAC = zeros(nperm, nChan, nTrial);
            
            parfor ip = 1:nperm
                
                % circular time-shift surrogates for power, all chan×trial at once
                % first dimension = time; trailing dims = chan x trial
                power_perm = surr_timeshift_vec(power_samples, [], use_mex);    % [time x chan x trial]
                
                % MVL for this permutation
                mvl_perm_3d = mvl_kernel(phase_samples, power_perm, use_mex);   % [1 x nChan x nTrial]
                permutedPAC(ip, :, :) = reshape(mvl_perm_3d, [nChan, nTrial]);
                
            end
            
            % Null distribution stats
            mu_null = squeeze(mean(permutedPAC, 1));        % [nChan x nTrial]
            sd_null = squeeze(std(permutedPAC, 0, 1));      % [nChan x nTrial]
            
            % Guard against zero std
            zero_sd = (sd_null == 0);
            sd_null(zero_sd) = NaN;
            
            % zMVL
            zMVL = (obsMVL - mu_null) ./ sd_null;           % [nChan x nTrial]
            
            % Assign to output: chan x loF x hiF x trial
            PAC(:, f1, f2, :) = reshape(zMVL, [nChan, 1, 1, nTrial]);
            
        else
            
            % Raw MVL
            PAC(:, f1, f2, :) = reshape(obsMVL, [nChan, 1, 1, nTrial]);
            
        end
        
    end
    
    % Update progress bar
    prog = 100 * (f1 / loF);
    fprintf(1, '\b\b\b\b%3.0f%%', prog);
    
end

fprintf(1, '\n');

end












