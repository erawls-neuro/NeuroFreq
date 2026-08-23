function PAC = nf_mipac_within(phas, pow, nperm, nbin, use_mex)
% NF_MIPAC_WITHIN    Phase-amplitude coupling via Tort Modulation Index (MI)
%
% GENERAL
% -------
% Averages phase-amplitude coupling using the Tort 2008 modulation index,
% optionally z-scored via permutation.
%
% MI is computed over time samples for each channel × low-freq ×
% high-freq × trial combination.
%
% USAGE
% -----
%   PAC = nf_mipac(phas, pow)
%   PAC = nf_mipac(phas, pow, nperm)
%   PAC = nf_mipac(phas, pow, nperm, nbin)
%   PAC = nf_mipac(phas, pow, nperm, nbin, use_mex)
%
% OUTPUT:
%   PAC  - MI or zMI values, size:
%          channels x loFreqs x hiFreqs x trials
%
% INPUTS:
%   phas    - phase matrix, size:
%             channels x loFreqs x times x trials  (radians, real)
%   pow     - power matrix, size:
%             channels x hiFreqs x times x trials (real)
%   nperm   - number of permutations for z-calculation
%             if 0 or omitted, returns raw MI
%   nbin    - number of phase bins (default 18)
%   use_mex - optional logical (default = true)
%             Controls usage of MI and time-shift MEX kernels.
%

% Basic checks and defaults
if nargin < 2 || isempty(phas) || isempty(pow)
    error('nf_mipac_within: at least phase and power matrices are required as input');
end
if nargin < 3 || isempty(nperm)
    nperm = 0;
    disp('nf_mipac_within: nperm = 0, returning raw MI (no permutation z-scoring)');
end
if nargin < 4 || isempty(nbin)
    nbin = 18;
    disp('nf_mipac_within: using 18 phase bins (default)');
end
if nargin < 5 || isempty(use_mex)
    use_mex = true;
end

if ~isscalar(nperm) || nperm < 0 || nperm ~= round(nperm)
    error('nf_mipac_within: nperm must be a non-negative integer scalar');
end
if ~isscalar(nbin) || nbin < 1 || nbin ~= round(nbin)
    error('nf_mipac_within: nbin must be a positive integer scalar');
end
if ~isscalar(use_mex) || (~islogical(use_mex) && ~ismember(use_mex, [0 1]))
    error('nf_mipac_within: use_mex must be a logical scalar');
end
if ~isreal(phas)
    error('nf_mipac_within: phas must be real-valued angles in radians');
end
if ~isreal(pow)
    error('nf_mipac_within: pow must be real-valued');
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
    error('nf_mipac_within: dimension mismatch - channels/times/trials must match for phase and power');
end

% checks done - you passed!

nChan  = nChan1;
nTimes = nTimes1;
nTrial = nTrial1;

% Preallocate output: chan x loF x hiF x trial
PAC = zeros(nChan, loF, hiF, nTrial);

% Progress bar
prog = 0;
fprintf(1, 'MI progress: %3d%%', prog);

% Main computation
for f1 = 1:loF
    
    % Phase for this low frequency: [chan x 1 x time x trial]
    phase_lo = phas(:, f1, :, :);
    
    % Reshape to [time x chan x trial] for MI kernel
    phase_tmp     = permute(phase_lo, [3, 1, 4, 2]);
    phase_samples = reshape(phase_tmp, [nTimes, nChan, nTrial]);
    
    for f2 = 1:hiF
        
        % Power for this high frequency: [chan x 1 x time x trial]
        power_hi = pow(:, f2, :, :);
        
        % [time x chan x trial]
        power_tmp     = permute(power_hi, [3, 1, 4, 2]);
        power_samples = reshape(power_tmp, [nTimes, nChan, nTrial]);
        
        % Observed MI using MEX-aware kernel
        mi_3d = mi_kernel(phase_samples, power_samples, nbin, use_mex);   % [1 x nChan x nTrial]
        obsMI = reshape(mi_3d, [nChan, nTrial]);                          % [nChan x nTrial]
        
        if nperm > 0
            
            permutedMI = zeros(nperm, nChan, nTrial);
            
            parfor ip = 1:nperm
                
                % circular time-shift surrogates for power, all chan×trial at once
                % first dimension = time; trailing dims = chan x trial
                power_perm = surr_timeshift_vec(power_samples, [], use_mex);    % [time x chan x trial]
                
                % MI for this permutation
                mi_perm_3d = mi_kernel(phase_samples, power_perm, nbin, use_mex);  % [1 x nChan x nTrial]
                permutedMI(ip, :, :) = reshape(mi_perm_3d, [nChan, nTrial]);
                
            end
            
            % Null distribution stats
            mu_null = squeeze(mean(permutedMI, 1));           % [nChan x nTrial]
            sd_null = squeeze(std(permutedMI, 0, 1));         % [nChan x nTrial]
            
            % Guard against zero std
            zero_sd = (sd_null == 0);
            sd_null(zero_sd) = NaN;
            
            % zMI
            zMI = (obsMI - mu_null) ./ sd_null;               % [nChan x nTrial]
            
            % Assign to output: chan x loF x hiF x trial
            PAC(:, f1, f2, :) = reshape(zMI, [nChan, 1, 1, nTrial]);
            
        else
            
            % Raw MI
            PAC(:, f1, f2, :) = reshape(obsMI, [nChan, 1, 1, nTrial]);
            
        end
        
    end
    
    % Update progress bar
    prog = 100 * (f1 / loF);
    fprintf(1, '\b\b\b\b%3.0f%%', prog);
    
end

fprintf(1, '\n');

end









'