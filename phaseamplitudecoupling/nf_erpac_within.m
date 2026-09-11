function PAC = nf_erpac_within(phas, pow, nperm, use_mex)
% NF_ERPAC_WITHIN    Across-trial time-resolved phase-amplitude coupling (ERPAC)
%
% GENERAL
% -------
% Across-trial time-resolved phase-amplitude coupling (ERPAC), using
% circular-linear correlations as in Voytek et al. 2013.
%
% ERPAC is computed over trials for each time point, and for each
% channel × low-frequency × high-frequency combination.
%
% USAGE
% -----
%   PAC = nf_erpac(phas, pow)
%   PAC = nf_erpac(phas, pow, nperm)
%   PAC = nf_erpac(phas, pow, nperm, use_mex)
%
% OUTPUT:
%   PAC  - ERPAC or zERPAC values with size:
%          channels x loFreqs x hiFreqs x times
%
% INPUTS:
%   phas    - phase matrix, size:
%             channels x loFreqs x times x trials  (radians, real)
%   pow     - power matrix, size:
%             channels x hiFreqs x times x trials (real)
%   nperm   - number of permutations for z-calculation
%             if 0 or omitted, returns raw ERPAC (no z-scoring)
%   use_mex - optional logical (default = true)
%             If true and circ_corrcl_mex is available, MEX is used.
%

% Input checks and defaults
if nargin < 2
    error('nf_erpac_within: at least phase and power matrices are required as input');
end
if nargin < 3 || isempty(nperm)
    nperm = 0;
    disp('nf_erpac_within: nperm = 0, returning raw ERPAC (no permutation z-scoring)');
end
if nargin < 4 || isempty(use_mex)
    use_mex = true;
end

if ~isscalar(nperm) || nperm < 0 || nperm ~= round(nperm)
    error('nf_erpac_within: nperm must be a non-negative integer scalar');
end
if ~isscalar(use_mex) || (~islogical(use_mex) && ~ismember(use_mex, [0 1]))
    error('nf_erpac_within: use_mex must be a logical scalar');
end
if ~isreal(phas)
    error('nf_erpac_within: phas must be real-valued angles in radians');
end
if ~isreal(pow)
    error('nf_erpac_within: pow must be real-valued');
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
    error('nf_erpac_within: dimension mismatch - channels/times/trials must match for phase and power');
end

% checks done - you passed!

nChan  = nChan1;
nTimes = nTimes1;
nTrial = nTrial1;

% Preallocate output: chan x loF x hiF x time
PAC = zeros(nChan, loF, hiF, nTimes);

% Initialize progress bar
prog = 0;
fprintf(1, 'ERPAC progress: %3d%%', prog);

% Main computation
for ch = 1:nChan
    
    % Power for this channel, all hiF at once
    % pow(ch, :, :, :) : 1 x hiF x time x trials
    % After permute/reshape: [nTrial x nTimes x hiF]
    pow_ch = permute(pow(ch, :, :, :), [4, 3, 2, 1]);
    pow_ch = reshape(pow_ch, [nTrial, nTimes, hiF]);
    
    for f1 = 1:loF
        
        % Phase for this channel and low frequency
        % phas(ch, f1, :, :) → [1 x 1 x nTimes x nTrial]
        % After permute/reshape: [nTrial x nTimes]
        phase_cf = permute(phas(ch, f1, :, :), [4, 3, 2, 1]);
        phase_cf = reshape(phase_cf, [nTrial, nTimes]);
        
        % Replicate phase across hiF dimension to match power
        % phase_rep: [nTrial x nTimes x hiF]
        phase_rep = repmat(phase_cf, [1, 1, hiF]);
        
        % Observed ERPAC: single call over all time x hiF
        % clcorr_kernel expects first dim = samples (trials)
        % Returns [1 x nTimes x hiF]
        [obsPAC_3d, ~] = clcorr_kernel(phase_rep, pow_ch, use_mex);
        obsPAC = reshape(obsPAC_3d, [nTimes, hiF]);   % [nTimes x hiF]
        
        if nperm > 0
            
            % Permutation null distribution
            permutedPAC = zeros(nperm, nTimes, hiF);
            
            parfor ip = 1:nperm
                
                % Permute trials identically across time and hiF
                perm_idx = randperm(nTrial);
                pow_perm = pow_ch(perm_idx, :, :); %#ok
                
                % Null ERPAC for this permutation
                [tmp_3d, ~] = clcorr_kernel(phase_rep, pow_perm, use_mex);   % [1 x nTimes x hiF]
                tmp = reshape(tmp_3d, [nTimes, hiF]);
                permutedPAC(ip, :, :) = tmp;
                
            end
            
            % Null distribution stats
            mu_null_3d = mean(permutedPAC, 1);   % [1 x nTimes x hiF]
            sd_null_3d = std(permutedPAC, 0, 1); % [1 x nTimes x hiF]
            mu_null = reshape(mu_null_3d, [nTimes, hiF]);
            sd_null = reshape(sd_null_3d, [nTimes, hiF]);
            
            % Guard against zero std
            zero_sd = (sd_null == 0);
            sd_null(zero_sd) = NaN;
            
            % z-score ERPAC
            zPAC = (obsPAC - mu_null) ./ sd_null;   % [nTimes x hiF]
            
            % Assign to output: PAC(ch, f1, :, :) = hiF x time
            PAC(ch, f1, :, :) = permute(zPAC, [2, 1]);
            
        else
            
            % Raw ERPAC (no z-scoring)
            PAC(ch, f1, :, :) = permute(obsPAC, [2, 1]);
            
        end
        
    end
    
    % Update progress bar
    prog = 100 * (ch / nChan);
    fprintf(1, '\b\b\b\b%3.0f%%', prog);
    
end

fprintf(1, '\n');

end










