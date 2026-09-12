function PAC = nf_erpac_cross(phas, pow, nperm, phaseCh, ampCh, use_mex)
% NF_ERPAC_CROSS    Cross-channel ERPAC (phase from phaseCh, amplitude from ampCh).
%
%   phas    : [nChan_phase_total x nLoF x nTimes x nTrial]
%   pow     : [nChan_amp_total   x nHiF x nTimes x nTrial]
%   nperm   : number of trial permutations for z-scoring (0 = raw ERPAC)
%   phaseCh : vector of indices into phas channels (default: all)
%   ampCh   : vector of indices into pow channels  (default: all)
%   use_mex : logical, whether to use MEX kernels (default: true)
%
%   Output:
%   PAC : [nPhaseChan x nAmpChan x nLoF x nHiF x nTimes]
%         If nperm == 0  -> raw ERPAC (rho)
%         If nperm  > 0  -> z-scored ERPAC

if nargin < 3 || isempty(nperm)
    nperm = 0;
end

if nargin < 4 || isempty(phaseCh)
    phaseCh = [];
end

if nargin < 5 || isempty(ampCh)
    ampCh = [];
end

if nargin < 6 || isempty(use_mex)
    use_mex = true;
end

if ndims(phas) ~= 4
    error('nf_erpac_cross: phas must be [nChan x nLoF x nTimes x nTrial].');
end

if ndims(pow) ~= 4
    error('nf_erpac_cross: pow must be [nChan x nHiF x nTimes x nTrial].');
end

[nChanPhaseTotal, nLoF, nTimes1, nTrial1] = size(phas);
[nChanAmpTotal,   nHiF, nTimes2, nTrial2] = size(pow);

if nTimes1 ~= nTimes2
    error('nf_erpac_cross: phas and pow must share the same nTimes.');
end

if nTrial1 ~= nTrial2
    error('nf_erpac_cross: phas and pow must share the same nTrial.');
end

if isempty(phaseCh)
    phaseCh = 1:nChanPhaseTotal;
end

if isempty(ampCh)
    ampCh = 1:nChanAmpTotal;
end

nPhaseChan = numel(phaseCh);
nAmpChan   = numel(ampCh);
nTimes     = nTimes1;
nTrial     = nTrial1;

PAC = zeros(nPhaseChan, nAmpChan, nLoF, nHiF, nTimes);

% Precompute amplitude side: flatten all amp channels and hi-freqs
pow_sel = pow(ampCh, :, :, :);
pow_all = permute(pow_sel, [4, 3, 2, 1]);
nVars   = nTimes * nHiF * nAmpChan;
pow_long = reshape(pow_all, [nTrial, nVars]);

for iPhase = 1:nPhaseChan
    chP = phaseCh(iPhase);
    
    for f1 = 1:nLoF
        phase_cf = squeeze(phas(chP, f1, :, :));
        phase_cf = reshape(phase_cf, [nTimes, nTrial]);
        phase_cf = phase_cf.';
        phase_rep = repmat(phase_cf, [1, 1, nHiF, nAmpChan]);
        phase_long = reshape(phase_rep, [nTrial, nVars]);
        
        rho_flat = clcorr_kernel(phase_long, pow_long, use_mex);
        
        if nperm > 0
            sum_null  = zeros(1, nVars);
            sum2_null = zeros(1, nVars);
            
            for ip = 1:nperm
                perm_idx = randperm(nTrial);
                pow_perm = pow_long(perm_idx, :);
                rho_perm = clcorr_kernel(phase_long, pow_perm, use_mex);
                
                sum_null  = sum_null  + rho_perm;
                sum2_null = sum2_null + rho_perm .^ 2;
            end
            
            mu = sum_null ./ nperm;
            
            if nperm > 1
                var = (sum2_null - nperm .* (mu .^ 2)) ./ (nperm - 1);
                var(var < 0) = 0;
                sd = sqrt(var);
            else
                sd = zeros(1, nVars);
            end
            
            sd(sd == 0) = NaN;
            out_flat = (rho_flat - mu) ./ sd;
        else
            out_flat = rho_flat;
        end
        
        out_3d = reshape(out_flat, [nTimes, nHiF, nAmpChan]);   % [time x hiF x ampChan]
        out_3d = permute(out_3d, [3, 2, 1]);                    % [ampChan x hiF x time]
        PAC(iPhase, :, f1, :, :) = out_3d;
        
    end
end
end







