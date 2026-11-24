function PAC = nf_mipac_cross(phas, pow, nperm, nbin, phaseCh, ampCh, use_mex)
% NF_MIPAC_CROSS    Cross-channel MI-based PAC (phase from phaseCh, amplitude from ampCh).
%
%   phas    : [nChan_phase_total x nLoF x nTimes x nTrial]
%   pow     : [nChan_amp_total   x nHiF x nTimes x nTrial]
%   nperm   : number of time-shift surrogates (0 = raw MI)
%   nbin    : phase bin count (default: 18)
%   phaseCh : vector of indices into phas channels (default: all)
%   ampCh   : vector of indices into pow channels  (default: all)
%   use_mex : logical, use MI MEX and timeshift MEX if available
%
%   Output:
%   PAC : [nPhaseChan x nAmpChan x nLoF x nHiF x nTrial]
%         If nperm == 0  -> raw MI
%         If nperm  > 0  -> z-scored MI

if nargin < 3 || isempty(nperm)
    nperm = 0;
end

if nargin < 4 || isempty(nbin)
    nbin = 18;
end

if nargin < 5 || isempty(phaseCh)
    phaseCh = [];
end

if nargin < 6 || isempty(ampCh)
    ampCh = [];
end

if nargin < 7 || isempty(use_mex)
    use_mex = true;
end

validateattributes(nperm, {'numeric'}, {'real','scalar','finite','integer','nonnegative'});
validateattributes(nbin, {'numeric'}, {'real','scalar','finite','integer','>=',2});
validateattributes(use_mex, {'numeric','logical'}, {'real','scalar','binary'});

if ndims(phas) > 4 || isempty(phas) || ~isreal(phas)
    error('nf_mipac_cross: phas must be [nChan x nLoF x nTimes x nTrial].');
end

if ndims(pow) > 4 || isempty(pow) || ~isreal(pow)
    error('nf_mipac_cross: pow must be [nChan x nHiF x nTimes x nTrial].');
end

[nChanPhaseTotal, nLoF, nTimes1, nTrial1] = size(phas);
[nChanAmpTotal,   nHiF, nTimes2, nTrial2] = size(pow);

if nTimes1 ~= nTimes2
    error('nf_mipac_cross: phas and pow must share the same nTimes.');
end

if nTrial1 ~= nTrial2
    error('nf_mipac_cross: phas and pow must share the same nTrial.');
end

if isempty(phaseCh)
    phaseCh = 1:nChanPhaseTotal;
end

if isempty(ampCh)
    ampCh = 1:nChanAmpTotal;
end

validateattributes(phaseCh, {'numeric'}, {'real','vector','integer','positive','<=',nChanPhaseTotal});
validateattributes(ampCh, {'numeric'}, {'real','vector','integer','positive','<=',nChanAmpTotal});
nPhaseChan = numel(phaseCh);
nAmpChan   = numel(ampCh);
nTimes     = nTimes1;
nTrial     = nTrial1;

PAC = zeros(nPhaseChan, nAmpChan, nLoF, nHiF, nTrial);

for iPhase = 1:nPhaseChan
    chP = phaseCh(iPhase);
    
    for iAmp = 1:nAmpChan
        chA = ampCh(iAmp);
        
        for f1 = 1:nLoF
            phase_ts = squeeze(phas(chP, f1, :, :));
            phase_ts = reshape(phase_ts, [nTimes, nTrial]);
            phase_samples = reshape(phase_ts, [nTimes, 1, nTrial]);
            
            for f2 = 1:nHiF
                pow_ts = squeeze(pow(chA, f2, :, :));
                pow_ts = reshape(pow_ts, [nTimes, nTrial]);
                power_samples = reshape(pow_ts, [nTimes, 1, nTrial]);
                
                mi_obs_3d = hnf_mi_kernel(phase_samples, power_samples, nbin, use_mex);
                mi_obs = reshape(mi_obs_3d, [1, nTrial]);
                
                out_vals = mi_obs;
                if nperm > 0
                    out_vals = hnf_pac_permutation_z(phase_samples, power_samples, mi_obs, ...
                        nperm, 'mi', nbin, use_mex, true);
                end

                PAC(iPhase, iAmp, f1, f2, :) = reshape(out_vals, [1, 1, 1, 1, nTrial]);
            end
        end
    end
end
end








