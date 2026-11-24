function loss = hnf_knee_loss(freqs,powerSpectrum,logFreqs,initial,cache)
% Choose the objective once per fit so each simplex evaluation does only
% the spectral math. The binary and M versions use the same log-domain SSE.
approved = ~isfield(cache,'approvedKneeKernel') || ...
    (~isempty(cache.approvedKneeKernel) && strcmp(which('hnf_knee_loss_mex'),cache.approvedKneeKernel));
if approved && exist('hnf_mex_available','file')==2 && hnf_mex_available('hnf_knee_loss_mex')
    try
        probe = hnf_knee_loss_mex(initial,logFreqs,powerSpectrum); %#ok<NASGU>
        loss = @(p) hnf_knee_loss_mex(p,logFreqs,powerSpectrum);
        return
    catch err
        if ~hnf_mex_load_error(err), rethrow(err); end
    end
end
loss = @(p) hnf_error_expo_function(p,freqs,powerSpectrum,logFreqs);
end
