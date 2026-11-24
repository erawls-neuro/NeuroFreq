function method = hnf_adjust_backend(requested,announce)
% The official MADE extractor dereferences MATLAB LinearModel diagnostics.
% Octave statistics fitlm has a different interface. Select the NF
% developmental ADJUST implementation and identify its actual provenance.
if nargin<2, announce=true; end
method = requested;
if exist('OCTAVE_VERSION','builtin') && strcmp(requested,'adjustedadjust')
    method = 'nfadjustedadjust';
    if ~announce, return; end
    warning('NeuroFreq:OctaveAdjustedADJUST', ...
        ['MADE adjusted_ADJUST requires MATLAB model diagnostics. Octave uses ' ...
         'NeuroFreq hnf_ADJUST (nfadjustedadjust). Its spectral rescue and ' ...
         'component decisions can differ from the original MADE classifier.']);
end
end
