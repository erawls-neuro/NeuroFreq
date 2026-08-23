function parm_spec = hnf_specparam_fit(spec, freqs, options)
% HNF_SPECPARAM_FIT Fit one spectrum using prevalidated options.
%
% This parser-free entry point returns the same full result structure as
% nf_specparam. Repeated TF callers should use
% hnf_specparam_fit_compact to avoid constructing fields they immediately
% discard.

[apFit, periodicFit, apData, periodicData, rSquare, mae, apParms, peakParms] = ...
    hnf_specparam_fit_compact(spec, freqs, options);

outputSize = size(spec);
apFit = reshape(apFit, outputSize);
periodicFit = reshape(periodicFit, outputSize);
apData = reshape(apData, outputSize);
periodicData = reshape(periodicData, outputSize);
modelFit = apFit + periodicFit;

parm_spec = struct();
parm_spec.options = options;
parm_spec.aperiodicParms = apParms;
parm_spec.peakParms = peakParms;
parm_spec.f = freqs;
parm_spec.data = spec;
parm_spec.modelFit = modelFit;
parm_spec.aPeriodicData = apData;
parm_spec.aPeriodicFit = apFit;
parm_spec.PeriodicData = periodicData;
parm_spec.PeriodicFit = periodicFit;
parm_spec.MAE = mae;
parm_spec.RSquare = rSquare;

end
