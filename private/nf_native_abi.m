function abi = nf_native_abi()
% Octave .mex files and MATLAB MEX files are distinct ABIs. Include the
% interpreter version and architecture even on hosts sharing a network drive.
if exist('OCTAVE_VERSION','builtin')
    abi = ['Octave/' version '/' feval('__octave_config_info__','api_version')];
else
    abi = ['MATLAB/' version '/R2018a'];
end
abi = [abi '/' computer '/' mexext];
end
