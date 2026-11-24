function file = nf_runtime_file(root, forWriting)
% Keep package and MEX selections separate when both runtimes share a checkout.
% Existing configurations keep working until the first Octave installation.
if nargin < 2, forWriting = false; end
cache = fullfile(root,'dependencies','_nf_local');
file = fullfile(cache,'runtime.mat');
if exist('OCTAVE_VERSION','builtin')
    octaveFile = fullfile(cache,'runtime-octave.mat');
    if forWriting || isfile(octaveFile), file = octaveFile; end
end
end
