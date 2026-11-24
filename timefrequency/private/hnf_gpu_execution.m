function [useGPU, device, unavailable] = hnf_gpu_execution(input, environment)
% Select a backend without initializing the GPU for the default host path.
if ~(ischar(environment) || (isstring(environment) && isscalar(environment))) || ...
        ~any(strcmpi(environment,{'auto','cpu','gpu'}))
    error('NeuroFreq:InvalidExecutionEnvironment','Use auto, cpu or gpu.');
end
resident = isa(input,'gpuArray');
device = [];
unavailable = '';
useGPU = strcmpi(environment,'gpu') || (resident && strcmpi(environment,'auto'));
if useGPU && exist('OCTAVE_VERSION','builtin')
    useGPU = false;
    unavailable = 'Octave uses the CPU FFT backend.';
elseif useGPU
    try
        device = gpuDevice();
    catch err
        if resident, rethrow(err); end
        useGPU = false;
        unavailable = err.message;
    end
end
end
