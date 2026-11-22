function result = hnf_mex_available(name,policy)
% Test a compiled NF kernel before dispatch. Direct function calls respect
% MATLAB's normal path. Within nf_run only validated binaries are eligible.
% nf_run scopes this policy with onCleanup, including nested calls/errors.
persistent active
if nargin==2
    if ~strcmp(name,'scope'), error('NeuroFreq:NativeScope','Expected scope.'); end
    result = active;
    active = policy;
    return
end
if isstruct(active)
    result = isfield(active,name) && strcmp(which(name),active.(name));
else
    result = exist(name,'file')==3;
end
end
