function mi = hnf_mi_kernel(phas, pow, nbin, use_mex)

if nargin < 4 || isempty(use_mex)
    use_mex = true;
end

phas_d = double(phas);
pow_d  = double(pow);

if use_mex && hnf_mex_available('hnf_mi_mex')
    try
        mi = hnf_mi_mex(phas_d, pow_d, nbin);
    catch failure
        if ~hnf_mex_load_error(failure), rethrow(failure); end
        warning('mi_kernel: unable to use mex version of MI; using native MATLAB version (slower).');
        mi = hnf_mi_vec(phas_d, pow_d, nbin);
    end
else
    mi = hnf_mi_vec(phas_d, pow_d, nbin);
end

end





