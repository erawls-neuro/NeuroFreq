function mvl = hnf_mvl_kernel(phas, pow, use_mex)

if nargin < 3 || isempty(use_mex)
    use_mex = true;
end

phas_d = double(phas);
pow_d  = double(pow);

if use_mex && hnf_mex_available('hnf_mvl_mex')
    try
        mvl = hnf_mvl_mex(phas_d, pow_d);
    catch failure
        if ~hnf_mex_load_error(failure), rethrow(failure); end
        warning('mvl_kernel: unable to use mex version of MVL; using native MATLAB version (slower).');
        mvl = hnf_mvl_vec(phas_d, pow_d);
    end
else
    mvl = hnf_mvl_vec(phas_d, pow_d);
end

end






