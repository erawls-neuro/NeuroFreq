function [rho, pval] = hnf_clcorr_kernel(alpha, x, use_mex)

if nargin < 3 || isempty(use_mex)
    use_mex = true;
end

alpha_d = double(alpha);
x_d     = double(x);

if use_mex && hnf_mex_available('hnf_circ_corrcl_mex')
    try
        [rho, pval] = hnf_circ_corrcl_mex(alpha_d, x_d);
    catch failure
        if ~hnf_mex_load_error(failure), rethrow(failure); end
        warning('clcorr_kernel: unable to use mex version of circ_corrcl; using native MATLAB version (slower).');
        [rho, pval] = hnf_circ_corrcl_vec(alpha_d, x_d);
    end
else
    [rho, pval] = hnf_circ_corrcl_vec(alpha_d, x_d);
end

end




