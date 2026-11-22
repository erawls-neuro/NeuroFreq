function mi = hnf_mi_kernel(phas, pow, nbin, use_mex)

if nargin < 4 || isempty(use_mex)
    use_mex = true;
end

phas_d = double(phas);
pow_d  = double(pow);

if use_mex && exist('hnf_mi_mex', 'file') == 3
    try
        mi = hnf_mi_mex(phas_d, pow_d, nbin);
    catch
        warning('mi_kernel: unable to use mex version of MI; using native MATLAB version (slower).');
        mi = hnf_mi_vec(phas_d, pow_d, nbin);
    end
else
    mi = hnf_mi_vec(phas_d, pow_d, nbin);
end

end





