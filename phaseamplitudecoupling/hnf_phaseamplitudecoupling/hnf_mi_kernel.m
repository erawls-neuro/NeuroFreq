function mi = hnf_mi_kernel(phas, pow, nbin, use_mex)

if nargin < 4 || isempty(use_mex)
    use_mex = true;
end

phas_d = double(phas);
pow_d  = double(pow);

if use_mex && exist('mi_mex', 'file') == 3
    try
        mi = mi_mex(phas_d, pow_d, nbin);
    catch
        warning('mi_kernel: unable to use mex version of MI; using native MATLAB version (slower).');
        mi = mi_vec(phas_d, pow_d, nbin);
    end
else
    mi = mi_vec(phas_d, pow_d, nbin);
end

end





