function mvl = hnf_mvl_kernel(phas, pow, use_mex)

if nargin < 3 || isempty(use_mex)
    use_mex = true;
end

phas_d = double(phas);
pow_d  = double(pow);

if use_mex && exist('mvl_mex', 'file') == 3
    try
        mvl = mvl_mex(phas_d, pow_d);
    catch
        warning('mvl_kernel: unable to use mex version of MVL; using native MATLAB version (slower).');
        mvl = mvl_vec(phas_d, pow_d);
    end
else
    mvl = mvl_vec(phas_d, pow_d);
end

end






