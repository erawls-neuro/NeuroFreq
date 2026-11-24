function [power,phase] = hnf_powphase(coefficients)
% Double power permits squared single coefficients outside the single range.
% For larger host arrays the compiled kernel produces both fields in one pass.
if ~isa(coefficients,'gpuArray') && numel(coefficients)>=4096 && ...
        exist('hnf_mex_available','file')==2 && hnf_mex_available('hnf_powphase_mex')
    try
        [power,phase] = hnf_powphase_mex(coefficients);
        return
    catch err
        if ~hnf_mex_load_error(err), rethrow(err); end
    end
end
power = abs(double(coefficients)).^2;
phase = angle(coefficients);
end
