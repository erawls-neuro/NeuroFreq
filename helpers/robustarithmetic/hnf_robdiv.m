function Z = hnf_robdiv(A, B)
% HNF_ROBDIV  Element-wise A./B with zero for invalid operands or zero B.
% The native loop fuses division and masking for large full arrays. MATLAB
% handles other types and broadcasting along individual array dimensions.
if max(numel(A),numel(B)) >= 65536 && (isa(A,'single') || isa(A,'double')) && ...
        isa(B,class(A)) && isreal(A) && isreal(B) && ~issparse(A) && ~issparse(B) && ...
        (isscalar(A) || isscalar(B) || isequal(size(A),size(B))) && ...
        exist('hnf_mex_available','file') == 2 && hnf_mex_available('hnf_robdiv_mex')
    try
        Z = hnf_robdiv_mex(A,B);
        return
    catch err
        if ~hnf_mex_load_error(err), rethrow(err); end
    end
end
Z = A ./ B;
invalid = ~isfinite(A) | ~isfinite(B) | B == 0;
Z(invalid) = 0;
end
