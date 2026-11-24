function hnf_aux_real(X)
% Full CPU real single/double; fail before accidental device/sparse expansion.
if ~(isa(X,'double') || isa(X,'single')) || ~isreal(X) || issparse(X) || ...
        isa(X,'gpuArray') || isa(X,'distributed') || ~isnumeric(X)
    error('hnf:Type', 'Input must be a full CPU real single or double array.');
end
end
