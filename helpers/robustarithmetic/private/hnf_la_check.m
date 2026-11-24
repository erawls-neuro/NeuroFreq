function hnf_la_check(A, allowMissing)
% Internal: base MATLAB dense CPU matrices only; never gather/densify.
if nargin < 2, allowMissing = false; end
if ~(isa(A,'double') || isa(A,'single')) || ~ismatrix(A) || issparse(A) ...
        || isa(A,'gpuArray') || isa(A,'distributed') || isa(A,'tall')
    error('hnf:UnsupportedInput','Expected a full CPU single/double matrix.');
end
if ~allowMissing && any(~isfinite(A(:)))
    error('hnf:NonfiniteInput','Matrix contains NaN or Inf.');
end
end
