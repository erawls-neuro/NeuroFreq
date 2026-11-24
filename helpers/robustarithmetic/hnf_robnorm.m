function [y, info] = hnf_robnorm(x, dim, opts)
% HNF_ROBNORM  Scaled Euclidean vector norm along one dimension.
%   Y = HNF_ROBNORM(X,DIM,OPTS), DIM defaults to first nonsingleton axis.
%   This computes vector 2-norms, not matrix spectral norms.
%   Options: nanPolicy='propagate' (default)|'omit',
%   inputPolicy='propagate'|'fill'|'error', fillValue=0,
%   overflow='error' (default)|'saturate'|'allow'.
%   Scaling real/imaginary components avoids intermediate |X|^2 overflow
%   and underflow. A mathematically unrepresentable norm obeys overflow.
%   An input Inf produces Inf; NaN dominates unless omitted. Empty or
%   all-omitted slices have norm zero. Full CPU single/double supported.
id = 'hnf_robnorm';
if nargin < 2, dim = []; end
if nargin < 3, opts = []; end
% String flags follow the native missing-data convention: only NaNs are
% omitted by 'omitnan', and an input Inf retains its mathematical meaning.
% The default policy can use the same fused real-array reduction. Resolve
% exceptional slices with the policy implementation below to retain its
% distinct missing-data and output-overflow messages.
if isempty(opts) && nargout<=1 && numel(x)>=65536 && isreal(x) && ...
        isfloat(x) && ~issparse(x) && ~isa(x,'gpuArray') && ...
        exist('hnf_mex_available','file')==2 && hnf_mex_available('hnf_scaled_length_mex')
    y=hnf_scaled_length(x,dim,'includenan',false);
    if all(isfinite(y(:))), return; end
end
if ischar(opts)
    if nargout > 1
        [y,info] = hnf_scaled_length(x,dim,opts,false);
    else
        y = hnf_scaled_length(x,dim,opts,false);
    end
    return
end
elem_float(x,id); dim = elem_dim(x,dim,id);
o = elem_options(struct('nanPolicy','propagate','inputPolicy','propagate', ...
    'fillValue',0,'overflow','error','underflow','allow'),opts,id);
o.nanPolicy = elem_choice(o.nanPolicy,{'propagate','omit'},id);
if ~strcmp(o.underflow,'allow')
    error([id ':Options'],'A Euclidean norm cannot underflow below a nonzero component.');
end
if nargout > 1, inputCount = nnz(~isfinite(x)); end
x = elem_input(x,o,id);
if strcmp(o.nanPolicy,'omit'), x(isnan(x)) = 0; end
if size(x,dim) == 0
    y = sum(real(x),dim);
else
    [scale,root,nanSlice,infSlice] = elem_normparts(x,dim);
    y = scale.*root;
    y(infSlice & ~nanSlice) = Inf; y(nanSlice) = NaN;
    y = elem_range(y,~nanSlice & ~infSlice,false(size(y)),o,id);
end
if nargout > 1
    info = struct('method','scaled-euclidean','dimension',dim, ...
        'nonfiniteInput',inputCount,'nonfiniteOutput',nnz(~isfinite(y)));
end
end
