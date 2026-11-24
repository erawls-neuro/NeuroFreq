function [Y, info] = hnf_robrms(X, dim, opts)
% HNF_ROBRMS  RMS without squaring unscaled real/complex input.
%   Y = hnf_robrms(X,dim,opts); dimension defaults to first nonsingleton.
%   opts.nonfinite='error' | 'omit' | 'propagate' (NaN for affected slices).
%   opts.overflow='error' | 'allow' | 'saturate'. Empty/all-omitted => NaN.
%   Dense CPU single/double only. Complex RMS is sqrt(mean(abs(X).^2)).
if nargin < 2, dim = []; end
if nargin < 3, opts = []; end
% String flags follow the native missing-data convention: only NaNs are
% omitted by 'omitnan', and an input Inf retains its mathematical meaning.
% The default policy can use the same fused real-array reduction. Resolve
% exceptional slices with the policy implementation below to retain its
% distinct missing-data and output-overflow messages.
if isempty(opts) && nargout<=1 && numel(X)>=65536 && isreal(X) && ...
        isfloat(X) && ~issparse(X) && ~isa(X,'gpuArray') && ...
        exist('hnf_mex_available','file')==2 && hnf_mex_available('hnf_scaled_length_mex')
    Y=hnf_scaled_length(X,dim,'includenan',true);
    if all(isfinite(Y(:))), return; end
end
if ischar(opts)
    if nargout > 1
        [Y,info] = hnf_scaled_length(X,dim,opts,true);
    else
        Y = hnf_scaled_length(X,dim,opts,true);
    end
    return
end
if ~(isa(X,'double') || isa(X,'single')) || issparse(X) || ...
        isa(X,'gpuArray') || isa(X,'distributed') || ~isnumeric(X)
    error('hnf:Type','Input must be full CPU single or double.');
end
dim = hnf_aux_dim(X,dim);
o = hnf_aux_options(opts,struct('nonfinite','error','overflow','error'));
if ~any(strcmp(o.nonfinite,{'error','omit','propagate'})) || ...
        ~any(strcmp(o.overflow,{'error','allow','saturate'}))
    error('hnf:Options','Invalid RMS options.');
end
valid = isfinite(X);
if strcmp(o.nonfinite,'error') && any(~valid(:))
    error('hnf:Nonfinite','RMS input must be finite.');
end
X(~valid) = 0; count = sum(valid,dim);
if size(X,dim)==0
    Y = sum(real(X),dim); Y(:) = NaN;
else
    scale = max(max(abs(real(X)),[],dim),max(abs(imag(X)),[],dim));
    scale(scale==0) = 1;
    T = X./scale;
    % Average before restoring physical units: norm(X) may overflow when RMS fits.
    q = sum(real(T).^2+imag(T).^2,dim)./max(count,1);
    Y = scale.*sqrt(q); Y(count==0) = NaN;
end
overflow = isinf(Y);
if any(overflow(:))
    if strcmp(o.overflow,'error'), error('hnf:Overflow','RMS exceeds the output precision range.'); end
    if strcmp(o.overflow,'saturate'), Y(overflow) = realmax(class(X)); end
end
if strcmp(o.nonfinite,'propagate'), Y(any(~valid,dim)) = NaN; end
if nargout>1, info = struct('count',count,'overflow',overflow); end
end
