function [y, info] = hnf_robclip(x, limits, opts)
% HNF_ROBCLIP  Explicit real amplitude clipping; no automatic threshold choice.
%   Y = HNF_ROBCLIP(X,[LOW HIGH],OPTS), LOW<=HIGH (infinite bounds allowed).
%   opts.nonfinite='error' (default) or 'propagate'. No imputation/omission.
%   opts.nonfinite='omit' is rejected because clipping preserves array size.
%   Bounds must be representable in X's class; NaN bounds are invalid.
%   INFO records the number of finite entries clipped and invalid entries.
if nargin < 3, opts = struct; end
opts = stat_options(opts,struct('nonfinite','error'));
if strcmp(opts.nonfinite,'omit')
    error('hnf:statistics:OptionValue','Clipping supports nonfinite error or propagate.');
end
if ~any(strcmp(class(x),{'single','double'})) || ~isreal(x) || issparse(x)
    error('hnf:statistics:Input','X must be real dense single or double.');
end
if ~isnumeric(limits) || ~isreal(limits) || numel(limits) ~= 2 || ...
        any(isnan(limits(:))) || limits(1) > limits(2)
    error('hnf:statistics:Limits','limits must be [LOW HIGH], LOW<=HIGH, without NaN.');
end
bounds = cast(limits,'like',x);
if any(isfinite(limits(:)) & ~isfinite(bounds(:)))
    error('hnf:statistics:Limits','A finite clipping bound exceeds the input class.');
end
finite = isfinite(x);
if strcmp(opts.nonfinite,'error') && any(~finite(:))
    error('hnf:statistics:Nonfinite','X contains NaN or Inf.');
end
y = max(bounds(1),min(bounds(2),x));
y(~finite) = x(~finite);
if nargout > 1
    info = struct('clippedCount',nnz(finite & (x < bounds(1) | x > bounds(2))), ...
        'invalidCount',nnz(~finite),'limits',bounds,'nonfinite',opts.nonfinite);
end
end
