function [y, info] = hnf_robsqrt(x, opts)
% HNF_ROBSQRT  Square root with explicit handling of negative values.
%   Options: domainPolicy='nan'|'clip'|'complex'|'error' (default 'nan'),
%   negativeTolerance=0, inputPolicy='propagate'|'fill'|'error', fillValue=0.
%   Finite real values in [-negativeTolerance,0) become zero. Larger
%   negatives become NaN, all become zero in clip mode, or throw in error
%   mode. Tolerance is absolute and has X's units. Complex mode preserves
%   MATLAB's principal square root and requires no tolerance or clipping.
id = 'hnf_robsqrt';
if nargin < 2, opts = []; end
elem_float(x,id);
o = elem_options(struct('domainPolicy','nan','negativeTolerance',0, ...
    'inputPolicy','propagate','fillValue',0),opts,id);
o.domainPolicy = elem_choice(o.domainPolicy,{'nan','clip','complex','error'},id);
t = o.negativeTolerance;
if ~isnumeric(t) || ~isscalar(t) || ~isreal(t) || ~isfinite(t) || t < 0
    error([id ':Options'],'negativeTolerance must be a finite nonnegative scalar.');
end
if nargout > 1, inputCount = nnz(~isfinite(x)); end
x = elem_input(x,o,id);
if strcmp(o.domainPolicy,'complex')
    y = sqrt(complex(x)); count = 0;
else
    if ~isreal(x), error([id ':Domain'],'Complex input requires complex domainPolicy.'); end
    tiny = isfinite(x) & x < 0 & x >= -t;
    x(tiny) = 0; bad = x < 0;
    if strcmp(o.domainPolicy,'error') && any(bad(:))
        error([id ':Domain'],'Negative input exceeds negativeTolerance.');
    elseif strcmp(o.domainPolicy,'clip')
        finiteBad = bad & isfinite(x); x(finiteBad) = 0; x(bad & ~finiteBad) = NaN;
    else
        x(bad) = NaN;
    end
    y = sqrt(x);
    if nargout > 1, count = nnz(tiny | bad); end
end
if nargout > 1
    info = struct('modifiedDomain',count,'nonfiniteInput',inputCount, ...
        'nonfiniteOutput',nnz(~isfinite(y)),'domainPolicy',o.domainPolicy);
end
end
