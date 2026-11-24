function [P, logP, info] = hnf_robnormtail(Z, opts)
% HNF_ROBNORMTAIL  Stable standard-normal tail probability and optional log.
%   P = hnf_robnormtail(Z) returns the two-sided probability Pr(|N|>=|Z|).
%   P = hnf_robnormtail(Z,struct('sides',1)) returns Pr(N>=Z), preserving
%   Z's sign. Here N has a standard normal distribution. Default sides=2.
%   [P,logP] also returns the natural logarithm of that probability, using
%   ERFCX to retain finite log tails when P underflows to zero (e.g. Z=40).
%   [P,logP,info] returns sides and underflow/range counts without printing.
%
%   Z must be a finite real full CPU single/double array. Shape and class
%   are preserved. No Statistics Toolbox is required. This is a NORMAL
%   tail helper, not Student-t, a normal inverse CDF, or an inference model.
%   Probabilities below the representable range can be zero. If the true
%   log probability itself exceeds that range in magnitude, logP is -Inf.
if nargin < 2, opts = []; end
hnf_aux_real(Z);
o = hnf_aux_options(opts,struct('sides',2));
if ~isnumeric(o.sides) || ~isreal(o.sides) || ~isscalar(o.sides) || ...
        ~isfinite(o.sides) || ~any(o.sides == [1 2])
    error('hnf:Options','sides must be 1 (signed upper tail) or 2 (two-sided).');
end
if any(~isfinite(Z(:)))
    error('hnf:Nonfinite','Normal test statistics must be finite.');
end
x = abs(Z) ./ cast(sqrt(2),'like',Z);
P = erfc(x);
if o.sides == 1
    P = P .* cast(0.5,'like',Z);
    negative = Z < 0;
    P(negative) = 1 - P(negative);
end
if nargout > 1
    logP = log(erfcx(x)) - x .* x;
    % ERFCX rounds to one very near zero; use the small complement there.
    small = x < cast(0.5,'like',Z);
    logP(small) = log1p(-erf(x(small)));
    if o.sides == 1
        logP = logP - cast(log(2),'like',Z);
        % For a negative Z the upper tail is close to one. LOG1P preserves
        % its small negative log without subtracting its probability from 1.
        logP(negative) = log1p(-exp(logP(negative)));
    end
end
if nargout > 2
    info = struct('sides',o.sides,'probabilityUnderflowCount',nnz(P == 0), ...
        'logOutOfRangeCount',nnz(isinf(logP)), ...
        'roundedToOneCount',nnz(P == 1 & Z ~= 0));
end
end
