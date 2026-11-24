function mi = hnf_mi_vec(phas, weights, nbin)
% HNF_MI_VEC  Tort's normalized KL divergence, independently per column.
% The first dimension contains time samples; trailing dimensions identify
% channels/trials. An unobserved bin has no estimable mean amplitude. This,
% nonfinite samples, negative weights or zero total weight yield NaN for
% that series. A sampled bin with mean zero contributes 0*log(0) = 0.
validateattributes(nbin, {'numeric'}, {'real','scalar','finite','integer','>=',2});
if ~isequal(size(phas),size(weights)) || ~isreal(phas) || ~isreal(weights) || isempty(phas)
    error('NeuroFreq:MIInput', 'Phase and weights must be real, nonempty arrays of equal size.');
end
shape = size(phas);
nSamples = shape(1);
shape(1) = 1;
mi = nan(shape);
if nbin > nSamples, return; end  % At least one phase bin must be empty.
phase = reshape(double(phas), nSamples, []);
weight = reshape(double(weights), nSamples, []);
good = all(isfinite(phase) & isfinite(weight) & weight >= 0, 1) & max(weight,[],1) > 0;
if ~any(good), return; end
phase = phase(:,good);
weight = weight(:,good);
% MI is scale invariant. Rescaling prevents overflow in bin sums.
weight = bsxfun(@rdivide, weight, max(weight,[],1));
% +pi and -pi represent the same direction in the half-open interval.
% For highly unwrapped angles, modulo with a rounded 2*pi loses the
% direction. libm's sine/cosine use high-precision argument reduction.
large = abs(phase) > 2^20;
phase(large) = atan2(sin(phase(large)),cos(phase(large)));
width = 2*pi/double(nbin);
bins = min(double(nbin), 1 + floor(mod(phase + pi, 2*pi) ./ width));
nVars = size(phase,2);
% MATLAB linear indexing varies the bin first, then the series.
key = bsxfun(@plus,bins,double(nbin)*(0:nVars-1));
count = reshape(accumarray(key(:),1,[nbin*nVars,1]), nbin, nVars);
total = reshape(accumarray(key(:),weight(:),[nbin*nVars,1]), nbin, nVars);
means = total ./ count;
p = bsxfun(@rdivide, means, sum(means,1));
terms = zeros(size(p));
positive = p > 0;
terms(positive) = p(positive) .* log(p(positive) .* nbin);
values = min(1,max(0,sum(terms,1) ./ log(double(nbin))));
values(any(count == 0,1)) = NaN;
mi(good) = values;
end
