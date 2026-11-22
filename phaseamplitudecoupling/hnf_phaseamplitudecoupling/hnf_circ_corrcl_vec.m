function [rho, pval] = hnf_circ_corrcl_vec(alpha, x)
% HNF_CIRC_CORRCL_VEC Circular-linear multiple correlation across axis 1.
% Regress the centered linear variable onto centered sin(alpha), cos(alpha).
% Orthogonalize the second predictor so the two squared projections sum
% to R^2 without dividing by (1-r_cs^2) near a singular design.
if ~isequal(size(alpha),size(x)) || ~isreal(alpha) || ~isreal(x) || isempty(alpha)
    error('NeuroFreq:CircInput','Phase and linear observations must be real nonempty arrays of equal size.');
end
shape = size(alpha);
samples = shape(1);
shape(1) = 1;
rho = nan(shape);
pval = nan(shape);
if samples < 3, return; end
A = reshape(double(alpha),samples,[]);
X = reshape(double(x),samples,[]);
S = sin(A);
C = cos(A);
% Subtract the first observation before centering. Constant inputs then
% have exactly zero variance even when a sum of identical floats rounds.
scale = max(abs(X),[],1);
scale(scale == 0) = 1;
X = bsxfun(@rdivide,X,scale);
X = center(X);
S = center(S);
C = center(C);
sn = sqrt(sum(S.^2,1));
cn = sqrt(sum(C.^2,1));
% Start with the larger phase predictor for a stable rank decision.
swap = cn > sn;
tmp = S(:,swap); S(:,swap) = C(:,swap); C(:,swap) = tmp;
firstNorm = max(sn,cn);
U = bsxfun(@rdivide, S, firstNorm);
residual = C - bsxfun(@times, U, sum(U.*C,1));
secondNorm = sqrt(sum(residual.^2,1));
V = bsxfun(@rdivide, residual, secondNorm);
xn = sqrt(sum(X.^2,1));
Y = bsxfun(@rdivide,X,xn);
r = hypot(sum(Y.*U,1),sum(Y.*V,1));
% Rank-two and nonconstant series are required. Retain NaNs explicitly.
valid = all(isfinite(A) & isfinite(X),1) & xn > 0 & ...
    firstNorm > 0 & secondNorm > samples*eps*firstNorm;
rho(valid) = min(1,r(valid));
% CircStat's large-sample chi-square(2) approximation, not a permutation p.
pval(valid) = exp(-0.5 * samples * rho(valid).^2);
end

function X = center(X)
X = bsxfun(@minus,X,X(1,:));
X = bsxfun(@minus,X,mean(X,1));
end
