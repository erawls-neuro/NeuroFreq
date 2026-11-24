function [rho, pval] = hnf_circ_corrcl_vec(alpha, x, orders)
% HNF_CIRC_CORRCL_VEC Circular-linear multiple correlation across axis 1.
% Regress the centered linear variable onto centered sin(alpha), cos(alpha).
% With optional trial orders [samples, permutations], x is [samples,
% times*blocks], alpha is [samples,times], and rho is [permutations,
% times*blocks]. All blocks share phase and the same trial permutation.
permuted = nargin > 2;
if isempty(alpha) || isempty(x) || ~isreal(alpha) || ~isreal(x)
    error('NeuroFreq:CircInput','Phase and linear observations must be nonempty and real.');
end
shape = size(x); samples = shape(1);
if permuted
    if ~ismatrix(alpha) || ~ismatrix(x) || size(alpha,1) ~= samples || ...
            mod(size(x,2),size(alpha,2)) ~= 0 || size(orders,1) ~= samples
        error('NeuroFreq:CircInput','Permuted observations must share trials and contain whole time blocks.');
    end
    shape = [size(orders,2),size(x,2)];
else
    if ~isequal(size(alpha),size(x))
        error('NeuroFreq:CircInput','Phase and linear observations must have equal size.');
    end
    shape(1) = 1;
end
rho = nan(shape);
if nargout > 1, pval = nan(shape); end
if samples < 3, return; end
A = reshape(double(alpha),samples,[]);
X = reshape(double(x),samples,[]);
% Scaling is independent of trial order. Centering below follows the order
% of the observations, preserving the single-permutation arithmetic.
scale = max(abs(X),[],1);
scale(scale == 0) = 1;
X = bsxfun(@rdivide,X,scale);
S = center(sin(A));
C = center(cos(A));
sn = sqrt(sum(S.^2,1));
cn = sqrt(sum(C.^2,1));
swap = cn > sn;
tmp = S(:,swap); S(:,swap) = C(:,swap); C(:,swap) = tmp;
firstNorm = max(sn,cn);
U = bsxfun(@rdivide,S,firstNorm);
residual = C - bsxfun(@times,U,sum(U.*C,1));
secondNorm = sqrt(sum(residual.^2,1));
V = bsxfun(@rdivide,residual,secondNorm);
phaseValid = all(isfinite(A),1) & firstNorm > 0 & secondNorm > samples*eps*firstNorm;
blocks = size(X,2)/size(A,2);
if blocks > 1
    U = repmat(U,1,blocks); V = repmat(V,1,blocks);
    phaseValid = repmat(phaseValid,1,blocks);
end
if permuted
    for p = 1:size(orders,2)
        rho(p,:) = project(X(orders(:,p),:),U,V,phaseValid);
    end
else
    rho = reshape(project(X,U,V,phaseValid),shape);
end
if nargout > 1
    % CircStat's large-sample chi-square(2) approximation.
    pval = exp(-0.5*samples*rho.^2);
end
end

function rho = project(X,U,V,phaseValid)
X = center(X);
xn = sqrt(sum(X.^2,1));
Y = bsxfun(@rdivide,X,xn);
r = hypot(sum(Y.*U,1),sum(Y.*V,1));
valid = phaseValid & all(isfinite(X),1) & xn > 0;
rho = nan(1,size(X,2));
rho(valid) = min(1,r(valid));
end

function X = center(X)
% Subtract the first observation before centering. Constant inputs then
% have exactly zero variance even when a sum of identical floats rounds.
X = bsxfun(@minus,X,X(1,:));
X = bsxfun(@minus,X,mean(X,1));
end
