function Results = hnf_ept_TFCE_glm(data1, X, chanlocs, varargin)
% OLS with one spatial maximum distribution per design column. The design
% is used as supplied (include a ones column if an intercept is wanted).
% For each test, preserve the fitted contribution of the other columns
% and permute the reduced-model residuals across independent observations
% (Freedman-Lane). A constant column instead requires symmetric errors
% and independently sign-flips each observation's entire residual map.
[eh, count, grid, neighbours, parallel] = hnf_tfce_options(data1, chanlocs, varargin{:});
Y = double(data1(:,:));
X = double(X);
[n,k] = size(X);
if n ~= size(Y,1) || ~isreal(Y) || ~isreal(X) || any(~isfinite(Y(:))) || any(~isfinite(X(:)))
    error('NeuroFreq:InvalidGLM', 'Supply finite real observations with matching design rows.');
end
% Fit in column-equilibrated units. The permutation design, t statistics
% and tested null are invariant to these positive changes of units.
[X,xScale] = hnf_stats_scale(X);
[Y,yScale] = hnf_stats_scale(Y);
df = n - k;
if df <= 0 || rank(X) ~= k
    error('NeuroFreq:InvalidGLM', 'Coefficient tests require an estimable full-rank design and residual degrees of freedom.');
end
% QR avoids squaring the condition number as X''X would.
[Q,R] = qr(X,0);
inverse = R \ Q';
B = inverse * Y;
seFactor = hnf_robnorm(inverse,2,'includenan');
rmse = hnf_robrms(Y-X*B,1,'includenan') * sqrt(n/df);
T = B ./ (seFactor .* rmse);
dims = size(data1);
shape = [k dims(2:end)];
enhanced = zeros(size(T));
p = zeros(size(T));
null = zeros(count,k);
for j = 1:k
    obs = enhance(T(j,:), dims, neighbours, eh, grid);
    enhanced(j,:) = obs(:)';
    nuisance = X(:,[1:j-1 j+1:k]);
    if isempty(nuisance)
        fitted = zeros(size(Y));
    else
        fitted = nuisance * (nuisance \ Y);
    end
    residual = Y - fitted;
    flip = all(X(:,j) == X(1,j));
    if parallel
        parfor b = 1:count
            null(b,j) = draw(fitted,residual,flip,X,inverse,j,seFactor(j),df,dims,neighbours,eh,grid);
        end
    else
        for b = 1:count
            null(b,j) = draw(fitted,residual,flip,X,inverse,j,seFactor(j),df,dims,neighbours,eh,grid);
        end
    end
    p(j,:) = hnf_tfce_pvalues(enhanced(j,:), null(:,j));
end
B = hnf_stats_restore_beta(B,yScale,xScale);
Results = struct('Beta',reshape(B,shape),'Obs',reshape(T,shape), ...
    'TFCE_Obs',reshape(enhanced,shape),'maxTFCE',sort(null,1), ...
    'P_Values',reshape(p,shape),'df',df);
end

function maxVal = draw(fitted,residual,flip,X,inverse,j,seFactor,df,dims,neighbours,eh,grid)
n = size(X,1);
if flip
    y = fitted + residual .* (2*(rand(n,1)>=0.5)-1);
else
    y = fitted + residual(randperm(n),:);
end
b = inverse*y;
rmse = hnf_robrms(y-X*b,1,'includenan') * sqrt(n/df);
map = enhance(b(j,:) ./ (seFactor*rmse), dims, neighbours, eh, grid);
maxVal = max(abs(map(:)));
end

function map = enhance(vector,dims,neighbours,eh,grid)
map = squeeze(reshape(vector,[1 dims(2:end)]));
map = hnf_tfce_transformation(map,neighbours,eh,grid);
end
