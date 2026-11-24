function test_hnf_statistics
% TEST_HNF_STATISTICS  Base-MATLAB assertions; invoke after addpath(module).
% These cases check estimator contracts and difficult numerical regimes.
x = [1 4 7; 2 5 8; 3 6 9];
assertClose(hnf_robmean(x),mean(x),1e-14);
assertClose(hnf_robstd(x),std(x),1e-14);
assertClose(hnf_robvar(x),var(x),1e-14);
assertClose(hnf_robstd(x,2,struct('normalization',1)),std(x,1,2),1e-14);
assertClose(hnf_robmean(reshape(1:24,[2 3 4]),3),mean(reshape(1:24,[2 3 4]),3),1e-14);
assert(isequal(size(hnf_robmean(x,4)),size(x)));
assert(isa(hnf_robmean(single(x)),'single'));
assert(isa(hnf_robstd(single(x)),'single'));
assert(hnf_robstd(2) == 0);
assert(all(isnan(hnf_robmean(zeros(0,3),1))));
assert(all(isnan(hnf_robstd(zeros(0,3),1))));
assert(isequal(size(hnf_robmean(zeros(2,0),1)),[1 0]));

% Finite means and population scales even when unscaled sums/ranges overflow.
big = realmax;
assert(hnf_robmean([big big],2) == big);
assertClose(hnf_robmean([-big big],2),0,0);
assertClose(hnf_robmean([1e16 1 -1e16],2),1/3,1e-15);
[cz,cm] = hnf_robzscore([1e16 1 -1e16],2);
assertClose(cm,1/3,1e-15); assert(cz(2) > 0);
assertClose(hnf_robstd([-big big],2,struct('normalization',1))/big,1,2e-15);
assertThrows(@() hnf_robstd([-big big],2),'hnf:statistics:Overflow');
[s,si] = hnf_robstd([-big big],2,struct('overflow','saturate'));
assert(s == big && si.overflow);
assertThrows(@() hnf_robvar([-1e200 1e200],2),'hnf:statistics:Overflow');
assertClose(hnf_robstd(1e12+[0 1 2],2),1,1e-14);
assertClose(hnf_robstd([-1e-200 1e-200],2,struct('normalization',1))/1e-200,1,1e-14);

% Ordered estimators, explicit interpolation, robust degeneracy, bounded IRLS.
assertClose(hnf_robquantile([0 10 20 30],[0 .25 .5 .75 1],2),[0 7.5 15 22.5 30],1e-14);
assertClose(hnf_robquantile([-big big],[0 .5 1],2),[-big 0 big],0);
q = hnf_robquantile(x,[.25 .75],1);
assert(isequal(size(q),[2 3])); assertClose(q,[1.5 4.5 7.5;2.5 5.5 8.5],1e-14);
assert(hnf_robmad([1 2 3],2) == 1);
% An extreme first outlier must not erase ordered central values.
ordered = [1e308 1e-100 2e-100];
assert(hnf_robmean(ordered,2,struct('method','median')) == 2e-100);
assert(hnf_robquantile(ordered,.5,2) == 2e-100);
assert(hnf_robmad(ordered,2) == 1e-100);
assert(hnf_robmean([-1 0 1 1e308],2,struct('method','median')) == .5);
tinyCenter = [1e308 -1e-100 1e-100 2e-100 3e-100];
tm = hnf_robmean(tinyCenter,2,struct('method','trim','trim',.2));
assert(abs(tm/1e-100-2) < 1e-14);
[hm,hi] = hnf_robmean(ordered,2,struct('method','huber'));
assert(isfinite(hm) && hm > 1e-100 && hm < 10e-100 && hi.converged);
[~,ml,ms] = hnf_robzscore([-1 0 1 1e308],2,struct('method','mad'));
assert(ml == .5 && abs(ms-1.482602218505602) < 1e-14);
assertClose(hnf_robstd([1 2 3],2,struct('method','mad')),1.482602218505602,1e-14);
assert(hnf_robmean([0 0 0 0 100],2,struct('method','median')) == 0);
assert(hnf_robmean([0 0 0 0 100],2,struct('method','trim','trim',.2)) == 0);
assertClose(hnf_robmean([0 0 0 0 100],2,struct('method','winsor','trim',.2)),4,1e-13);
[mu,mi] = hnf_robmean([0 0 0 0 100],2,struct('method','huber'));
assert(mu == 0 && mi.zeroScale && mi.iterations == 0);
u = [-3 -2 -1 0 1 2 100];
[mu,mi] = hnf_robmean(u,2,struct('method','huber','maxiter',1,'tol',0));
assert(isfinite(mu) && mi.iterations <= 1);
[a,ai] = hnf_robmean(u,2,struct('method','huber','tol',1e-10));
[b,bi] = hnf_robmean(u+1e12,2,struct('method','huber','tol',1e-10));
assert(abs((b-1e12)-a) < .001 && ai.converged && bi.converged);

% Missingness is a policy. Omission removes both NaN and Inf explicitly.
assertThrows(@() hnf_robmean([1 NaN],2),'hnf:statistics:Nonfinite');
[mu,mi] = hnf_robmean([1 NaN Inf 3],2,struct('nonfinite','omit'));
assert(mu == 2 && mi.count == 2 && mi.invalid);
assert(isnan(hnf_robmean([NaN Inf],2,struct('nonfinite','omit'))));
assert(isnan(hnf_robstd([NaN Inf],2,struct('nonfinite','omit'))));
assert(isnan(hnf_robmean([1 NaN],2,struct('nonfinite','propagate'))));
assert(all(isnan(hnf_robquantile([NaN Inf],[0 .5 1],2,struct('nonfinite','omit')))));
[z,mu,s,zi] = hnf_robzscore([1 NaN 3],2,struct('nonfinite','omit'));
assert(isnan(z(2)) && mu == 2 && zi.count == 2);
assertClose(z([1 3]),[-1 1]/sqrt(2),1e-14); assertClose(s,sqrt(2),1e-14);
assert(isequal(hnf_robzscore([2 2],2),[0 0]));
assert(all(isnan(hnf_robzscore([2 2],2,struct('constant','nan')))));
assertThrows(@() hnf_robzscore([2 2],2,struct('constant','error')),'hnf:statistics:ZeroScale');
z = hnf_robzscore([-big big],2); assertClose(z,[-1 1]/sqrt(2),1e-14);

% Complex classical statistics retain imaginary values and conjugation.
cx = [1+2i 2-i; 3-1i -1+3i; -2+4i 4+2i];
assertClose(hnf_robmean(cx),mean(cx),1e-14);
assertClose(hnf_robstd(cx),std(cx),1e-14);
expected = bsxfun(@minus,cx,mean(cx)); expected = expected'*expected/2;
[c,mu,ci] = hnf_robcov(cx);
assertClose(c,expected,1e-13); assertClose(mu,mean(cx),1e-14);
assertClose(c,c',1e-14); assert(ci.psdByConstruction);
r = hnf_robcorr(cx); sd = sqrt(real(diag(expected)));
assertClose(r,bsxfun(@rdivide,bsxfun(@rdivide,expected,sd),sd.'),1e-13);
assertClose(hnf_robmean([complex(big,big) complex(big,big)],2),complex(big,big),0);
assertThrows(@() hnf_robmean(cx,1,struct('method','median')),'hnf:statistics:Input');
assertThrows(@() hnf_robcov(cx,struct('method','winsor')),'hnf:statistics:ComplexOrder');

% Scaled covariance products, PSD construction, and correlation without covariance.
y = [1 2;2 1;3 5;4 4];
[c,mu,ci] = hnf_robcov(y); assertClose(c,cov(y),1e-13); assertClose(mu,mean(y),1e-14);
assert(ci.count == 4 && ci.invalidRowCount == 0);
assertClose(hnf_robcorr(y),corrcoef(y),1e-13);
assertClose(hnf_robcov(y*1e-150)/1e-300,cov(y),1e-12);
huge = [-big -big; big big];
assertClose(hnf_robcorr(huge),ones(2),1e-14);
assertThrows(@() hnf_robcov(huge),'hnf:statistics:Overflow');
[c,~,ci] = hnf_robcov([y;NaN 9],struct('nonfinite','omit'));
assertClose(c,cov(y),1e-13); assert(ci.count == 4 && ci.invalidRowCount == 1);
assert(all(all(isnan(hnf_robcov([y;NaN 9],struct('nonfinite','propagate'))))));
assert(all(all(isnan(hnf_robcov(zeros(0,2))))));
assert(isequal(hnf_robcov([2 3]),zeros(2)));
assert(all(all(isnan(hnf_robcorr(ones(3,2))))));
assert(isequal(hnf_robcorr(ones(3,2),struct('constant','zero')),zeros(2)));
c = hnf_robcov(y,struct('shrinkage',1,'loading',.5));
assertClose(c,diag(diag(cov(y)))+.5*eye(2),1e-13);
c = hnf_robcov([y y(:,1)],struct('psdrepair',true,'eigenfloor',.01));
ev = eig(c); assert(min(ev) > 0 && min(ev) >= .009999*max(ev));
[w,~,wi] = hnf_robcov([zeros(9,2);100 100],struct('method','winsor','trim',.2));
assert(isequal(w,zeros(2)) && strcmp(wi.method,'winsor'));

% Explicit clipping and option validation.
[v,vi] = hnf_robclip([-3 0 8],[-1 2]); assert(isequal(v,[-1 0 2]) && vi.clippedCount == 2);
v = hnf_robclip([NaN Inf -Inf 5],[0 1],struct('nonfinite','propagate'));
assert(isnan(v(1)) && v(2) == Inf && v(3) == -Inf && v(4) == 1);
assertThrows(@() hnf_robmean(x,1,struct('typo',1)),'hnf:statistics:Option');
assertThrows(@() hnf_robquantile(x,1.1),'hnf:statistics:Probability');
assertThrows(@() hnf_robclip(x,[2 1]),'hnf:statistics:Limits');
end

function assertClose(actual, expected, tolerance)
assert(isequal(size(actual),size(expected)));
if tolerance == 0
    assert(isequaln(actual,expected));
else
    delta = abs(actual-expected);
    assert(all(delta(:) <= tolerance*max(1,max(abs(expected(:))))));
end
end

function assertThrows(action, identifier)
caught = false;
try
    action();
catch exception
    caught = strcmp(exception.identifier,identifier);
end
assert(caught);
end
