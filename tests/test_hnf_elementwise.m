function result = test_hnf_elementwise
% TEST_HNF_ELEMENTWISE  Numerical edge cases and observable policy contracts.
% Run in MATLAB from any current directory. Throws on the first failure.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; restorePath = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'helpers','robustarithmetic'));
R = realmax; r = realmin; tiny = realmin.*eps;

% Division: broadcasting, domain policy, representable complex extremes.
y = hnf_robdiv([1;2],[1 0 2]);
assert(isequal(y,[1 0 .5;2 0 1]));
assert(isequaln(hnf_robdiv([1;2],[1 0 2],struct()),[1 NaN .5;2 NaN 1]));
assert(hnf_robdiv(NaN,2)==0 && isnan(hnf_robdiv(NaN,2,struct())));
assert(hnf_robdiv(NaN,2,struct('inputPolicy','fill','fillValue',4)) == 2);
assert(hnf_robdiv(1,0,struct('denominatorPolicy','native')) == Inf);
assert(hnf_robdiv(1,-1e-20,struct('denominatorPolicy','floor','denominatorFloor',1e-6)) == -1e6);
q = hnf_robdiv(1,(3+4i).*1e-300,struct('denominatorPolicy','floor','denominatorFloor',1e-6));
assert(abs(q-(.6-.8i).*1e6) < 1e-8);
assert(abs(hnf_robdiv(complex(R,R),complex(R,R))-1) < 8*eps);
assert(abs(hnf_robdiv(complex(R,R),complex(R,-R))-1i) < 8*eps);
q = hnf_robdiv(complex(1,1e-320),r);
assert(abs(imag(q)./(1e-320./r)-1) < 8*eps);
q = hnf_robdiv(complex(R,0),complex(2,tiny));
assert(abs(imag(q)./(-(R./4).*tiny)-1) < 8*eps);
[q,divisionInfo] = hnf_robdiv([1+2i 3-4i],[2-3i 4+5i]);
assert(isequal(q,[1+2i 3-4i]./[2-3i 4+5i]) && divisionInfo.rescuedComplexQuotients == 0);
assert(isinf(hnf_robdiv(R,eps)));
elem_assertError(@() hnf_robdiv(R,eps,struct()),'hnf_robdiv:Overflow');
assert(hnf_robdiv(R,eps,struct('overflow','saturate')) == R);
elem_assertError(@() hnf_robdiv(tiny,2,struct('underflow','error')),'hnf_robdiv:Underflow');
elem_assertError(@() hnf_robdiv(1e300,single([1 2]),struct()),'hnf_robdiv:Type');
elem_assertError(@() hnf_robdiv(single([1 2]),1e-300,struct()),'hnf_robdiv:Type');

% Domain boundaries must retain missingness and distinguish clipping.
assert(isequaln(hnf_roblog([0 -1 NaN 1]),[-Inf NaN NaN 0]));
assert(abs(hnf_roblog(100,struct('base',10))-2) < 8*eps);
assert(hnf_roblog(0,struct('domainPolicy','floor','floorValue',1)) == 0);
z = hnf_roblog(complex(R,R),struct('domainPolicy','complex'));
assert(isfinite(z) && abs(real(z)-(log(R)+log(2)/2)) < 8*eps(log(R)));
assert(abs(imag(z)-pi/4) < 8*eps);
assert(isequaln(hnf_robsqrt([4 -eps -1 NaN],struct('negativeTolerance',eps)),[2 0 NaN NaN]));
assert(hnf_robsqrt(-1,struct('domainPolicy','complex')) == 1i);
assert(hnf_robsqrt(-3,struct('domainPolicy','clip')) == 0);
elem_assertError(@() hnf_robsqrt(-1,struct('domainPolicy','error')),'hnf_robsqrt:Domain');

% Exp: honest overflow and explicit saturation/floor policies.
assert(max(abs(hnf_robexp([-1 0 1])-exp([-1 0 1]))) == 0);
elem_assertError(@() hnf_robexp(1000),'hnf_robexp:Overflow');
assert(hnf_robexp(1000,struct('overflow','saturate')) == R);
assert(hnf_robexp(-1000,struct('underflow','floor')) == tiny);
assert(isnan(hnf_robexp(NaN)) && hnf_robexp(-Inf) == 0);
z = hnf_robexp(complex(log(R)+.2,pi/4));
assert(isfinite(z) && real(z) > R/2 && imag(z) > R/2);

% Norm and unit-vector routines distinguish representable from overflowed
% true norms, and normalize even when that norm itself cannot be stored.
assert(abs(hnf_robnorm([3e200 4e200],2)/5e200-1) < 8*eps);
assert(abs(hnf_robnorm([3e-200 4e-200],2)/5e-200-1) < 8*eps);
assert(hnf_robnorm(tiny) == tiny);
elem_assertError(@() hnf_robnorm([R R],2),'hnf_robnorm:Overflow');
assert(hnf_robnorm([R R],2,struct('overflow','saturate')) == R);
assert(isnan(hnf_robnorm([Inf NaN],2)));
assert(hnf_robnorm([Inf NaN],2,struct('nanPolicy','omit')) == Inf);
assert(isequal(hnf_robnorm(zeros(0,3),1),zeros(1,3)));
u = hnf_robunit([R R],2); assert(max(abs(u-1/sqrt(2))) < 8*eps);
u = hnf_robunit(complex(R,R)); assert(abs(u-(1+1i)/sqrt(2)) < 8*eps);
assert(hnf_robunit(tiny) == 1);
assert(all(isnan(hnf_robunit([0 0],2))));
assert(isequal(hnf_robunit([0 0],2,struct('zeroPolicy','zero')),[0 0]));
assert(all(isnan(hnf_robunit([1 Inf],2))));

% Selective compensation rescues cancellation and intermediate overflow.
assert(abs(hnf_robsum([1e16 1 -1e16],2)-1) < 8*eps);
assert(hnf_robsum([R R -R],2) == R);
assert(hnf_robsum([R R -R -R],2) == 0);
assert(isequal(hnf_robsum([1e16 2;1 3;-1e16 4],1),[1 9]));
assert(isequal(hnf_robsum([1e16 1 -1e16;2 3 4],2),[1;9]));
assert(hnf_robsum([1 NaN 2],2,struct('nanPolicy','omit')) == 3);
assert(isnan(hnf_robsum([1 NaN 2],2)));
assert(isnan(hnf_robsum([Inf -Inf],2)));
elem_assertError(@() hnf_robsum([R R],2),'hnf_robsum:Overflow');

% Log-sum-exp retains tiny positive log sums as well as huge logits.
assert(abs(hnf_roblogsumexp([1000 1000],2)-(1000+log(2))) < 8*eps(1000));
assert(abs(hnf_roblogsumexp([0 -40],2)/log1p(exp(-40))-1) < 8*eps);
assert(hnf_roblogsumexp([-Inf -Inf],2) == -Inf);
assert(hnf_roblogsumexp([Inf 2],2) == Inf);
assert(isnan(hnf_roblogsumexp([NaN Inf],2)));
assert(hnf_roblogsumexp([NaN 2],2,struct('nanPolicy','omit')) == 2);
assert(isequal(hnf_roblogsumexp(zeros(0,3),1),[-Inf -Inf -Inf]));
x = reshape(1:24,2,3,4); x = double(x);
difference = hnf_roblogsumexp(x,2)-log(sum(exp(x),2));
assert(max(abs(difference(:))) < 1e-12);
assert(isequal(hnf_roblogsumexp(x,5),x));

% Softmax checks both stable temperature ordering and explicit conventions.
p = hnf_robsoftmax([1000 1001 1002],2);
assert(abs(sum(p)-1) < 8*eps && all(p > 0));
p = hnf_robsoftmax([-R R],2,struct('temperature',R));
assert(max(abs(p-[exp(-2) 1]./(1+exp(-2)))) < 8*eps);
assert(all(isnan(hnf_robsoftmax([Inf 1],2))));
assert(isequal(hnf_robsoftmax([Inf Inf 0],2,struct('infinityPolicy','equal')),[.5 .5 0]));
assert(isequal(hnf_robsoftmax([-Inf -Inf],2,struct('degeneratePolicy','uniform')),[.5 .5]));
assert(all(isnan(hnf_robsoftmax([NaN 1],2))));
elem_assertError(@() hnf_robsoftmax([0 -1000],2,struct('underflow','error')),'hnf_robsoftmax:Underflow');

% Power logging must precede squaring for extreme finite coefficients.
assert(hnf_robpower(3+4i) == 25);
assert(hnf_robabs2(3+4i) == 25);
elem_assertError(@() hnf_robpower(R),'hnf_robpower:Overflow');
assert(abs(hnf_robpower(complex(R,R),struct('output','log'))-(2*log(R)+log(2))) < 8*eps(2*log(R)));
assert(isfinite(hnf_robpower(tiny,struct('output','db'))));
assert(hnf_robpower(0,struct('output','log')) == -Inf);
assert(hnf_robpower(0,struct('output','log','floorValue',1)) == 0);
assert(isnan(hnf_robpower(NaN,struct('output','log','floorValue',1))));
v = sqrt(tiny).*sqrt(.3); assert(hnf_robpower(complex(v,v)) == tiny);

% Class, unsupported input, strict options, and diagnostic contracts.
assert(isa(hnf_robdiv(single([1 2]),single(2)),'single'));
assert(isa(hnf_robnorm(single([3 4]),2),'single'));
assert(isa(hnf_robsoftmax(single([1 2]),2),'single'));
assert(isa(hnf_roblogsumexp(single([1 2]),2),'single'));
assert(isa(hnf_robpower(single(3+4i)),'single'));
elem_assertError(@() hnf_robnorm(int16(1)),'hnf_robnorm:Type');
elem_assertError(@() hnf_robsum(sparse(eye(2)),1),'hnf_robsum:Type');
elem_assertError(@() hnf_robexp(1,struct('overfow','saturate')),'hnf_robexp:Options');
[~,d] = hnf_robdiv([1 NaN],0);
assert(d.smallDenominator == 2 && d.nonfiniteInput == 1);
result = struct('passed',true,'suite','elementwise');
end

function elem_assertError(f,expected)
try
    f();
catch exception
    assert(strcmp(exception.identifier,expected), ...
        'Expected %s, got %s.',expected,exception.identifier);
    return;
end
error('test_hnf_elementwise:MissingError','Expected %s.',expected);
end
