function test_hnf_auxiliary
% Known-answer checks for probability, Fisher, RMS, and streaming helpers.
P = hnf_robprob([realmax realmax],2);
assert(isequal(P,[.5 .5]));
assert(abs(hnf_robentropy(P,2,struct('base',2))-1)<1e-14);
assert(isequal(hnf_robentropy([1 0]),0));
assert(isequal(hnf_robprob(zeros(2,3),2,struct('zero','uniform')),ones(2,3)/3));
Pzero = hnf_robprob(zeros(2,3),2,struct('zero','nan'));
assert(all(isnan(Pzero(:))));
assert(isequal(hnf_robprob([0 0],2,struct('pseudocount',1)),[.5 .5]));
must_error(@() hnf_robprob([0 0]),'hnf:ZeroMass');
must_error(@() hnf_robprob([-1 2]),'hnf:Domain');
must_error(@() hnf_robentropy([.2 .2]),'hnf:ProbabilitySum');
Z = hnf_robfisher([-.9 0 .9]);
assert(max(abs(Z-atanh([-.9 0 .9])))<1e-14);
must_error(@() hnf_robfisher(1),'hnf:Boundary');
assert(isfinite(hnf_robfisher(1,struct('boundary','clip'))));
assert(isinf(hnf_robfisher(1,struct('boundary','allow'))));
assert(hnf_robrms([realmax realmax])==realmax);
assert(abs(hnf_robrms([3+4i 0])-5/sqrt(2))<1e-14);
assert(hnf_robrms([NaN 3],2,struct('nonfinite','omit'))==3);
X = [1 4; 2 7; 3 10; 4 13];
[s,mu,sd] = hnf_robmoments(X,[]);
assert(max(abs(mu-mean(X,1)))<1e-14);
assert(max(abs(sd-std(X,0,1)))<1e-14);
s1 = hnf_robmoments(X(1:2,:),[]);
[s1,m1,d1] = hnf_robmoments(X(3:4,:),s1);
assert(max(abs(m1-mu))<1e-14 && max(abs(d1-sd))<1e-14);
[~,m2,d2] = hnf_robmoments([],s1);
assert(isequal(m1,m2) && isequal(d1,d2));
assert(isequal(s.count,[4 4]));
[~,mu,sd] = hnf_robmoments([realmax;realmax],[]);
assert(mu==realmax && sd==0);
[~,mu,sd] = hnf_robmoments([1e-300;3e-300],[]);
assert(abs(mu/2e-300-1)<1e-14 && abs(sd/(sqrt(2)*1e-300)-1)<1e-14);
[~,mu,sd] = hnf_robmoments(1e16+[0;2;4],[]);
assert(mu==1e16+2 && abs(sd-2)<1e-14);
sExtreme = hnf_robmoments(-1e308,[]);
[~,mu,sd] = hnf_robmoments([-1e308;repmat(1e308,30,1)],sExtreme);
assert(isfinite(mu) && abs(mu/8.75e307-1)<1e-14 && isfinite(sd));
sOffset = hnf_robmoments(1e16,[]);
[~,mu,sd] = hnf_robmoments(1e16+[2;4],sOffset);
assert(mu==1e16+2 && abs(sd-2)<1e-14);
[~,mu,sd] = hnf_robmoments([-realmax;realmax],[],struct('normalization',1));
assert(mu==0 && sd==realmax);
[~,mu,sd] = hnf_robmoments([1 NaN;3 2],[],struct('nonfinite','omit'));
assert(isequal(mu,[2 2]) && abs(sd(1)-sqrt(2))<1e-14 && isnan(sd(2)));
[s,mu,sd] = hnf_robmoments(single(X),[]);
assert(isa(mu,'single') && isa(sd,'single') && isa(s.count,'double'));
assert(hnf_robrank([1 2;2 4;3 6])==1);
assert(hnf_robrank(diag([1e308 1e308]))==2);
assert(hnf_robrank(zeros(4,3))==0);
assert(hnf_robrank(diag([1 1e-8]),struct('absoluteTolerance',1e-7))==1);
Xc = [1+2i 3; 2 1-1i; -1+1i 2; 3-2i -1; 0 2i];
[C,mu] = hnf_robcov(Xc); W = hnf_robwhiten(C);
Yc = (Xc-mu)*W';
assert(norm(hnf_robcov(Yc)-eye(2),'fro')<1e-12);
end
function must_error(f,id)
try, f(); catch ME, assert(strcmp(ME.identifier,id)); return; end
error('test:ExpectedError','Expected error was not raised.');
end
