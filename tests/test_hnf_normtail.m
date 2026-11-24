function test_hnf_normtail
% Standard-normal known answers and tails beyond probability range.
% Reference values are upper-tail integrals of exp(-t^2/2)/sqrt(2*pi).
Z = [0 1.96 8 40];
[P,L,info] = hnf_robnormtail(Z,struct('sides',1));
reference = [0.5 0.0249978951482204358 6.22096057427178412e-16];
assert(max(abs(P(1:3)./reference-1)) < 5e-13);
assert(P(4)==0 && abs(L(4)-(-804.608442013753788))<2e-12);
assert(info.probabilityUnderflowCount==1 && info.logOutOfRangeCount==0);
assert(max(abs(L(1:3)-log(reference)))<5e-13);
[P2,L2] = hnf_robnormtail(Z);
assert(P2(1)==1 && max(abs(P2(1:3)./(2.*reference)-1))<5e-13);
assert(abs(L2(4)-(-804.608442013753788+log(2)))<2e-12);
[PN,LN] = hnf_robnormtail(-Z,struct('sides',1));
assert(max(abs(P(1:3)+PN(1:3)-1))<2e-16);
assert(abs(LN(3)/(-6.22096057427178609e-16)-1)<5e-13);
assert(isequal(P2,hnf_robnormtail(-Z)));
[Pbig,Lbig,ibig] = hnf_robnormtail(realmax);
assert(Pbig==0 && isinf(Lbig) && Lbig<0 && ibig.logOutOfRangeCount==1);
[Ps,Ls] = hnf_robnormtail(single([0 1.96 8 40]));
assert(isa(Ps,'single') && isa(Ls,'single') && isequal(size(Ps),[1 4]));
assert(abs(double(Ps(2))/(2*reference(2))-1)<3e-6);
assert(isfinite(Ls(4)) && abs(double(Ls(4))-double(L2(4)))<2e-4);
[~,Lnear] = hnf_robnormtail(1e-200);
assert(abs(Lnear/(-sqrt(2/pi)*1e-200)-1)<2e-14);
[~,LnearSingle] = hnf_robnormtail(single(1e-20));
assert(abs(double(LnearSingle)/(-sqrt(2/pi)*1e-20)-1)<2e-6);
E=hnf_robnormtail(zeros(0,3)); assert(isequal(size(E),[0 3]));
must_error(@() hnf_robnormtail(NaN),'hnf:Nonfinite');
must_error(@() hnf_robnormtail(Inf),'hnf:Nonfinite');
must_error(@() hnf_robnormtail(1i),'hnf:Type');
must_error(@() hnf_robnormtail(sparse(1)),'hnf:Type');
must_error(@() hnf_robnormtail(0,struct('sides',3)),'hnf:Options');
must_error(@() hnf_robnormtail(0,struct('unknown',1)),'hnf:Options');
end
function must_error(f,id)
try, f(); catch ME, assert(strcmp(ME.identifier,id)); return; end
error('test:ExpectedError','Expected error was not raised.');
end
