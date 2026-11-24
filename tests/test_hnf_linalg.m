function result = test_hnf_linalg()
% Focused deterministic MATLAB assertions; no toolboxes, random data or output.
% Run after addpath(fullfile(packageRoot,'helpers','robustarithmetic')).
A=[1 0;1 1;1 2;1 3]; y=[1;3;5;7];
[b,rc,p,inf]=hnf_robols(A,y);
assert(norm(b-[1;2])<1e-12 && rc>0 && ~p && strcmp(inf.method,'qr'));
Ac=A+1i*[0 1;1 0;0 -1;2 1]; bc=[2+1i;-3+2i];
assert(norm(hnf_robols(Ac,Ac*bc)-bc)<1e-11);
Ar=[1 2;2 4;3 6]; yr=[1;2;3];
[b,~,p]=hnf_robols(Ar,yr);
assert(p && norm(b-[.2;.4])<1e-12);
Aw=[1 0 1;0 1 1]; yw=[2;3];
assert(norm(hnf_robols(Aw,yw)-pinv(Aw)*yw)<1e-12);
for M={A,Ar,Aw,Ac,Ac.'}
    X=M{1}; P=hnf_robpinv(X);
    assert(norm(X*P*X-X,'fro')<1e-10*max(1,norm(X,'fro')));
    assert(norm(P*X*P-P,'fro')<1e-10*max(1,norm(P,'fro')));
    assert(norm(X*P-(X*P)','fro')<1e-10);
    assert(norm(P*X-(P*X)','fro')<1e-10);
end
P=hnf_robpinv(diag([1 1e-5]),1e-4);
assert(isequal(P,diag([1 0])));
assert(isequal(hnf_robpinv(zeros(3,2)),zeros(2,3)));
assert(isequal(size(hnf_robpinv(zeros(0,3))),[3 0]));
assert(isequal(size(hnf_robols(zeros(0,3),zeros(0,2))),[3 2]));
b=hnf_robols([1;NaN;1],[1;50;3],struct('missing','omitrows'));
assert(abs(b-2)<1e-12);
b=hnf_robols([1;1;1],[1;2;100],struct('weights',[1;1;0]));
assert(abs(b-1.5)<1e-12);
b=hnf_robols([realmin;realmax],[1;0],struct('weights',[1;0]));
assert(abs(b*realmin-1)<1e-12);
b=hnf_robols(1,[1e308 1e-308]);
assert(b(1)==1e308 && abs(b(2)/1e-308-1)<1e-12);
b=hnf_robols(1e308,1e308);
assert(abs(b-1)<1e-12);
tiny=realmin*eps;
b=hnf_robols(2,[5*tiny 1]);
assert(b(1)==2*tiny && b(2)==.5); % round-to-nearest-even subnormal
P=hnf_robpinv(complex(realmax,realmax));
assert(isfinite(P) && real(P)>0 && imag(P)<0);
base=[ones(21,1),(0:20)']; truth=[2;.25]; yy=base*truth; yy(11)=yy(11)+200;
[b,~,~,rinfo]=hnf_robols(base,yy,struct('robust','huber','maxIterations',100));
assert(norm(b-truth)<1e-4 && all(rinfo.converged) && rinfo.weights(11)<1e-4);
[b,~,~,rinfo]=hnf_robols(base,yy,struct('robust','bisquare','maxIterations',100));
assert(norm(b-truth)<1e-4 && all(rinfo.converged));
[br,~,~,ri]=hnf_robols(A,y,struct('ridge',.1));
an=A/max(abs(A(:))); yn=y/max(abs(y));
expected=[an;sqrt(.1)*eye(2)]\[yn;zeros(2,1)]*max(abs(y))/max(abs(A(:)));
assert(norm(br-expected)<1e-12 && strcmp(ri.method,'svd-ridge'));
C=[2 .5;.5 1];
[L,ci,Cused]=hnf_robchol(C);
assert(norm(L*L'-C,'fro')<1e-12 && ci.shift==0 && isequal(Cused,C));
[L,ci,Cused]=hnf_robchol([1 1;1 1]);
assert(ci.repaired && ci.shift>0 && norm(L*L'-Cused,'fro')<1e-12);
[L,ci,Cused]=hnf_robchol(C,.125);
assert(abs(ci.shift-.125)<eps && norm(Cused-(C+.125*eye(2)),'fro')<1e-12);
[L,ci,Cused]=hnf_robchol(diag([1 -1]),struct('maxRelativeShift',2));
assert(ci.shift>1 && norm(L*L'-Cused,'fro')<1e-12);
[L,ci]=hnf_robchol(zeros(2));
assert(all(isfinite(L(:))) && ci.shift>0);
[X,si]=hnf_robsolve(C,[1;3],struct('structure','positive'));
assert(norm(C*X-[1;3])<1e-12 && strcmp(si.method,'chol'));
[X,si]=hnf_robsolve([1 1;1 1],[2;2],struct('structure','positive'));
assert(norm(X-[1;1])<1e-12 && si.usedPseudoinverse);
[V,d]=hnf_robeig(C);
assert(norm(V*diag(d)*V'-C,'fro')<1e-12);
[V,d,ei]=hnf_robeig(diag([-1 2]),struct('psd','clip'));
assert(isequal(d,[0;2]) && ei.modifiedEigenvalues==1 && norm(V'*V-eye(2))<1e-12);
Cc=[2 1i;-1i 2]; [V,d]=hnf_robeig(Cc);
assert(norm(V*diag(d)*V'-Cc,'fro')<1e-12);
[W,wi,Cused]=hnf_robwhiten(Cc);
assert(norm(W*Cused*W'-eye(2),'fro')<1e-12 && wi.rank==2);
[W,wi,Cused]=hnf_robwhiten(diag([2 0]),struct('method','pca'));
assert(isequal(size(W),[1 2]) && wi.rank==1 && abs(W*Cused*W'-1)<1e-12);
[W,wi,Cused]=hnf_robwhiten(diag([2 0]));
assert(norm(W*Cused*W'-diag([1 0]),'fro')<1e-12 && wi.rank==1);
[W,wi,Cused]=hnf_robwhiten(zeros(2),struct('mode','floor','relativeTolerance',1e-6));
assert(norm(W*Cused*W'-eye(2),'fro')<1e-12 && wi.rank==2);
[v,s]=hnf_roblogdet(diag([1e200 1e200]));
assert(abs(v-400*log(10))<1e-10 && s==1);
[v,s]=hnf_roblogdet(diag([1e308 1e-308]));
assert(abs(v)<1e-10 && s==1);
[v,s]=hnf_roblogdet(diag([1e308 1e-308]),struct('structure','positive'));
assert(abs(v)<1e-10 && s==1);
[v,s]=hnf_roblogdet([1e308 1e-308;1e308 0]);
assert(abs(v)<1e-10 && s==-1);
[v,s]=hnf_roblogdet([0 1;1 0]); assert(abs(v)<eps && s==-1);
[v,s]=hnf_roblogdet(diag([1i 2])); assert(abs(v-log(2))<1e-12 && abs(s-1i)<1e-12);
[v,s]=hnf_roblogdet(zeros(2)); assert(isinf(v) && v<0 && s==0);
[v,s]=hnf_roblogdet(zeros(0)); assert(v==0 && s==1);
[v,s]=hnf_roblogdet(C,struct('structure','positive'));
assert(abs(v-log(det(C)))<1e-12 && s==1);
assert(isa(hnf_robols(single(A),single(y)),'single'));
assert(isa(hnf_robpinv(single(A)),'single'));
assert(isa(hnf_robchol(single(C)),'single'));
assert(isa(hnf_robwhiten(single(C)),'single'));
assert(isa(hnf_robsolve(single(C),single([1;3]),struct('structure','positive')),'single'));
musterror(@()hnf_robols([1;NaN],[1;2]),'hnf:NonfiniteInput');
assert(isequal(hnf_robpinv(speye(2)),speye(2)));
musterror(@()hnf_robpinv(sparse(eye(2)),struct()),'hnf:UnsupportedInput');
musterror(@()hnf_robpinv(single(eye(2)),struct('ridge',1e-100)),'hnf:InvalidOption');
musterror(@()hnf_robchol(diag([1 -1]),struct()),'hnf:ShiftLimit');
musterror(@()hnf_robchol([1 1;0 1],struct()),'hnf:NotHermitian');
musterror(@()hnf_robwhiten(diag([1 -1])),'hnf:NotPositiveSemidefinite');
musterror(@()hnf_robpinv(realmin*eps,0),'hnf:UnrepresentableOutput');
musterror(@()hnf_robols(A,y,struct('rigde',.1)),'hnf:UnknownOption');
musterror(@()hnf_robols(A,y+1i,struct('robust','huber')),'hnf:ComplexRobust');
before=warning; transcript=evalc('hnf_robols(A,y); hnf_robpinv(Ar); hnf_robchol([1 1;1 1]);');
after=warning;
assert(isempty(strtrim(transcript)) && isequal(before,after));
result=struct('passed',true,'suite','hnf_linalg');
end

function musterror(f,id)
try
    f();
catch e
    assert(strcmp(e.identifier,id),'Expected %s, got %s.',id,e.identifier);
    return;
end
error('test_hnf_linalg:ExpectedError','Expected %s.',id);
end
