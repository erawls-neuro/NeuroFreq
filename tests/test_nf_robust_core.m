function test_nf_robust_core()
% Scale tests exercise the public functions through the scoped entry point.
nf_run(@check_core);
end

function check_core()
R=realmax; tiny=realmin*eps;
assert(isequaln(hnf_robdiv([NaN;Inf;2],[3,0]),[0 0;0 0;2/3 0]));
assert(isinf(hnf_robdiv(R,tiny)));
assert(hnf_robdiv(complex(R,R),complex(R,-R))==1i);
assert(hnf_robdiv(1,0,struct('denominatorPolicy','floor','denominatorFloor',.5))==2);
for name={'single','double'}
    cls=name{1}; maxval=realmax(cls); minval=realmin(cls);
    A=diag([maxval,minval]); L=hnf_robchol(A);
    assert(all(isfinite(L(:))) && all(diag(L)>0));
    [L,ci]=hnf_robchol(repmat(maxval,2));
    assert(all(isfinite(L(:))) && ci.repaired);
    X=cast([1,0;0,1],cls).*maxval;
    P=hnf_robpinv(X);
    assert(norm(double(X*P)-eye(2)) < 16*eps(cls));
    z=cast([1;0.5],cls).*maxval;
    [b,~,p]=hnf_robols(z,z);
    assert(abs(b-1)<16*eps(cls) && ~p);
    assert(hnf_robrms([maxval maxval],2,'includenan')==maxval);
    assert(hnf_robrms([minval minval],2,'includenan')==minval);
end
assert(isequal(hnf_robpinv([NaN 0;0 2]),diag([0 .5])));
P=hnf_robpinv(diag([1e308,1e-200]),0);
assert(P(1,1)>0 && abs(P(2,2)/1e200-1)<eps);
[~,~,p]=hnf_robols(diag([1,1e-10]),[2;2]); assert(p);
assert(abs(hnf_robmean([R,1e-200,-R])/1e-200-1/3)<eps);
assert(abs(hnf_robstd(complex([R R],[-1e-200 1e-200]),1,2)/1e-200-1)<eps);
assert(isequaln(hnf_robmean([Inf NaN;1 1],1,'omitnan'),[Inf 1]));
assert(isequaln(hnf_robstd([Inf NaN;1 1],0,1,'omitnan'),[NaN 0]));
assert(isequaln(hnf_robrms([NaN NaN Inf;NaN 3 NaN],1,'omitnan'),[NaN 3 Inf]));
assert(isequaln(hnf_robnorm([NaN NaN Inf;NaN 3 NaN],1,'omitnan'),[0 3 Inf]));
assert(isinf(hnf_robnorm([R R],2,'includenan')));
assert(isequal(size(hnf_robrms(zeros(0,3),1,'omitnan')),[1 3]));
x=reshape(mod(0:131071,101)-50,128,64,16);
x(:,1,2)=NaN; x(1,1,2)=Inf;
for dim=[1 2 3 5]
    for omit=[false true]
        flag='includenan'; if omit, flag='omitnan'; end
        a=hnf_robrms(x,dim,flag);
        [ref,~]=hnf_robrms(x,dim,flag); % request metadata to force .m path
        assert(isequal(size(a),size(ref)) && isequal(isnan(a),isnan(ref)) ...
            && isequal(isinf(a),isinf(ref)));
        finite=isfinite(ref);
        assert(all(abs(a(finite)-ref(finite))<=8*eps.*max(1,ref(finite))));
    end
end
end
