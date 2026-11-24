function test_nf_compensated_sum()
root=fileparts(fileparts(mfilename('fullpath')));
old=path; restore=onCleanup(@()path(old)); %#ok<NASGU>
addpath(fullfile(root,'helpers','robustarithmetic'));
for cls={'single','double'}
    x = realmax(cls{1}); tiny = realmin(cls{1})*eps(cls{1});
    a = [x x tiny -x -x];
    assert(hnf_robsum(a,2)==tiny);
    assert(hnf_robsum(-a,2)==-tiny);
    assert(hnf_robsum(a+1i*a,2)==complex(tiny,tiny));
    assert(hnf_robsum([x -x 1 eps(cls{1})/2 tiny],2)==1+eps(cls{1}));
    assert(hnf_robsum([x -x 1 eps(cls{1})/2 -tiny],2)==1);
    assert(hnf_robsum([x -x tiny tiny],2)==2*tiny);
end
fprintf('Finite cancellation, signed subnormal and tie-rounding checks passed.\n');
end
