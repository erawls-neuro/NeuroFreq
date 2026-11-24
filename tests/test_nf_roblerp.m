function test_nf_roblerp()
% Stress the same interpolation used in type-7 quantiles and TF resampling.
root=fileparts(fileparts(mfilename('fullpath')));
saved=path; cleanup=onCleanup(@()path(saved)); %#ok<NASGU>
addpath(root,fullfile(root,'helpers','robustarithmetic'));
x=realmax; q=[0;.25;.5;.75;1];
y=hnf_roblerp([-x 3],[x 5],q);
assert(max(abs(y(:,1)./x-[-1;-.5;0;.5;1])) <= eps);
assert(isequal(y(:,2),[3;3.5;4;4.5;5]));
assert(isequal(hnf_roblerp(1e300,1,1),1));
assert(isequal(hnf_roblerp(complex(-x,x),complex(x,-x),.5),0));
tf=struct('times',[0 1],'freqs',[-x x],'Fs',1,'nsensor',1,'power',[-x x;x -x]);
z=nf_run('nf_resample',tf,[],[-x 0 x],false);
assert(isequal(z.power,[-x x;0 0;x -x]));
% At an exact knot the other endpoint's NaN has no influence in the helper.
assert(hnf_roblerp(2,NaN,0)==2 && hnf_roblerp(NaN,4,1)==4);
fprintf('Robust interpolation and public resampling checks passed.\n');
end
