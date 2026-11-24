function test_nf_statistics_scaling()
% Unit changes should preserve observed and randomized effects and restore
% physical-unit beta. Use paired seeds so the same null draws are compared.
rng(918);
n = 24;
Y = randn(n,3);
X = [ones(n,1), randn(n,1)];
f = [1e-220,2,1e220];
signs = 2*(rand(n,1)>.5)-1;
order = randperm(n);
for native = [false,true]
    if native && exist('hnf_permuted_t_mex','file') ~= 3, continue; end
    t = hnf_permuted_t_map(Y,signs,false);
    u = hnf_permuted_t_map(Y.*f,signs,native);
    assert(max(abs(t-u)) < 5e-13);
    t = hnf_permuted_t_map(Y,order,false,9);
    u = hnf_permuted_t_map(Y.*f,order,native,9);
    assert(max(abs(t-u)) < 5e-13);
    a = [10;10;10;1e-220;2e-220;3e-220];
    u = hnf_permuted_t_map(a,1:6,native,3);
    assert(abs(u/1e221-sqrt(3)) < 2e-13);
    a = ([0;1;2;3;5;8;13;21]+2^51);
    t0 = hnf_permuted_t_map(a-2^51,1:8,false,4);
    u = hnf_permuted_t_map(a,1:8,native,4);
    assert(abs(t0-u) < 1e-12);
end
rng(12); base = nf_glm_TFCE(Y,X,[],[],[],3,0);
rng(12); scaled = nf_glm_TFCE(Y.*f,X.*[1 1e40],[],[],[],3,0);
beta = scaled.beta .* [1;1e40] ./ f;
assert(max(abs(beta(:)-base.beta(:))) < 1e-12);
assert(max(abs(scaled.tstat(:)-base.tstat(:))) < 1e-12);
assert(isequal(scaled.p_vals,base.p_vals));
for method={'pearson','spearman','kendall'}
    a = hnf_stats_corr(Y,X(:,2),method{1});
    b = hnf_stats_corr(Y.*f,X(:,2)*1e-230,method{1});
    assert(max(abs(a-b)) < 5e-13);
end
% Matrix products across multiple interaction factors can return to range.
if exist('table','file') || exist('table','class')
    a = [1;2;3]*1e200; b=a; c=[1;2;3]*1e-200;
    T = table(a,b,c,'VariableNames',{'a','b','c'});
    [D,~] = nf_makefactorialdesign(T,'maxOrder',3);
    assert(max(abs(D(:,end)./1e200-[1;8;27])) < 1e-12);
end
fprintf('PASS scale-invariant correlation, OLS, Student and Welch maps\n');
end
