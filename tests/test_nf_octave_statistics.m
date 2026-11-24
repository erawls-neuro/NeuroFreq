function test_nf_octave_statistics()
% Run after nf_run, with datatypes/Tablicious and the chosen LME backend.
% Numeric paths below also run on base Octave without any table package.
x = (1:10)';
y = [3.1;5.2;6.9;9.1;11;12.9;15.2;17.1;19;5];
b = hnf_robustfit(x,y);
% Default bisquare reference from GNU statistics robustfit, verified there
% against MATLAB. Last response is a gross outlier.
assert(max(abs(b-[1.08223788958791;1.9947584179332])) < 2e-8);
trial = struct('power',reshape(y,[1 1 1 10]),'phase',reshape(y,[1 1 1 10]), ...
    'times',0,'freqs',10);
rho = corrcoef(x,y);
expected = [rho(1,2),0.56363636363636,0.64444444444444,b(2)*std(x)/std(y)];
methods = {'pearson','spearman','kendall','robust'};
for k = 1:numel(methods)
    r = nf_stregress(trial,x,'method',methods{k},'permutations',0,'progress',0);
    assert(abs(r.power.corrcoef-expected(k)) < 2e-8);
    assert(isfinite(r.phase.corrcoef));
    r = nf_stregress(trial,x,'method',methods{k},'permutations',3,'perm_seed',3,'progress',0);
    assert(isfinite(r.power.corrcoef));
end
% Tied pairs: two concordant, one tie on each variable and two further
% concordant pairs. Tau-b denominator is sqrt(5*5).
assert(abs(hnf_stats_corr([1;1;2;3],[1;2;2;3],'kendall')-4/5)<1e-12);

% Independently profiled Gaussian block covariance Vg = Zg*Psi*Zg'+s2*I.
% Log determinant, GLS coefficients and their standard errors were computed
% outside both fitters. Unbalanced groups distinguish full from diagonal.
counts = repmat([7;9;10],4,1);
subject = repelem((1:12)',counts);
i = (1:numel(subject))';
x = sin(i*1.39) + .1*subject;
y = 1 + (.7+.4*cos(subject)+.2*sin(subject)).*x + .7*cos(subject) + ...
    .5*sin(subject) + .19*sin(i*4.417) + .09*cos(i*.651);
design = table(categorical(subject),x,'VariableNames',{'subject','x'});
[~,design2] = nf_makefactorialdesign(design,'keepSubjectInDesignTab',true);
assert(isequal(design2.Properties.VariableNames,{'subject','x'}));
assert(isequal(hnf_subject_codes(categorical({'b';'a';'b'})),[2;1;2]));
ref = [0.67114785 17.3886610; 0.67950327 6.9318012; 0.67745211 7.2687165];
slopes = {[],{'x'},{'x'}};
methods = {'REML','REML','ML'};
for k = 1:3
    r = nf_lme_TFCE([y 2*y],design,[],[],1:2,1,0,'randomSlopes',slopes{k}, ...
        'uncorrelated',k~=3,'fitMethod',methods{k});
    assert(max(abs(r.beta(:) - ref(k,1)*[1;2])) < 5e-5);
    assert(max(abs(r.tstat(:) - ref(k,2))) < .002);
    assert(all(r.p_vals(:) >= .5 & r.p_vals(:) <= 1));
    fprintf('%s, %s: slope %.7f, t %.6f\n',r.fitBackend,methods{k},r.beta(1),r.tstat(1));
end
fprintf('PASS single-trial inference, ranks, and 3 Gaussian mixed-model specifications\n');
end
