function bench_nf_permutations(root,bin)
% Full public one-sample workload, with exact control of the sign-flip RNG.
old=rng; savedPath=path;
if exist('OCTAVE_VERSION','builtin'), pkg load statistics; end
addpath(root,genpath(fullfile(root,'statistics')),genpath(fullfile(root,'helpers')), ...
    genpath(fullfile(root,'phaseamplitudecoupling')),bin);
pol=hnf_mex_available('scope',struct()); guard=onCleanup(@()cleanup(pol,old,savedPath)); %#ok<NASGU>
rng(2401);
% Fixed seed reuses exactly the same per-subject signs in either run.
x=randn(40,24,96); x(:,8:11,40:49)=x(:,8:11,40:49)+.5;
x(:,16:19,54:62)=x(:,16:19,54:62)-.6;
p=struct();
for name={'hnf_permuted_t_mex','hnf_tfce_tail_mex'}
    p.(name{1})=which(name{1});
end
for pass=1:2
    if pass==2,p.hnf_tfce_clusters_mex=which('hnf_tfce_clusters_mex'); end
    hnf_mex_available('scope',p); rng(1960);
    clock=tic; evalc('res=nf_onesample_TFCE(x,0,[],1:24,1:96,1000,0);'); t=toc(clock);
    fprintf('1000 permutations pass %d: %.6fs\n',pass,t);
    if pass==1,ret=res;next=rand(1,8);
    else,assert(isequaln(ret.tstat,res.tstat)&&isequaln(ret.p_vals,res.p_vals)&&isequal(next,rand(1,8)));end
end
rng(1960); clock=tic; evalc('res=nf_onesample_TFCE(x,0,[],1:24,1:96,10000,0);'); t=toc(clock);
fprintf('10000 permutations native only %.6fs; significant %d, min corrected p %.12g\n',t,sum(res.p_vals(:)<=.05),min(res.p_vals(:)));
end
function cleanup(p,r,oldPath),hnf_mex_available('scope',p);rng(r);path(oldPath);end
