function elapsed = bench_nf_qc(root,bin)
% Public QC, exact fixture and per-kernel dispatch policy. Run after nf_run.
% Supply a directory with current-runtime builds of hnf_qc_window_mex and
% hnf_qc_percentiles_mex (or use build_nf_bench_kernels).
orig=path; state=rng; exits=onCleanup(@() restore(orig,state)); %#ok<NASGU>
addpath(bin,genpath(fullfile(root,'helpers')),genpath(fullfile(root,'preprocessing')), ...
    fullfile(root,'phaseamplitudecoupling','hnf_phaseamplitudecoupling'));
assert(exist('hnf_qc_window_mex','file')==3 && exist('hnf_qc_percentiles_mex','file')==3);
old=hnf_mex_available('scope',struct()); scope=onCleanup(@() hnf_mex_available('scope',old)); %#ok<NASGU>
rng(21,'twister'); e.data=randn(16,512,16);e.nbchan=16;e.pnts=512;e.trials=16;
e.srate=256;e.xmin=0;e.xmax=511/256;
e.chanlocs=struct('labels',arrayfun(@num2str,1:16,'UniformOutput',false));
post=e;post.data(2,:,:)=NaN;post.data(8,47,2)=Inf;
policies={struct(),struct('hnf_qc_percentiles_mex',which('hnf_qc_percentiles_mex'))};
policies{3}=policies{2};policies{3}.hnf_qc_window_mex=which('hnf_qc_window_mex');
times=zeros(3,6);
for rep=0:6
    order=mod((0:2)+rep,3)+1;
    for k=order
        hnf_mex_available('scope',policies{k});
        clock=tic; q=nf_eegquality(e,post,'plot',false); seconds=toc(clock);
        if rep>0, times(k,rep)=seconds; end
        if rep==0 && k==1, reference=q; end
        for label={'preclean','postclean'}
            a=reference.metrics.(label{1});b=q.metrics.(label{1});
            assert(isequal(a.artifact.flagBits,b.artifact.flagBits));
            assert(isequaln(a.psd.channelDb,b.psd.channelDb));
            for field=fieldnames(a.artifact.measurements)'
                x=a.artifact.measurements.(field{1});y=b.artifact.measurements.(field{1});
                assert(isequal(isnan(x),isnan(y)));
                delta=abs(x(isfinite(x))-y(isfinite(x)));
                assert(isempty(delta) || max(delta)<1e-11);
            end
        end
    end
end
elapsed=median(times,2);
fprintf('QC seconds: M-files %.4f, percentile MEX %.4f, both MEX %.4f; fused-to-prior speedup %.2fx.\n', ...
    elapsed(1),elapsed(2),elapsed(3),elapsed(2)/elapsed(3));
end

function restore(paths,state)
rng(state);path(paths);
end
