function build_nf_bench_kernels(root,bin)
% Build an isolated numerical comparison, leaving normal installation alone.
if ~exist(bin,'dir'),mkdir(bin);end
for directory={'phaseamplitudecoupling','timefrequency'}
    files=dir(fullfile(root,directory{1},['hnf_' directory{1}],'hnf_*_mex.c'));
    for k=1:numel(files)
        source=fullfile(files(k).folder,files(k).name);
        [~,name]=fileparts(source);
        if exist('OCTAVE_VERSION','builtin')
            [msg,code]=mkoctfile('--mex','-O',source,'-o',fullfile(bin,[name '.' mexext]));
            assert(code==0,'%s',msg);
        else
            mex('-O','-R2018a','-outdir',bin,source);
        end
    end
end
end
