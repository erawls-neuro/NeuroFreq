function bench_nf_native(root,bin)
% Compare the same revised public functions under matched native policies.
% bin contains ABI-local binaries built by build_nf_bench_kernels(root,bin).
% Timings exclude fixture construction and numerical comparisons.
p=path; r=rng;
if exist('OCTAVE_VERSION','builtin')
    pkg load signal
    pkg load statistics
    previousThreads=fftw('threads'); fftw('threads',1);
else
    previousThreads=[];
end
for directory={'preprocessing','timefrequency','setutilities','phaseamplitudecoupling', ...
        'spectralparameterization','statistics','helpers'}
    addpath(genpath(fullfile(root,directory{1})));
end
addpath(root,bin);
oldPolicy=hnf_mex_available('scope',struct());
guard=onCleanup(@()restore(p,r,oldPolicy,previousThreads)); %#ok<NASGU>
src=dir(fullfile(bin,['*.' mexext])); native=struct();
for i=1:numel(src),[~,n]=fileparts(src(i).name);native.(n)=fullfile(bin,src(i).name);end
new={'hnf_tfce_clusters_mex','hnf_compsum_mex','hnf_knee_loss_mex','hnf_qc_window_mex'};
shipped=rmfield(native,intersect(fieldnames(native),new)); off=struct();
assert(isfield(native,'hnf_tfce_clusters_mex'),'Build current kernels first.');
rng(3406);
s=load(fullfile(root,'demo/timefrequency/eeglab_sample.set'),'-mat');
a=double(s.data(1:16,:,1:32)); fs=double(s.srate);
tf=struct('freqs',1:24,'times',(0:191)/256,'nsensor',8,'ntrls',8,'Fs',256);
tf.power=rand(8,24,192,8); tf.phase=2*pi*rand(size(tf.power))-pi;
tf.event=struct('latency',{3,81,160}); tf.scale='linear';
gp=struct('freqs',[4,5,20,21],'times',(0:255)/128,'Fs',128,'nsensor',2,'ntrls',12);
gp.phase=2*pi*rand(2,4,256,12)-pi; gp.power=(1+.3*cos(gp.phase)+rand(size(gp.phase))).^2;
f=linspace(1,60,64); spl=struct('freqs',f,'nsensor',24,'ntrls',8,'times',1:64,'scale','log10');
spl.power=single(randn(24,64,64,8)/10 + reshape(3-1.1*log10(f),1,64,1,1));
draw=randn(40,24,96); vec=randn(128,256,32); b=randn(size(vec)); b(1:193:end)=0;
cases={
'nf_prepdata sample16x384x32' ,@() nf_prepdata(a,1,1,1,0,0),off,shipped;
'nf_stft sample16x384x32', @() nf_stft(a,fs,.25,80,1),off,shipped;
'nf_filterhilbert sample8x384x12',@() nf_filterhilbert(a(1:8,:,1:12),fs,4:4:32,4,3),off,shipped;
'nf_resample (8x24x192x8) to 37x241', @() nf_resample(tf,linspace(tf.times(1),tf.times(end),241),linspace(1,24,37),false),off,shipped;
'nf_pac cross-MVL 2x4x256x12 P25',@() nf_pac(gp,[4 5],[20 21],'method','MVL','mode','cross','nperm',25,'weighting','amplitude'),off,shipped;
'nf_tfspecparam 24x64x64x8',@() nf_tfspecparam(spl,'maxPeaks',0,'FoD','data'),off,shipped;
'hnf_robdiv 128x256x32', @() hnf_robdiv(vec,b),off,shipped;
'nf_onesample_TFCE 40x24x96 P20 (shipped kernels -> +clusters)',@() nf_onesample_TFCE(draw,0,[],1:24,1:96,20,0),shipped,native;
'nf_onesample_TFCE 40x24x96 P20 (all off -> all on)',@() nf_onesample_TFCE(draw,0,[],1:24,1:96,20,0),off,native;
};
fprintf('Runtime %s; sample %dx%dx%d Fs %.1f\n',version,size(a),fs);
for j=1:size(cases,1)
    measure(cases{j,:},3);
end
for dims={[64 256],[32 24 80]}
    sz=dims{1}; x=randn(sz); n=sz(1); row=(1:n)'; ch=[row,mod(row,n)+1,mod(row-2,n)+1];
    measure(['hnf_tfce_integrate ',mat2str(sz)],@() hnf_tfce_integrate(x,ch,[.66 2]),shipped,native,3);
end
% Larger 64-channel volume. Known exact population-level integrand for two
% equal-height opposite-signed components with no neighbours across halves.
sz=[64 32 256]; n=sz(1); row=(1:n)'; ch=[row,mod(row,n)+1,mod(row-2,n)+1];
x=ones(sz); x(:,:,129:end)=-1; h=((.5:49.5)/50);
expected=sign(x)*(numel(x)/2)^.66*(sum(h.^2)/50);
hnf_mex_available('scope',native);
clock=tic; result=hnf_tfce_integrate(x,ch,[.66 2]); seconds=toc(clock);
assert(max(abs(result(:)-expected(:)))<1e-10);
fprintf('LARGE plateau 64x32x256 voxels=%d %.6fs abs_error=%.4g\n',numel(x),seconds,max(abs(result(:)-expected(:))));
end

function measure(label,fn,control,candidate,n)
times=zeros(2,n);
for rep=0:n
    for mode=1+mod((0:1)+rep,2)
        if mode==1, hnf_mex_available('scope',control); else, hnf_mex_available('scope',candidate); end
        rng(1960); clock=tic; evalc('value=fn();'); elapsed=toc(clock);
        if rep>0, times(mode,rep)=elapsed; end
        if mode==1,before=compare_fields(value);else,after=compare_fields(value);end
    end
    assert(isequal(size(before),size(after)) && isequal(isnan(before),isnan(after)) && isequal(isinf(before),isinf(after)));
    good=isfinite(before)&isfinite(after); difference=abs(before(good)-after(good)); scale=max(1,max(abs(before(good))));
    err=max([0;difference(:)]); assert(err<1e-7*scale,label);
end
med=median(times,2);fprintf('%s   %.7f -> %.7f   %.3fx   max_abs %.8g; raw %s\n',label,med(1),med(2),med(1)/med(2),err,mat2str(times,5));
end

function v=compare_fields(result)
if ~isstruct(result), v=result(:); return; end
if isfield(result,'power_ap')
    v=[double(result.power_ap(:)); double(result.power_osc(:)); ...
       result.specparam.ap_offset(:); result.specparam.ap_exponent(:)];
elseif isfield(result,'phase')
    v=[double(result.power(:)); real(exp(1i*result.phase(:))); imag(exp(1i*result.phase(:)))];
elseif isfield(result,'PAC')
    v=result.PAC(:);
else
    v=[result.tstat(:);result.p_vals(:)];
end
end

function restore(p,r,policy,threads)
hnf_mex_available('scope',policy); rng(r);
if ~isempty(threads),fftw('threads',threads);end
path(p);
end
