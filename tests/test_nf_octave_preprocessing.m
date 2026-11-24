function test_nf_octave_preprocessing()
% Run after nf_run. Table tests need tablicious or datatypes; graphs need a toolkit.
% Test the source contracts behind the Octave specializations, without test shims.
assert(exist('OCTAVE_VERSION','builtin')~=0);
[p,x,w,h]=hnf_mara_peaks([-2 1 1 -1 2 -2],1:6);
assert(isequal(p,[1 2]) && isequal(x,[2 5]) && isequal(h,[2 4]));
assert(max(abs(w-[1+1/3+1/2 2/3+1/2]))<1e-12);
csv=[tempname '.csv'];
fid=fopen(csv,'wb');
fprintf(fid,'trials.thisN,key,participant,rt\r\n0,"left, ""up""",007,.25\r\n1,"two\nlines",007,\r\n');
fclose(fid); cleaner=onCleanup(@() delete(csv)); %#ok<NASGU>
[t,info]=nf_readpsychopy(csv,'Verbose',false);
assert(isequal(info.sourceRows,[1;2]) && isequal(t.participant,{'007';'007'}));
assert(isequal(t.key,{'left, "up"';sprintf('two\nlines')}) && isnan(t.rt(2)));
b=nf_makebehavior(t,-1);
assert(isnan(b(1).rt_1back) && b(2).rt_1back==.25);
% Compare the same recordings through a range where raw squares overflow
% or underflow. QC power density is in dB while PCA RMS is in microvolts.
state=rng; restore=onCleanup(@() rng(state)); %#ok<NASGU>
rng(417,'twister'); e.data=randn(3,256,2); e.nbchan=3;e.pnts=256;e.trials=2;
e.srate=128;e.xmin=0;e.xmax=255/128;e.chanlocs=struct('labels',{'Fz','Cz','Oz'});
reference=nf_eegquality(e,e,'plot',false);
for gain=[1e-190,1e190]
    x=e;x.data=x.data*gain;
    q=nf_eegquality(x,x,'plot',false);
    a=reference.metrics.preclean.artifact.measurements;
    z=q.metrics.preclean.artifact.measurements;
    assert(max(abs(z.rmsMicrovolts(:)/gain-a.rmsMicrovolts(:)))<1e-10);
    assert(max(abs(z.spatialCorrelation(:)-a.spatialCorrelation(:)))<1e-10);
    assert(max(abs(z.muscleRatioDb(:)-a.muscleRatioDb(:)))<1e-8);
    psd=q.metrics.preclean.psd.channelDb-reference.metrics.preclean.psd.channelDb;
    assert(max(abs(psd(:)-20*log10(gain)))<1e-8);
    assert(q.metrics.pca.available && abs(q.metrics.pca.candidatePreRms/gain-reference.metrics.pca.candidatePreRms)<1e-10);
end
fprintf('PASS signed peak widths/prominence, quoted CSV, metadata/lag alignment, QC gain and covariance invariance.\n');
end
