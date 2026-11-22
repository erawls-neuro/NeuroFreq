function test_nf_pac()
% TEST_NF_PAC  Scalar reference calculations for the six PAC estimators.
% Run addpath('tests'); test_nf_pac from the NeuroFreq root. If compiled
% hnf_*_mex functions are on the path, both native and MATLAB paths are run.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; oldRng = rng;
cleanup = onCleanup(@() restore(oldPath,oldRng)); %#ok<NASGU>
addpath(genpath(fullfile(root,'phaseamplitudecoupling')));
addpath(fullfile(root,'setutilities','hnf_setutilities'),fullfile(root,'timefrequency','hnf_timefrequency'));
names = {'hnf_mvl_mex','hnf_mi_mex','hnf_circ_corrcl_mex','hnf_surr_timeshift_mex'};
native = all(cellfun(@(f) exist(f,'file') == 3,names));
rng(33,'twister');
[f,t,channels,trials] = deal([4 6 22 30 44],(0:71)/128,2,9);
TF = struct('freqs',f,'times',t,'Fs',128,'nsensor',channels,'ntrls',trials, ...
    'chanlocs',struct('labels',{'A','B'}),'scale','linear');
TF.phase = 2*pi*rand(channels,numel(f),numel(t),trials)-pi;
amp = 1.3+0.7*cos(TF.phase(:,1,:,:))+0.3*rand(size(TF.phase));
TF.power = amp.^2;
checkKernels(native);
for useMex = unique([false native])
    for m = {'MVL','MI','ERPAC'}
        method = m{1};
        args = {'method',method,'nbin',4,'nperm',0,'use_mex',logical(useMex)};
        pac = nf_pac(TF,[4 6],[22 44],args{:},'weighting','amplitude');
        cross = nf_pac(TF,[4 6],[22 44],args{:},'weighting','amplitude', ...
            'mode','cross','phaseCh',[2 1],'ampCh',[1 2]);
        for ic=1:channels
            for jc=1:channels
                for lo=1:2
                    for hi=1:3
                        phase = reshape(TF.phase(ic,lo,:,:),numel(t),trials);
                        weight = reshape(amp(jc,hi+2,:,:),numel(t),trials);
                        expected = reference(method,phase,weight,4);
                        got = reshape(cross.PAC(3-ic,jc,lo,hi,:),1,[]);
                        checkClose(got,expected);
                        if ic == jc
                            checkClose(reshape(pac.PAC(ic,lo,hi,:),1,[]),expected);
                        end
                    end
                end
            end
        end
        checkNull(method,TF,useMex);
    end
end
% Historical power weighting and envelope weighting are separate estimands.
one = TF; one.nsensor = 1; one.chanlocs = one.chanlocs(1);
one.phase = squeeze(TF.phase(1,:,:,:)); one.power = squeeze(TF.power(1,:,:,:));
raw = nf_pac(one,[4 4],[22 22],'nperm',0,'mode','cross','toi',[t(9) t(23)]);
phase = reshape(one.phase(1,9:23,:),15,trials);
power = reshape(one.power(3,9:23,:),15,trials);
checkClose(reshape(raw.PAC,1,[]),reference('MVL',phase,power,4));
assert(isequal(raw.times,t(9:23)) && strcmp(raw.weighting,'power'));
one.scale = 'log10'; one.power = log10(one.power);
converted = nf_pac(one,[4 4],[22 22],'nperm',0,'mode','cross','toi',[t(9) t(23)]);
checkClose(converted.PAC,raw.PAC);
singleTrial = TF;
singleTrial.power = singleTrial.power(:,:,:,1);
singleTrial.phase = singleTrial.phase(:,:,:,1); singleTrial.ntrls = 1;
for m={'MVL','MI','ERPAC'}
    value = nf_pac(singleTrial,[4 4],[22 22],'nperm',0,'mode','cross','nbin',4,'method',m{1});
    if strcmp(m{1},'ERPAC'), assert(all(isnan(value.PAC(:)))); else, assert(all(isfinite(value.PAC(:)))); end
end
% A demodulated carrier has constant coefficient phase for a pure tone.
% Local phase must rotate in time before testing amplitude coupling.
carrier=TF; carrier.freqs=[4 22]; carrier.nsensor=1; carrier.chanlocs=[];
carrier.ntrls=1; carrier.times=(0:127)/128;
carrier.phase=zeros(2,128); carrier.power=ones(2,128);
carrier.power(2,:)=(2+cos(2*pi*4*carrier.times)).^2;
carrier.phaseReference='epoch-start-carrier'; carrier.phaseReferenceSeconds=carrier.times;
coupling=nf_pac(carrier,[4 4],[22 22],'method','MVL','nperm',0,'weighting','amplitude');
checkClose(coupling.PAC,0.5);
fprintf('PAC tests passed (compiled MEX backends exercised: %d).\n',native);
end

function checkKernels(native)
% Uniform sampled bin means give zero. A single positive bin gives one.
a = repmat([-2.6;-1;0.7;2],4,1);
amp = [ones(16,1),double(a<-2),2+cos(a)];
phase = repmat(a,1,3);
mi = hnf_mi_vec(phase,amp,4);
checkClose(mi,[0 1 scalarMI(a,amp(:,3),4)]);
assert(isnan(hnf_mi_vec(a,zeros(size(a)),4)));
assert(isnan(hnf_mi_vec(zeros(size(a)),ones(size(a)),4)));  % no phase coverage
edge = [-pi;pi;-1;0.3;2.1];
checkClose(hnf_mi_vec(edge,[1;3;2;2;2],4),0);
bad = amp; bad(3,2)=NaN;
miBad = hnf_mi_vec(phase,bad,4); assert(isnan(miBad(2))); checkClose(miBad([1 3]),mi([1 3]));
[r,p] = hnf_circ_corrcl_vec(a,cos(a));
checkClose([r p],[1 exp(-numel(a)/2)]);
assert(isnan(hnf_circ_corrcl_vec(a,0.1*ones(size(a)))));
assert(isnan(hnf_circ_corrcl_vec(repmat([0;1],8,1),amp(:,3))));
checkClose(hnf_mvl_vec(zeros(4,1),ones(4,1)*1e307)/1e307,1);
x = reshape(1:36,12,3);
rng(5); shifts = randi([2,11],1,3);
expected = x;
for j=1:3, expected(:,j)=circshift(x(:,j),shifts(j)); end
rng(5); assert(isequal(hnf_surr_timeshift_vec(x,2,false),expected));
if native
    checkClose(hnf_mi_mex(phase,amp,4),mi);
    assert(isnan(hnf_mi_mex(zeros(size(a)),ones(size(a)),4)));
    checkClose(hnf_mi_mex(edge,[1;3;2;2;2],4),0);
    checkClose(hnf_mi_mex(phase,bad,4),miBad);
    checkClose(hnf_mvl_mex(phase,amp),hnf_mvl_vec(phase,amp));
    [r1,p1]=hnf_circ_corrcl_mex(a,cos(a)); checkClose([r1,p1],[r,p]);
    assert(isequal(hnf_surr_timeshift_mex(x,shifts),expected));
    rng(5); assert(isequal(hnf_surr_timeshift_vec(x,2,true),expected));
end
end

function checkNull(method,TF,useMex)
% Compute a seeded null with circshift, randperm and independent statistics.
phase = reshape(TF.phase(1,1,:,:),numel(TF.times),TF.ntrls);
weight = reshape(TF.power(1,3,:,:),numel(TF.times),TF.ntrls);
observed = reference(method,phase,weight,4);
n = 9; null = zeros(n,numel(observed)); rng(91,'twister');
for i=1:n
    if strcmp(method,'ERPAC')
        surrogate = weight(:,randperm(TF.ntrls));
    else
        shifts = randi([1,size(weight,1)-1],1,TF.ntrls);
        surrogate = weight;
        for j=1:TF.ntrls, surrogate(:,j)=circshift(weight(:,j),shifts(j)); end
    end
    null(i,:) = reference(method,phase,surrogate,4);
end
expected = (observed-mean(null,1))./std(null,0,1);
one = TF; one.nsensor = 1; one.chanlocs = one.chanlocs(1);
one.phase = squeeze(TF.phase(1,:,:,:)); one.power = squeeze(TF.power(1,:,:,:));
for mode = {'within','cross'}
    rng(91,'twister');
    got = nf_pac(one,[4 4],[22 22],'method',method,'nperm',n,'nbin',4, ...
        'use_mex',logical(useMex),'mode',mode{1});
    checkClose(reshape(got.PAC,1,[]),expected);
end
end

function result=reference(method, phase, weights, bins)
% Independent scalar formulas. ERPAC is sqrt(R^2) from a linear fit with
% intercept and the two periodic predictors, across trials at each time.
if strcmp(method,'ERPAC')
    result = nan(1,size(phase,1));
    for t=1:numel(result)
        x=weights(t,:)'; angles=phase(t,:)';
        design=[ones(size(x)),cos(angles),sin(angles)];
        fit=design*(design\x);
        result(t)=sqrt(max(0,1-sum((x-fit).^2)/sum((x-mean(x)).^2)));
    end
else
    result = nan(1,size(phase,2));
    for j=1:numel(result)
        if strcmp(method,'MI')
            result(j)=scalarMI(phase(:,j),weights(:,j),bins);
        else
            result(j)=abs(sum(weights(:,j).*exp(1i*phase(:,j)))/size(phase,1));
        end
    end
end
end

function result=scalarMI(phase,weights,bins)
angles=mod(phase+pi,2*pi)-pi; edges=linspace(-pi,pi,bins+1); p=zeros(1,bins);
for j=1:bins
    mask=angles>=edges(j) & angles<edges(j+1);
    if ~any(mask), result=NaN; return; end
    p(j)=mean(weights(mask));
end
p=p/sum(p); positive=p>0;
result=1+sum(p(positive).*log(p(positive)))/log(bins);
end

function checkClose(actual,expected)
assert(isequal(size(actual),size(expected)) && isequal(isnan(actual),isnan(expected)));
keep = isfinite(expected);
assert(all(abs(actual(keep)-expected(keep))<2e-10));
end

function restore(p,r)
path(p); rng(r);
end
