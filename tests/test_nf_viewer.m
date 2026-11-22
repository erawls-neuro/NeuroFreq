function test_nf_viewer()
% App package integrity and numerical data mapping. In desktop MATLAB this
% also opens the shipped app and executes callbacks after nf_run has returned.
root = fileparts(fileparts(mfilename('fullpath')));
original = path;
restore = onCleanup(@() path(original)); %#ok<NASGU>
addpath(root,fullfile(root,'private'));
TF.freqs = [4 2];
TF.times = [-.1 0 .1];
TF.Fs = 10;
TF.chanlocs = struct('labels',{'Fz','Cz'},'theta',{0,90},'radius',{.2,.4});
TF.power = reshape(1:48,[2 2 3 4]);
TF.phase = TF.power*.05;
TF.itpc = TF.power/50;
map = nf_viewer_data(TF,'power');
assert(isequal(map.ctf,TF.power(:,[2 1],:,:)) && isequal(map.t,[-100 0 100]));
TF.times = map.t;
assert(isequal(nf_viewer_data(TF,'itpc').t,map.t));
TF.times = [-.1 0 .1];
TF.aggtype = 'avg'; TF.trial_condindex = [1 2 1 2];
grouped = nf_viewer_data(TF);
assert(isequal(grouped.ctf(:,:,:,1),mean(TF.power(:,[2 1],:,[1 3]),4)));
TF.phase(:,:,:,1) = pi-.05; TF.phase(:,:,:,3) = -pi+.05;
circular = nf_viewer_data(TF,'phase');
assert(max(abs(abs(reshape(circular.ctf(:,:,:,1),1,[]))-pi)) < 1e-12);
TF = rmfield(TF,{'aggtype','trial_condindex'});
TF.power_evoked = 2*TF.power; TF.power_induced = 3*TF.power;
TF.power_ap = TF.power/2; TF.power_osc = TF.power/3;
cpm = nf_viewer_data(TF,'cpm'); erp = nf_viewer_data(TF,'erp');
spr = nf_viewer_data(TF,'sprint');
assert(cpm.trls==8 && erp.trls==8 && spr.trls==8);
assert(isequal(cpm.ctf(:,:,:,1:4),erp.ctf(:,:,:,5:8)));
assert(isequal(erp.ctf(:,:,:,1:4),TF.power_induced(:,[2 1],:,:)));
assert(isequal(spr.ctf(:,:,:,1:4),TF.power_ap(:,[2 1],:,:)));

% GLM/LME maps place coefficients first. Single-trial regressions keep the
% observation axis last. The p-value overlay must follow the identical map.
TF.fixedNames = {'one','two','three','four'};
TF.tstat = permute(TF.power,[4 1 2 3]); TF.p_vals = TF.tstat/100;
glm = nf_viewer_data(TF,'stat');
assert(isequal(glm.ctf,map.ctf) && isequal(glm.mask,glm.ctf<5));
TF.primary_measure = 'power'; TF.tstat = TF.power; TF.p_vals = TF.power/100;
reg = nf_viewer_data(TF,'stat');
assert(isequal(reg.ctf,glm.ctf) && isequal(reg.mask,glm.mask));
% Duplicate coordinates, wrapping at +/-pi, and a rank-deficient montage.
[~,~,z] = nf_viewer_scalp([0 0 1 0],[0 0 0 1], ...
    [pi-.05 -pi+.05 pi-.02 -pi+.02],0:.5:1,0:.5:1,true);
assert(abs(abs(z(1,1))-pi)<1e-10);
[~,~,missing] = nf_viewer_scalp([0 1],[0 1],[2 3],0:1,0:1,false);
assert(all(isnan(missing(:))));

% The editable representation and the runtime XML both have to change when
% a callback is edited. Otherwise opening and saving the app reverts the fix.
temp = tempname;
mkdir(temp);
clean = onCleanup(@() rmdir(temp,'s')); %#ok<NASGU>
unzip(fullfile(root,'nf_viewerapp.mlapp'),temp);
runtime = fileread(fullfile(temp,'matlab','document.xml'));
saved = load(fullfile(temp,'appdesigner','appModel.mat'),'code');
code = saved.code;
section = strjoin(reshape(code.EditableSectionCode,1,[]),sprintf('\n'));
assert(~isempty(strfind(runtime,section))); %#ok<STREMP>
callbacks = [reshape(code.Callbacks,1,[]) code.StartupCallback];
for cb = callbacks
    body = strjoin(reshape(cb.Code,1,[]),sprintf('\n'));
    signature = regexp(runtime,['function ' cb.Name '\([^\n]*\)\n'],'end','once');
    assert(~isempty(signature) && strncmp(runtime(signature+1:end),body,numel(body)));
end

if ~exist('OCTAVE_VERSION','builtin') && usejava('desktop')
    % Only the public root is needed during and after construction.
    rmpath(fullfile(root,'private'));
    TF = rmfield(TF,{'tstat','p_vals','fixedNames','primary_measure'});
    before = path;
    app = nf_run('nf_viewerapp',TF,'power');
    closer = onCleanup(@() delete(app)); %#ok<NASGU>
    assert(strcmp(path,before));
    app.GoToGlobalMax.ButtonPushedFcn(app.GoToGlobalMax,[]);
    assert(app.tf_pickelectrode.Value==2 && str2double(app.tf_trial2plot.Value)==4);
    app.tf_trialM.ButtonPushedFcn(app.tf_trialM,[]);
    assert(str2double(app.tf_trial2plot.Value)==3);
    fprintf('App Designer launch and callback smoke checks passed.\n');
else
    fprintf('Skipped MATLAB desktop GUI launch.\n');
end
fprintf('Viewer numerical and editable/runtime checks passed.\n');
end
