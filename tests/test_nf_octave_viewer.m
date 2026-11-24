function test_nf_octave_viewer()
% Exercise scoped dispatch and graphics callbacks after the temporary path closes.
root = fileparts(fileparts(mfilename('fullpath')));
before = path;
cleanup = onCleanup(@() path(before));
addpath(root);
TF.times = [-.2 .05 .3]; TF.freqs = [11 5];
TF.chanlocs = struct('labels',{'Fz','C4','C3','EOG'}, ...
    'theta',{0,90,-90,[]},'radius',{.4,.3,.3,[]});
TF.power = reshape(1:48,4,2,3,2);
TF.power(4,1,3,2) = NaN;
% No desktop is needed for creating figure objects and invoking callbacks.
if exist('OCTAVE_VERSION','builtin') && ~isempty(available_graphics_toolkits())
    visible = get(0,'DefaultFigureVisible');
    restore = onCleanup(@() set(0,'DefaultFigureVisible',visible));
    set(0,'DefaultFigureVisible','off');
    initialPath = path;
    fig = nf_run('nf_viewerapp',TF,'power');
    closer = onCleanup(@() close(fig));
    assert(strcmp(path,initialPath));
    callback = get(fig,'KeyPressFcn');
    callback(fig,struct('Key','m','Modifier',{{}}));
    assert(isequal(getappdata(fig,'NeuroFreqSelection'),[4 1 3 2]));
    callback(fig,struct('Key','n','Modifier',{{'shift'}}));
    assert(isequal(getappdata(fig,'NeuroFreqSelection'),[4 2 1 1]));
    markers = findall(fig,'Tag','NeuroFreqSensor');
    assert(numel(markers)==3);
    pick = get(markers(1),'ButtonDownFcn'); sensor = get(markers(1),'UserData');
    pick(markers(1),[]);
    selection = getappdata(fig,'NeuroFreqSelection');
    assert(selection(1)==sensor);
    drawnow();
end
% The one-sided discrete Wigner–Ville distribution conserves analytic energy
% and places an integer-bin complex sinusoid at its known frequency.
Fs=128; n=63; bin=9;
x=exp(2i*pi*bin*(0:n-1)/(2*n));
[d,f,t] = nf_run(@() hnf_wvd(x,Fs));
assert(isequal(size(d),[n n]) && isequal(t,(0:n-1)/Fs));
assert(max(abs(sum(d,1)-abs(x).^2))<1e-12);
[peak,idx]=max(d(:,(n+1)/2));
assert(abs(peak-1)<1e-12 && abs(f(idx)-bin*Fs/(2*n))<1e-12);
high = nf_run(@() hnf_wvd(1e154*x,Fs));
assert(all(isfinite(high(:))) && abs(high(idx,(n+1)/2)/1e308-1)<1e-12);
fprintf('Octave viewer and demo numerical checks passed.\n');
end
