function fig = nf_viewerapp(TF,mode,palette)
% NF_VIEWERAPP  Explore a NeuroFreq TF structure in a regular Octave figure.
%   fig = nf_run('nf_viewerapp',TF,'phase');
%   fig = nf_viewerapp(TF,mode,colormap);
% Modes and the optional colormap have the same meanings as in nf_viewerapp.
% Click the TF map to pick a time/frequency bin, or a site to pick a sensor.
% Arrows step through time/frequency. A/D select a sensor, Z/X an observation.
% M/N find the global maximum/minimum; Shift+M/N restrict the search to the
% current sensor. L toggles the log-frequency axis, C toggles the crosshairs.
% In a noninteractive graphics toolkit only the plot is available.
if nargin < 2, mode = []; end
if nargin < 3 || isempty(palette), palette = 'viridis'; end
data = nf_viewer_data(TF,mode);
if ischar(palette) && strcmpi(palette,'parula') && isempty(which('parula'))
    warning('NeuroFreq:OctaveColormap','Octave has no parula colormap. Using viridis.');
    palette = 'viridis';
end
[~,f] = min(abs(data.f-6));
[~,t] = min(abs(data.t-300));
% Selection follows the source array's [sensor,frequency,time,observation] order.
selection = [1 f t 1];
count = [numel(data.chanlocs) numel(data.f) numel(data.t) data.trls];
showCrosshairs = true;
frequencyScale = 'linear';
if data.angular
    limits = [-pi pi];
else
    % Fix one scale across sensors/observations so their colors are comparable.
    finite = data.ctf(isfinite(data.ctf));
    bound = 1;
    if ~isempty(finite)
        bound = max(abs(prctile(finite,[1 99])));
        if bound==0, bound = 1; end
    end
    limits = [-bound bound];
    clear finite
end
fig = figure('Name','NeuroFreq TF viewer','NumberTitle','off');
try
    tf = subplot(2,2,1,'Parent',fig); topo = subplot(2,2,2,'Parent',fig);
    time = subplot(2,2,3,'Parent',fig); freq = subplot(2,2,4,'Parent',fig);
    if ischar(palette)
        colors = feval(palette,256);
    else
        colors = palette;
    end
    set(fig,'Colormap',colors,'KeyPressFcn',@navigate);
    redraw();
catch err
    close(fig);
    rethrow(err)
end

    function navigate(~,event)
        switch lower(event.Key)
            case 'rightarrow', selection(3) = selection(3)+1;
            case 'leftarrow', selection(3) = selection(3)-1;
            case 'uparrow', selection(2) = selection(2)+1;
            case 'downarrow', selection(2) = selection(2)-1;
            case 'a', selection(1) = selection(1)-1;
            case 'd', selection(1) = selection(1)+1;
            case 'z', selection(4) = selection(4)-1;
            case 'x', selection(4) = selection(4)+1;
            case 'c', showCrosshairs = ~showCrosshairs;
            case 'l'
                if any(data.f<=0), return; end
                if strcmp(frequencyScale,'linear'), frequencyScale = 'log';
                else, frequencyScale = 'linear'; end
            case {'m','n'}
                sensorOnly = isfield(event,'Modifier') && any(strcmpi(event.Modifier,'shift'));
                values = data.ctf;
                if sensorOnly, values = values(selection(1),:,:,:); end
                finite = find(isfinite(values));
                if isempty(finite), return; end
                if strcmpi(event.Key,'m'), [~,idx] = max(values(finite));
                else, [~,idx] = min(values(finite)); end
                [s,fbin,tbin,obs] = ind2sub(size(values),finite(idx));
                if sensorOnly, s = selection(1); end
                selection = [s fbin tbin obs];
            otherwise, return
        end
        selection = min(count,max(1,selection));
        redraw();
    end

    function pickBin(~,~)
        p = get(tf,'CurrentPoint');
        [~,selection(2)] = min(abs(data.f-p(1,2)));
        [~,selection(3)] = min(abs(data.t-p(1,1)));
        redraw();
    end

    function pickSite(source,~)
        selection(1) = get(source,'UserData');
        redraw();
    end

    function redraw()
        s = selection(1); fi = selection(2); ti = selection(3); o = selection(4);
        plane = reshape(data.ctf(s,:,:,o),count(2),count(3));
        cla(tf); hold(tf,'on');
        if count(2)>1 && count(3)>1
            surface(tf,data.t,data.f,zeros(size(plane)),plane,'EdgeColor','none', ...
                'FaceColor','interp','ButtonDownFcn',@pickBin); view(tf,2);
            if isfield(data,'mask')
                mask = reshape(data.mask(s,:,:,o),size(plane));
                if any(mask(:)) && ~all(mask(:))
                    contour(tf,data.t,data.f,double(mask),[.5 .5],'k','HitTest','off');
                end
            end
        else
            [ts,fs] = meshgrid(data.t,data.f);
            scatter(tf,ts(:),fs(:),32,plane(:),'filled','ButtonDownFcn',@pickBin);
        end
        axis(tf,'tight');
        if showCrosshairs
            line(tf,[data.t(ti) data.t(ti)],ylim(tf),'Color','k','HitTest','off');
            line(tf,xlim(tf),[data.f(fi) data.f(fi)],'Color','k','HitTest','off');
        end
        set(tf,'ButtonDownFcn',@pickBin,'YScale',frequencyScale);
        title(tf,sprintf('%s, observation %d/%d',data.chanlocs(s).labels,o,count(4)),'Interpreter','none');
        xlabel(tf,'Time (ms)'); ylabel(tf,'Frequency (Hz)'); colorbar(tf);
        if data.angular
            % Collapse phase on the unit circle and leave cancellation undefined.
            p = exp(1i*plane);
            series = mean(p,1,'omitnan'); spectrum = mean(p,2,'omitnan');
            series(abs(series)<8*eps(class(series))) = NaN;
            spectrum(abs(spectrum)<8*eps(class(spectrum))) = NaN;
            series = angle(series); spectrum = angle(spectrum);
        else
            series = mean(plane,1,'omitnan'); spectrum = mean(plane,2,'omitnan');
        end
        plot(time,data.t,series); xlabel(time,'Time (ms)'); ylabel(time,data.measure);
        plot(freq,data.f,spectrum); xlabel(freq,'Frequency (Hz)'); ylabel(freq,data.measure);
        scalp = data.ctf(:,fi,ti,o);
        drawScalp(topo,data.chanlocs,scalp,data.pillow,data.angular,s,@pickSite);
        title(topo,sprintf('%g ms, %g Hz, %s = %.4g',data.t(ti),data.f(fi),data.measure,scalp(s)));
        set(freq,'XScale',frequencyScale);
        set([tf topo],'CLim',limits);
        setappdata(fig,'NeuroFreqSelection',selection);
    end
end

function drawScalp(ax,locs,values,pillow,angular,selected,pick)
% EEGLAB polar coordinates use theta=0 at the nose and +90 to the right.
if pillow, first = 'X'; second = 'Y'; else, first = 'theta'; second = 'radius'; end
x = nan(size(values)); y = x;
for c = 1:numel(values)
    if ~isfield(locs,first) || ~isfield(locs,second), continue; end
    a = locs(c).(first); b = locs(c).(second);
    if ~isnumeric(a) || ~isscalar(a) || ~isnumeric(b) || ~isscalar(b), continue; end
    if pillow, x(c) = b; y(c) = a;
    else, x(c) = b*sind(a); y(c) = b*cosd(a); end
end
valid = find(isfinite(x) & isfinite(y) & isfinite(values));
cla(ax); hold(ax,'on');
if isempty(valid)
    text(ax,.5,.5,'No mapped sensors at this bin','Units','normalized','HorizontalAlignment','center');
else
    if numel(valid)>=3 && rank([ones(numel(valid),1) x(valid) y(valid)])==3
        if pillow
            xi = linspace(min(x(valid)),max(x(valid)),67);
            yi = linspace(min(y(valid)),max(y(valid)),67);
        else
            radius = max([.5; hypot(x(valid),y(valid))]);
            xi = linspace(-radius,radius,67); yi = xi;
        end
        [X,Y,Z] = nf_viewer_scalp(x(valid),y(valid),values(valid),xi,yi,angular);
        if ~pillow, Z(hypot(X,Y)>radius) = NaN; end
        height = zeros(size(Z)); height(~isfinite(Z)) = NaN;
        surface(ax,X,Y,height,Z,'EdgeColor','none','FaceColor','interp'); view(ax,2);
    end
    for c = valid(:).'
        % Display the original sensor index even with incomplete coordinates.
        dot = scatter(ax,x(c),y(c),36,values(c),'filled','MarkerEdgeColor','k');
        set(dot,'UserData',c,'ButtonDownFcn',pick,'Tag','NeuroFreqSensor');
        if c==selected, set(dot,'MarkerEdgeColor','m','LineWidth',2); end
    end
    colorbar(ax);
end
axis(ax,'equal'); axis(ax,'off');
end
