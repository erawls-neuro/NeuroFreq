function data = nf_viewer_data(TF, mode)
% Convert NF's stored maps to [sensor frequency time displayed observation].
% Only the selected measure is pooled. A single-trial set retains every trial.
if nargin < 2 || isempty(mode)
    if isfield(TF,'tstat'), mode = 'stat'; else, mode = 'power'; end
end
mode = lower(char(mode));
data.angular = strcmp(mode,'phase');
data.measure = mode;
data.t = double(TF.times(:).');
validateattributes(data.t,{'numeric'},{'real','finite','nonempty','increasing'});
% Current NF reports seconds. Older datasets using milliseconds can be
% identified by their time step relative to the declared sampling rate.
isMs = numel(data.t)>1 && isfield(TF,'Fs') && isscalar(TF.Fs) && TF.Fs>0 && ...
    abs(log(median(diff(data.t))*TF.Fs/1000)) < abs(log(median(diff(data.t))*TF.Fs));
if ~isMs, data.t = 1000*data.t; end
[data.f, order] = sort(double(TF.freqs(:).'));
validateattributes(data.f,{'numeric'},{'real','finite','nonempty','increasing'});
data.chanlocs = [];
if isfield(TF,'chanlocs'), data.chanlocs = TF.chanlocs; end
if isempty(data.chanlocs)
    if isfield(TF,'nsensor'), n = TF.nsensor; else, n = 1; end
    data.chanlocs = repmat(struct('labels',''),1,n);
end
for c = 1:numel(data.chanlocs)
    if ~isfield(data.chanlocs,'labels') || isempty(data.chanlocs(c).labels)
        data.chanlocs(c).labels = sprintf('Sensor %d',c);
    end
end
data.pillow = isfield(TF,'custommontage') && strcmp(TF.custommontage,'pillow64');
dims = [numel(data.chanlocs), numel(data.f), numel(data.t)];
switch mode
    case {'power','phase','itpc'}
        fields = {TF.(mode)};
    case 'shap'
        fields = {TF.SHAP};
    case 'sprint'
        if isfield(TF,'power_ap')
            fields = {TF.power_ap, TF.power_osc};
        else
            fields = {TF.SPRiNT.ap_power, TF.SPRiNT.osc_power};
        end
    case 'cpm'
        if isfield(TF,'power_evoked')
            fields = {TF.power_evoked, TF.power_induced};
        else
            fields = {TF.cpm.PLmean, TF.cpm.NPLmean};
        end
    case 'erp'
        if isfield(TF,'power_induced')
            fields = {TF.power_induced, TF.power_evoked};
        else
            fields = {TF.erprem.erprempow, TF.erprem.erppow};
        end
    case 'stat'
        % GLM/LME TFCE numeric results store coefficient first. Single-trial
        % regression exposes [C,F,T,K], even when its method was TFCE.
        if ~isequal(size(TF.tstat),size(TF.p_vals))
            error('nf_viewerapp:Statistics','t-stat and p-value maps must have identical axes.');
        end
        coefficientFirst = isfield(TF,'fixedNames') && ~isfield(TF,'primary_measure') && ...
            ~isfield(TF,'corrcoef') && ~isfield(TF,'tfce_tstat');
        if coefficientFirst
            k = numel(TF.fixedNames);
            stats = permute(grid(TF.tstat,[k dims]),[2 3 4 1]);
            mask = permute(grid(TF.p_vals,[k dims]),[2 3 4 1]);
        else
            stats = grid(TF.tstat,[dims numel(TF.tstat)/prod(dims)]);
            mask = grid(TF.p_vals,[dims numel(TF.tstat)/prod(dims)]);
        end
        fields = {stats};
        data.statflag = true;
        data.mask = mask(:,order,:,:) < .05;
    otherwise
        error('nf_viewerapp:Mode','Unknown viewer mode %s.',mode);
end
for k = 1:numel(fields)
    map = fields{k};
    if isfield(TF,'aggtype') && strcmpi(TF.aggtype,'avg') && ~strcmp(mode,'stat')
        map = grid(map,[dims numel(TF.trial_condindex)]);
        indices = unique(TF.trial_condindex,'stable');
        grouped = cell(1,numel(indices));
        for c = 1:numel(indices)
            grouped{c} = pool(map(:,:,:,TF.trial_condindex == indices(c)),4,data.angular);
        end
        map = cat(4,grouped{:});
    elseif isfield(TF,'nsubjects') && ~isfield(TF,'aggtype') && ~strcmp(mode,'stat')
        nDisplayed = numel(map)/(TF.nsubjects*prod(dims));
        map = grid(map,[TF.nsubjects dims nDisplayed]);
        map = pool(map,1,data.angular);
        map = reshape(map,[dims nDisplayed]);
    else
        map = grid(map,[dims numel(map)/prod(dims)]);
    end
    fields{k} = map(:,order,:,:);
end
data.ctf = cat(4,fields{:});
data.trls = size(data.ctf,4);
end

function map = grid(map,dims)
% Squeezed axes may be restored, but transposed maps of equal size must
% never silently reassign values or p-values to another sensor/bin.
if ~isnumeric(map) || ~isreal(map) || any(~isfinite(dims) | dims<1 | dims~=fix(dims))
    error('nf_viewerapp:Map','Expected a real numeric map on the declared TF grid.');
end
actual = size(map);
if ~isequal(actual(actual~=1),dims(dims~=1))
    error('nf_viewerapp:Map','Map of size %s is incompatible with axes %s.',mat2str(actual),mat2str(dims));
end
map = reshape(map,dims);
end

function result = pool(data,dim,angular)
if angular
    result = mean(exp(1i*data),dim,'omitnan');
    undefined = abs(result) < 8*eps(class(result));
    result = angle(result);
    result(undefined) = NaN;
else
    result = mean(data,dim,'omitnan');
end
end
