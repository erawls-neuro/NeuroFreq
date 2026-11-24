function TF = nf_resample(TF,tVec,fVec,verbose)
% NF_RESAMPLE Interpolate TF measurements to a new time or frequency grid.
%   TF = nf_resample(TF,tVec,fVec) leaves an axis alone when its vector is
%   empty. tVec uses the same time units as TF.times. Fs is always in Hz.
%   Phase angles are unwrapped along the interpolated axis and then wrapped
%   back to radians. This interpolation does not apply an anti-alias filter.
if nargin < 1 || isempty(TF), error('nf_resample:TFRequired','Supply a TF set.'); end
if nargin < 2, tVec = []; end
if nargin < 3, fVec = []; end
if nargin < 4 || isempty(verbose), verbose = true; end
oldTimes = double(TF.times(:));
oldFreqs = double(TF.freqs(:));
targetTimes = oldTimes;
targetFreqs = oldFreqs;
if ~isempty(tVec)
    targetTimes = double(tVec(:));
    dt = diff(targetTimes);
    if isempty(dt) || any(~isfinite(dt)) || any(dt<=0) || ...
            any(abs(dt-mean(dt)) > 1e-6*mean(dt))
        error('nf_resample:UniformRequired','Supply at least two uniformly increasing times.');
    end
    if numel(oldTimes)<2, error('nf_resample:TimeAxis','Need two original time samples.'); end
    % Infer seconds vs milliseconds from the original axis and sample rate.
    unitsPerSecond = 1;
    if abs(log(median(diff(oldTimes))*TF.Fs/1000)) < abs(log(median(diff(oldTimes))*TF.Fs))
        unitsPerSecond = 1000;
    end
    TF.Fs = unitsPerSecond / mean(dt);
end
if ~isempty(fVec)
    targetFreqs = double(fVec(:));
    if any(~isfinite(targetFreqs)) || any(diff(targetFreqs)<=0)
        error('nf_resample:FrequencyAxis','Frequencies must be finite and strictly increasing.');
    end
end
if verbose, fprintf('[nf_resample]: %d time points, %d frequencies.\n',numel(targetTimes),numel(targetFreqs)); end
nCh = TF.nsensor;
nF = numel(oldFreqs);
nT = numel(oldTimes);
nObs = numel(TF.power)/(nCh*nF*nT);
fields = fieldnames(TF);
for k = 1:numel(fields)
    name = fields{k};
    if ~isnumeric(TF.(name)) || isempty(regexpi(name,'^(power|phase|wpli|itpc|complex)($|_)','once'))
        continue
    end
    X = TF.(name);
    [X,compact] = hnf_tf_grid(X,nCh,nF,nT,nObs,name);
    isAngle = ~isempty(regexpi(name,'^phase($|_)','once')) && isempty(strfind(lower(name), 'itpc')) && isempty(strfind(lower(name), 'itc'));
    if ~isempty(tVec), X = interpolate(X,oldTimes,targetTimes,3,isAngle); end
    if ~isempty(fVec), X = interpolate(X,oldFreqs,targetFreqs,2,isAngle); end
    if compact, X = squeeze(X); end
    TF.(name) = X;
end
if ~isempty(tVec)
    TF = hnf_event_grid(TF,oldTimes,targetTimes,1:nObs,[]);
    if isfield(TF,'phaseReference') && strcmp(TF.phaseReference,'epoch-start-carrier')
        TF.phaseReferenceSeconds = interp1(oldTimes,TF.phaseReferenceSeconds,double(tVec),'linear',NaN);
    end
    TF.times = double(tVec);
end
if ~isempty(fVec), TF.freqs = double(fVec); end
TF.power = double(TF.power);
end

function out = interpolate(X,axis,query,dim,isAngle)
order = [dim,setdiff(1:4,dim,'stable')];
X = permute(X,order);
dims = [size(X,2),size(X,3),size(X,4)];
X = reshape(X,numel(axis),[]);
if numel(axis)==1
    if ~isequal(query,axis), error('nf_resample:Singleton','Cannot interpolate a singleton axis.'); end
    values = X;
else
    if isAngle, X = unwrap(X,[],1); end
    native = numel(X)>=2048 && isa(X,'double') && isreal(X) && ~issparse(X) && ...
        exist('hnf_grid_interp_mex','file')==3 && ...
        (~exist('hnf_mex_available','file') || hnf_mex_available('hnf_grid_interp_mex')) && ...
        all(isfinite(axis)) && all(diff(axis)>0) && all(isfinite(X(:)));
    if native
        try
            values = hnf_grid_interp_mex(axis,query,X);
        catch err
            if ~exist('hnf_mex_load_error','file') || ~hnf_mex_load_error(err), rethrow(err); end
            native = false;
        end
    end
    if ~native
        % interp1 may form a slope that overflows for finite endpoints on
        % a narrow interval. Compute each fraction before applying it.
        if isfloat(X) && ~issparse(X) && all(isfinite(axis)) && all(diff(axis)>0) && all(isfinite(X(:)))
            [~,bin] = histc(query,axis);
            valid = bin > 0;
            lo = min(bin(valid),numel(axis)-1);
            a = axis(lo); b = axis(lo+1); q = query(valid);
            width = b-a;
            w = (q-a)./width;
            huge = ~isfinite(width);
            w(huge) = (.5*q(huge)-.5*a(huge))./(.5*b(huge)-.5*a(huge));
            values = nan(numel(query),size(X,2),'like',X);
            % Bound the temporary gathered endpoints on multi-trial datasets.
            block = max(1,floor(131072/max(1,numel(lo))));
            for first = 1:block:size(X,2)
                cols = first:min(size(X,2),first+block-1);
                values(valid,cols) = hnf_roblerp(X(lo,cols),X(lo+1,cols),w);
            end
        else
            values = interp1(axis,X,query,'linear',NaN);
        end
    end
    if isAngle, values = angle(exp(1i*values)); end
end
out = ipermute(reshape(values,[numel(query),dims]),order);
end
