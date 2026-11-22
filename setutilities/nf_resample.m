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
    isAngle = ~isempty(regexpi(name,'^phase($|_)','once')) && ~contains(lower(name),'itpc') && ~contains(lower(name),'itc');
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
    values = interp1(axis,X,query,'linear',NaN);
    if isAngle, values = angle(exp(1i*values)); end
end
out = ipermute(reshape(values,[numel(query),dims]),order);
end
