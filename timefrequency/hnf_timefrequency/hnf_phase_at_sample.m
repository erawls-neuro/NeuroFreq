function phase = hnf_phase_at_sample(TF,phase,freqIndices,timeIndices)
% Express linear transform coefficients as oscillator phase at sample time.
% Pass a phase array [C,F,T,R] (or [F,T,R] for one sensor). If it contains
% a frequency/time subset, supply the matching indices into the original
% TF axes. The original sensor/trial dimensions and storage shape are kept.
if nargin < 2, phase = TF.phase; end
if nargin < 3, freqIndices = 1:numel(TF.freqs); end
if nargin < 4, timeIndices = 1:numel(TF.times); end
if isfield(TF,'phaseReference')
    reference = TF.phaseReference;
elseif isfield(TF,'method') && any(strcmpi(TF.method,{'stransform','stft','RID-Rihaczek'}))
    error('NeuroFreq:PhaseProvenance', ...
        'This historical %s output lacks phase-reference information. Recompute the transform.',TF.method);
else
    return % Ordinary analytic phase in an older or user-created TF structure.
end
switch reference
    case 'instantaneous'
        return
    case 'rihaczek'
        error('NeuroFreq:BilinearPhase','Rihaczek distribution phase is not a linear oscillator phase.');
    case 'epoch-start-carrier'
        offset = TF.phaseReferenceSeconds;
        if numel(offset) ~= numel(TF.times)
            error('NeuroFreq:PhaseGrid','Stockwell phase offsets must match the time grid.');
        end
        offset = reshape(offset(timeIndices),1,[]);
    case 'frame-start-carrier'
        offset = TF.phaseReferenceSeconds;
        validateattributes(offset,{'numeric'},{'real','finite','scalar'});
    otherwise
        error('NeuroFreq:PhaseReference','Unknown phase reference "%s".',reference);
end
nf = numel(freqIndices);
nt = numel(timeIndices);
if isfield(TF,'ntrls') && isscalar(TF.ntrls)
    nObs = TF.ntrls;
elseif size(phase,2) == nf && size(phase,3) == nt
    nObs = size(phase,4);
else
    nObs = size(phase,3); % one sensor, sensor dimension omitted
end
nChan = numel(phase)/(nf*nt*nObs);
if ~isfinite(nChan) || nChan ~= fix(nChan) || nChan < 1
    error('NeuroFreq:PhaseDimensions','Phase size disagrees with the selected axes and trial count.');
end
originalSize = size(phase);
phase = reshape(phase,nChan,nf,nt,nObs);
cycles = reshape(TF.freqs(freqIndices),[],1)*offset;
advance = reshape(2*pi*mod(cycles,1),1,nf,[]);
phase = reshape(mod(bsxfun(@plus,phase,advance)+pi,2*pi)-pi,originalSize);
end
