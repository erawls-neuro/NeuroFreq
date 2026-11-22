function [X,compact] = hnf_tf_grid(X,nChan,nFreq,nTime,nObs,field)
% Recover NF's logical [sensor frequency time observation] axes. MATLAB
% drops trailing singleton dimensions; older single-sensor transforms also
% squeeze their output. Never infer permutations from equal element counts.
dims = [nChan nFreq nTime nObs];
stored = size(X);
stored(end+1:4) = 1;
compact = ~isequal(stored,dims);
if compact
    reduced = dims(dims~=1);
    reduced(end+1:2) = 1;
    sensorCompact = nChan==1 && isequal(stored,[nFreq nTime nObs 1]);
    if ~sensorCompact && ~isequal(size(X),reduced)
        error('NeuroFreq:TFGrid','%s has size %s; expected %s (or its singleton-squeezed layout).', ...
            field,mat2str(size(X)),mat2str(dims));
    end
end
X = reshape(X,dims);
end
