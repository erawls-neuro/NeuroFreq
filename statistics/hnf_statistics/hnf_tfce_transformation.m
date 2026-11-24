function enhanced = hnf_tfce_transformation(stat, ChN, E_H, flag_ft)
% EPT expects channel x sample or channel x frequency x time maps.
% On a regular frequency-time grid, adjacent frequency bins (including
% the bin itself) replace the electrode-neighbour rows. This counts each
% bin once and applies the requested [extent height] in either geometry.
if ~isreal(stat) || any(~isfinite(stat(:)))
    error('NeuroFreq:UndefinedTFCE', ...
        'TFCE requires finite real statistics. Inspect constant features, inputs and model fits.');
end
if ~any(stat(:)), enhanced = zeros(size(stat)); return; end
if ~isreal(E_H) || numel(E_H) ~= 2 || any(~isfinite(E_H)) || any(E_H < 0)
    error('NeuroFreq:InvalidTFCEExponents', 'Supply two finite nonnegative TFCE exponents.');
end
if flag_ft || isequal(ChN, 1)
    if ~ismatrix(stat)
        error('NeuroFreq:MissingChanlocs', 'A channel-frequency-time volume requires channel locations.');
    end
    row = (1:size(stat,1))';
    ChN = [row, max(1,row-1), min(row(end),row+1)];
elseif isvector(stat) && size(ChN,1) == numel(stat)
    stat = stat(:); % Channels only, one sample per electrode.
end
if size(ChN,1) ~= size(stat,1) || any(~isfinite(ChN(:))) || ...
        any(ChN(:) < 0 | ChN(:) > size(stat,1) | fix(ChN(:)) ~= ChN(:))
    error('NeuroFreq:NeighbourMismatch', 'The neighbour matrix must describe every channel in the statistic.');
end
% Native and interpreted integrators use the same signed cluster definition.
% A positive cluster must never change the seed later used for a negative
% cluster (the external EPT C code reuses that index while scanning edges).
if ndims(stat) > 3
    error('NeuroFreq:InvalidTFCEMap', 'Expected a channel, frequency or time map with at most three axes.');
end
enhanced = hnf_tfce_integrate(full(double(stat)),full(double(ChN)),full(double(E_H)));
end
