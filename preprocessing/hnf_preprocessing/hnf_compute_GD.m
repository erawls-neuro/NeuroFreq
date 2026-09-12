function res = hnf_compute_GD(topografie, canali, num_componenti)
% HNF_COMPUTE_GD  Computes the Generic Discontinuity spatial feature.
%
% Usage:
%   res = compute_GD_feat(topografie, canali, num_componenti);
%
% Inputs:
%   topografie     - Real IC topographies, components x channels.
%   canali         - Channel locations, in the SAME order as map columns.
%                    Each channel requires finite scalar X, Y, Z values.
%   num_componenti - Number of components to evaluate.
%
% Output:
%   res            - One GDSF value per component (row vector).
%
% Internal test modification, 2026-09-10:
%   Implement Mognon et al. (2011), p. 233, Eq. 8 using all other channels,
%   include every electrode in the maximum, and exclude self by index.
%   Preserve exp(-distance) in the supplied coordinate units and the original
%   mean(weight .* topography), i.e. divide by N-1, NOT sum(weight).
%   This changes the feature relative to the released ten-neighbor code.
%
nChan = numel(canali);
if nChan < 2
    error('compute_GD_feat:TooFewChannels', ...
        'GDSF requires at least two channel locations.');
end
if size(topografie, 2) ~= nChan
    error('compute_GD_feat:ChannelMismatch', ...
        ['Topography columns and channel locations must match. ' ...
        'Apply any channel subset to both, in the same order.']);
end
% Reject missing coordinates before concatenation can discard empty entries.
coordinates = [{canali.X}, {canali.Y}, {canali.Z}];
if any(~cellfun(@(v) isnumeric(v) && isreal(v) && isscalar(v) && ...
        isfinite(v), coordinates))
    error('compute_GD_feat:InvalidCoordinates', ...
        'Each channel must have finite, real scalar X, Y, Z coordinates.');
end
xpos = [canali.X]; ypos = [canali.Y]; zpos = [canali.Z];
pos = [xpos', ypos', zpos'];
res = zeros(1, num_componenti);
for ic = 1:num_componenti
    aux = zeros(1, nChan);
    for el = 1:nChan
        P = pos(el, :);
        d = pos - repmat(P, nChan, 1);
        dist = sqrt(sum(d .* d, 2));
        % All channels except this electrode; valid even with tied positions.
        repchas = [1:el-1, el+1:nChan];
        weightchas = exp(-dist(repchas));
        % Eq. 8 averages weighted map values over the N-1 other electrodes.
        aux(el) = abs(topografie(ic, el) - ...
            mean(weightchas .* topografie(ic, repchas)'));
    end
    res(ic) = max(aux);
end
end
