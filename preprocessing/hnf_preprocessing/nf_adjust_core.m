function [rejected, classification] = nf_adjust_core(EEG, profile, varargin)
% NF_ADJUST_CORE  Shared native ADJUST-family component classifier.
%
% This private core implements the spatial and temporal intentions of
% ADJUST and adjusted-ADJUST without their square-ICA, fixed-montage,
% numerical, plotting, or report-file assumptions. Public calls are made
% through nfADJUST for the adult profile. NeuroFreq's pipeline uses
% hnf_ADJUST for its separate developmental nfadjustedadjust route.
%
% SPDX-License-Identifier: GPL-2.0-or-later

profile = nf_profile_local(profile);
options = nf_options_local(profile, varargin{:});
decomposition = nf_decomposition_local(EEG);
geometry = nf_geometry_local( ...
    EEG, decomposition.icaChannels, options.referenceChanlocs);

maps = double(EEG.icawinv.');
mapScale = sqrt(sum(maps .* maps, 2));
if any(~isfinite(mapScale)) || any(mapScale <= 0)
    error('nf_adjust_core:InvalidComponentMaps', ...
        'Every ICA component must have a finite nonzero scalp map.');
end

mapsL2 = maps ./ mapScale;
mapMean = mean(maps, 2);
mapStandardDeviation = std(maps, 0, 2);
mapZScoreValid = isfinite(mapStandardDeviation) & ...
    mapStandardDeviation > 0;
mapsZ = nan(size(maps));
if any(mapZScoreValid)
    mapsZ(mapZScoreValid, :) = ...
        (maps(mapZScoreValid, :) - mapMean(mapZScoreValid)) ./ ...
        mapStandardDeviation(mapZScoreValid);
end
activations = decomposition.activations .* mapScale;
needKurtosis = strcmp(profile, 'adult');
temporal = nf_temporal_features_local( ...
    activations, EEG, options.epochLength, needKurtosis);
features = nf_features_local(mapsL2, mapsZ, geometry, temporal);
features.mapZScoreValid = mapZScoreValid;

if strcmp(profile, 'adult')
    [rawByClass, evaluatedByClass, statusByClass, thresholds] = ...
        nf_adult_decisions_local(features, geometry);
else
    [rawByClass, evaluatedByClass, statusByClass, thresholds] = ...
        nf_developmental_decisions_local(features, geometry);
end

rawCandidateMask = rawByClass.horizontalEyeMovement | ...
    rawByClass.verticalEyeMovement | ...
    rawByClass.blink | ...
    rawByClass.discontinuity;

spread = nf_spread_rescue_local( ...
    mapsZ, rawByClass.horizontalEyeMovement, rawByClass.blink, ...
    geometry, options.spatialRescue);

postSpatialByClass = rawByClass;
postSpatialByClass.horizontalEyeMovement = ...
    rawByClass.horizontalEyeMovement & ~spread.horizontalMask;
postSpatialByClass.blink = ...
    rawByClass.blink & ~spread.blinkMask;
postSpatialCandidateMask = ...
    postSpatialByClass.horizontalEyeMovement | ...
    postSpatialByClass.verticalEyeMovement | ...
    postSpatialByClass.blink | ...
    postSpatialByClass.discontinuity;
periodic = nf_periodic_rescue_local( ...
    activations, EEG, postSpatialCandidateMask, options);

finalByClass = postSpatialByClass;
classNames = fieldnames(finalByClass);
for classIndex = 1:numel(classNames)
    className = classNames{classIndex};
    finalByClass.(className) = ...
        postSpatialByClass.(className) & ...
        ~periodic.rescueMask;
end
rejectedMask = finalByClass.horizontalEyeMovement | ...
    finalByClass.verticalEyeMovement | ...
    finalByClass.blink | finalByClass.discontinuity;
rescueMask = rawCandidateMask & ~rejectedMask;
rejected = find(rejectedMask).';

if options.emitWarnings == true
    nf_warn_unavailable_local(evaluatedByClass, statusByClass);
    nf_warn_spread_unavailable_local(spread);
    nf_warn_periodic_unavailable_local(periodic);
end

classification = nf_result_local( ...
    profile, options, decomposition, features, thresholds, geometry, ...
    rawByClass, evaluatedByClass, statusByClass, rawCandidateMask, ...
    spread, postSpatialByClass, postSpatialCandidateMask, periodic, ...
    finalByClass, rescueMask, rejectedMask, rejected);
end

function profile = nf_profile_local(profile)

if ~(ischar(profile) || ...
        (isstring(profile) && isscalar(profile)))
    error('nf_adjust_core:InvalidProfile', ...
        'profile must be ''adult'' or ''developmental''.');
end

profile = lower(strtrim(char(profile)));
if ~ismember(profile, {'adult', 'developmental'})
    error('nf_adjust_core:InvalidProfile', ...
        'profile must be ''adult'' or ''developmental''.');
end
end

function options = nf_options_local(profile, varargin)

if strcmp(profile, 'adult')
    defaultEpochLength = 5;
    defaultPeriodicRescue = false;
    defaultSpatialRescue = false;
else
    defaultEpochLength = 1;
    defaultPeriodicRescue = true;
    defaultSpatialRescue = true;
end

parser = inputParser;
parser.FunctionName = 'nf_adjust_core';
addParameter(parser, 'epochLength', defaultEpochLength, ...
    @nf_positive_scalar_local);
addParameter(parser, 'periodicRescue', defaultPeriodicRescue, ...
    @nf_logical_scalar_local);
addParameter(parser, 'spatialRescue', defaultSpatialRescue, ...
    @nf_logical_scalar_local);
addParameter(parser, 'rescueBand', [5 15], ...
    @nf_increasing_pair_local);
% Forty hertz avoids fitting the MADE preprocessing roll-off near 45 Hz.
addParameter(parser, 'spectrumRange', [1 40], ...
    @nf_increasing_pair_local);
addParameter(parser, 'narrowPeakWidth', [0.5 12], ...
    @nf_increasing_pair_local);
addParameter(parser, 'plateauPeakWidth', [8 24], ...
    @nf_increasing_pair_local);
addParameter(parser, 'referenceChanlocs', [], ...
    @nf_chanlocs_or_empty_local);
addParameter(parser, 'emitWarnings', true, ...
    @nf_logical_scalar_local);
varargin = nf_option_arguments_local(varargin);
parse(parser, varargin{:});

options = parser.Results;
options.epochLength = double(options.epochLength);
options.periodicRescue = logical(options.periodicRescue);
options.spatialRescue = logical(options.spatialRescue);
options.emitWarnings = logical(options.emitWarnings);
options.rescueBand = reshape(double(options.rescueBand), 1, 2);
options.spectrumRange = reshape(double(options.spectrumRange), 1, 2);
options.narrowPeakWidth = ...
    reshape(double(options.narrowPeakWidth), 1, 2);
options.plateauPeakWidth = ...
    reshape(double(options.plateauPeakWidth), 1, 2);

if options.rescueBand(1) < options.spectrumRange(1) || ...
        options.rescueBand(2) > options.spectrumRange(2)
    error('nf_adjust_core:InvalidPeriodicRange', ...
        'rescueBand must be contained within spectrumRange.');
end
end

function args = nf_option_arguments_local(args)

if numel(args) ~= 1 || ~isstruct(args{1})
    return
end

optionStructure = args{1};
if ~isscalar(optionStructure)
    error('nf_adjust_core:InvalidOptions', ...
        'An options structure must be scalar.');
end

names = fieldnames(optionStructure);
args = cell(1, 2 .* numel(names));
for index = 1:numel(names)
    args{2 .* index - 1} = names{index};
    args{2 .* index} = optionStructure.(names{index});
end
end

function valid = nf_positive_scalar_local(value)

valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
    isfinite(value) && value > 0;
end

function valid = nf_logical_scalar_local(value)

valid = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
    isfinite(value) && ismember(value, [0 1]);
end

function valid = nf_increasing_pair_local(value)

valid = isnumeric(value) && isreal(value) && numel(value) == 2 && ...
    all(isfinite(value)) && value(1) >= 0 && value(2) > value(1);
end

function valid = nf_chanlocs_or_empty_local(value)

valid = isempty(value) || isstruct(value);
end

function decomposition = nf_decomposition_local(EEG)

required = {'data', 'srate', 'chanlocs', 'icaweights', 'icasphere', ...
    'icawinv', 'icachansind'};
if ~isstruct(EEG) || numel(EEG) ~= 1
    error('nf_adjust_core:InvalidEEG', ...
        'EEG must be one EEGLAB dataset structure.');
end

for index = 1:numel(required)
    if ~isfield(EEG, required{index})
        error('nf_adjust_core:InvalidEEG', ...
            'EEG is missing required field %s.', required{index});
    end
end

if ~isnumeric(EEG.data) || ~isreal(EEG.data) || isempty(EEG.data) || ...
        any(~isfinite(EEG.data(:)))
    error('nf_adjust_core:InvalidEEG', ...
        'EEG.data must contain finite real numeric data.');
end

if ~nf_positive_scalar_local(EEG.srate)
    error('nf_adjust_core:InvalidSamplingRate', ...
        'EEG.srate must be a positive finite scalar.');
end

if ~isstruct(EEG.chanlocs) || ...
        numel(EEG.chanlocs) ~= size(EEG.data, 1)
    error('nf_adjust_core:InvalidChannelLocations', ...
        'EEG.chanlocs must contain one entry per data channel.');
end

if ~isnumeric(EEG.icachansind) || ...
        ~isreal(EEG.icachansind) || ...
        ~isvector(EEG.icachansind)
    error('nf_adjust_core:InvalidICAChannels', ...
        'EEG.icachansind must be a real numeric index vector.');
end

icaChannels = double(EEG.icachansind(:).');
if isempty(icaChannels) || any(~isfinite(icaChannels)) || ...
        any(icaChannels < 1) || any(icaChannels ~= round(icaChannels)) || ...
        any(icaChannels > size(EEG.data, 1)) || ...
        numel(unique(icaChannels)) ~= numel(icaChannels)
    error('nf_adjust_core:InvalidICAChannels', ...
        'EEG.icachansind must contain unique valid channel indices.');
end

if ~isnumeric(EEG.icaweights) || ~isreal(EEG.icaweights) || ...
        ~isnumeric(EEG.icasphere) || ~isreal(EEG.icasphere) || ...
        ~isnumeric(EEG.icawinv) || ~isreal(EEG.icawinv) || ...
        ~ismatrix(EEG.icaweights) || ...
        ~ismatrix(EEG.icasphere) || ~ismatrix(EEG.icawinv)
    error('nf_adjust_core:InvalidICA', ...
        'ICA weights, sphere, and inverse weights must be real numeric matrices.');
end

weights = double(EEG.icaweights);
sphere = double(EEG.icasphere);
inverseWeights = double(EEG.icawinv);
if isempty(weights) || isempty(sphere) || isempty(inverseWeights) || ...
        any(~isfinite(weights(:))) || any(~isfinite(sphere(:))) || ...
        any(~isfinite(inverseWeights(:)))
    error('nf_adjust_core:InvalidICA', ...
        'ICA weights, sphere, and inverse weights must be finite.');
end

nComponents = size(weights, 1);
nMapChannels = numel(icaChannels);
if size(weights, 2) ~= size(sphere, 1) || ...
        size(sphere, 2) ~= nMapChannels || ...
        size(inverseWeights, 1) ~= nMapChannels || ...
        size(inverseWeights, 2) ~= nComponents || ...
        nComponents > size(sphere, 1)
    error('nf_adjust_core:InvalidICADimensions', ...
        ['ICA matrices are dimensionally inconsistent with ' ...
        'EEG.icachansind.']);
end

data = reshape( ...
    double(EEG.data(icaChannels, :, :)), nMapChannels, []);
activations = (weights * sphere) * data;
if size(activations, 1) ~= nComponents || ...
        size(activations, 2) ~= size(data, 2) || ...
        any(~isfinite(activations(:)))
    error('nf_adjust_core:InvalidActivations', ...
        'ICA activations could not be reconstructed consistently.');
end

decomposition = struct();
decomposition.icaChannels = icaChannels;
decomposition.nMapChannels = nMapChannels;
decomposition.nComponents = nComponents;
decomposition.activations = activations;
decomposition.sampleCount = size(activations, 2);
end

function geometry = nf_geometry_local(EEG, icaChannels, referenceChanlocs)

mapChanlocs = EEG.chanlocs(icaChannels);
if isempty(referenceChanlocs)
    referenceChanlocs = EEG.chanlocs;
end

atlas = nf_adjust_role_atlas();
[referenceRoleUnit, referenceSpatialUnit, referenceValid] = ...
    nf_location_units_local(referenceChanlocs);
[mapRoleUnit, mapSpatialUnit, mapValid] = ...
    nf_location_units_local(mapChanlocs);

referenceLabels = nf_labels_local(referenceChanlocs);
mapLabels = nf_labels_local(mapChanlocs);
nf_validate_unique_labels_local(referenceLabels, 'reference');
nf_validate_unique_labels_local(mapLabels, 'ICA-map');

nReference = numel(referenceChanlocs);
nMapChannels = numel(mapChanlocs);
referenceFlags = zeros(nReference, 1, 'uint32');
referenceBestAngle = nan(nReference, 1);
for index = 1:nReference
    if referenceValid(index) == true
        [referenceFlags(index), referenceBestAngle(index)] = ...
            nf_nearest_atlas_local(referenceRoleUnit(index, :), atlas);
    end
end

mapFlags = zeros(nMapChannels, 1, 'uint32');
mapBestAngle = nan(nMapChannels, 1);
mapStatus = repmat({''}, nMapChannels, 1);
for index = 1:nMapChannels
    label = mapLabels{index};
    match = [];
    if ~isempty(label)
        match = find(strcmpi(referenceLabels,label));
    end

    if numel(match) == 1 && referenceValid(match) == true
        mapFlags(index) = referenceFlags(match);
        mapBestAngle(index) = referenceBestAngle(match);
        hadMapCoordinates = mapValid(index);
        mapSpatialUnit(index, :) = referenceSpatialUnit(match, :);
        mapValid(index) = true;
        if hadMapCoordinates == false
            mapStatus{index} = ...
                'referenceLabelInheritedCoordinates';
        else
            mapStatus{index} = 'referenceLabel';
        end
    elseif mapValid(index) == true
        [mapFlags(index), mapBestAngle(index)] = ...
            nf_nearest_atlas_local(mapRoleUnit(index, :), atlas);
        mapStatus{index} = 'directAtlas';
    else
        mapStatus{index} = 'missingCoordinates';
    end
end

adjust = struct();
adjust.front = bitget(mapFlags, 1).' ~= 0;
adjust.posterior = bitget(mapFlags, 2).' ~= 0;
adjust.leftEye = bitget(mapFlags, 3).' ~= 0;
adjust.rightEye = bitget(mapFlags, 4).' ~= 0;

beall = struct();
beall.horizontalLeftEye = bitget(mapFlags, 5).' ~= 0;
beall.horizontalRightEye = bitget(mapFlags, 6).' ~= 0;
beall.blinkLeftEye = bitget(mapFlags, 7).' ~= 0;
beall.blinkRightEye = bitget(mapFlags, 8).' ~= 0;
beall.blinkCenterEye = bitget(mapFlags, 9).' ~= 0;
beall.control = bitget(mapFlags, 10).' ~= 0;
beall.posteriorLeft = bitget(mapFlags, 11).' ~= 0;
beall.posteriorRight = bitget(mapFlags, 12).' ~= 0;
beall.blinkCenterFallback = bitget(mapFlags, 17).' ~= 0;
beall.blinkCenterFallbackUsed = false;
if ~any(beall.blinkCenterEye) && any(beall.blinkCenterFallback)
    beall.blinkCenterEye = beall.blinkCenterFallback;
    beall.blinkCenterFallbackUsed = true;
end

spread = struct();
spread.outerRing = bitget(mapFlags, 13).' ~= 0;
spread.horizontalEye = bitget(mapFlags, 14).' ~= 0;
spread.blinkEye = bitget(mapFlags, 15).' ~= 0;
spread.higherThreshold = bitget(mapFlags, 16).' ~= 0;

[gdsfWeights, spreadNearMatrix] = ...
    nf_spatial_graphs_local(mapSpatialUnit, mapValid);

branchStatus = struct();
branchStatus.adjustHorizontal = ...
    any(adjust.leftEye) && any(adjust.rightEye);
branchStatus.adjustVertical = ...
    branchStatus.adjustHorizontal && any(adjust.front) && ...
    any(adjust.posterior);
branchStatus.adjustBlink = branchStatus.adjustVertical;
branchStatus.beallHorizontal = ...
    any(beall.horizontalLeftEye) && ...
    any(beall.horizontalRightEye) && ...
    any(beall.control) && any(beall.posteriorLeft) && ...
    any(beall.posteriorRight);
branchStatus.beallBlink = ...
    any(beall.blinkLeftEye) && any(beall.blinkRightEye) && ...
    any(beall.blinkCenterEye) && any(beall.control) && ...
    any(beall.posteriorLeft) && any(beall.posteriorRight);
branchStatus.gdsf = sum(mapValid) >= 2;

geometry = struct();
geometry.source = atlas.source;
geometry.coordinateFrame = 'EEGLAB X-nose/Y-left/Z';
geometry.mapIndices = icaChannels(:).';
geometry.nMapChannels = nMapChannels;
geometry.roleFlags = mapFlags;
geometry.roleNames = atlas.roleNames;
geometry.referenceRoleFlags = referenceFlags;
geometry.bestAtlasAngle = mapBestAngle;
geometry.referenceBestAtlasAngle = referenceBestAngle;
geometry.bestAtlasAngleUnits = 'radians';
geometry.mappingStatus = mapStatus;
geometry.referenceLabels = referenceLabels;
geometry.mapLabels = mapLabels;
geometry.mapCoordinateValid = mapValid;
geometry.referenceCoordinateValid = referenceValid;
geometry.mapUnitCoordinates = mapSpatialUnit;
geometry.referenceUnitCoordinates = referenceSpatialUnit;
geometry.adjust = adjust;
geometry.beall = beall;
geometry.spread = spread;
geometry.gdsfWeights = gdsfWeights;
geometry.spreadNearMatrix = spreadNearMatrix;
geometry.branchStatus = branchStatus;
geometry.gdsfNeighborCount = min(10, max(sum(mapValid) - 1, 0));
geometry.spreadChordThreshold = 0.55;
end

function [roleUnit, spatialUnit, valid] = ...
    nf_location_units_local(chanlocs)

nChannels = numel(chanlocs);
xyz = nan(nChannels, 3);
roleUnit = nan(nChannels, 3);
valid = false(nChannels, 1);
for index = 1:nChannels
    [coordinate, direction, available] = ...
        nf_one_location_local(chanlocs(index));
    if available == true
        xyz(index, :) = coordinate;
        roleUnit(index, :) = direction;
        valid(index) = true;
    end
end

spatialUnit = nan(nChannels, 3);
if any(valid)
    center = nf_sphere_center_local(xyz(valid, :));
    centered = xyz(valid, :) - center;
    lengths = sqrt(sum(centered .* centered, 2));
    usable = isfinite(lengths) & lengths > 0;
    validIndices = find(valid);
    invalidIndices = validIndices(~usable);
    valid(invalidIndices) = false;
    validIndices = validIndices(usable);
    spatialUnit(validIndices, :) = ...
        centered(usable, :) ./ lengths(usable);
end
end

function [coordinate, direction, available] = ...
    nf_one_location_local(chanloc)

coordinate = [NaN NaN NaN];
direction = [NaN NaN NaN];
available = false;

hasXYZ = isfield(chanloc, 'X') && isfield(chanloc, 'Y') && ...
    isfield(chanloc, 'Z');
if hasXYZ
    candidate = [chanloc.X chanloc.Y chanloc.Z];
    hasXYZ = isnumeric(candidate) && isreal(candidate) && ...
        numel(candidate) == 3 && all(isfinite(candidate)) && ...
        norm(double(candidate)) > 0;
end

hasTopography = isfield(chanloc, 'theta') && ...
    isfield(chanloc, 'radius');
if hasTopography
    theta = chanloc.theta;
    radius = chanloc.radius;
    hasTopography = isnumeric(theta) && isreal(theta) && ...
        isscalar(theta) && isfinite(theta) && ...
        isnumeric(radius) && isreal(radius) && ...
        isscalar(radius) && isfinite(radius) && radius >= 0;
end

if hasXYZ
    direction = double(candidate) ./ norm(double(candidate));
elseif hasTopography
    azimuth = double(theta) .* pi ./ 180;
    colatitude = double(radius) .* pi;
    % EEGLAB's 2-D polar theta has the opposite sign of its Cartesian
    % Y axis: negative theta points left, while positive Cartesian Y
    % points left. Keep both representations in the canonical
    % X-nose/Y-left/Z frame.
    direction = [ ...
        sin(colatitude) .* cos(azimuth), ...
        -sin(colatitude) .* sin(azimuth), ...
        cos(colatitude)];
end

if hasXYZ
    coordinate = double(candidate);
else
    if hasTopography
        coordinate = direction;
    end
end

available = all(isfinite(coordinate)) && ...
    all(isfinite(direction)) && norm(direction) > 0;
end

function center = nf_sphere_center_local(xyz)

center = [0 0 0];
if size(xyz, 1) < 4
    return
end

design = [2 .* xyz ones(size(xyz, 1), 1)];
target = sum(xyz .* xyz, 2);
if rank(design) < 4
    return
end

solution = design \ target;
candidate = solution(1:3).';
if all(isfinite(candidate))
    center = candidate;
end
end

function labels = nf_labels_local(chanlocs)

labels = repmat({''}, numel(chanlocs), 1);
for index = 1:numel(chanlocs)
    if isfield(chanlocs(index), 'labels')
        label = chanlocs(index).labels;
        if ischar(label)
            labels{index} = strtrim(label);
        elseif isstring(label) && isscalar(label)
            labels{index} = strtrim(char(label));
        end
    end
end
end

function nf_validate_unique_labels_local(labels, labelType)

nonempty = ~cellfun(@isempty, labels);
lowerLabels = cellfun(@lower,labels(nonempty),'UniformOutput',false);
if numel(unique(lowerLabels)) ~= numel(lowerLabels)
    error('nf_adjust_core:DuplicateChannelLabels', ...
        '%s channel labels must be unique, ignoring case.', labelType);
end
end

function [flags, angle] = nf_nearest_atlas_local(direction, atlas)

similarity = direction * atlas.coordinates.';
similarity = min(max(similarity, -1), 1);
[bestSimilarity, index] = max(similarity);
flags = atlas.flags(index);
angle = acos(bestSimilarity);
end

function [weights, nearMatrix] = ...
    nf_spatial_graphs_local(unitCoordinates, valid)

nChannels = size(unitCoordinates, 1);
weights = zeros(nChannels, nChannels);
nearMatrix = false(nChannels, nChannels);
validIndices = find(valid);
nValid = numel(validIndices);
if nValid < 2
    return
end

coordinates = unitCoordinates(validIndices, :);
similarity = coordinates * coordinates.';
similarity = min(max(similarity, -1), 1);
angularDistance = acos(similarity);
chordDistance = sqrt(max(0, 2 - 2 .* similarity));
near = chordDistance < 0.55;
near(1:nValid + 1:end) = false;
nearMatrix(validIndices, validIndices) = near;

neighborCount = min(10, nValid - 1);
for center = 1:nValid
    distance = angularDistance(center, :);
    distance(center) = Inf;
    [sortedDistance, order] = sort(distance, 'ascend');
    neighbors = order(1:neighborCount);
    neighborDistance = sortedDistance(1:neighborCount);
    neighborWeight = exp(-neighborDistance);
    neighborWeight = neighborWeight ./ sum(neighborWeight);
    weights(validIndices(center), validIndices(neighbors)) = ...
        neighborWeight;
end
end

function temporal = nf_temporal_features_local( ...
    activations, EEG, epochLength, needKurtosis)

[windowStarts, windowSamples, buffer] = ...
    nf_temporal_windows_local(EEG, epochLength);
nComponents = size(activations, 1);
nWindows = numel(windowStarts);
maximumEpochVariance = nan(nComponents, 1);
temporalKurtosis = nan(nComponents, 1);
buffer.trimFraction = 0.01;
buffer.trimCount = 0;

if nWindows == 0 || windowSamples < 2
    temporal = struct();
    temporal.maximumEpochVariance = maximumEpochVariance;
    temporal.temporalKurtosis = temporalKurtosis;
    temporal.buffer = buffer;
    return
end

indices = (0:windowSamples - 1).' + windowStarts;
trimCount = 0;
if nWindows > 100
    trimCount = floor(0.01 .* nWindows);
end

for component = 1:nComponents
    componentActivation = activations(component, :);
    chunks = reshape( ...
        componentActivation(indices), windowSamples, nWindows);
    centered = chunks - mean(chunks, 1);
    epochVariance = sum(centered .* centered, 1) ./ ...
        (windowSamples - 1);
    sortedVariance = sort(epochVariance, 'ascend');
    retainedVariance = ...
        sortedVariance(1:nWindows - trimCount);
    varianceMean = mean(retainedVariance);
    if isfinite(varianceMean) && varianceMean > 0
        maximumEpochVariance(component) = ...
            retainedVariance(end) ./ varianceMean;
    end

    if needKurtosis == true
        secondMoment = mean(centered .* centered, 1);
        fourthMoment = mean((centered .* centered) .^ 2, 1);
        epochKurtosis = fourthMoment ./ (secondMoment .^ 2);
        finiteKurtosis = epochKurtosis(isfinite(epochKurtosis));
        if ~isempty(finiteKurtosis)
            finiteKurtosis = sort(finiteKurtosis, 'ascend');
            kurtosisTrim = 0;
            if numel(finiteKurtosis) > 100
                kurtosisTrim = floor(0.01 .* numel(finiteKurtosis));
            end
            retainedKurtosis = finiteKurtosis( ...
                1:numel(finiteKurtosis) - kurtosisTrim);
            temporalKurtosis(component) = mean(retainedKurtosis);
        end
    end
end

buffer.trimFraction = 0.01;
buffer.trimCount = trimCount;
temporal = struct();
temporal.maximumEpochVariance = maximumEpochVariance;
temporal.temporalKurtosis = temporalKurtosis;
temporal.buffer = buffer;
end

function [validStarts, windowSamples, buffer] = ...
    nf_temporal_windows_local(EEG, epochLength)

totalSamples = numel(EEG.data(1, :, :));
trials = size(EEG.data, 3);
pnts = size(EEG.data, 2);
if trials > 1
    windowSamples = pnts;
    starts = 1 + (0:(trials - 1)) .* pnts;
    tailSamples = 0;
    bufferMode = 'inputEpochs';
else
    windowSamples = round(epochLength .* double(EEG.srate));
    if windowSamples < 2
        error('nf_adjust_core:InvalidEpochLength', ...
            'epochLength must contain at least two samples.');
    end
    starts = 1:windowSamples:(totalSamples - windowSamples + 1);
    tailSamples = mod(totalSamples, windowSamples);
    bufferMode = 'fixedContinuousChunks';
end

cuts = nf_data_cuts_local(EEG);
[validStarts, crossingDropped] = ...
    nf_drop_crossing_windows_local(starts, windowSamples, cuts);

buffer = struct();
buffer.windowSamples = windowSamples;
buffer.windowSeconds = windowSamples ./ double(EEG.srate);
buffer.mode = bufferMode;
buffer.chunkCount = numel(validStarts);
buffer.crossingChunksDropped = crossingDropped;
buffer.tailSamplesDropped = tailSamples;
buffer.cutLatencies = cuts;
end

function cuts = nf_data_cuts_local(EEG)

pnts = size(EEG.data, 2);
trials = size(EEG.data, 3);
cuts = (1:(trials - 1)) .* pnts + 0.5;
if ~isfield(EEG, 'event') || isempty(EEG.event)
    cuts = unique(double(cuts));
    return
end

for index = 1:numel(EEG.event)
    if ~isfield(EEG.event(index), 'type') || ...
            ~isfield(EEG.event(index), 'latency')
        continue
    end
    eventType = EEG.event(index).type;
    isBoundary = false;
    if ischar(eventType)
        isBoundary = strcmpi(strtrim(eventType), 'boundary');
    elseif isstring(eventType) && isscalar(eventType)
        isBoundary = strcmpi(strtrim(char(eventType)), 'boundary');
    end
    latency = EEG.event(index).latency;
    if isBoundary == true && isnumeric(latency) && ...
            isreal(latency) && isscalar(latency) && ...
            isfinite(latency)
        cuts(end + 1) = double(latency);
    end
end

cuts = unique(double(cuts));
end

function [validStarts, dropped] = ...
    nf_drop_crossing_windows_local(starts, windowSamples, cuts)

valid = true(size(starts));
leftEdge = starts - 0.5;
rightEdge = starts + windowSamples - 0.5;
for index = 1:numel(cuts)
    crosses = cuts(index) > leftEdge & cuts(index) < rightEdge;
    valid = valid & ~crosses;
end

dropped = sum(~valid);
validStarts = starts(valid);
end

function features = nf_features_local(mapsL2, mapsZ, geometry, temporal)

features = struct();
features.leftEyeMean = ...
    nf_mask_mean_local(mapsL2, geometry.adjust.leftEye, false);
features.rightEyeMean = ...
    nf_mask_mean_local(mapsL2, geometry.adjust.rightEye, false);
features.frontMean = ...
    nf_mask_mean_local(mapsL2, geometry.adjust.front, false);
features.posteriorMean = ...
    nf_mask_mean_local(mapsL2, geometry.adjust.posterior, false);
features.spatialEyeDifference = ...
    abs(features.leftEyeMean - features.rightEyeMean);
features.spatialAverageDifference = ...
    abs(features.frontMean) - abs(features.posteriorMean);
features.frontVariance = ...
    nf_mask_variance_local(mapsL2, geometry.adjust.front);
features.posteriorVariance = ...
    nf_mask_variance_local(mapsL2, geometry.adjust.posterior);
features.spatialVarianceDifference = ...
    features.frontVariance - features.posteriorVariance;

if geometry.branchStatus.gdsf == true
    neighborMean = mapsL2 * geometry.gdsfWeights.';
    difference = abs(mapsL2 - neighborMean);
    features.discontinuitySpatial = ...
        max(difference(:, geometry.mapCoordinateValid), [], 2);
else
    features.discontinuitySpatial = nan(size(mapsL2, 1), 1);
end

features.maximumEpochVariance = temporal.maximumEpochVariance;
features.temporalKurtosis = temporal.temporalKurtosis;
features.temporalBuffer = temporal.buffer;
features.beall = nf_beall_features_local(mapsZ, geometry);
end

function values = nf_mask_mean_local(data, mask, useAbsolute)

if ~any(mask)
    values = nan(size(data, 1), 1);
    return
end

selected = data(:, mask);
if useAbsolute == true
    selected = abs(selected);
end
values = mean(selected, 2);
end

function values = nf_mask_variance_local(data, mask)

if ~any(mask)
    values = nan(size(data, 1), 1);
    return
end

values = var(data(:, mask), 0, 2);
end

function beall = nf_beall_features_local(mapsZ, geometry)

beall = struct();
beall.horizontalLeftEye = nf_mask_mean_local( ...
    mapsZ, geometry.beall.horizontalLeftEye, true);
beall.horizontalRightEye = nf_mask_mean_local( ...
    mapsZ, geometry.beall.horizontalRightEye, true);
beall.blinkLeftEye = nf_mask_mean_local( ...
    mapsZ, geometry.beall.blinkLeftEye, true);
beall.blinkRightEye = nf_mask_mean_local( ...
    mapsZ, geometry.beall.blinkRightEye, true);
beall.blinkCenterEye = nf_mask_mean_local( ...
    mapsZ, geometry.beall.blinkCenterEye, false);
beall.control = nf_mask_mean_local( ...
    mapsZ, geometry.beall.control, true);
beall.posteriorLeft = nf_mask_mean_local( ...
    mapsZ, geometry.beall.posteriorLeft, true);
beall.posteriorRight = nf_mask_mean_local( ...
    mapsZ, geometry.beall.posteriorRight, true);

horizontalFrontMask = geometry.beall.horizontalLeftEye | ...
    geometry.beall.horizontalRightEye;
blinkFrontMask = geometry.beall.blinkLeftEye | ...
    geometry.beall.blinkRightEye;
posteriorMask = geometry.beall.posteriorLeft | ...
    geometry.beall.posteriorRight;
beall.horizontalFrontVariance = ...
    nf_mask_variance_local(mapsZ, horizontalFrontMask);
beall.blinkFrontVariance = ...
    nf_mask_variance_local(mapsZ, blinkFrontMask);
beall.posteriorVariance = ...
    nf_mask_variance_local(mapsZ, posteriorMask);
beall.absolutePosteriorVariance = ...
    nf_mask_variance_local(abs(mapsZ), posteriorMask);
end

function [raw, evaluated, status, thresholds] = ...
    nf_adult_decisions_local(features, geometry)

nComponents = numel(features.maximumEpochVariance);
raw = nf_empty_class_masks_local(nComponents);
evaluated = nf_empty_class_status_local(false);
status = nf_empty_class_text_local('notEvaluated');

thresholds = struct();
thresholds.spatialEyeDifference = ...
    nf_em_threshold_local(features.spatialEyeDifference);
thresholds.spatialAverageDifference = ...
    nf_em_threshold_local(features.spatialAverageDifference);
thresholds.temporalKurtosis = ...
    nf_em_threshold_local(features.temporalKurtosis);
thresholds.discontinuitySpatial = ...
    nf_em_threshold_local(features.discontinuitySpatial);
thresholds.maximumEpochVariance = ...
    nf_em_threshold_local(features.maximumEpochVariance);

evaluated.horizontalEyeMovement = ...
    geometry.branchStatus.adjustHorizontal && ...
    thresholds.spatialEyeDifference.valid && ...
    thresholds.maximumEpochVariance.valid;
if evaluated.horizontalEyeMovement == true
    raw.horizontalEyeMovement = ...
        features.spatialEyeDifference >= ...
        thresholds.spatialEyeDifference.value & ...
        features.leftEyeMean .* features.rightEyeMean < 0 & ...
        features.maximumEpochVariance >= ...
        thresholds.maximumEpochVariance.value;
    status.horizontalEyeMovement = 'evaluated';
else
    status.horizontalEyeMovement = nf_class_reason_local( ...
        geometry.branchStatus.adjustHorizontal, ...
        [thresholds.spatialEyeDifference.valid ...
        thresholds.maximumEpochVariance.valid]);
end

evaluated.verticalEyeMovement = ...
    geometry.branchStatus.adjustVertical && ...
    thresholds.spatialAverageDifference.valid && ...
    thresholds.maximumEpochVariance.valid;
if evaluated.verticalEyeMovement == true
    raw.verticalEyeMovement = ...
        features.spatialAverageDifference >= ...
        thresholds.spatialAverageDifference.value & ...
        features.leftEyeMean .* features.rightEyeMean > 0 & ...
        features.spatialVarianceDifference > 0 & ...
        features.maximumEpochVariance >= ...
        thresholds.maximumEpochVariance.value;
    status.verticalEyeMovement = 'evaluated';
else
    status.verticalEyeMovement = nf_class_reason_local( ...
        geometry.branchStatus.adjustVertical, ...
        [thresholds.spatialAverageDifference.valid ...
        thresholds.maximumEpochVariance.valid]);
end

evaluated.blink = geometry.branchStatus.adjustBlink && ...
    thresholds.spatialAverageDifference.valid && ...
    thresholds.temporalKurtosis.valid;
if evaluated.blink == true
    raw.blink = ...
        features.spatialAverageDifference >= ...
        thresholds.spatialAverageDifference.value & ...
        features.leftEyeMean .* features.rightEyeMean > 0 & ...
        features.temporalKurtosis >= ...
        thresholds.temporalKurtosis.value & ...
        features.spatialVarianceDifference > 0;
    status.blink = 'evaluated';
else
    status.blink = nf_class_reason_local( ...
        geometry.branchStatus.adjustBlink, ...
        [thresholds.spatialAverageDifference.valid ...
        thresholds.temporalKurtosis.valid]);
end

evaluated.discontinuity = geometry.branchStatus.gdsf && ...
    thresholds.discontinuitySpatial.valid && ...
    thresholds.maximumEpochVariance.valid;
if evaluated.discontinuity == true
    raw.discontinuity = ...
        features.discontinuitySpatial >= ...
        thresholds.discontinuitySpatial.value & ...
        features.maximumEpochVariance >= ...
        thresholds.maximumEpochVariance.value;
    status.discontinuity = 'evaluated';
else
    status.discontinuity = nf_class_reason_local( ...
        geometry.branchStatus.gdsf, ...
        [thresholds.discontinuitySpatial.valid ...
        thresholds.maximumEpochVariance.valid]);
end
end

function [raw, evaluated, status, thresholds] = ...
    nf_developmental_decisions_local(features, geometry)

nComponents = numel(features.maximumEpochVariance);
raw = nf_empty_class_masks_local(nComponents);
evaluated = nf_empty_class_status_local(false);
status = nf_empty_class_text_local('notEvaluated');

thresholds = struct();
thresholds.spatialEyeDifference = ...
    nf_unused_threshold_local('unusedByDevelopmentalProfile');
thresholds.spatialAverageDifference = ...
    nf_em_threshold_local(features.spatialAverageDifference);
thresholds.temporalKurtosis = ...
    nf_unused_threshold_local('unusedByDevelopmentalProfile');
thresholds.discontinuitySpatial = ...
    nf_em_threshold_local(features.discontinuitySpatial);
thresholds.maximumEpochVariance = ...
    nf_em_threshold_local(features.maximumEpochVariance);

evaluated.horizontalEyeMovement = ...
    geometry.branchStatus.beallHorizontal;
if evaluated.horizontalEyeMovement == true
    controlMean = mean([ ...
        abs(features.beall.posteriorLeft), ...
        abs(features.beall.posteriorRight), ...
        abs(features.beall.control)], 2);
    raw.horizontalEyeMovement = ...
        (features.beall.horizontalLeftEye > 2 | ...
        features.beall.horizontalRightEye > 2) & ...
        controlMean < 1 & ...
        (features.beall.posteriorVariance < 0.15 | ...
        features.beall.absolutePosteriorVariance < 0.075) & ...
        features.beall.horizontalFrontVariance > 2.5;
    status.horizontalEyeMovement = 'evaluated';
else
    status.horizontalEyeMovement = 'missingRequiredGeometry';
end

evaluated.blink = geometry.branchStatus.beallBlink;
if evaluated.blink == true
    controlMean = mean([ ...
        abs(features.beall.posteriorLeft), ...
        abs(features.beall.posteriorRight), ...
        abs(features.beall.control)], 2);
    raw.blink = ...
        (features.beall.blinkLeftEye > 2 | ...
        abs(features.beall.blinkCenterEye) > 2.5 | ...
        features.beall.blinkRightEye > 2) & ...
        controlMean < 1 & ...
        (features.beall.posteriorVariance < 0.15 | ...
        features.beall.absolutePosteriorVariance < 0.075) & ...
        features.beall.blinkFrontVariance > ...
        features.beall.posteriorVariance;
    status.blink = 'evaluated';
else
    status.blink = 'missingRequiredGeometry';
end

evaluated.verticalEyeMovement = ...
    geometry.branchStatus.adjustVertical && ...
    thresholds.spatialAverageDifference.valid && ...
    thresholds.maximumEpochVariance.valid;
if evaluated.verticalEyeMovement == true
    raw.verticalEyeMovement = ...
        features.spatialAverageDifference >= ...
        thresholds.spatialAverageDifference.value & ...
        features.leftEyeMean .* features.rightEyeMean > 0 & ...
        features.spatialVarianceDifference > 0 & ...
        features.maximumEpochVariance >= ...
        thresholds.maximumEpochVariance.value;
    status.verticalEyeMovement = 'evaluated';
else
    status.verticalEyeMovement = nf_class_reason_local( ...
        geometry.branchStatus.adjustVertical, ...
        [thresholds.spatialAverageDifference.valid ...
        thresholds.maximumEpochVariance.valid]);
end

evaluated.discontinuity = geometry.branchStatus.gdsf && ...
    thresholds.discontinuitySpatial.valid && ...
    thresholds.maximumEpochVariance.valid;
if evaluated.discontinuity == true
    raw.discontinuity = ...
        features.discontinuitySpatial >= ...
        thresholds.discontinuitySpatial.value & ...
        features.maximumEpochVariance >= ...
        thresholds.maximumEpochVariance.value;
    status.discontinuity = 'evaluated';
else
    status.discontinuity = nf_class_reason_local( ...
        geometry.branchStatus.gdsf, ...
        [thresholds.discontinuitySpatial.valid ...
        thresholds.maximumEpochVariance.valid]);
end
end

function masks = nf_empty_class_masks_local(nComponents)

masks = struct();
masks.horizontalEyeMovement = false(nComponents, 1);
masks.verticalEyeMovement = false(nComponents, 1);
masks.blink = false(nComponents, 1);
masks.discontinuity = false(nComponents, 1);
end

function status = nf_empty_class_status_local(value)

status = struct();
status.horizontalEyeMovement = value;
status.verticalEyeMovement = value;
status.blink = value;
status.discontinuity = value;
end

function status = nf_empty_class_text_local(value)

status = struct();
status.horizontalEyeMovement = value;
status.verticalEyeMovement = value;
status.blink = value;
status.discontinuity = value;
end

function reason = nf_class_reason_local(geometryAvailable, thresholdValid)

if geometryAvailable == false
    reason = 'missingRequiredGeometry';
elseif ~all(thresholdValid)
    reason = 'thresholdNotIdentifiable';
else
    reason = 'notEvaluated';
end
end

function threshold = nf_unused_threshold_local(status)

threshold = struct();
threshold.value = Inf;
threshold.valid = false;
threshold.status = status;
threshold.converged = false;
threshold.iterations = 0;
threshold.sampleCount = 0;
threshold.lowMean = NaN;
threshold.highMean = NaN;
threshold.lowVariance = NaN;
threshold.highVariance = NaN;
threshold.lowPrior = NaN;
threshold.highPrior = NaN;
end

function threshold = nf_em_threshold_local(values)

threshold = nf_unused_threshold_local('notIdentifiable');
x = double(values(isfinite(values))).';
threshold.sampleCount = numel(x);
if numel(x) < 4
    return
end

minimum = min(x);
maximum = max(x);
rangeValue = maximum - minimum;
if ~isfinite(rangeValue) || rangeValue <= ...
        eps(max(1, max(abs(x))))
    return
end

z = (x - minimum) ./ rangeValue;
lowTraining = z(z < 0.495);
highTraining = z(z > 0.505);
if isempty(lowTraining) || isempty(highTraining)
    return
end

varianceFloor = 1e-12;
lowMean = mean(lowTraining);
highMean = mean(highTraining);
lowVariance = nf_initial_variance_local( ...
    lowTraining, varianceFloor);
highVariance = nf_initial_variance_local( ...
    highTraining, varianceFloor);
lowPrior = numel(lowTraining) ./ ...
    (numel(lowTraining) + numel(highTraining));
highPrior = 1 - lowPrior;

tolerance = 1e-6;
maximumIterations = 500;
converged = false;
iterations = 0;
for iteration = 1:maximumIterations
    iterations = iteration;
    old = [lowMean highMean lowVariance highVariance ...
        lowPrior highPrior];

    logLow = log(lowPrior) - 0.5 .* ...
        (log(2 .* pi .* lowVariance) + ...
        ((z - lowMean) .^ 2) ./ lowVariance);
    logHigh = log(highPrior) - 0.5 .* ...
        (log(2 .* pi .* highVariance) + ...
        ((z - highMean) .^ 2) ./ highVariance);
    stabilizer = max(logLow, logHigh);
    lowLikelihood = exp(logLow - stabilizer);
    highLikelihood = exp(logHigh - stabilizer);
    denominator = lowLikelihood + highLikelihood;
    lowResponsibility = lowLikelihood ./ denominator;
    highResponsibility = highLikelihood ./ denominator;

    lowCount = sum(lowResponsibility);
    highCount = sum(highResponsibility);
    if lowCount <= eps || highCount <= eps || ...
            any(~isfinite([lowCount highCount]))
        return
    end

    lowMean = sum(lowResponsibility .* z) ./ lowCount;
    highMean = sum(highResponsibility .* z) ./ highCount;
    lowVariance = sum(lowResponsibility .* ...
        ((z - lowMean) .^ 2)) ./ lowCount;
    highVariance = sum(highResponsibility .* ...
        ((z - highMean) .^ 2)) ./ highCount;
    lowVariance = max(lowVariance, varianceFloor);
    highVariance = max(highVariance, varianceFloor);
    lowPrior = max(lowCount ./ numel(z), eps);
    highPrior = max(highCount ./ numel(z), eps);
    priorSum = lowPrior + highPrior;
    lowPrior = lowPrior ./ priorSum;
    highPrior = highPrior ./ priorSum;

    current = [lowMean highMean lowVariance highVariance ...
        lowPrior highPrior];
    if any(~isfinite(current))
        return
    end
    if max(abs(current - old)) <= tolerance
        converged = true;
        break
    end
end

if highMean < lowMean
    temporaryMean = lowMean;
    lowMean = highMean;
    highMean = temporaryMean;
    temporaryVariance = lowVariance;
    lowVariance = highVariance;
    highVariance = temporaryVariance;
    temporaryPrior = lowPrior;
    lowPrior = highPrior;
    highPrior = temporaryPrior;
end

if highMean - lowMean <= sqrt(varianceFloor)
    return
end

threshold.converged = converged;
threshold.iterations = iterations;
threshold.lowMean = minimum + rangeValue .* lowMean;
threshold.highMean = minimum + rangeValue .* highMean;
threshold.lowVariance = rangeValue .^ 2 .* lowVariance;
threshold.highVariance = rangeValue .^ 2 .* highVariance;
threshold.lowPrior = lowPrior;
threshold.highPrior = highPrior;
[boundary, boundaryStatus] = nf_mixture_boundary_local( ...
    lowMean, highMean, lowVariance, highVariance, ...
    lowPrior, highPrior);
if ~isfinite(boundary)
    threshold.status = boundaryStatus;
    return
end

if converged == false
    threshold.status = 'maxIterations';
    return
end
threshold.value = minimum + rangeValue .* boundary;
threshold.valid = true;
threshold.status = boundaryStatus;
end

function value = nf_initial_variance_local(samples, floorValue)

if numel(samples) > 1
    value = var(samples, 0);
else
    value = floorValue;
end
value = max(value, floorValue);
end

function [boundary, status] = nf_mixture_boundary_local( ...
    lowMean, highMean, lowVariance, highVariance, lowPrior, highPrior)

coefficientA = 1 ./ (2 .* lowVariance) - ...
    1 ./ (2 .* highVariance);
coefficientB = highMean ./ highVariance - ...
    lowMean ./ lowVariance;
coefficientC = lowMean .^ 2 ./ (2 .* lowVariance) - ...
    highMean .^ 2 ./ (2 .* highVariance) + ...
    log(highPrior ./ sqrt(highVariance)) - ...
    log(lowPrior ./ sqrt(lowVariance));

boundary = NaN;
status = 'noPosteriorIntersection';
coefficientScale = max(1, max(abs([coefficientA ...
    coefficientB coefficientC])));
if abs(coefficientA) <= eps(coefficientScale) .* 16
    if abs(coefficientB) > eps(coefficientScale) .* 16
        boundary = -coefficientC ./ coefficientB;
        if boundary >= lowMean && boundary <= highMean
            status = 'posteriorIntersection';
        else
            boundary = NaN;
        end
    end
else
    discriminant = coefficientB .^ 2 - ...
        4 .* coefficientA .* coefficientC;
    if discriminant >= 0
        squareRoot = sqrt(discriminant);
        if coefficientB >= 0
            q = -0.5 .* (coefficientB + squareRoot);
        else
            q = -0.5 .* (coefficientB - squareRoot);
        end
        if abs(q) > eps(coefficientScale) .* 16
            roots = [q ./ coefficientA coefficientC ./ q];
        else
            roots = [ ...
                (-coefficientB + squareRoot) ./ ...
                (2 .* coefficientA), ...
                (-coefficientB - squareRoot) ./ ...
                (2 .* coefficientA)];
        end
        interior = roots(isfinite(roots) & ...
            roots >= lowMean & roots <= highMean);
        if ~isempty(interior)
            [~, index] = min(abs(interior - ...
                (lowMean + highMean) ./ 2));
            boundary = interior(index);
            status = 'posteriorIntersection';
        end
    end
end
end

function spread = nf_spread_rescue_local( ...
    mapsZ, horizontalMask, blinkMask, geometry, requested)

nComponents = size(mapsZ, 1);
horizontalRescue = false(nComponents, 1);
blinkRescue = false(nComponents, 1);
spread = struct();
spread.requested = logical(requested);
spread.horizontalCandidateMask = logical(horizontalMask);
spread.blinkCandidateMask = logical(blinkMask);
spread.candidateMask = horizontalMask | blinkMask;
spread.horizontalGeometryAvailable = false;
spread.blinkGeometryAvailable = false;
spread.horizontalStatus = 'notRequested';
spread.blinkStatus = 'notRequested';
if requested == false
    spread.status = 'notRequested';
    spread.horizontalMask = horizontalRescue;
    spread.blinkMask = blinkRescue;
    spread.rescueMask = horizontalRescue | blinkRescue;
    return
end

horizontalComponents = find(horizontalMask);
commonGeometry = sum(geometry.mapCoordinateValid) >= 3;
spread.horizontalGeometryAvailable = commonGeometry;
spread.blinkGeometryAvailable = commonGeometry;
if isempty(horizontalComponents)
    spread.horizontalStatus = 'noCandidates';
elseif spread.horizontalGeometryAvailable == false
    spread.horizontalStatus = 'notEvaluatedMissingGeometry';
else
    spread.horizontalStatus = 'evaluated';
    for index = 1:numel(horizontalComponents)
        component = horizontalComponents(index);
        horizontalRescue(component) = ...
            nf_one_spread_rescue_local( ...
            mapsZ(component, :), geometry.spread.outerRing, ...
            geometry.spread.horizontalEye, ...
            geometry.spread.higherThreshold, ...
            geometry.spreadNearMatrix, ...
            geometry.mapCoordinateValid);
    end
end

blinkComponents = find(blinkMask);
if isempty(blinkComponents)
    spread.blinkStatus = 'noCandidates';
elseif spread.blinkGeometryAvailable == false
    spread.blinkStatus = 'notEvaluatedMissingGeometry';
else
    spread.blinkStatus = 'evaluated';
    for index = 1:numel(blinkComponents)
        component = blinkComponents(index);
        blinkRescue(component) = nf_one_spread_rescue_local( ...
            mapsZ(component, :), geometry.spread.outerRing, ...
            geometry.spread.blinkEye, ...
            geometry.spread.higherThreshold, ...
            geometry.spreadNearMatrix, ...
            geometry.mapCoordinateValid);
    end
end

if isempty(horizontalComponents) && isempty(blinkComponents)
    spread.status = 'noCandidates';
elseif strcmp(spread.horizontalStatus, 'evaluated') && ...
        (strcmp(spread.blinkStatus, 'evaluated') || ...
        strcmp(spread.blinkStatus, 'noCandidates')) || ...
        strcmp(spread.blinkStatus, 'evaluated') && ...
        strcmp(spread.horizontalStatus, 'noCandidates')
    spread.status = 'evaluated';
elseif strcmp(spread.horizontalStatus, 'evaluated') || ...
        strcmp(spread.blinkStatus, 'evaluated')
    spread.status = 'partiallyEvaluated';
else
    spread.status = 'notEvaluatedMissingGeometry';
end
spread.horizontalMask = horizontalRescue;
spread.blinkMask = blinkRescue;
spread.rescueMask = horizontalRescue | blinkRescue;
end

function rescue = nf_one_spread_rescue_local( ...
    zMap, outerMask, eyeMask, higherThresholdMask, nearMatrix, validMask)

active = abs(zMap) > 2;
active(~validMask) = false;
active(outerMask | eyeMask) = false;
higherActive = active & higherThresholdMask;
active(higherActive & abs(zMap) < 2.5) = false;
indices = find(active);
if numel(indices) < 3
    rescue = false;
    return
end

connections = sum(nearMatrix(indices, indices), 2);
rescue = any(connections >= 2);
end

function periodic = nf_periodic_rescue_local( ...
    activations, EEG, candidateMask, options)

nComponents = size(activations, 1);
periodic = nf_empty_periodic_local(nComponents, options);
periodic.candidateMask = candidateMask;
periodic.componentIndices = find(candidateMask).';
if options.periodicRescue == false
    periodic.status = 'notRequested';
    return
end
if ~any(candidateMask)
    periodic.status = 'noCandidates';
    return
end
requiredHelpers = {'hnf_specparam_fit_compact', ...
    'hnf_specparam_frequency_cache', 'hnf_robust_ap_fit', ...
    'hnf_simple_ap_fit', 'hnf_specparam_fixed_ap_fit', ...
    'hnf_gen_aperiodic', 'hnf_fit_peaks', 'hnf_fit_peak_guess', ...
    'hnf_specparam_peak_model', 'hnf_gaussian', 'hnf_error_model', ...
    'hnf_error_expo_nk_function'};
availableHelpers = false(size(requiredHelpers));
for helper = 1:numel(requiredHelpers)
    availableHelpers(helper) = ...
        exist(requiredHelpers{helper}, 'file') > 0;
end
if ~all(availableHelpers)
    missingHelpers = requiredHelpers(~availableHelpers);
    error('nf_adjust_core:MissingSpecparam', ...
        ['Native periodic rescue is missing required specparam ' ...
        'helpers: %s.'], strjoin(missingHelpers, ', '));
end

samplingRate = double(EEG.srate);
windowSamples = round(samplingRate);
if windowSamples < 2
    periodic.status = 'insufficientFrequencySupport';
    return
end

totalSamples = size(activations, 2);
trialCount = size(EEG.data, 3);
samplesPerTrial = size(EEG.data, 2);
if trialCount > 1
    withinTrialStarts = ...
        1:windowSamples:(samplesPerTrial - windowSamples + 1);
    trialOffsets = (0:trialCount - 1) .* samplesPerTrial;
    starts = reshape( ...
        withinTrialStarts(:) + trialOffsets, 1, []);
    tailSamples = trialCount .* mod(samplesPerTrial, windowSamples);
    bufferMode = 'perTrialNonoverlappingOneSecond';
else
    starts = 1:windowSamples:(totalSamples - windowSamples + 1);
    tailSamples = mod(totalSamples, windowSamples);
    bufferMode = 'globalNonoverlappingOneSecond';
end
cuts = nf_data_cuts_local(EEG);
[starts, crossingDropped] = ...
    nf_drop_crossing_windows_local(starts, windowSamples, cuts);
periodic.buffer.windowSamples = windowSamples;
periodic.buffer.windowSeconds = windowSamples ./ samplingRate;
periodic.buffer.mode = bufferMode;
periodic.buffer.chunkCount = numel(starts);
periodic.buffer.crossingChunksDropped = crossingDropped;
periodic.buffer.tailSamplesDropped = tailSamples;
periodic.buffer.cutLatencies = cuts;
if isempty(starts)
    periodic.status = 'noCompleteChunks';
    return
end

nPositive = floor(windowSamples ./ 2) + 1;
frequency = (0:nPositive - 1) .* samplingRate ./ windowSamples;
nyquist = samplingRate ./ 2;
fitMaximum = min(options.spectrumRange(2), nyquist);
fitMask = frequency >= options.spectrumRange(1) & ...
    frequency <= fitMaximum & frequency > 0;
fitFrequency = frequency(fitMask);
if numel(fitFrequency) < 4 || ...
        fitFrequency(1) > options.rescueBand(1) || ...
        fitFrequency(end) < options.rescueBand(2)
    periodic.status = 'insufficientFrequencySupport';
    return
end

candidateIndices = find(candidateMask);
nCandidates = numel(candidateIndices);
linearPower = nan(nCandidates, numel(fitFrequency));
logPower = nan(nCandidates, numel(fitFrequency));
indices = (0:windowSamples - 1).' + starts;
if windowSamples == 1
    window = 1;
else
    sample = (0:windowSamples - 1).';
    window = 0.54 - 0.46 .* cos( ...
        2 .* pi .* sample ./ (windowSamples - 1));
end
windowEnergy = sum(window .* window);

for candidate = 1:nCandidates
    component = candidateIndices(candidate);
    componentActivation = activations(component, :);
    chunks = reshape( ...
        componentActivation(indices), windowSamples, []);
    chunks = chunks .* window;
    coefficients = fft(chunks, windowSamples, 1);
    power = abs(coefficients(1:nPositive, :)) .^ 2 ./ ...
        (samplingRate .* windowEnergy);
    if rem(windowSamples, 2) == 0
        if nPositive > 2
            power(2:end - 1, :) = 2 .* power(2:end - 1, :);
        end
    else
        if nPositive > 1
            power(2:end, :) = 2 .* power(2:end, :);
        end
    end
    meanPower = mean(power, 2).';
    meanPower = meanPower(fitMask);
    powerScale = max(meanPower);
    powerFloor = max(realmin('double'), ...
        powerScale .* eps('double'));
    meanPower = max(meanPower, powerFloor);
    linearPower(candidate, :) = meanPower;
    logPower(candidate, :) = log10(meanPower);
end

periodic.frequencyHz = fitFrequency;
periodic.linearPower = linearPower;
periodic.log10Power = logPower;
frequencyCache = hnf_specparam_frequency_cache(fitFrequency);

for candidate = 1:nCandidates
    component = candidateIndices(candidate);
    spectrum = logPower(candidate, :);
    aperiodicGuess = -(spectrum(end) - spectrum(1)) ./ ...
        log10(fitFrequency(end) ./ fitFrequency(1));
    periodic.aperiodicExponentGuess(component) = ...
        aperiodicGuess;
    narrowOptions = nf_specparam_options_local( ...
        options.narrowPeakWidth, 6, 2, aperiodicGuess);
    broadOptions = nf_specparam_options_local( ...
        options.plateauPeakWidth, 1, 1, aperiodicGuess);
    try
        [~, ~, ~, ~, narrowR2, narrowMae, narrowAP, ...
            narrowPeaks] = hnf_specparam_fit_compact( ...
            spectrum, fitFrequency, narrowOptions, ...
            frequencyCache, true);
        narrowEvidence = nf_peak_evidence_local( ...
            narrowPeaks, options.rescueBand, ...
            options.narrowPeakWidth, fitFrequency);
        periodic.narrowMask(component) = narrowEvidence;
        periodic.narrowPeakParameters{component} = narrowPeaks;
        periodic.narrowPeakBandwidthHz{component} = ...
            nf_peak_bandwidth_local(narrowPeaks);
        periodic.narrowAperiodicParameters{component} = narrowAP;
        periodic.narrowRSquare(component) = narrowR2;
        periodic.narrowMAE(component) = narrowMae;
        periodic.narrowFitStatus{component} = 'evaluated';
    catch narrowError
        periodic.narrowFitStatus{component} = ...
            ['failed:' narrowError.identifier];
    end
    try
        [~, ~, ~, ~, broadR2, broadMae, broadAP, ...
            broadPeaks] = hnf_specparam_fit_compact( ...
            spectrum, fitFrequency, broadOptions, ...
            frequencyCache, true);
        broadEvidence = nf_peak_evidence_local( ...
            broadPeaks, options.rescueBand, ...
            options.plateauPeakWidth, fitFrequency);
        periodic.broadMask(component) = broadEvidence;
        periodic.broadPeakParameters{component} = broadPeaks;
        periodic.broadPeakBandwidthHz{component} = ...
            nf_peak_bandwidth_local(broadPeaks);
        periodic.broadAperiodicParameters{component} = broadAP;
        periodic.broadRSquare(component) = broadR2;
        periodic.broadMAE(component) = broadMae;
        periodic.broadFitStatus{component} = 'evaluated';
    catch broadError
        periodic.broadFitStatus{component} = ...
            ['failed:' broadError.identifier];
    end
    narrowEvaluated = strcmp( ...
        periodic.narrowFitStatus{component}, 'evaluated');
    broadEvaluated = strcmp( ...
        periodic.broadFitStatus{component}, 'evaluated');
    if narrowEvaluated && broadEvaluated
        periodic.fitStatus{component} = 'evaluated';
    elseif narrowEvaluated || broadEvaluated
        periodic.fitStatus{component} = 'partiallyEvaluated';
    else
        periodic.fitStatus{component} = 'failed';
    end
end

periodic.rescueMask = candidateMask & ...
    (periodic.narrowMask | periodic.broadMask);
candidateStatus = periodic.fitStatus(candidateIndices);
if all(strcmp(candidateStatus, 'evaluated'))
    periodic.status = 'evaluated';
elseif any(strcmp(candidateStatus, 'evaluated')) || ...
        any(strcmp(candidateStatus, 'partiallyEvaluated'))
    periodic.status = 'partiallyEvaluated';
else
    periodic.status = 'fitFailed';
end
end

function periodic = nf_empty_periodic_local(nComponents, options)

periodic = struct();
periodic.requested = options.periodicRescue;
periodic.status = 'notEvaluated';
periodic.candidateMask = false(nComponents, 1);
periodic.componentIndices = [];
periodic.narrowMask = false(nComponents, 1);
periodic.broadMask = false(nComponents, 1);
periodic.rescueMask = false(nComponents, 1);
periodic.frequencyHz = [];
periodic.linearPower = [];
periodic.log10Power = [];
periodic.aperiodicExponentGuess = nan(nComponents, 1);
periodic.narrowPeakParameters = cell(nComponents, 1);
periodic.broadPeakParameters = cell(nComponents, 1);
periodic.narrowPeakBandwidthHz = cell(nComponents, 1);
periodic.broadPeakBandwidthHz = cell(nComponents, 1);
periodic.narrowAperiodicParameters = cell(nComponents, 1);
periodic.broadAperiodicParameters = cell(nComponents, 1);
periodic.narrowRSquare = nan(nComponents, 1);
periodic.broadRSquare = nan(nComponents, 1);
periodic.narrowMAE = nan(nComponents, 1);
periodic.broadMAE = nan(nComponents, 1);
periodic.fitStatus = repmat({'notCandidate'}, nComponents, 1);
periodic.narrowFitStatus = repmat({'notCandidate'}, nComponents, 1);
periodic.broadFitStatus = repmat({'notCandidate'}, nComponents, 1);
periodic.buffer = struct();
periodic.options = struct();
periodic.options.spectrumRange = options.spectrumRange;
periodic.options.spectrumRangeRationale = ...
    ['NF defaults to 1-40 Hz to avoid the MADE preprocessing ' ...
    'roll-off near 45 Hz; the released MADE fit extended to 50 Hz'];
periodic.options.rescueBand = options.rescueBand;
periodic.options.narrowPeakWidth = options.narrowPeakWidth;
periodic.options.broadPeakWidth = options.plateauPeakWidth;
periodic.options.narrowMaxPeaks = 6;
periodic.options.broadMaxPeaks = 1;
periodic.options.narrowPeakThreshold = 2;
periodic.options.broadPeakThreshold = 1;
periodic.options.minPeakHeight = 0;
periodic.options.aperiodicMode = 'fixed';
periodic.options.peakType = 'gaussian';
periodic.options.threshAfter = 1;
periodic.options.optim = 0;
periodic.options.aperiodicGuessFormula = ...
    '-(logPowerLast-logPowerFirst)/log10(fmax/fmin)';
periodic.options.removeMean = false;
periodic.options.window = 'symmetricHamming';
periodic.options.averageDomain = 'linearPower';
periodic.options.logTransform = 'log10AfterAverage';
end

function options = nf_specparam_options_local( ...
    peakWidth, maxPeaks, peakThreshold, aperiodicGuess)

options = struct();
options.peakWidthLims = peakWidth;
options.maxPeaks = maxPeaks;
options.minPeakHeight = 0;
options.aPeriodicMode = 'fixed';
options.peakThreshold = peakThreshold;
options.peakType = 'gaussian';
options.threshAfter = 1;
options.optim = 0;
options.ag = aperiodicGuess;
options.plt = 0;
end

function evidence = nf_peak_evidence_local( ...
    peaks, rescueBand, widthLimits, frequency)

evidence = false;
if isempty(peaks) || size(peaks, 2) ~= 3
    return
end

center = peaks(:, 1);
height = peaks(:, 2);
sigma = peaks(:, 3);
bandwidth = 2 .* sigma;
valid = all(isfinite(peaks), 2);
valid = valid & center >= rescueBand(1) & center <= rescueBand(2);
valid = valid & height >= 0;
valid = valid & bandwidth >= widthLimits(1) & ...
    bandwidth <= widthLimits(2);
valid = valid & center - sigma > frequency(1) & ...
    center + sigma < frequency(end);
evidence = any(valid);
end

function bandwidth = nf_peak_bandwidth_local(peaks)

if isempty(peaks) || size(peaks, 2) ~= 3
    bandwidth = zeros(0, 1);
else
    bandwidth = 2 .* peaks(:, 3);
end
end

function nf_warn_unavailable_local(evaluated, status)

names = fieldnames(evaluated);
for index = 1:numel(names)
    name = names{index};
    if evaluated.(name) == false
        warning('nf_adjust_core:ClassNotEvaluated', ...
            '%s was not evaluated: %s.', name, status.(name));
    end
end
end

function nf_warn_spread_unavailable_local(spread)

missing = {};
if strcmp(spread.horizontalStatus, ...
        'notEvaluatedMissingGeometry')
    missing{end + 1} = 'horizontal';
end
if strcmp(spread.blinkStatus, ...
        'notEvaluatedMissingGeometry')
    missing{end + 1} = 'blink';
end
if ~isempty(missing)
    warning('nf_adjust_core:SpreadRescueNotEvaluated', ...
        ['Spatial-spread rescue was not evaluated for %s ' ...
        'candidate components because the montage lacks required ' ...
        'adaptive role or neighbor support. Classification continues.'], ...
        strjoin(missing, ' and '));
end
end

function nf_warn_periodic_unavailable_local(periodic)

warningStatuses = {'noCompleteChunks', ...
    'insufficientFrequencySupport', 'fitFailed', ...
    'partiallyEvaluated'};
if periodic.requested && ismember(periodic.status, warningStatuses)
    warning('nf_adjust_core:PeriodicRescueIncomplete', ...
        ['Periodic rescue was requested but its status is %s. ' ...
        'Components without valid periodic evidence are not rescued.'], ...
        periodic.status);
end
end

function classification = nf_result_local( ...
    profile, options, decomposition, features, thresholds, geometry, ...
    rawByClass, evaluatedByClass, statusByClass, rawCandidateMask, ...
    spread, postSpatialByClass, postSpatialCandidateMask, periodic, ...
    finalByClass, rescueMask, rejectedMask, rejected)

classification = struct();
classification.schemaVersion = '2.0.0';
if strcmp(profile, 'adult')
    classification.name = 'nfADJUST';
    classification.rule = ...
        'native ADJUST feature conjunctions; no rescue by default';
else
    classification.name = 'nfadjustedADJUST';
    classification.rule = ...
        ['native Beall horizontal/blink plus ADJUST vertical/discontinuity; ' ...
        'class-specific spread rescue followed by union-level periodic rescue'];
end
classification.profile = profile;
classification.classNames = {'horizontalEyeMovement', ...
    'verticalEyeMovement', 'blink', 'discontinuity'};
classification.evaluatedByClass = evaluatedByClass;
classification.statusByClass = statusByClass;
classification.rawCandidateMasksByClass = rawByClass;
classification.rawMasksByClass = rawByClass;
classification.rawCandidatesByClass = struct();
classification.rawCandidatesByClass.horizontalEyeMovement = ...
    find(rawByClass.horizontalEyeMovement).';
classification.rawCandidatesByClass.verticalEyeMovement = ...
    find(rawByClass.verticalEyeMovement).';
classification.rawCandidatesByClass.blink = ...
    find(rawByClass.blink).';
classification.rawCandidatesByClass.discontinuity = ...
    find(rawByClass.discontinuity).';
classification.rawCandidateMask = rawCandidateMask;
classification.rawCandidates = find(rawCandidateMask).';
classification.spatialRescue = spread;
classification.postSpatialCandidateMasksByClass = ...
    postSpatialByClass;
classification.postSpatialCandidateMask = postSpatialCandidateMask;
classification.postSpatialCandidates = ...
    find(postSpatialCandidateMask).';
classification.periodicEvidence = periodic;
classification.rescueMask = rescueMask;
classification.rescued = find(rescueMask).';
classification.rescueEvidenceMask = ...
    spread.rescueMask | periodic.rescueMask;
classification.rescueRule = ...
    ['spread rescue is eye-class-specific; periodic rescue applies ' ...
    'to the post-spread class union; rescueMask marks components ' ...
    'fully removed from the final union'];
classification.rejectedMask = rejectedMask;
classification.rejected = rejected;
classification.finalCandidateMasksByClass = finalByClass;
classification.finalMasksByClass = finalByClass;
classification.finalCandidatesByClass = struct();
for classIndex = 1:numel(classification.classNames)
    className = classification.classNames{classIndex};
    classification.finalCandidatesByClass.(className) = ...
        find(finalByClass.(className)).';
end
classification.horizontalEyeMovement = ...
    classification.finalCandidatesByClass.horizontalEyeMovement;
classification.verticalEyeMovement = ...
    classification.finalCandidatesByClass.verticalEyeMovement;
classification.blink = ...
    classification.finalCandidatesByClass.blink;
classification.discontinuity = ...
    classification.finalCandidatesByClass.discontinuity;
classification.finalCandidateMask = rejectedMask;
classification.finalCandidates = rejected;
classification.features = features;
classification.thresholds = thresholds;
classification.geometry = geometry;
classification.componentCount = decomposition.nComponents;
classification.parameters = struct();
classification.parameters.epochLength = options.epochLength;
classification.parameters.temporalEpochRule = ...
    '3-D trials are authoritative; epochLength applies to continuous data';
classification.parameters.periodicRescue = options.periodicRescue;
classification.parameters.spatialRescue = options.spatialRescue;
classification.parameters.rescueBand = options.rescueBand;
classification.parameters.spectrumRange = options.spectrumRange;
classification.parameters.narrowPeakWidth = options.narrowPeakWidth;
classification.parameters.plateauPeakWidth = options.plateauPeakWidth;
classification.parameters.referenceChanlocs = struct();
classification.parameters.referenceChanlocs.supplied = ...
    ~isempty(options.referenceChanlocs);
classification.parameters.referenceChanlocs.channelCount = ...
    numel(geometry.referenceLabels);
classification.parameters.referenceChanlocs.labels = ...
    geometry.referenceLabels;
classification.parameters.emitWarnings = options.emitWarnings;
classification.parameters.trimFraction = 0.01;
classification.parameters.gdsfNeighborCount = ...
    geometry.gdsfNeighborCount;
classification.parameters.emTolerance = 1e-6;
classification.parameters.emMaxIterations = 500;
end
