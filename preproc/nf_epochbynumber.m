function [EEG, indices, com, info] = nf_epochbynumber(EEG, epochs, xmin, xmax)
% NF_EPOCHBYNUMBER  Select existing EEGLAB epochs without rebuilding events.
%
% [EEG, INDICES, COM, INFO] = NF_EPOCHBYNUMBER(EEG, EPOCHS, XMIN, XMAX)
% preserves the requested epoch order, crops to the requested time window,
% and keeps EEG.etc.behavior aligned when that metadata is present.

nf_validate_epoched_eeg(EEG);
originalTrials = EEG.trials;

if nargin < 2 || isempty(epochs)
    epochs = 1:originalTrials;
end
if nargin < 3 || isempty(xmin)
    xmin = EEG.xmin;
end
if nargin < 4 || isempty(xmax)
    xmax = EEG.xmax;
end

if islogical(epochs)
    if ~isvector(epochs) || numel(epochs) ~= originalTrials
        error('nf_epochbynumber:InvalidMask', ...
            'A logical epoch mask must contain EEG.trials elements.');
    end
    epochs = find(epochs);
end

if ~isnumeric(epochs) || ~isvector(epochs) || isempty(epochs) || ...
        any(~isfinite(epochs)) || any(epochs ~= round(epochs)) || ...
        any(epochs < 1) || any(epochs > originalTrials)
    error('nf_epochbynumber:InvalidEpochs', ...
        'epochs must contain finite integer indices within EEG.trials.');
end
epochs = unique(epochs(:)', 'stable');

if ~isnumeric(xmin) || ~isscalar(xmin) || ~isfinite(xmin) || ...
        ~isnumeric(xmax) || ~isscalar(xmax) || ~isfinite(xmax) || ...
        xmin >= xmax || xmin < EEG.xmin || xmax > EEG.xmax
    error('nf_epochbynumber:InvalidLimits', ...
        'XMIN and XMAX must be finite limits inside the existing epoch.');
end

behavior = [];
hasBehavior = isfield(EEG, 'etc') && isfield(EEG.etc, 'behavior');
if hasBehavior
    behavior = nf_subset_trial_metadata(EEG.etc.behavior, ...
        epochs, originalTrials, 'EEG.etc.behavior');
end

epochIds = 1:originalTrials;
if isfield(EEG, 'etc') && isfield(EEG.etc, 'nf_epoch_ids')
    epochIds = nf_subset_trial_metadata(EEG.etc.nf_epoch_ids, ...
        epochs, originalTrials, 'EEG.etc.nf_epoch_ids');
else
    epochIds = epochIds(epochs);
end

EEG = pop_select(EEG, 'trial', epochs, 'sorttrial', 'off', ...
    'time', [xmin xmax]);
EEG = eeg_checkset(EEG);

if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
if hasBehavior
    EEG.etc.behavior = behavior;
end
EEG.etc.nf_epoch_ids = epochIds;

indices = epochs;
info = struct();
info.schemaVersion = '2.0.0';
info.originalTrials = originalTrials;
info.selectedIndices = indices;
info.originalEpochIds = epochIds;
info.requestedTimeLimits = [xmin xmax];
info.actualTimeLimits = [EEG.xmin EEG.xmax];
info.outputTrials = EEG.trials;
info.outputPnts = EEG.pnts;

com = sprintf('EEG = nf_epochbynumber(EEG, [%s], %.12g, %.12g);', ...
    num2str(indices), xmin, xmax);

end

function value = nf_subset_trial_metadata(value, indices, originalTrials, fieldName)
if istable(value)
    if height(value) ~= originalTrials
        error('nf_epochbynumber:MetadataMismatch', ...
            '%s has %d rows but EEG has %d trials.', ...
            fieldName, height(value), originalTrials);
    end
    value = value(indices, :);
    return
end

if isvector(value) && numel(value) == originalTrials
    value = value(indices);
    return
end

if ~isempty(value) && size(value, 1) == originalTrials
    dimensions = repmat({':'}, 1, ndims(value));
    dimensions{1} = indices;
    value = value(dimensions{:});
    return
end

if isempty(value) && originalTrials == 0
    return
end

error('nf_epochbynumber:MetadataMismatch', ...
    ['%s must have one row per trial or be a vector with one element ' ...
    'per trial.'], fieldName);
end

function nf_validate_epoched_eeg(EEG)
if ~isstruct(EEG) || numel(EEG) ~= 1 || ~isfield(EEG, 'data') || ...
        ~isfield(EEG, 'trials') || ~isfield(EEG, 'pnts') || ...
        ~isfield(EEG, 'xmin') || ~isfield(EEG, 'xmax')
    error('nf_epochbynumber:InvalidEEG', ...
        'EEG must be one valid EEGLAB dataset structure.');
end
hasEpochStructure = isfield(EEG, 'epoch') && ~isempty(EEG.epoch);
if EEG.trials < 1 || (EEG.trials == 1 && ~hasEpochStructure)
    error('nf_epochbynumber:ContinuousData', ...
        'The input must contain existing EEGLAB epochs.');
end
if ~isnumeric(EEG.data) || ~isreal(EEG.data) || ...
        size(EEG.data, 2) ~= EEG.pnts || size(EEG.data, 3) ~= EEG.trials
    error('nf_epochbynumber:InvalidEEG', ...
        'EEG.data dimensions do not match EEG.pnts and EEG.trials.');
end
end
