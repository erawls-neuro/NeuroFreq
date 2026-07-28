function [EEG, info] = nf_filter(EEG, lowpass, highpass, notch, targetRate)
% NF_FILTER  Filtering and resampling for raw M/EEG.
%
% [EEG, INFO] = NF_FILTER(EEG, LOWPASS, HIGHPASS, NOTCH, TARGETRATE)
%
% The positional order is intentionally LOWPASS then HIGHPASS for
% compatibility with the released NeuroFreq helper. Use 0 or [] to disable
% an operation. NOTCH uses a +/-2-Hz stop band. A notch wholly above the
% retained low-pass band is skipped and recorded as redundant.

nf_validate_eeg(EEG);

if nargin < 2 || isempty(lowpass)
    lowpass = 45;
end
if nargin < 3 || isempty(highpass)
    highpass = 0.3;
end
if nargin < 4 || isempty(notch)
    notch = 60;
end
if nargin < 5 || isempty(targetRate)
    targetRate = 250;
end

nf_validate_settings(EEG, lowpass, highpass, notch, targetRate);

info = struct();
info.schemaVersion = '2.0.0';
info.started = datestr(now, 30); %#ok<TNOW1,DATST>
info.inputSrate = EEG.srate;
info.inputPnts = EEG.pnts;
info.inputTrials = EEG.trials;
info.requested.lowpassHz = lowpass;
info.requested.highpassHz = highpass;
info.requested.notchHz = notch;
info.requested.notchHalfWidthHz = 2;
info.requested.targetRateHz = targetRate;
info.applied.highpass = false;
info.applied.notch = false;
info.applied.lowpass = false;
info.applied.resample = false;
info.skipped = {};

if highpass > 0
    EEG = pop_eegfiltnew(EEG, 'locutoff', highpass);
    EEG = eeg_checkset(EEG);
    info.applied.highpass = true;
end

notchIsRedundant = lowpass > 0 && lowpass <= notch - 2;
if notch > 0 && ~notchIsRedundant
    EEG = pop_eegfiltnew(EEG, 'locutoff', notch - 2, ...
        'hicutoff', notch + 2, 'revfilt', 1);
    EEG = eeg_checkset(EEG);
    info.applied.notch = true;
elseif notch > 0
    info.skipped{end + 1} = ...
        'Notch skipped because the low-pass cutoff is below the stop band.';
end

if lowpass > 0
    EEG = pop_eegfiltnew(EEG, 'hicutoff', lowpass);
    EEG = eeg_checkset(EEG);
    info.applied.lowpass = true;
end

if targetRate > 0 && targetRate ~= EEG.srate
    EEG = pop_resample(EEG, targetRate);
    EEG = eeg_checkset(EEG);
    info.applied.resample = true;
end

info.finished = datestr(now, 30); %#ok<TNOW1,DATST>
info.outputSrate = EEG.srate;
info.outputPnts = EEG.pnts;
info.outputTrials = EEG.trials;
info.eeglabHistoryLength = 0;
if isfield(EEG, 'history') && ~isempty(EEG.history)
    info.eeglabHistoryLength = numel(EEG.history);
end

if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
EEG.etc.nf_filter = info;
if ~isfield(EEG.etc, 'nf_filter_history')
    EEG.etc.nf_filter_history = {};
elseif ~iscell(EEG.etc.nf_filter_history)
    EEG.etc.nf_filter_history = {EEG.etc.nf_filter_history};
end
EEG.etc.nf_filter_history{end + 1} = info;

end

function nf_validate_settings(EEG, lowpass, highpass, notch, targetRate)
values = {lowpass, highpass, notch, targetRate};
for valueIndex = 1:numel(values)
    value = values{valueIndex};
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
            ~isfinite(value) || value < 0
        error('nf_filter:InvalidSettings', ...
            'Filter settings must be finite, nonnegative numeric scalars.');
    end
end
if lowpass > 0 && highpass > 0 && highpass >= lowpass
    error('nf_filter:InvalidPassband', ...
        'highpass must be below lowpass. The legacy call order is lowpass, highpass.');
end

inputNyquist = EEG.srate / 2;
if lowpass >= inputNyquist
    error('nf_filter:InvalidLowpass', ...
        'lowpass must be below the input Nyquist frequency.');
end
if highpass >= inputNyquist
    error('nf_filter:InvalidHighpass', ...
        'highpass must be below the input Nyquist frequency.');
end
if notch > 0 && notch <= 2
    error('nf_filter:InvalidNotch', ...
        'notch must exceed its 2-Hz half-width.');
end
if notch > 0 && notch + 2 >= inputNyquist && ...
        ~(lowpass > 0 && lowpass <= notch - 2)
    error('nf_filter:InvalidNotch', ...
        'The notch stop band must be below the input Nyquist frequency.');
end
if targetRate > 0 && lowpass > 0 && lowpass >= targetRate / 2
    error('nf_filter:InvalidResampleRate', ...
        'lowpass must be below the target Nyquist frequency.');
end
if targetRate > 0 && highpass > 0 && highpass >= targetRate / 2
    error('nf_filter:InvalidResampleRate', ...
        'highpass must be below the target Nyquist frequency.');
end
if EEG.trials > 1 && targetRate > 0 && targetRate ~= EEG.srate
    error('nf_filter:EpochedResampling', ...
        ['Resampling epoched data can damage event provenance. Resample the ' ...
        'continuous dataset before epoching.']);
end
end

function nf_validate_eeg(EEG)
if ~isstruct(EEG) || numel(EEG) ~= 1 || ~isfield(EEG, 'data') || ...
        ~isfield(EEG, 'srate') || ~isfield(EEG, 'nbchan') || ...
        ~isfield(EEG, 'pnts') || ~isfield(EEG, 'trials')
    error('nf_filter:InvalidEEG', 'EEG must be one valid EEGLAB dataset structure.');
end
if ~isnumeric(EEG.nbchan) || ~isscalar(EEG.nbchan) || ...
        ~isfinite(EEG.nbchan) || EEG.nbchan < 1 || ...
        EEG.nbchan ~= round(EEG.nbchan) || ...
        ~isnumeric(EEG.pnts) || ~isscalar(EEG.pnts) || ...
        ~isfinite(EEG.pnts) || EEG.pnts < 2 || EEG.pnts ~= round(EEG.pnts) || ...
        ~isnumeric(EEG.trials) || ~isscalar(EEG.trials) || ...
        ~isfinite(EEG.trials) || EEG.trials < 1 || ...
        EEG.trials ~= round(EEG.trials) || ...
        ~isnumeric(EEG.data) || ~isreal(EEG.data) || isempty(EEG.data) || ...
        EEG.nbchan ~= size(EEG.data, 1) || EEG.pnts ~= size(EEG.data, 2) || ...
        EEG.trials ~= size(EEG.data, 3) || any(~isfinite(EEG.data(:))) || ...
        ~isnumeric(EEG.srate) || ~isreal(EEG.srate) || ...
        ~isscalar(EEG.srate) || EEG.srate <= 0 || ~isfinite(EEG.srate)
    error('nf_filter:InvalidEEG', 'EEG data or sampling rate is invalid.');
end
end
