function TF = nf_tfspecparam(TF, varargin)
% NF_TFSPECPARAM    Computes a time-resolved spectral parameterization
%
% GENERAL
% -------
% Time resolved spectral parameterization (Wilson et al., 2022) rewritten
% for use with TF structures output by nf_tftransform.m. Computes a
% time-resolved spectral parameterization using specparam.m.
%
% A nonfinite spectrum at any sensor or time removes the affected trial
% or condition across the output. Errors in the solver or configuration
% interrupt the calculation. An entirely invalid dataset raises an error.
%
% Accepts TF transforms calculated using any method. Automatically detects
% if power values are on a linear scale and converts to log scale.
%
% OUTPUT
% ------
% TF.specparam fields (all log10 scale)
%   Aperiodic:
%     TF.specparam.rsquare     : C×T×R
%     TF.specparam.mae         : C×T×R
%     TF.specparam.ap_offset   : C×T×R
%     TF.specparam.ap_knee     : C×T×R (NaN if 'fixed')
%     TF.specparam.ap_exponent : C×T×R
%   Peaks (ALL peaks per slice, parfor-safe numeric arrays; padded with NaN):
%     TF.specparam.pk_n        : C×T×R    (#peaks detected)
%     TF.specparam.pk_f        : C×K×T×R  (Hz)        K = maxPeaks
%     TF.specparam.pk_aLog10   : C×K×T×R  (log10 power height)
%     TF.specparam.pk_sd       : C×K×T×R  (Hz; Gaussian SD or Cauchy HWHM)
%
% INPUT
% -----
%   TF.power : C×F×T[×R] power
%   TF.scale : 'linear'|'log10'  (if 'linear', converted to log10)
%   TF.freqs : F×1 (Hz)
%
% Name-Value
% ----------
%   'peakWidthLims' : [min max] Hz (default [0.5 12])
%   'maxPeaks'      : max #peaks per slice; 0 disables peak fitting
%                     (default numel(TF.freqs))
%   'minPeakHeight' : minimum peak height (log10 power units if scale->log10) (default 0)
%   'aPeriodicMode' : 'fixed'|'knee' (default 'fixed')
%   'peakThreshold' : noise SD threshold for detection (default 2.0)
%   'peakType'      : 'gaussian'|'cauchy' (default 'gaussian')
%   'threshAfter'   : re-threshold when no optim (default 1)
%   'optim'         : use Optimization Toolbox if available (default 1)
%   'FoD'           : 'fit' returns model fit; 'data' returns modeled data (default 'data')
%


p = inputParser;

valid2Scalar = @(x) isnumeric(x) && isreal(x) && numel(x) == 2 && ...
    all(isfinite(x(:))) && x(1) > 0 && x(2) >= x(1);
validThreshold = @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x >= 0;
expectedPeakType = @(x) any(validatestring(x, {'gaussian', 'cauchy'}));
expectedAPModes = @(x) any(validatestring(x, {'knee', 'fixed'}));
expectedFoD = @(x) any(validatestring(x, {'fit', 'data'}));

addRequired(p, 'TF');

addParameter(p, 'peakWidthLims', [0.5 12], valid2Scalar);
addParameter(p, 'maxPeaks', [], @nf_valid_max_peaks_local);
addParameter(p, 'minPeakHeight', 0, validThreshold);
addParameter(p, 'aPeriodicMode', 'fixed', expectedAPModes);
addParameter(p, 'peakThreshold', 2.0, validThreshold);
addParameter(p, 'peakType', 'gaussian', expectedPeakType);
addParameter(p, 'threshAfter', 1, @(x) isscalar(x) && ismember(x, [0 1]));
addParameter(p, 'optim', 1, @(x) isscalar(x) && ismember(x, [0 1]));
addParameter(p, 'FoD', 'data', expectedFoD);

parse(p, TF, varargin{:});
opt = rmfield(p.Results, 'TF');
clear p
opt.aPeriodicMode = validatestring(opt.aPeriodicMode, {'fixed','knee'});
opt.peakType = validatestring(opt.peakType, {'gaussian','cauchy'});
opt.FoD = validatestring(opt.FoD, {'fit','data'});

if ~isfield(TF, 'freqs')
    error('nf_tfspecparam:FreqsRequired','TF.freqs is required.');
end

if ~isfield(TF, 'scale')
    error('nf_tfspecparam:ScaleRequired','TF.scale is required.');
end

if ~isfield(TF, 'nsensor')
    error('nf_tfspecparam:NSensorRequired','TF.nsensor is required.');
end

if ~isfield(TF, 'ntrls')
    error('nf_tfspecparam:NTrialsRequired','TF.ntrls is required.');
end

if isnumeric(TF.nsensor) == false
    error('nf_tfspecparam:BadSensorCount', ...
        'TF.nsensor must be a numeric scalar count.');
end

if isscalar(TF.nsensor) == false
    error('nf_tfspecparam:BadSensorCount', ...
        'TF.nsensor must be a numeric scalar count.');
end

if isnan(TF.nsensor) == true || isinf(TF.nsensor) == true
    error('nf_tfspecparam:BadSensorCount', ...
        'TF.nsensor must be a positive integer count.');
end

if TF.nsensor < 1 || TF.nsensor ~= fix(TF.nsensor)
    error('nf_tfspecparam:BadSensorCount', ...
        'TF.nsensor must be a positive integer count.');
end

frequencyCache = hnf_specparam_frequency_cache(TF.freqs);
if strcmp(opt.aPeriodicMode,'knee') && numel(TF.freqs) < 3
    error('nf_tfspecparam:TooFewFrequencies','A three-parameter knee fit needs at least three frequencies.');
end
freqs = TF.freqs(:);
freqRow = double(freqs).';

if isempty(opt.maxPeaks)
    opt.maxPeaks = numel(freqs);
end

fitOptions = nf_compact_options_local(opt);

if ~(strcmpi(TF.scale, 'linear') || strcmpi(TF.scale, 'log10'))
    error('nf_tfspecparam:TFScaleWrong','TF.scale must be ''linear'' or ''log10''.');
end

powFields = nf_find_power_fields_local(TF);

if ~any(strcmpi(powFields,'power'))
    error('nf_tfspecparam:MissingPower','A numeric TF.power field is required.');
end
if strcmpi(TF.scale, 'linear')
    for i = 1:numel(powFields)
        name = powFields{i};
        X = TF.(name);
        if ~isreal(X)
            error('nf_tfspecparam:ComplexPower', 'Power field %s must be real.', name);
        end
        if ~isfloat(X), X = double(X); end
        X(X == 0) = eps; % Preserve the existing floor for a zero power bin.
        X(X < 0) = NaN; % Signed energy has no real log10-power representation.
        TF.(name) = log10(X);
    end
    TF.scale = 'log10';
end

if isfield(TF, 'conds')
    nSeries = TF.conds;
else
    nSeries = TF.ntrls;
end

if isnumeric(nSeries) == false
    error('nf_tfspecparam:BadSeriesCount', ...
        'TF.ntrls or TF.conds must be a numeric scalar count.');
end

if isscalar(nSeries) == false
    error('nf_tfspecparam:BadSeriesCount', ...
        'TF.ntrls or TF.conds must be a numeric scalar count.');
end

if isnan(nSeries) == true || isinf(nSeries) == true
    error('nf_tfspecparam:BadSeriesCount', ...
        'TF.ntrls or TF.conds must be a positive integer count.');
end

if nSeries < 1 || nSeries ~= fix(nSeries)
    error('nf_tfspecparam:BadSeriesCount', ...
        'TF.ntrls or TF.conds must be a positive integer count.');
end

badTrialsAll = false(1, nSeries);

disp('[nf_tfspecparam]: beginning time-frequency spectral parameterization');

for i = 1:numel(powFields)

    srcName = powFields{i};

    if nf_is_param_power_field_local(srcName) == 1
        continue
    end

    surf = TF.(srcName);

    [ap_log, osc_log, specparam, badTrialsThis] = ...
        nf_tfspecparam_core_local( ...
        surf, freqs, opt, strcmpi(srcName, 'power'), ...
        TF.nsensor, nSeries, frequencyCache, fitOptions);
    % The one-sensor transform uses F x T x R. Return derived surfaces in
    % the same layout as their input, including singleton sensor/time axes.
    ap_log = reshape(ap_log,size(surf));
    osc_log = reshape(osc_log,size(surf));
    clear surf

    TF.([srcName '_osc']) = osc_log;
    TF.([srcName '_ap']) = ap_log;

    if isfield(TF,'conds')
        if numel(badTrialsThis) ~= TF.conds
            error('nf_tfspecparam:BadCondMaskMisMatch','Bad-trial mask length mismatch: got %d, expected %d.', numel(badTrialsThis), TF.conds);
        end
    elseif numel(badTrialsThis) ~= TF.ntrls
        error('nf_tfspecparam:BadTrialMaskMisMatch','Bad-trial mask length mismatch: got %d, expected %d.', numel(badTrialsThis), TF.ntrls);
    end

    badTrialsAll = badTrialsAll | badTrialsThis(:).';

    if strcmpi(srcName, 'power')
        TF.specparam = specparam;
        TF.specparam.method = 'nf_tfspecparam';
        TF.specparam.options = opt;
        if strcmp(opt.peakType,'gaussian')
            TF.specparam.peak_width_definition = 'Gaussian standard deviation (Hz)';
        else
            TF.specparam.peak_width_definition = 'Cauchy half-width at half-maximum (Hz)';
        end
        TF.specparam.power_osc_field = 'power_osc';
        TF.specparam.power_ap_field = 'power_ap';
    end

end

if any(badTrialsAll)
    keepIdx = ~badTrialsAll;
    removedIdx = find(badTrialsAll);

    if ~any(keepIdx)
        error('nf_tfspecparam:NoValidSpectra', ...
            'Every trial or condition contained a failed spectrum. Inspect the frequency axis, power and scale.');
    end
    TF = nf_remove_trials_local(TF, keepIdx);

    if isfield(TF, 'specparam') && isstruct(TF.specparam)
        TF.specparam.removed_trials = removedIdx(:);
        TF.specparam.n_removed_trials = numel(removedIdx);
        TF.specparam.removal_reason = 'nonfinite power in at least one fitted spectrum';
    end
end

TF = nf_powerfront_local(TF);

disp('[nf_tfspecparam]: time-frequency spectral parameterization complete');

end

function [ap_log, osc_log, specparam, badTrials] = ...
    nf_tfspecparam_core_local( ...
    surf, freqs, opt, keepspecparam, nSensors, nSeries, ...
    frequencyCache, fitOptions)

F = numel(freqs);
if ~isreal(surf)
    error('nf_tfspecparam:ComplexPower','Power fields must be real-valued.');
end
if ~isfloat(surf), surf = double(surf); end

surf = nf_canonical_power_local( ...
    surf, F, nSensors, nSeries);

[C, Fstored, T, R] = size4d_local(surf);

if Fstored ~= F
    error('nf_tfspecparam:BadPowerSecondDim', ...
        ['Second dimension of power (F=%d) differs from ' ...
        'numel(freqs) (%d).'], Fstored, F);
end

K = opt.maxPeaks;
storeParameters = keepspecparam == 1;

if K > 200
    warning('nf_tfspecparam:PeaksGonnaBBig', ...
        ['maxPeaks is %d. Peak arrays will be large: C×K×T×R. ' ...
        'Consider reducing maxPeaks.'], K);
end

nSpectra = C * T * R;

spectra = reshape( ...
    permute(surf, [2 1 3 4]), ...
    [F nSpectra]);

dataPrototype = surf(1);
clear surf
freqRow = frequencyCache.freqs;
returnFit = strcmpi(opt.FoD, 'fit');
fixedMode = strcmpi(opt.aPeriodicMode, 'fixed');

if fixedMode == true && K == 0
    fixedPredictor = -frequencyCache.log10Freqs(:);
    fixedPredictorCentered = fixedPredictor - mean(fixedPredictor);
    fixedPredictorSumSquares = sum( ...
        fixedPredictorCentered .* fixedPredictorCentered);

    validFixedPredictor = all(isfinite(fixedPredictor));
    validFixedPredictor = validFixedPredictor ...
        && all(freqRow > 0) ...
        && fixedPredictorSumSquares > 0;

    if validFixedPredictor == true
        [ap_log, osc_log, specparam, badTrials] = ...
            nf_fixed_no_peaks_matrix_local( ...
            spectra, fixedPredictor, ...
            fixedPredictorCentered, ...
            fixedPredictorSumSquares, ...
            C, T, R, dataPrototype, ...
            returnFit, storeParameters);
        return
    end
end

apFlat = nan(F, nSpectra, 'like', dataPrototype);
oscFlat = nan(F, nSpectra, 'like', dataPrototype);

badSlice = false(1, nSpectra);

if storeParameters == true
    rsqFlat = nan(1, nSpectra);
    maeFlat = nan(1, nSpectra);

    apOffsetFlat = nan(1, nSpectra);
    apKneeFlat = nan(1, nSpectra);
    apExponentFlat = nan(1, nSpectra);

    pkNFlat = nan(1, nSpectra);
    pkFFlat = nan(K, nSpectra, 'like', dataPrototype);
    pkAFlat = nan(K, nSpectra, 'like', dataPrototype);
    pkSdFlat = nan(K, nSpectra, 'like', dataPrototype);

    parfor sliceIndex = 1:nSpectra
        apSlice = nan(F, 1, 'like', dataPrototype);
        oscSlice = nan(F, 1, 'like', dataPrototype);

        rsqSlice = NaN;
        maeSlice = NaN;

        apOffsetSlice = NaN;
        apKneeSlice = NaN;
        apExponentSlice = NaN;

        pkNSlice = NaN;
        pkFSlice = nan(K, 1, 'like', dataPrototype);
        pkASlice = nan(K, 1, 'like', dataPrototype);
        pkSdSlice = nan(K, 1, 'like', dataPrototype);

        dataX = spectra(:, sliceIndex).';

        sliceFailed = any(~isfinite(dataX));
        if ~sliceFailed
            [apFit, periodicFit, apData, periodicData, ...
                rSquare, maeValue, apParms, peakParms] = ...
                hnf_specparam_fit_compact( ...
                dataX, freqRow, fitOptions, frequencyCache, true);

            if returnFit == true
                apSlice = apFit(:);
                oscSlice = periodicFit(:);
            else
                apSlice = apData(:);
                oscSlice = periodicData(:);
            end

            rsqSlice = rSquare;
            maeSlice = maeValue;

            apVector = apParms(:).';

            apOffsetSlice = apVector(1);

            if fixedMode == true
                apKneeSlice = NaN;
                apExponentSlice = apVector(2);
            else
                apKneeSlice = apVector(2);
                apExponentSlice = apVector(3);
            end

            nPeaks = size(peakParms, 1);
            pkNSlice = nPeaks;

            if nPeaks > 0
                nUse = min(nPeaks, K);

                if nUse > 0
                    pkFSlice(1:nUse) = peakParms(1:nUse, 1);
                    pkASlice(1:nUse) = peakParms(1:nUse, 2);
                    pkSdSlice(1:nUse) = peakParms(1:nUse, 3);
                end
            end
        end

        apFlat(:, sliceIndex) = apSlice;
        oscFlat(:, sliceIndex) = oscSlice;

        rsqFlat(sliceIndex) = rsqSlice;
        maeFlat(sliceIndex) = maeSlice;

        apOffsetFlat(sliceIndex) = apOffsetSlice;
        apKneeFlat(sliceIndex) = apKneeSlice;
        apExponentFlat(sliceIndex) = apExponentSlice;

        pkNFlat(sliceIndex) = pkNSlice;
        pkFFlat(:, sliceIndex) = pkFSlice;
        pkAFlat(:, sliceIndex) = pkASlice;
        pkSdFlat(:, sliceIndex) = pkSdSlice;

        badSlice(sliceIndex) = sliceFailed;
    end
else
    parfor sliceIndex = 1:nSpectra
        apSlice = nan(F, 1, 'like', dataPrototype);
        oscSlice = nan(F, 1, 'like', dataPrototype);

        dataX = spectra(:, sliceIndex).';

        sliceFailed = any(~isfinite(dataX));
        if ~sliceFailed
            [apFit, periodicFit, apData, periodicData] = ...
                hnf_specparam_fit_compact( ...
                dataX, freqRow, fitOptions, frequencyCache, false);

            if returnFit == true
                apSlice = apFit(:);
                oscSlice = periodicFit(:);
            else
                apSlice = apData(:);
                oscSlice = periodicData(:);
            end
        end

        apFlat(:, sliceIndex) = apSlice;
        oscFlat(:, sliceIndex) = oscSlice;
        badSlice(sliceIndex) = sliceFailed;
    end
end

clear spectra

ap_log = permute( ...
    reshape(apFlat, [F C T R]), ...
    [2 1 3 4]);
clear apFlat
osc_log = permute( ...
    reshape(oscFlat, [F C T R]), ...
    [2 1 3 4]);
clear oscFlat

badBySeries = reshape(badSlice, [C * T R]);
badTrials = any(badBySeries, 1);

if storeParameters == true
    rsq = reshape(rsqFlat, [C T R]);
    mae = reshape(maeFlat, [C T R]);

    ap_offset = reshape(apOffsetFlat, [C T R]);
    ap_knee = reshape(apKneeFlat, [C T R]);
    ap_exponent = reshape(apExponentFlat, [C T R]);

    pk_n = reshape(pkNFlat, [C T R]);
    pk_f = double(permute( ...
        reshape(pkFFlat, [K C T R]), ...
        [2 1 3 4]));
    clear pkFFlat
    pk_aLog10 = double(permute( ...
        reshape(pkAFlat, [K C T R]), ...
        [2 1 3 4]));
    clear pkAFlat
    pk_sd = double(permute( ...
        reshape(pkSdFlat, [K C T R]), ...
        [2 1 3 4]));
    clear pkSdFlat

    specparam = struct();
    specparam.rsquare = rsq;
    specparam.mae = mae;
    specparam.ap_offset = ap_offset;
    specparam.ap_knee = ap_knee;
    specparam.ap_exponent = ap_exponent;
    specparam.pk_n = pk_n;
    specparam.pk_f = pk_f;
    specparam.pk_aLog10 = pk_aLog10;
    specparam.pk_sd = pk_sd;
else
    specparam = struct();
end

end

function [ap_log, osc_log, specparam, badTrials] = ...
    nf_fixed_no_peaks_matrix_local( ...
    spectra, predictor, predictorCentered, ...
    predictorSumSquares, C, T, R, dataPrototype, ...
    returnFit, storeParameters)

[F, ~] = size(spectra);
useNative = numel(spectra) >= 65536 && (isa(spectra,'single') || isa(spectra,'double')) && isreal(spectra) && ...
    ~issparse(spectra) && exist('hnf_mex_available','file') == 2 && ...
    hnf_mex_available('hnf_fixed_spectra_mex');
if useNative
    args = {spectra, predictor, predictorCentered, mean(predictor), ...
        predictorSumSquares, double(returnFit)};
    try
        if storeParameters
            [apFlat, oscFlat, badSlice, summary] = hnf_fixed_spectra_mex(args{:});
            offsetFlat = summary(1,:);
            exponentFlat = summary(2,:);
            maeFlat = summary(3,:);
            rsqFlat = summary(4,:);
            pkNFlat = zeros(size(badSlice));
            pkNFlat(badSlice) = NaN;
        else
            [apFlat, oscFlat, badSlice] = hnf_fixed_spectra_mex(args{:});
        end
    catch err
        if ~hnf_mex_load_error(err), rethrow(err); end
        useNative = false;
    end
end
if ~useNative
    [apFlat, oscFlat, badSlice, offsetFlat, exponentFlat, maeFlat, rsqFlat, pkNFlat] = ...
        nf_fixed_arrays_matlab_local(spectra, predictor, predictorCentered, ...
        predictorSumSquares, dataPrototype, returnFit, storeParameters);
end
clear spectra args

ap_log = permute( ...
    reshape(apFlat, [F C T R]), ...
    [2 1 3 4]);
clear apFlat
osc_log = permute( ...
    reshape(oscFlat, [F C T R]), ...
    [2 1 3 4]);
clear oscFlat

badBySeries = reshape(badSlice, [C * T R]);
badTrials = any(badBySeries, 1);

if storeParameters == true
    specparam = struct();
    specparam.rsquare = reshape(rsqFlat, [C T R]);
    specparam.mae = reshape(maeFlat, [C T R]);
    specparam.ap_offset = reshape(offsetFlat, [C T R]);
    specparam.ap_knee = nan(C, T, R);
    specparam.ap_exponent = reshape(exponentFlat, [C T R]);
    specparam.pk_n = reshape(pkNFlat, [C T R]);
    specparam.pk_f = nan(C, 0, T, R);
    specparam.pk_aLog10 = nan(C, 0, T, R);
    specparam.pk_sd = nan(C, 0, T, R);
else
    specparam = struct();
end

end

function [apFlat, oscFlat, badSlice, offsetFlat, exponentFlat, maeFlat, rsqFlat, pkNFlat] = ...
    nf_fixed_arrays_matlab_local(spectra, predictor, predictorCentered, ...
    predictorSumSquares, dataPrototype, returnFit, storeParameters)
% Block the MATLAB implementation to bound temporary storage on large TF data.
[F, nSpectra] = size(spectra);

apFlat = nan(F, nSpectra, 'like', dataPrototype);
oscFlat = nan(F, nSpectra, 'like', dataPrototype);

badSlice = any(isfinite(spectra) == false, 1);
goodSlice = badSlice == false;

if storeParameters == true
    offsetFlat = nan(1, nSpectra);
    exponentFlat = nan(1, nSpectra);
    maeFlat = nan(1, nSpectra);
    rsqFlat = nan(1, nSpectra);
    pkNFlat = nan(1, nSpectra);
else
    offsetFlat = [];
    exponentFlat = [];
    maeFlat = [];
    rsqFlat = [];
    pkNFlat = [];
end

if any(goodSlice) == true
    goodIndices = find(goodSlice);
    maxBlockElements = 4000000;
    blockColumns = max(1, floor(maxBlockElements ./ F));
    predictorMean = mean(predictor);

    for blockStart = 1:blockColumns:numel(goodIndices)
        blockStop = min( ...
            blockStart + blockColumns - 1, ...
            numel(goodIndices));
        blockIndices = goodIndices(blockStart:blockStop);

        blockSpectra = double(spectra(:, blockIndices));
        responseMean = mean(blockSpectra, 1);
        responseCentered = blockSpectra - responseMean;

        predictorResponseProduct = ...
            predictorCentered.' * responseCentered;
        exponent = predictorResponseProduct ./ predictorSumSquares;
        offset = responseMean - exponent .* predictorMean;

        apFit = predictor .* exponent + offset;

        needResidual = storeParameters == true || returnFit == false;

        if needResidual == true
            residual = blockSpectra - apFit;
        else
            residual = [];
        end

        if returnFit == true
            apFlat(:, blockIndices) = apFit;
            oscFlat(:, blockIndices) = 0;
        else
            apFlat(:, blockIndices) = blockSpectra;
            oscFlat(:, blockIndices) = residual;
        end

        if storeParameters == true
            offsetFlat(blockIndices) = offset;
            exponentFlat(blockIndices) = exponent;
            maeFlat(blockIndices) = mean(abs(residual), 1);
            pkNFlat(blockIndices) = 0;

            observedSumSquares = sum( ...
                responseCentered .* responseCentered, 1);
            fittedSumSquares = ...
                exponent .* exponent .* predictorSumSquares;
            denominator = observedSumSquares .* fittedSumSquares;
            numerator = exponent .* predictorResponseProduct;

            validRSquare = denominator > 0;
            rSquare = nan(1, numel(blockIndices));
            rSquare(validRSquare) = ...
                numerator(validRSquare) ...
                ./ sqrt(denominator(validRSquare));
            rSquare(validRSquare) = ...
                rSquare(validRSquare) .* rSquare(validRSquare);

            rsqFlat(blockIndices) = rSquare;
        end
    end
end

end

function fitOptions = nf_compact_options_local(opt)

fitOptions = struct();
fitOptions.peakWidthLims = opt.peakWidthLims;
fitOptions.maxPeaks = opt.maxPeaks;
fitOptions.minPeakHeight = opt.minPeakHeight;
fitOptions.aPeriodicMode = opt.aPeriodicMode;
fitOptions.peakThreshold = opt.peakThreshold;
fitOptions.peakType = opt.peakType;
fitOptions.threshAfter = opt.threshAfter;
fitOptions.optim = opt.optim;
fitOptions.ag = [];
fitOptions.plt = 0;

end

function surf = nf_canonical_power_local( ...
    surf, nFreqs, nSensors, nSeries)

nDimensions = ndims(surf);

if nDimensions == 2

    isSingleSensorTime = ...
        nSensors == 1 && ...
        size(surf, 1) == nFreqs && ...
        nSeries == 1;

    isSensorSpectrum = ...
        size(surf, 1) == nSensors && ...
        size(surf, 2) == nFreqs && ...
        nSeries == 1;

    if isSingleSensorTime == true
        surf = reshape( ...
            surf, ...
            [1 nFreqs size(surf, 2) 1]);
    elseif isSensorSpectrum == true
        surf = reshape( ...
            surf, ...
            [nSensors nFreqs 1 1]);
    else
        error('nf_tfspecparam:BadPowerFreqLayout', ...
            ['A 2D power field must be F×T for one sensor or ' ...
            'C×F for one time point.']);
    end

elseif nDimensions == 3

    isSingleSensorTrials = ...
        nSensors == 1 && ...
        size(surf, 1) == nFreqs && ...
        size(surf, 3) == nSeries;

    isSensorTime = ...
        size(surf, 1) == nSensors && ...
        size(surf, 2) == nFreqs && ...
        nSeries == 1;

    if isSingleSensorTrials == true
        surf = reshape( ...
            surf, ...
            [1 nFreqs size(surf, 2) nSeries]);
    elseif isSensorTime == true
        surf = reshape( ...
            surf, ...
            [nSensors nFreqs size(surf, 3) 1]);
    else
        error('nf_tfspecparam:BadPowerFreqLayout', ...
            ['A 3D power field must be F×T×R for one sensor or ' ...
            'C×F×T for one trial or condition.']);
    end

elseif nDimensions == 4

    if size(surf, 1) ~= nSensors
        error('nf_tfspecparam:BadPowerSensorDim', ...
            ['First power dimension (%d) differs from TF.nsensor ' ...
            '(%d).'], size(surf, 1), nSensors);
    end

    if size(surf, 2) ~= nFreqs
        error('nf_tfspecparam:BadPowerSecondDim', ...
            ['Second power dimension (%d) differs from numel(TF.freqs) ' ...
            '(%d).'], size(surf, 2), nFreqs);
    end

    if size(surf, 4) ~= nSeries
        error('nf_tfspecparam:BadPowerSeriesDim', ...
            ['Fourth power dimension (%d) differs from the expected ' ...
            'trial or condition count (%d).'], size(surf, 4), nSeries);
    end

else
    error('nf_tfspecparam:BadPowerFieldDim', ...
        'Power field must have two, three, or four dimensions.');
end

end

function tf = nf_valid_max_peaks_local(value)

if isempty(value) == true
    tf = true;
    return
end

tf = false;

if isnumeric(value) == false
    return
end

if isreal(value) == false
    return
end

if isscalar(value) == false
    return
end

if isnan(value) == true || isinf(value) == true
    return
end

if value < 0
    return
end

if value ~= fix(value)
    return
end

tf = true;

end

function TF = nf_remove_trials_local(TF, keepIdx)

% Use the known series axis for each schema. A time or frequency axis
% can have exactly the same length as the trial count.
nOld = numel(keepIdx);
if isfield(TF,'epoch') || isfield(TF,'event')
    if ~isfield(TF,'times') || isempty(TF.times)
        error('nf_tfspecparam:NoEventGrid','TF.times is required to remap rejected-trial events.');
    end
    TF = hnf_event_grid(TF, TF.times, TF.times, find(keepIdx), []);
end
names = fieldnames(TF);
for i = 1:numel(names)
    name = names{i};
    val = TF.(name);
    if ismember(name,{'epoch','event'}), continue; end
    if strcmp(name,'behavior') || strcmp(name,'trialinfo') || ...
            startsWith(name,'trial_') || strcmp(name,'trlerp')
        if istable(val) && height(val)==nOld
            TF.(name) = val(keepIdx,:);
        elseif isvector(val) && numel(val)==nOld
            TF.(name) = val(keepIdx);
        elseif ismatrix(val) && size(val,1)==nOld
            TF.(name) = val(keepIdx,:);
        end
    elseif isnumeric(val) || islogical(val)
        if ndims(val)>=3 && size(val,ndims(val))==nOld
            TF.(name) = last_axis(val,keepIdx);
        end
    end
end
sp = TF.specparam;
metrics = {'ap_offset','ap_knee','ap_exponent','pk_n','mae','rsquare', ...
    'pk_f','pk_sd','pk_aLog10'};
for i = 1:numel(metrics)
    name = metrics{i};
    sp.(name) = last_axis(sp.(name),keepIdx);
end
TF.specparam = sp;
if isfield(TF,'conds')
    TF.conds = nnz(keepIdx);
    if isfield(TF,'normbase') && isfield(TF.normbase,'conditions')
        TF.normbase.conditions = subset_condition_fields(TF.normbase.conditions, keepIdx);
    end
    if isfield(TF,'cpm')
        % Numeric moment/validity maps use the same last condition axis as
        % power. Recompute summed diagnostics from the retained conditions.
        TF.cpm = subset_cpm(TF.cpm,keepIdx);
    end
else
    TF.ntrls = nnz(keepIdx);
end
end

function out = last_axis(val,keepIdx)
index = repmat({':'},1,ndims(val));
index{end} = keepIdx;
out = val(index{:});
end

function s = subset_condition_fields(s,keepIdx)
for name = {'values','trialCounts'}
    k = name{1};
    if isfield(s,k), s.(k) = s.(k)(keepIdx); end
end
end

function cpm = subset_cpm(cpm,keepIdx)
if isfield(cpm,'conditions')
    cpm.conditions = subset_condition_fields(cpm.conditions,keepIdx);
end
groups = {'total','phaseLocked','nonPhaseLocked','raw','validity'};
for k = 1:numel(groups)
    g = groups{k};
    if ~isfield(cpm,g), continue; end
    names = fieldnames(cpm.(g));
    for j = 1:numel(names)
        name = names{j};
        x = cpm.(g).(name);
        if (isnumeric(x) || islogical(x)) && size(x,ndims(x))==numel(keepIdx)
            cpm.(g).(name) = last_axis(x,keepIdx);
        end
    end
end
if isfield(cpm,'leakageRatio'), cpm.leakageRatio = last_axis(cpm.leakageRatio,keepIdx); end
if isfield(cpm,'validity') && isfield(cpm.validity,'conditionSummaries')
    summaries = cpm.validity.conditionSummaries(keepIdx);
    cpm.validity.conditionSummaries = summaries;
    cpm.validity.summary = sum_condition_counts(summaries);
end
if isfield(cpm,'baselineDiagnosticsByCondition')
    summaries = cpm.baselineDiagnosticsByCondition(keepIdx);
    cpm.baselineDiagnosticsByCondition = summaries;
    cpm.baselineDiagnostics = sum_condition_counts(summaries);
end
end

function summary = sum_condition_counts(items)
summary = struct();
names = fieldnames(items);
for i = 1:numel(names)
    name = names{i};
    summary.(name) = sum([items.(name)]);
end
end

function tf = nf_is_param_power_field_local(name)

tf = 0;

if endsWith(lower(name), '_ap')
    tf = 1;
end

if endsWith(lower(name), '_osc')
    tf = 1;
end

end

function powFields = nf_find_power_fields_local(TF)

fn = fieldnames(TF);

powFields = {};

for i = 1:numel(fn)

    name = fn{i};
    val = TF.(name);

    if isnumeric(val)
        if contains(lower(name), 'power')
            powFields{end + 1, 1} = name; %#ok
        end
    end

end

idx = find(strcmpi(powFields, 'power'));

if ~isempty(idx)
    powFields(idx) = [];
    powFields = [{'power'}; powFields];
end

end

function TF = nf_powerfront_local(TF)

fn = fieldnames(TF);

pow = {};
other = {};

for i = 1:numel(fn)

    name = fn{i};
    val = TF.(name);

    if isnumeric(val)
        if contains(lower(name), 'power')
            pow{end + 1, 1} = name; %#ok
        else
            other{end + 1, 1} = name; %#ok
        end
    else
        other{end + 1, 1} = name; %#ok
    end

end

idx = find(strcmpi(pow, 'power'));

if ~isempty(idx)
    pow(idx) = [];
    pow = [{'power'}; pow];
end

order = [pow; other];

TF = orderfields(TF, order);

end

function [C, F, T, R] = size4d_local(X)

sz = size(X);

C = sz(1);
F = sz(2);

T = 1;
R = 1;

if numel(sz) >= 3
    T = sz(3);
end

if numel(sz) >= 4
    R = sz(4);
end

end
