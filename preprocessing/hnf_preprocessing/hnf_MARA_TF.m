function features = hnf_MARA_TF(EEG)
% HNF_MARA_TF  Rescue ICs with spectral bumps/plateaus for developmental preprocessing
% features = list of components with bumps
%
% Eric Rawls modified 9/10/2026.
% Recompute activations from the ICA training channels in their original order.
icaChannels = EEG.icachansind(:).';
if isempty(icaChannels)
    icaChannels = 1:size(EEG.data,1);
end
% Flatten samples and epochs after selecting the training channels.
% EEG.nbchan is not necessarily the number of channels used for ICA.
data = reshape(EEG.data(icaChannels,:,:), numel(icaChannels), []);
% resample to 100 Hz with anti aliasing
data = resample(data, 100, EEG.srate, 'Dimension',2);
fs = 100;
%compute ICA activation and standardize variance to 1
icacomps = (EEG.icaweights * EEG.icasphere * data)';
icacomps = icacomps./repmat(std(icacomps,0,1),length(icacomps(:,1)),1);
icacomps = icacomps';
%creating a matrix to help us check for peaks in data
change=[0 linspace(-1,1,50)];
components_with_bumps = [];
components_with_bumps_2 = [];
components_with_plats = [];
components_with_plats_2 = [];
pxx_matrix = zeros(length(icacomps(:,1)),51);
resids_pxx_matrix = zeros(length(icacomps(:,1)),50);
%run through all the components
for ic=1:length(icacomps(:,1))
    % One-second Welch windows and a 1-Hz frequency grid.
    % Requires the actual sample rate fs to be an integer >= 100 Hz.
    [pxx,freq] = pwelch(icacomps(ic,:), fs, 0, fs, fs, 'onesided');
    % Retain 0–50 Hz and log-transform the linear PSD once.
    logpxx = log10(pxx(1:51)).';
    freq = freq(1:51).';
    % Keep the dB spectrum for the subsequent raw-spectrum peak checks.
    % The original factor of 50 adds only a constant; it is unnecessary.
    pxx = 10*logpxx;
    pxx_matrix(ic,1:51) = pxx;
    % Fit 1–50 Hz, excluding DC; suppress plotting.
    fooofd = nf_specparam(logpxx(2:51), freq(2:51));
    % Residual above the aperiodic baseline, including oscillatory peaks.
    r = fooofd.PeriodicData(:);
    % Ordinary-regression leverage for a fixed aperiodic model:
    % log10(power) = offset - exponent*log10(frequency).
    x = log10(fooofd.f(:));
    x = x - mean(x);
    nFreq = numel(r);
    h = 1/nFreq + x.^2 / sum(x.^2);
    % MADE-style residual scaling, retaining peaks in the sum of squares.
    % Two fitted aperiodic parameters: offset and exponent.
    scale2 = sum(r.^2) / (nFreq - 2);
    % Approximate standardized residuals; no additional dB conversion.
    resids_pxx = (r ./ sqrt(scale2 * (1-h))).';
    resids_pxx_matrix(ic,1:50) = resids_pxx;
    % Checking for bumps between 5 and 15 Hz
    % look for peaks with the residual values
    [pks, ~, widths, proms] = findpeaks(resids_pxx(4:16),4:16); %looks for peaks 5-15 Hz (1 less and more so any bumps at 5 or 15 are caught)
    for numPeaks = 1:length(pks) %in case we have mult. peaks
        if (proms(numPeaks) > 0.3 && pks(numPeaks)>1 && widths(numPeaks) > 0.9) %conservative
            components_with_bumps = [components_with_bumps ic]; %#ok <AGROW>
        end
    end
    % look for peaks with the NON-residual values
    [pks, ~, ~, proms] = findpeaks(pxx(5:17),5:17); %looks for peaks 5-15 Hz (1 less and more so any bumps at 5 or 15 are caught)
    for numPeaks = 1:length(pks) %in case we have mult. peaks
        if (proms(numPeaks) > 0.15) %conservative
            components_with_bumps_2 = [components_with_bumps_2 ic]; %#ok <AGROW>
        end
    end
    % Checking for plateaus between 5 and 15 Hz for first 7 components
    if ic < 8
        [pks, ~, ~, proms] = findpeaks(resids_pxx(4:16),4:16);
        for numPeaks = 1:length(pks) %in case we have mult. peaks
            if (proms(numPeaks) > 0.15 && pks(numPeaks)>1) %conservative
                components_with_plats = [components_with_plats ic]; %#ok <AGROW>
            end
        end
        [pks, ~, ~, proms] = findpeaks(pxx(5:17)+(10*change(5:17)),5:17);
        for numPeaks = 1:length(pks) %in case we have mult. peaks
            if (proms(numPeaks) > 0.05) %conservative
                components_with_plats_2 = [components_with_plats_2 ic]; %#ok <AGROW>
            end
        end
    end
end
% initialize these for later
features = []; % this will be our output matrix (our list of components with bumps)
features_2 = [];
% combine the two components with bumps matrices
components_with_bumps = unique(components_with_bumps);
components_with_bumps_2 = unique(components_with_bumps_2);
for check=1:length(components_with_bumps)
    if sum(components_with_bumps(check) == components_with_bumps_2)
        % if they both have a component, add it to the features matrix
        features = [features components_with_bumps(check)]; %#ok <AGROW>
    end
end
% combine the two components with plateaus matrices
components_with_plats = unique(components_with_plats);
components_with_plats_2 = unique(components_with_plats_2);
for check=1:length(components_with_plats)
    if sum(components_with_plats(check) == components_with_plats_2)
        features_2 = [features_2 components_with_plats(check)]; %#ok <AGROW>
    end
end
% combine the components with bumps matrix and the plateau matrix
features = unique([features features_2]); %output matrix (list of components with bumps)
end
