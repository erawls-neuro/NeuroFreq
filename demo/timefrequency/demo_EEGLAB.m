function demo_EEGLAB()
% DEMO_EEGLAB  Compare power and phase consistency on the bundled EEG epochs.
% Run nf_install, then call this demo with its own folder on the MATLAB path.
dataFile = fullfile(fileparts(mfilename('fullpath')), 'eeglab_sample.set');
nf_run(@compareEEG,dataFile);
end

function compareEEG(dataFile)
if exist('OCTAVE_VERSION','builtin') && (isempty(which('pop_loadset')) || isempty(which('pop_select')))
    % This particular sample .set embeds its numeric data as MAT variables.
    % Keep only the fields the demo passes to nf_prepdata and nf_tftransform.
    EEG = load('-mat',dataFile,'data','pnts','trials','times','xmin','xmax', ...
        'srate','chanlocs','event','epoch','etc');
    keep = ~ismember({EEG.chanlocs.labels},{'EOG1','EOG2'});
    EEG.data = EEG.data(keep,:,:);
    EEG.chanlocs = EEG.chanlocs(keep);
    EEG.nbchan = nnz(keep);
else
    EEG = pop_loadset('filename',dataFile);
    EEG = pop_select(EEG,'nochannel',{'EOG1','EOG2'});
end
EEG = nf_prepdata(EEG);
pz = find(strcmpi({EEG.chanlocs.labels},'Pz'),1);
assert(~isempty(pz),'The bundled EEG montage must contain Pz.');
methods = {'stft','filterhilbert','demodulation','dcwt','cwt','stransform', ...
    'ridbinomial','ridbornjordan','ridrihaczek'};
powerFig = figure('Name','Alpha scalp power and Pz time-frequency surface');
phaseFig = figure('Name','Delta scalp ITPC and Pz phase consistency');
for k = 1:numel(methods)
    if strcmp(methods{k},'cwt') && ~exist('OCTAVE_VERSION','builtin') && isempty(which('cwt'))
        fprintf('Skipping CWT: MATLAB Wavelet Toolbox is required.\n'); continue
    end
    TF = nf_tftransform(EEG,'method',methods{k},'freqs',1:30,'times',[-0.5 1]);
    if startsWith(methods{k},'rid')
        % These are signed quadratic distributions. Average the signed
        % coefficients directly without a floor or logarithmic baseline.
        power = mean(TF.power,4);
    else
        average = nf_normbase(TF,'BaselineMethod','none');
        power = average.power;
    end
    drawPair(powerFig,k,TF,power,pz,TF.freqs >= 8 & TF.freqs <= 12,methods{k});
    if ~startsWith(methods{k},'rid') && isfield(TF,'phase')
        % nf_normbase stores the mean unit-vector magnitude in .itpc.
        % Rihaczek distribution phase is not a linear oscillator phase.
        drawPair(phaseFig,k,TF,average.itpc,pz,TF.freqs > 0 & TF.freqs < 5,methods{k});
    end
end
end

function drawPair(fig,index,TF,measure,pz,band,label)
% Select the destination without forcing a hidden figure onto the desktop.
set(groot,'CurrentFigure',fig);
topo = reshape(mean(mean(measure(:,band,:),2),3),1,[]);
subplot(3,9,3*index-2);
if exist('OCTAVE_VERSION','builtin') && isempty(which('topoplot'))
    % Display the same measured band values at the actual sensor positions.
    locs = TF.chanlocs;
    scatter([locs.radius].*sind([locs.theta]),[locs.radius].*cosd([locs.theta]), ...
        34,topo,'filled'); axis equal off;
else
    topoplot(topo,TF.chanlocs,'numcontour',3);
end
subplot(3,9,[3*index-1 3*index]);
contourf(TF.times,TF.freqs,squeeze(measure(pz,:,:)),40,'linestyle','none');
axis xy; title(label); xlabel('Time (s)'); ylabel('Hz'); colorbar;
end
