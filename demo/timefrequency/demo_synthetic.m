function demo_synthetic()
% DEMO_SYNTHETIC  Compare transforms on the bundled clean and noisy signals.
dataFile = fullfile(fileparts(mfilename('fullpath')), 'synthetic_dataset.csv');
nf_run(@compareSignals,dataFile);
end

function compareSignals(dataFile)
% The first CSV row contains labels. Use the measured time grid and read
% the two amplitude columns by position, independent of MATLAB release.
% textscan preserves the supplied numeric values and accepts spaces in headers.
fid = fopen(dataFile,'r');
if fid < 0, error('NeuroFreq:DemoInput','Cannot open %s.',dataFile); end
closeFile = onCleanup(@() fclose(fid));
columns = textscan(fid,'%f%f%f','Delimiter',',','HeaderLines',1,'CollectOutput',true);
clear closeFile
values = columns{1};
t = values(:,1).';
Fs = 1/median(diff(t));
data = nf_prepdata(values(:,2:3).');
foi = 0.25:0.25:150;
methods = {'STFT, 0.5 s','STFT, 5 s','Filter-Hilbert','Complex demodulation', ...
    'CWT','Morlet DCWT','Stockwell'};
linearFig = figure('Name','Linear transforms, clean and noisy signals');
for ch = 1:2
    signal = data(ch,:);
    set(groot,'CurrentFigure',linearFig); subplot(8,4,2*ch-1:2*ch); plot(t,real(signal));
    if ch == 1, title('Clean signal'); else, title('Noise-corrupted signal'); end
    for k = 1:numel(methods)
        switch k
            case 1, TF = nf_stft(signal,Fs,0.5,80,0.1);
            case 2, TF = nf_stft(signal,Fs,5,95,0.1);
            case 3, TF = nf_filterhilbert(signal,Fs,foi(2:end-2),0.5);
            case 4, TF = nf_demodulation(signal,Fs,foi);
            case 5
                if ~exist('OCTAVE_VERSION','builtin') && isempty(which('cwt')), fprintf('Skipping CWT: Wavelet Toolbox required.\n'); continue; end
                TF = nf_cwt(signal,Fs);
            case 6, TF = nf_dcwt(signal,Fs,foi);
            case 7, TF = nf_stransform(signal,Fs);
        end
        showTF(linearFig,8,k+1,ch,TF.times+t(1),TF.freqs,TF.power,methods{k},[t(1) t(end)],[15 20]);
        clear TF
    end
end

% Quadratic transforms form time-by-lag arrays. Use the same five-second
% excerpt in all four comparisons so the demo can run on a laptop.
use = t>=15 & t<=20; clipTime=t(use);
names = {'Wigner-Ville','Binomial RID','Born-Jordan RID','RID-Rihaczek'};
quadraticFig=figure('Name','Quadratic distributions, clean and noisy signals');
for ch=1:2
    clip = data(ch,use);
    set(groot,'CurrentFigure',quadraticFig); subplot(5,4,2*ch-1:2*ch); plot(clipTime,real(clip));
    if ch == 1, title('Clean signal'); else, title('Noise-corrupted signal'); end
    for k=1:numel(names)
        switch k
            case 1
                if exist('OCTAVE_VERSION','builtin')
                    [power,freq,time] = hnf_wvd(clip,Fs);
                elseif ~isempty(which('wvd'))
                    [power,freq,time] = wvd(clip,Fs);
                else
                    fprintf('Skipping Wigner-Ville: wvd required.\n'); continue
                end
            case 2, TF=nf_ridbinomial(clip,Fs,0.5,0.25);
            case 3, TF=nf_ridbornjordan(clip,Fs,0.5,0.25);
            case 4, TF=nf_ridrihaczek(clip,Fs,0.5);
        end
        if k>1, power=TF.power; freq=TF.freqs; time=TF.times; end
        showTF(quadraticFig,5,k+1,ch,time+clipTime(1),freq,power,names{k},[15 20],[15 20]);
        clear TF power
    end
end
end

function showTF(fig,nRows,row,ch,time,freq,power,label,fullWindow,zoomWindow)
% Select bins once on the original axes before indexing the surface. Draw
% the actual frequency coordinates, including the nonuniform CWT grid.
set(groot,'CurrentFigure',fig);
for zoom = 0:1
    if zoom, window = zoomWindow; maxHz=3; else, window=fullWindow; maxHz=150; end
    ti = find(time>=window(1) & time<=window(2)); fi=find(freq>=0 & freq<=maxHz);
    if numel(ti)<2 || numel(fi)<2, continue; end
    % Limit rendered vertices only. The calculations use every sample.
    ti=ti(round(linspace(1,numel(ti),min(numel(ti),640))));
    fi=fi(round(linspace(1,numel(fi),min(numel(fi),300))));
    subplot(nRows,4,(row-1)*4+2*ch-1+zoom);
    surface(time(ti),freq(fi),power(fi,ti),'EdgeColor','none'); view(2);
    axis tight; title(label); xlabel('Time (s)'); ylabel('Hz');
end
end
