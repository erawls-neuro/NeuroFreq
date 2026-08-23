%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%sgTF DEMO: time-frequency transforms with different methods%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
data = readtable('synthetic_dataset.csv');
dataClean = data.CleanSignalAmplitude';
dataNoisy = data.NoiseSignalAmplitude';
Fs = 500;
foi=.25:.25:150;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%CALL THE tf_fun DIRECTLY%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataClean = nf_prepdata(dataClean);
dataNoisy = nf_prepdata(dataNoisy);

%First Class: Linear Methods (clean)
TF1 = nf_stft(dataClean,Fs,.5,80,0.1);                 %brief t-spectrogram
TF1.freqs = TF1.freqs(TF1.freqs<=150);
TF1.power = TF1.power(TF1.freqs<=150,:);
TF2 = nf_stft(dataClean,Fs,5,95,0.1);                  %long t-spectrogram
TF2.freqs = TF2.freqs(TF2.freqs<=150);
TF2.power = TF2.power(TF2.freqs<=150,:);
TF3 = nf_filterhilbert(dataClean,Fs,foi(2:end-2),.5);  %filter-Hilbert
TF4 = nf_demodulation(dataClean,Fs,foi);               %complex demodulation
TF5 = nf_cwt(dataClean,Fs);                            %Continuous Wavelet
TF5.freqs = TF5.freqs(TF5.freqs<=150);
TF5.power = TF5.power(TF5.freqs<=150,:);
TF6 = nf_dcwt(dataClean,Fs,foi);                      %Morlet dcwt:
TF7 = nf_stransform(dataClean,Fs);                     %Stockwell Transform
TF7.freqs = TF7.freqs(TF7.freqs<=150);
TF7.power = TF7.power(TF7.freqs<=150,:);
%First Class: Linear Methods (corrupted)
TF8 = nf_stft(dataNoisy,Fs,0.5,80,0.1);                %brief t-spectrogram
TF8.freqs = TF8.freqs(TF8.freqs<=150);
TF8.power = TF8.power(TF8.freqs<=150,:);
TF9 = nf_stft(dataNoisy,Fs,5,95,0.1);                  %long t-spectrogram
TF9.freqs = TF9.freqs(TF9.freqs<=150);
TF9.power = TF9.power(TF9.freqs<=150,:);
TF10 = nf_filterhilbert(dataNoisy,Fs,foi(2:end-2),.5); %filter-Hilbert
TF11 = nf_demodulation(dataNoisy,Fs,foi);              %complex demodulation
TF12 = nf_cwt(dataNoisy,Fs);                           %Continuous Wavelet
TF12.freqs = TF12.freqs(TF12.freqs<=150);
TF12.power = TF12.power(TF12.freqs<=150,:);
TF13 = nf_dcwt(dataNoisy,Fs,foi);                   %Morlet dcwt:
TF14 = nf_stransform(dataNoisy,Fs);                    %Stockwell Transform
TF14.freqs = TF14.freqs(TF14.freqs<=150);
TF14.power = TF14.power(TF14.freqs<=150,:);

%plot the outputs
figure;
subplot(8,6,[1,2]); plot(TF3.times,dataClean); xlim('tight'); title('Synthetic Data (Clean)');
subplot(8,6,[4,5]); plot(TF3.times,dataNoisy); xlim('tight'); title('Synthetic Data (Noise Corrupted)');
subplot(8,6,[7,8]); imagesc(TF1.times,TF1.freqs,TF1.power); caxis([-abs(prctile(TF1.power(:),98)) abs(prctile(TF1.power(:),98))]); set(gca,'YDir','normal'); xlim([0 21]); title('STFT (brief window)');
subplot(8,6,[10,11]); imagesc(TF8.times,TF8.freqs,TF8.power); caxis([-abs(prctile(TF8.power(:),98)) abs(prctile(TF8.power(:),98))]); set(gca,'YDir','normal'); xlim([0 21]); title('STFT (brief window)');
subplot(8,6,[13,14]); imagesc(TF2.times,TF2.freqs,TF2.power); caxis([-abs(prctile(TF2.power(:),98)) abs(prctile(TF2.power(:),98))]); set(gca,'YDir','normal'); xlim([0 21]); title('STFT (long window)');
subplot(8,6,[16,17]); imagesc(TF9.times,TF9.freqs,TF9.power); caxis([-abs(prctile(TF9.power(:),98)) abs(prctile(TF9.power(:),98))]); set(gca,'YDir','normal'); xlim([0 21]); title('STFT (long window)');
subplot(8,6,[19,20]); imagesc(TF3.times,TF3.freqs,TF3.power); caxis([-abs(prctile(TF3.power(:),98)) abs(prctile(TF3.power(:),98))]); set(gca,'YDir','normal'); title('Filter-Hilbert');
subplot(8,6,[22,23]); imagesc(TF10.times,TF10.freqs,TF10.power); caxis([-abs(prctile(TF10.power(:),98)) abs(prctile(TF10.power(:),98))]); set(gca,'YDir','normal'); title('Filter-Hilbert');
subplot(8,6,[25,26]); imagesc(TF4.times,TF4.freqs,TF4.power); caxis([-abs(prctile(TF4.power(:),98)) abs(prctile(TF4.power(:),98))]); set(gca,'YDir','normal'); title('Complex Demodulation');
subplot(8,6,[28,29]); imagesc(TF11.times,TF11.freqs,TF11.power); caxis([-abs(prctile(TF11.power(:),98)) abs(prctile(TF11.power(:),98))]); set(gca,'YDir','normal'); title('Complex Demodulation');
subplot(8,6,[31,32]); imagesc(TF5.times,TF5.freqs,TF5.power); caxis([-abs(prctile(TF5.power(:),98)) abs(prctile(TF5.power(:),98))]); set(gca,'YDir','normal','YScale','log'); yticks([1,10,100]); title('CWT');
subplot(8,6,[34,35]); imagesc(TF12.times,TF12.freqs,TF12.power); caxis([-abs(prctile(TF12.power(:),98)) abs(prctile(TF12.power(:),98))]); set(gca,'YDir','normal','YScale','log'); yticks([1,10,100]); title('CWT');
subplot(8,6,[37,38]); imagesc(TF6.times,TF6.freqs,TF6.power); caxis([-abs(prctile(TF6.power(:),98)) abs(prctile(TF6.power(:),98))]); set(gca,'YDir','normal'); title('DCWT');
subplot(8,6,[40,41]); imagesc(TF13.times,TF13.freqs,TF13.power); caxis([-abs(prctile(TF13.power(:),98)) abs(prctile(TF13.power(:),98))]); set(gca,'YDir','normal'); title('DCWT');
subplot(8,6,[43,44]); imagesc(TF7.times,TF7.freqs,TF7.power); caxis([-abs(prctile(TF7.power(:),98)) abs(prctile(TF7.power(:),98))]); set(gca,'YDir','normal'); title('S-transform');
subplot(8,6,[46,47]); imagesc(TF14.times,TF14.freqs,TF14.power); caxis([-abs(prctile(TF14.power(:),98)) abs(prctile(TF14.power(:),98))]); set(gca,'YDir','normal'); title('S-transform');
%zoomed in (all but CWT)
subplot(8,6,9); imagesc(TF1.times(149:198),TF1.freqs(1:31),TF1.power(1:31,149:198)); caxis([-abs(prctile(TF1.power(:),98)) abs(prctile(TF1.power(:),98))]); set(gca,'YDir','normal');
subplot(8,6,12); imagesc(TF8.times(149:198),TF8.freqs(1:31),TF8.power(1:31,149:198)); caxis([-abs(prctile(TF8.power(:),98)) abs(prctile(TF8.power(:),98))]); set(gca,'YDir','normal');
subplot(8,6,15); imagesc(TF2.times(51:65),TF2.freqs(1:31),TF2.power(1:31,51:65)); caxis([-abs(prctile(TF2.power(:),98)) abs(prctile(TF2.power(:),98))]); set(gca,'YDir','normal');
subplot(8,6,18); imagesc(TF9.times(51:65),TF9.freqs(1:31),TF9.power(1:31,51:65)); caxis([-abs(prctile(TF9.power(:),98)) abs(prctile(TF9.power(:),98))]); set(gca,'YDir','normal');
subplot(8,6,21); imagesc(TF3.times(7501:10001),TF3.freqs(1:11),TF3.power(1:11,7501:10001)); caxis([-abs(prctile(TF3.power(:),98)) abs(prctile(TF3.power(:),98))]); set(gca,'YDir','normal'); 
subplot(8,6,24); imagesc(TF10.times(7501:10001),TF10.freqs(1:11),TF10.power(1:11,7501:10001)); caxis([-abs(prctile(TF10.power(:),98)) abs(prctile(TF10.power(:),98))]); set(gca,'YDir','normal');
subplot(8,6,27); imagesc(TF4.times(7501:10001),TF4.freqs(1:12),TF4.power(1:12,7501:10001)); caxis([-abs(prctile(TF4.power(:),98)) abs(prctile(TF4.power(:),98))]); set(gca,'YDir','normal');
subplot(8,6,30); imagesc(TF11.times(7501:10001),TF11.freqs(1:12),TF11.power(1:12,7501:10001)); caxis([-abs(prctile(TF11.power(:),98)) abs(prctile(TF11.power(:),98))]); set(gca,'YDir','normal');
subplot(8,6,39); imagesc(TF6.times(7501:10001),TF6.freqs(1:12),TF6.power(1:12,7501:10001)); caxis([-abs(prctile(TF6.power(:),98)) abs(prctile(TF6.power(:),98))]); set(gca,'YDir','normal'); 
subplot(8,6,42); imagesc(TF13.times(7501:10001),TF13.freqs(1:12),TF13.power(1:12,7501:10001)); caxis([-abs(prctile(TF13.power(:),98)) abs(prctile(TF13.power(:),98))]); set(gca,'YDir','normal'); 
subplot(8,6,45); imagesc(TF7.times(7501:10001),TF7.freqs(1:64),TF7.power(1:64,7501:10001)); caxis([-abs(prctile(TF7.power(:),98)) abs(prctile(TF7.power(:),98))]); set(gca,'YDir','normal'); 
subplot(8,6,48); imagesc(TF14.times(7501:10001),TF14.freqs(1:64),TF14.power(1:64,7501:10001)); caxis([-abs(prctile(TF14.power(:),98)) abs(prctile(TF14.power(:),98))]); set(gca,'YDir','normal');
colormap(flipud(brewermap([],'PuOr')))










%Second Class: Quadratic Distributions (clean)
TF15 = nf_ridbinomial(dataClean,Fs);               %Type-II Binomial RID
TF15.freqs = TF15.freqs(TF15.freqs<=150);
TF15.power = TF15.power(TF15.freqs<=150,:);
TF16 = nf_ridbornjordan(dataClean,Fs);             %Type-II Born-Jordan RID
TF16.freqs = TF16.freqs(TF16.freqs<=150);
TF16.power = TF16.power(TF16.freqs<=150,:);
TF17 = nf_ridrihaczek(dataClean,Fs);               %RID-Rihaczek TFD
TF17.freqs = TF17.freqs(TF17.freqs<=150);
TF17.power = TF17.power(TF17.freqs<=150,:);
%Second Class: Quadratic Distributions (corrupted)
TF18 = nf_ridbinomial(dataNoisy,Fs);               %Type-II Binomial RID
TF18.freqs = TF18.freqs(TF18.freqs<=150);
TF18.power = TF18.power(TF18.freqs<=150,:);
TF19 = nf_ridbornjordan(dataNoisy,Fs);             %Type-II Born-Jordan RID
TF19.freqs = TF19.freqs(TF19.freqs<=150);
TF19.power = TF19.power(TF19.freqs<=150,:);
TF20 = nf_ridrihaczek(dataNoisy,Fs);               %RID-Rihaczek TFD
TF20.freqs = TF20.freqs(TF20.freqs<=150);
TF20.power = TF20.power(TF20.freqs<=150,:);
%calculate an additional wigner-ville transform (not included)
[d1,f1,t1]=wvd(dataClean, Fs); %clean
d1(f1>150,:)=[];
f1(f1>150)=[];
[d2,f2,t2]=wvd(dataNoisy, Fs); %corrupted
d2(f2>150,:)=[];
f2(f2>150)=[];

%plot the quadratic outputs
figure;
subplot(5,6,[1,2]); plot(TF15.times,dataClean); xlim('tight'); title('Synthetic Data (Clean)');
subplot(5,6,[4,5]); plot(TF15.times,dataNoisy); xlim('tight'); title('Synthetic Data (Noise Corrupted)');
subplot(5,6,[7,8]); imagesc(t1,f1,d1); caxis([prctile(d1(:),2), prctile(d1(:),98)]); set(gca,'YDir','normal'); title('Wigner-Ville Distribution');
subplot(5,6,[10,11]); imagesc(t2,f2,d2); caxis([prctile(d2(:),2), prctile(d2(:),98)]); set(gca,'YDir','normal'); title('Wigner-Ville Distribution');
subplot(5,6,[13,14]); imagesc(TF15.times,TF15.freqs,TF15.power); caxis([-abs(prctile(TF15.power(:),98)) abs(prctile(TF15.power(:),98))]); set(gca,'YDir','normal'); title('Binomial RID');
subplot(5,6,[16,17]); imagesc(TF18.times,TF18.freqs,TF18.power); caxis([-abs(prctile(TF18.power(:),98)) abs(prctile(TF18.power(:),98))]); set(gca,'YDir','normal'); title('Binomial RID');
subplot(5,6,[19,20]); imagesc(TF16.times,TF16.freqs,TF16.power); caxis([-abs(prctile(TF16.power(:),98)) abs(prctile(TF16.power(:),98))]); set(gca,'YDir','normal'); title('Born-Jordan RID');
subplot(5,6,[22,23]); imagesc(TF19.times,TF19.freqs,TF19.power); caxis([-abs(prctile(TF19.power(:),98)) abs(prctile(TF19.power(:),98))]); set(gca,'YDir','normal'); title('Born-Jordan RID');
subplot(5,6,[25,26]); imagesc(TF17.times,TF17.freqs,TF17.power); caxis([-abs(prctile(TF17.power(:),98)) abs(prctile(TF17.power(:),98))]); set(gca,'YDir','normal'); title('RID-Rihaczek');
subplot(5,6,[28,29]); imagesc(TF20.times,TF20.freqs,TF20.power); caxis([-abs(prctile(TF20.power(:),98)) abs(prctile(TF20.power(:),98))]); set(gca,'YDir','normal'); title('RID-Rihaczek');
%zoomed in
subplot(5,6,9); imagesc(t1(15001:20001),f1(1:127),d1(1:127,15001:20001)); caxis([prctile(d1(:),2), prctile(d1(:),98)]); set(gca,'YDir','normal');
subplot(5,6,12); imagesc(t2(15001:20001),f2(1:127),d2(1:127,15001:20001)); caxis([prctile(d2(:),2), prctile(d2(:),98)]); set(gca,'YDir','normal');
subplot(5,6,15); imagesc(TF15.times(7501:10001),TF15.freqs(1:127),TF15.power(1:127,7501:10001)); caxis([-abs(prctile(TF15.power(:),98)) abs(prctile(TF15.power(:),98))]); set(gca,'YDir','normal');
subplot(5,6,18); imagesc(TF18.times(7501:10001),TF18.freqs(1:127),TF18.power(1:127,7501:10001)); caxis([-abs(prctile(TF18.power(:),98)) abs(prctile(TF18.power(:),98))]); set(gca,'YDir','normal');
subplot(5,6,21); imagesc(TF16.times(7501:10001),TF16.freqs(1:127),TF16.power(1:127,7501:10001)); caxis([-abs(prctile(TF16.power(:),98)) abs(prctile(TF16.power(:),98))]); set(gca,'YDir','normal');
subplot(5,6,24); imagesc(TF19.times(7501:10001),TF19.freqs(1:127),TF19.power(1:127,7501:10001)); caxis([-abs(prctile(TF19.power(:),98)) abs(prctile(TF19.power(:),98))]); set(gca,'YDir','normal');
subplot(5,6,27); imagesc(TF17.times(7501:10001),TF17.freqs(1:127),TF17.power(1:127,7501:10001)); caxis([-abs(prctile(TF17.power(:),98)) abs(prctile(TF17.power(:),98))]); set(gca,'YDir','normal'); 
subplot(5,6,30); imagesc(TF20.times(7501:10001),TF20.freqs(1:127),TF20.power(1:127,7501:10001)); caxis([-abs(prctile(TF20.power(:),98)) abs(prctile(TF20.power(:),98))]); set(gca,'YDir','normal');
%colormap
colormap(flipud(brewermap([],'PuOr')));








