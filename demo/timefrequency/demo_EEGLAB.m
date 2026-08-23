%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%NeuroFreq DEMO: time-frequency transforms with different methods%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
EEG = pop_loadset('eeglab_sample.set');
EEG = pop_select( EEG, 'nochannel', {'EOG1','EOG2'});
foi=1:1:30;
toi=[-.5 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%BEGIN BY PREPARING DATA%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG = nf_prepdata( EEG );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%USE tfUtility.m FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First Class: Linear Decompositions
TF1 = nf_tftransform(EEG,'freqs',foi,'times',toi,...
    'method', 'stft');               %matlab spectrogram (STFT)
TF2 = nf_tftransform(EEG,'freqs',foi,'times',toi,...
    'method','filterhilbert');      %filter-Hilbert
TF3 = nf_tftransform(EEG,'freqs',foi,'times',toi,...
    'method', 'demodulation');        %complex demodulation
TF4 = nf_tftransform(EEG,'freqs',foi,'times',toi,...
    'method',  'dcwt');           %Morlet discretized wavelet
TF5 = nf_tftransform(EEG,'freqs',foi,'times',toi,...
    'method',  'cwt');                  %Continuous Wavelet Transform
TF6 = nf_tftransform(EEG,'freqs',foi,'times',toi,...
    'method',  'stransform');           %Stockwell Transform
%Second Class: Quadratic Distributions
TF7 = nf_tftransform(EEG,'freqs',foi,'times',toi,...
    'method',  'ridbinomial');        %Type-II Binomial RID
TF8 = nf_tftransform(EEG,'freqs',foi,'times',toi,...
    'method',  'ridbornjordan');      %Type-II Born-Jordan RID
TF9 = nf_tftransform(EEG,'freqs',foi,'times',toi,...
    'method',  'ridrihaczek');        %RID-Rihaczek TFD



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%NOW AVERAGE THE TF (NO BASELINE CORRECTION)%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TFA1 = nf_avebase(TF1,'none');
TFA2 = nf_avebase(TF2,'none');
TFA3 = nf_avebase(TF3,'none');
TFA4 = nf_avebase(TF4,'none');
TFA5 = nf_avebase(TF5,'none');
TFA6 = nf_avebase(TF6,'none');
TFA7 = nf_avebase(TF7,'none');
TFA8 = nf_avebase(TF8,'none');
TFA9 = nf_avebase(TF9,'none');


%plot the outputs
%topos for alpha
topo1 = squeeze(mean(mean(mean(TFA1.power(:,(TFA1.freqs<12)+(TFA1.freqs>8)==2,:),2),4),3));
topo2 = squeeze(mean(mean(mean(TFA2.power(:,(TFA2.freqs<12)+(TFA2.freqs>8)==2,:),2),4),3));
topo3 = squeeze(mean(mean(mean(TFA3.power(:,(TFA3.freqs<12)+(TFA3.freqs>8)==2,:),2),4),3));
topo4 = squeeze(mean(mean(mean(TFA4.power(:,(TFA4.freqs<12)+(TFA4.freqs>8)==2,:),2),4),3));
topo5 = squeeze(mean(mean(mean(TFA5.power(:,(TFA5.freqs<12)+(TFA5.freqs>8)==2,:),2),4),3));
topo6 = squeeze(mean(mean(mean(TFA6.power(:,(TFA6.freqs<12)+(TFA6.freqs>8)==2,:),2),4),3));
topo7 = squeeze(mean(mean(mean(TFA7.power(:,(TFA7.freqs<12)+(TFA7.freqs>8)==2,:),2),4),3));
topo8 = squeeze(mean(mean(mean(TFA8.power(:,(TFA8.freqs<12)+(TFA8.freqs>8)==2,:),2),4),3));
topo9 = squeeze(mean(mean(mean(TFA9.power(:,(TFA9.freqs<12)+(TFA9.freqs>8)==2,:),2),4),3));
%surfaces for Pz
surf1 = squeeze(TFA1.power(20,:,:));
surf2 = squeeze(TFA2.power(20,:,:));
surf3 = squeeze(TFA3.power(20,:,:));
surf4 = squeeze(TFA4.power(20,:,:));
surf5 = squeeze(TFA5.power(20,:,:));
surf6 = squeeze(TFA6.power(20,:,:));
surf7 = squeeze(TFA7.power(20,:,:));
surf8 = squeeze(TFA8.power(20,:,:));
surf9 = squeeze(TFA9.power(20,:,:));
%make the plots
figure;
subplot(3,9,1); topoplot(topo1, TF1.chanlocs,'numcontour',3,'maplimits',[0 abs(prctile(TFA1.power(:),99))]);
subplot(3,9,[2,3]); contourf(TFA1.times,TFA1.freqs,surf1,100,'linestyle','none'); caxis([0 abs(prctile(TFA1.power(:),99))]); title('STFT');
subplot(3,9,4); topoplot(topo2, TF2.chanlocs,'numcontour',3,'maplimits',[0 abs(prctile(TFA2.power(:),99))]);
subplot(3,9,[5,6]); contourf(TFA2.times,TFA2.freqs,surf2,100,'linestyle','none'); caxis([0 abs(prctile(TFA2.power(:),99))]); title('Filter-Hilbert');
subplot(3,9,7); topoplot(topo3, TF3.chanlocs,'numcontour',3,'maplimits',[0 abs(prctile(TFA3.power(:),99))]);
subplot(3,9,[8,9]); contourf(TFA3.times,TFA3.freqs,surf3,100,'linestyle','none'); caxis([0 abs(prctile(TFA3.power(:),99))]); title('Complex Demodulation');
subplot(3,9,10); topoplot(topo4, TF4.chanlocs,'numcontour',3,'maplimits',[0 abs(prctile(TFA4.power(:),99))]);
subplot(3,9,[11,12]); contourf(TFA4.times,TFA4.freqs,surf4,100,'linestyle','none'); caxis([0 abs(prctile(TFA4.power(:),99))]); title('DCWT');
subplot(3,9,13); topoplot(topo5, TF5.chanlocs,'numcontour',3,'maplimits',[0 abs(prctile(TFA5.power(:),99))]);
subplot(3,9,[14,15]); contourf(TFA5.times,TFA5.freqs,surf5,100,'linestyle','none'); caxis([0 abs(prctile(TFA5.power(:),99))]); title('CWT');
subplot(3,9,16); topoplot(topo6, TF6.chanlocs,'numcontour',3,'maplimits',[0 abs(prctile(TFA6.power(:),99))]);
subplot(3,9,[17,18]); contourf(TFA6.times,TFA6.freqs,surf6,100,'linestyle','none'); caxis([0 abs(prctile(TFA6.power(:),99))]); title('S-transform');
subplot(3,9,19); topoplot(topo7, TF7.chanlocs,'numcontour',3,'maplimits',[0 abs(prctile(TFA7.power(:),99))]);
subplot(3,9,[20,21]); contourf(TFA7.times,TFA7.freqs,surf7,100,'linestyle','none'); caxis([0 abs(prctile(TFA7.power(:),99))]); title('Binomial RID');
subplot(3,9,22); topoplot(topo8, TF8.chanlocs,'numcontour',3,'maplimits',[0 abs(prctile(TFA8.power(:),99))]);
subplot(3,9,[23,24]); contourf(TFA8.times,TFA8.freqs,surf8,100,'linestyle','none'); caxis([0 abs(prctile(TFA8.power(:),99))]); title('Born-Jordan RID');
subplot(3,9,25); topoplot(topo9, TF9.chanlocs,'numcontour',3,'maplimits',[0 abs(prctile(TFA9.power(:),99))]);
subplot(3,9,[26,27]); contourf(TFA9.times,TFA9.freqs,surf9,100,'linestyle','none'); caxis([0 abs(prctile(TFA9.power(:),99))]); title('RID-Rihaczek');
colormap(brewermap([],'OrRd'));










%plot phase-locking outputs
%topos for delta
topo1 = squeeze(mean(mean(mean(TFA1.phase(:,TFA1.freqs<5,:),2),4),3));
topo2 = squeeze(mean(mean(mean(TFA2.phase(:,TFA2.freqs<5,:),2),4),3));
topo3 = squeeze(mean(mean(mean(TFA3.phase(:,TFA3.freqs<5,:),2),4),3));
topo4 = squeeze(mean(mean(mean(TFA4.phase(:,TFA4.freqs<5,:),2),4),3));
topo5 = squeeze(mean(mean(mean(TFA5.phase(:,TFA5.freqs<5,:),2),4),3));
topo6 = squeeze(mean(mean(mean(TFA6.phase(:,TFA6.freqs<5,:),2),4),3));
topo9 = squeeze(mean(mean(mean(TFA9.phase(:,TFA9.freqs<5,:),2),4),3));
%surfaces for Pz
surf1 = squeeze(TFA1.phase(3,:,:));
surf2 = squeeze(TFA2.phase(3,:,:));
surf3 = squeeze(TFA3.phase(3,:,:));
surf4 = squeeze(TFA4.phase(3,:,:));
surf5 = squeeze(TFA5.phase(3,:,:));
surf6 = squeeze(TFA6.phase(3,:,:));
surf9 = squeeze(TFA9.phase(3,:,:));
%make the plots
figure;
subplot(4,6,1); topoplot(topo1, TF1.chanlocs,'numcontour',3,'maplimits',[0 .6]);
subplot(4,6,[2,3]); contourf(TFA1.times,TFA1.freqs,surf1,100,'linestyle','none'); caxis([0 .6]); title('STFT');
subplot(4,6,4); topoplot(topo2, TF2.chanlocs,'numcontour',3,'maplimits',[0 .6]);
subplot(4,6,[5,6]); contourf(TFA2.times,TFA2.freqs,surf2,100,'linestyle','none'); caxis([0 .6]); title('Filter-Hilbert');
subplot(4,6,7); topoplot(topo3, TF3.chanlocs,'numcontour',3,'maplimits',[0 .6]);
subplot(4,6,[8,9]); contourf(TFA3.times,TFA3.freqs,surf3,100,'linestyle','none'); caxis([0 .6]); title('Complex Demodulation');
subplot(4,6,10); topoplot(topo4, TF4.chanlocs,'numcontour',3,'maplimits',[0 .6]);
subplot(4,6,[11,12]); contourf(TFA4.times,TFA4.freqs,surf4,100,'linestyle','none'); caxis([0 .6]); title('DCWT');
subplot(4,6,13); topoplot(topo5, TF5.chanlocs,'numcontour',3,'maplimits',[0 .6]);
subplot(4,6,[14,15]); contourf(TFA5.times,TFA5.freqs,surf5,100,'linestyle','none'); caxis([0 .6]); title('CWT');
subplot(4,6,16); topoplot(topo6, TF6.chanlocs,'numcontour',3,'maplimits',[0 .6]);
subplot(4,6,[17,18]); contourf(TFA6.times,TFA6.freqs,surf6,100,'linestyle','none'); caxis([0 .6]); title('S-transform');
subplot(4,6,19); topoplot(topo9, TF9.chanlocs,'numcontour',3,'maplimits',[0 .6]);
subplot(4,6,[20,21]); contourf(TFA9.times,TFA9.freqs,surf9,100,'linestyle','none'); caxis([0 .6]); title('RID-Rihaczek');
colormap(brewermap([],'OrRd'));


