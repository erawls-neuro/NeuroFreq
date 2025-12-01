function EEG = nf_icaclean2eeglab(inpath,montage,saveset)
% GENERAL
% -------
% ICAClean is a software environment running in MATLAB that is used for EEG
% preprocessing and analysis by the Sponheim lab at the Minneapolis VA.
% This function provides an interface to apply NeuroFreq analyses to
% ICAClean-formatted data by converting ICAClean to EEGLAB format. ICAClean
% files typically contain subject-level directories within a data.mat
% directory. The function can be used two ways. First, one can convert a
% single ICAClean file by selecting the subject directory within data.mat,
% or, one can convert multiple ICAClean files by selecting the parent
% data.mat directory. If saveset is selected, files will be saved to the 
% same directory the ecnt file is located.
%
% OUTPUT
% ------
% EEG - EEGLAB .set structure
%
% INPUT
% -----
% inpath - input path for file. this should be the path to either the 
%       subject-level directory or the data.mat directory. if blank, a gui 
%       pops up
% montage - montage (required because ICAClean does not
%       store electrode locations). Can be blank (gui appears),
%       EEGLAB-formatted chanlocs in memory, or a filename in which case
%       eeglab readlocs() function is applied to read the locations.
% saveset - save the formatted set in the same directory? 0 or 1
%
%

%if no ecnt provided, use gui
if nargin<1 || isempty(inpath)
    disp('select input filepath:');
    inpath=uigetdir('*.mat','Select Directory:');
    inecnt=dir([inpath '/*/*ecnt_af.mat']);
    if isempty(inecnt) %no files
        inecnt=dir([inpath '/*/*ecnt_af.mat']); %is it multiple files?
    end   
    if isempty(inecnt) %still no files
        error('no input ecnt_af files found. try a different directory');
    else %got some
        disp(['converting ' num2str(length(inecnt)) ' files to EEGLAB format']);
    end
end
%default saveset=0
if nargin<3 || isempty(saveset)
    disp('not saving output');
    saveset = 0;
end
%channel lookup if not supplied
if nargin<2 || isempty(montage)
    disp('no montage supplied - provide channel locations file');
    montage = pop_readlocs();
elseif exist('montage','var')
    disp('using variable in memory');
elseif isstring(montage)
    disp('montage supplied as filename - importing');
    montage = readlocs(montage);
else
    error('montage not found - please provide a valid montage file');
end

%enumerate over datasets
for i=1:length(inecnt)
    
    %get one file
    ecnt=load([inecnt(i).folder filesep inecnt(i).name]);
    ecnt=ecnt.ecnt_af;

    %check channel size
    if size(ecnt.data,1)~=length(montage)
        error('different channel locations between ecnt and montage?');
    end
    
    %error if all epochs are not same length
    if numel(unique(ecnt.ntbins))>1
        error('EEGLAB cannot handle epochs of different lengths');
    end

    %get data and events
    data = ecnt.data(ecnt.eegElecs,:);
    montage = montage(ecnt.eegElecs);
    event = ecnt.event;
    srate = ecnt.samplerate;

    %get sizes
    nChans = numel(ecnt.eegElecs);
    nTrials = numel(ecnt.ntbins);
    nTimes = ecnt.ntbins(1);

    %concatenate (last channel is events)
    importDat = [data;event];

    %reshape to 3D
    importDat = reshape(importDat,[nChans+1 nTimes nTrials]);

    %import
    EEG = pop_importdata('dataformat','matlab','nbchan',nChans+1,'data',...
        importDat,'srate',srate);
    EEG = pop_chanevent(EEG,size(importDat,1),'edge','leading','edgelen',0);
    EEG.chanlocs = montage;

    %add sweep to .etc field and finalize
    EEG.etc.behavior = ecnt.sweep;
    EEG = eeg_checkset(EEG);

    %save optional
    if saveset==1
        [~,b]=fileparts(inecnt(i).name);
        pop_saveset(EEG, [inpath filesep b '_eeglab.set']);
    end
    
end







