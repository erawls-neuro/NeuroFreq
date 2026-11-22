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
    disp('[nf_icaclean2eeglab]: select input filepath:');
    inpath = uigetdir(pwd, 'Select ICAClean directory');
    if isequal(inpath, 0), EEG = []; return; end
end
inecnt = [dir(fullfile(inpath,'*ecnt_af.mat')); dir(fullfile(inpath,'*','*ecnt_af.mat'))];
if isempty(inecnt)
    error('nf_icaclean2eeglab:NoInput','No ecnt_af files found in %s or its subject folders.',inpath);
end
%default saveset=0
if nargin<3 || isempty(saveset)
    disp('[nf_icaclean2eeglab]: not saving output');
    saveset = 0;
end
%channel lookup if not supplied
if nargin<2 || isempty(montage)
    disp('[nf_icaclean2eeglab]: no montage supplied - provide channel locations file');
    [file, folder] = uigetfile('*.*','Select channel locations');
    if isequal(file,0), EEG = []; return; end
    montage = readlocs(fullfile(folder,file));
elseif ischar(montage) || (isstring(montage) && isscalar(montage))
    montage = readlocs(char(montage));
elseif isstruct(montage)
    disp('[nf_icaclean2eeglab]: using supplied channel locations');
else
    error('nf_icaclean2eeglab:NoMontage',...
        'montage not found - please provide a valid montage file');
end

%enumerate over datasets
for i=1:length(inecnt)
    
    %get one file
    ecnt=load([inecnt(i).folder filesep inecnt(i).name]);
    ecnt=ecnt.ecnt_af;

    %check channel size
    if size(ecnt.data,1)~=length(montage)
        error('nf_icaclean2eeglab:BadChanlocs',...
            'different channel locations between ecnt and montage?');
    end
    
    %error if all epochs are not same length
    if numel(unique(ecnt.ntbins))>1
        error('nf_icaclean2eeglab:EpochLengthIssue',...
            'EEGLAB cannot handle epochs of different lengths');
    end

    %get data and events
    data = ecnt.data(ecnt.eegElecs,:);
    subjectMontage = montage(ecnt.eegElecs);
    event = reshape(ecnt.event,1,[]);
    srate = ecnt.samplerate;

    %get sizes
    nChans = size(data,1);
    nTrials = numel(ecnt.ntbins);
    nTimes = ecnt.ntbins(1);

    % Leading edges are measured separately within each epoch. ICAClean
    % stores the actual trigger code at the first active sample. Preserve
    % that code, including a trigger on the first sample of an epoch.
    data = reshape(data,[nChans nTimes nTrials]);
    EEG = pop_importdata('dataformat','matlab','nbchan',nChans,'data', ...
        data,'srate',srate);
    triggers = reshape(event,nTimes,nTrials);
    rises = abs(triggers) > [zeros(1,nTrials); abs(triggers(1:end-1,:))];
    onset = find(rises);
    EEG.event = struct('type',num2cell(triggers(onset)), ...
        'latency',num2cell(onset),'epoch',num2cell(ceil(onset/nTimes)));
    EEG = eeg_checkset(EEG,'makeur');
    EEG = eeg_checkset(EEG,'eventconsistency');
    EEG.chanlocs = subjectMontage;

    %add sweep to .etc field and finalize
    EEG.etc.behavior = ecnt.sweep;
    EEG = eeg_checkset(EEG);

    %save optional
    if saveset==1
        [~,b]=fileparts(inecnt(i).name);
        EEG = pop_saveset(EEG, 'filename', [b '_eeglab.set'], 'filepath', inecnt(i).folder);
    end

    disp('[nf_icaclean2eeglab]: dataset conversion complete');
    
end







