function blinks = hnf_beall_blink(cwb,topog,chanlocs,n)
% HNF_BEALL_BLINK   Detect blinks for developmental data
%
% Note: this function was rewritten in large part by Eric Rawls 9/2026.
% Check intended use case before applying.
%
% Author: Daniel J. Beall
% Email: DJBEALL1101@gmail.com
%
% Inputs:
%   cwb        - list of components with bumps
%   topog      - topographies vector
%   chanlocs   - EEG.chanlocs struct
%   n          - number of ICs
%   nchannels  - number of channels
% Outputs:
%   blinks     - list of detected blinks
%

blinks = [];
nchannels=length(chanlocs);

% Eric Rawls modified 9/10/2026.
% Standardize each component across all channels; allow non-square ICA.
zmatrix = zeros(n,nchannels);
for i = 1:n
    curr_vect = topog(i,1:nchannels).';
    z_vect = zscore(curr_vect);
    zmatrix(i,:) = z_vect.';
end

% Define scalp zones
% Find electrodes in Frontal Area (FA)
dimlefteyes=0; %number of electrodes
index1=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
    if (chanlocs(1,k).theta > -60 && chanlocs(1,k).theta < 0) && (chanlocs(1,k).radius>0.45) && (chanlocs(1,k).radius<0.60)
        dimlefteyes=dimlefteyes+1; %count electrodes
        index1(1,dimlefteyes)=k;
    end
end
dimrighteyes=0; %number of electrodes
index2=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
    if (chanlocs(1,k).theta < 60 && chanlocs(1,k).theta > 0) && (chanlocs(1,k).radius>0.45) && (chanlocs(1,k).radius<0.60)
        dimrighteyes=dimrighteyes+1; %count electrodes
        index2(1,dimrighteyes)=k;
    end
end
dimcentereyes=0; %number of electrodes
index3=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
    % Eric Rawls modified 9/10/2026: manuscript alignment, added radius
    % maximum of 0.60
    if(abs(chanlocs(1,k).theta) < 20) && (chanlocs(1,k).radius>0.45) && (chanlocs(1,k).radius<0.60)
        dimcentereyes=dimcentereyes+1; %count electrodes
        index3(1,dimcentereyes)=k;
    end
end
%add in code in case E5 and E10 are deleted in a 64-chan net
if dimcentereyes==0
    for k=1:nchannels
        if(abs(chanlocs(1,k).theta) == 0) && (chanlocs(1,k).radius>0.39)
            dimcentereyes=dimcentereyes+1; %count electrodes
            index3(1,dimcentereyes)=k;
        end
    end
end
% Find electrodes in everywhere else
dimcenter=0; %number of electrodes
indexc=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
    if( (abs(chanlocs(1,k).theta)) > 35 && abs(chanlocs(1,k).theta) < 109 && (chanlocs(1,k).radius<0.45) )
        dimcenter=dimcenter+1; %count electrodes
        indexc(1,dimcenter)=k;
    end
end
dimbackleft=0; %number of electrodes
indexbl=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
    if (chanlocs(1,k).theta <= -109) && (chanlocs(1,k).radius<0.55) %electrodes are in FA
        dimbackleft=dimbackleft+1; %count electrodes
        indexbl(1,dimbackleft)=k;
    end
end
dimbackright=0; %number of electrodes
indexbr=zeros(1,nchannels); %indexes of electrodes
for k=1:nchannels
    if (chanlocs(1,k).theta >= 109) && (chanlocs(1,k).radius<0.55) %electrodes are in FA
        dimbackright=dimbackright+1; %count electrodes
        indexbr(1,dimbackright)=k;
    end
end
% Skip this detector if any required region has no electrodes.
% Do not calculate regional means with a zero denominator.
if any([dimlefteyes dimrighteyes dimcentereyes ...
        dimcenter dimbackleft dimbackright] == 0)
    warning('hnf_beall_blink_detection:MissingRegion', ...
        'Blink detection skipped: at least one required scalp region is empty.');
    return
end
% z-scores of the front and back
blinksCount = 1;
for i=1:n % for each topography
    % create FA electrodes vector
    zfrontlefteyes=zeros(1,dimlefteyes);
    for h=1:dimlefteyes
        zfrontlefteyes(1,h)=zmatrix(i,index1(1,h));
    end
    zfrontrighteyes=zeros(1,dimrighteyes);
    for h=1:dimrighteyes
        zfrontrighteyes(1,h)=zmatrix(i,index2(1,h));
    end
    zcentereyes=zeros(1,dimcentereyes);
    for h=1:dimcentereyes
        zcentereyes(1,h)=zmatrix(i,index3(1,h));
    end
    % create other electrodes vectors
    zcenter=zeros(1,dimcenter);
    for h=1:dimcenter
        zcenter(1,h)=zmatrix(i,indexc(1,h));
    end
    zbackl=zeros(1,dimbackleft);
    for h=1:dimbackleft
        zbackl(1,h)=zmatrix(i,indexbl(1,h));
    end
    zbackr=zeros(1,dimbackright);
    for h=1:dimbackright
        zbackr(1,h)=zmatrix(i,indexbr(1,h));
    end
    % Calculates the mean of each z-score vector
    % Sums the absolute value of each element in the vector,
    % then divides that by the number of elements.
    le = 0;
    for j=1:dimlefteyes
        le = le + abs(zfrontlefteyes(1,j));
    end
    le = le/dimlefteyes;
    re = 0;
    for j=1:dimrighteyes
        re = re + abs(zfrontrighteyes(1,j));
    end
    re = re/dimrighteyes;
    ce = 0;
    for j=1:dimcentereyes
        ce = ce + zcentereyes(1,j); %removed abs() around zcentereyes b/c they should be the same sign (1/29/19)
    end
    ce = ce/dimcentereyes;
    ch = 0;
    for j=1:dimcenter
        ch = ch + abs(zcenter(1,j));
    end
    ch = ch/dimcenter;
    bl = 0;
    for j=1:dimbackleft
        bl = bl + abs(zbackl(1,j));
    end
    bl = bl/dimbackleft;
    br = 0;
    for j=1:dimbackright
        br = br + abs(zbackr(1,j));
    end
    br = br/dimbackright;
    % add var calc
    vP    = var([zbackl zbackr],0,2);
    vAbsP = var(abs([zbackl zbackr]),0,2);
    vEyes = var([zfrontlefteyes zfrontrighteyes],0,2);
    % if one of the means of the eye vectors pass a certain threshold and
    % all of the back are below a certain threshold, then we reject
    if (abs(le) > 2 || abs(ce) > 2.5 || abs(re) > 2) %if high activity at front
        if (mean([abs(bl) abs(br) abs(ch)]) < 1) % if low activity in back
            % Eric Rawls modified 9/10/2026 for alignment with manuscript
            % equations
            if (vP < .15 && vEyes > vP) || (vAbsP < .075 && vEyes > vAbsP)
                blinks(blinksCount) = i; %#ok <AGROW>
                blinksCount = blinksCount + 1;
            end
        end
    end
end

%STEPH EDIT TO CODE 1/18/19
%find all ICs with bumps btw 5 and 15 Hz
%cwb=MARAcode_alphapeak(EEG);
%loop through artifacted ICs and remove any ICs that have bumps
remove_indices=[];
for ii = 1:length(blinks)
    for jj = 1:length(cwb)%components_with_bumps)
        if blinks(ii) == cwb(jj)%components_with_bumps(jj)
            remove_indices = [remove_indices ii]; %#ok <AGROW>
        end
    end
end
blinks(remove_indices)=[];


%STEPH EDIT TO ADJUST CODE 2/8/19
if ~isempty(blinks)
    % eye_electrodes is entered into the func below to check
    % that the eye artifacts doesn't go too far into the head
    eye_electrodes=[];
    for e=1:nchannels
        if(chanlocs(1,e).radius>0.45 && chanlocs(1,e).radius<0.54)
            if(abs(chanlocs(1,e).theta)<17)
                eye_electrodes(e) = e; %#ok <AGROW>
            end
        end
    end
    eye_electrodes = nonzeros(eye_electrodes);
    sec = hnf_spatial_info_eyes(blinks,eye_electrodes,topog,chanlocs,n);
    remove_indices=[];
    for ii = 1:length(blinks)
        for jj = 1:length(sec)%blinks_to_keep
            if blinks(ii) == sec(jj)
                remove_indices = [remove_indices ii]; %#ok <AGROW>
            end
        end
    end
    blinks(remove_indices)=[];
end