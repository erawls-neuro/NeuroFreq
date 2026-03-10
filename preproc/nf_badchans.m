%
% Automatically detects and interpolates bad channels in EEG recordings.
% Works adaptively to reject no more than a specified number of channels.
%
function EEG = nf_badchans(EEG,maxbad,interp)

% set interp default
if nargin<3 || isempty(interp)
    interp = 1;
end
%set max default
if nargin<2 || isempty(maxbad)
    maxbad = floor(EEG.nbchan/10);
end
%copy EEG
nC = EEG.nbchan;
%average reference
EEG = pop_reref(EEG,[]);
%original channels
ogChan = EEG.chanlocs;
EEG.etc.ogchan = ogChan;
%flatline channels
EEG = clean_flatlines( EEG, 5);
%start iteration
thresh = 0.8;
%channel correlation for bads detection
try
    EEG1 = clean_channels( EEG, thresh );
    chanret = EEG1.nbchan;
catch
    chanret = 0;
end
%iterate and find an acceptable cutoff
while (nC-chanret)>maxbad
    thresh=thresh-.05;
    disp(['threshold=' num2str(thresh)]);
    try
        if thresh<0
            break
        end
        EEG1 = clean_channels( EEG, thresh );
        chanret = EEG1.nbchan;
    catch
        chanret = 0;
        continue
    end
end
EEG = clean_channels(EEG,thresh);
%record missing percent
EEG.etc.badchans = nC-EEG.nbchan;
if interp==1
    %interpolate missing channels
    EEG = pop_interp(EEG, ogChan, 'sphericalKang');
    %average reference again following interpolation
    EEG = pop_reref(EEG,[]);
else
    % record original
    EEG.etc.ogchan = ogChan;
end

end