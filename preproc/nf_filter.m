function EEG = nf_filter( EEG, lopass, hipass, notch, resample )

if nargin < 5 || isempty(resample)
    disp('resampling to 250 Hz (default)');
    resample = 250;
end
if nargin < 4 || isempty(notch)
    disp('notch filtering at 60 Hz (default)');
    notch = 60;
end
if nargin < 3 || isempty(hipass)
    disp('using hi pass filter at 0.3 Hz (default)');
    hipass = 0.3;
end
if nargin < 2 || isempty(lopass)
    disp('using low pass filter at 35 Hz (default)');
    lopass = 35;
end
if nargin < 1 || isempty(EEG)
    error('at least an EEG set is required!');
end

%notch filter
EEG = pop_eegfiltnew( EEG, 'locutoff', notch-2, 'hicutoff', notch+2, 'revfilt', 1);
%high-pass filter
EEG = pop_eegfiltnew( EEG, 'locutoff', hipass);
%low-pass filter
EEG = pop_eegfiltnew( EEG, 'hicutoff', lopass);
%downsample
EEG = pop_resample(EEG, resample);

end