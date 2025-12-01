function dOut = nf_prepdata( EEG, center, detrend, taper, analytic)
% NF_PREPDATA    Prepares data cube or EEGLAB .set for TF computation.
%
% GENERAL
% -------
%
% Prepares data matrices for TF computation. Preparation includes removing
% quadratic polynomial from time series, cnetering data, cosine-square 
% tapering the outer 5% (per end) intervals of each epoch, and makes the
% signal analytic to reduce misspecification. All optional.
%
% Function applies equivalently to EEGLAB .set structures in memory or to
% 1/2/3D data matrices (dimensions must be channel, time, trial if a data
% matrix is entered).
%
% OUTPUT
% ------
% EEG - EEGLAB .set struct, OR
% data - 1/2/3D prepared data matrix
%
% INPUT
% -----
% 1) EEG - EEGLAB .set struct, OR
%      data - 1/2/3D prepared data matrix
% 2) center - center data? 0 or 1
% 3) detrend - remove polynomial? 0 or 1
% 4) taper - taper outer 5% of data points? 0 or 1
% 5) analytic: make signal analytic? 0 or 1
%
%


if nargin<5 || isempty(analytic)
    disp('no argument for analytic - making signal analytic (default)');
    analytic=1;
end
if nargin<4 || isempty(taper)
    disp('no argument for taper - tapering edges (default)');
    taper=1;
end
if nargin<3 || isempty(detrend)
    disp('no argument for detrend - detrending signal (default)');
    detrend=1;
end
if nargin<2 || isempty(center)
    disp('no argument for center - centering signal (default)');
    center=1;
end
if nargin<1 || isempty(EEG)
    error('data is a required input');
end
%get dimensions
if isstruct( EEG )
    flag = 1;
    disp('detected structure input - getting info from fields');
    nChan = EEG.nbchan;
    nTimes = EEG.pnts;
    nTrls = EEG.trials;
    dEEG = EEG.data;
else
    flag = 0;
    disp('input is not a structure - determining dimensions');
    dEEG = EEG;
    nChan = size(dEEG,1);
    nTimes = size(dEEG,2);
    if ndims(dEEG)==3
        nTrls = size(dEEG,3);
    else
        nTrls = 1;
    end
end

%begin preparation
disp('Prepping data');
%progress
prog=1;
fprintf(1,'Data preparation progress: %3d%%\n',prog);
for eloc=1:nChan
    for trl=1:nTrls
        %one stream of data
        data = squeeze(dEEG(eloc,:,trl));
        if center==1
            %center it
            data = data-mean(data);
        end
        if detrend==1
            %remove quadratic trends
            ind = 0:nTimes-1;
            r = polyfit(ind,data,2);
            fit = polyval(r,ind);
            data = data-fit;
        end
        if taper==1
            %cosine square taper
            w = tukeywin(numel(data),.1);
            data = data(:).*w(:);
        end
        if analytic==1
            data = hilbert(data);
        end
        %put it back in dEEG
        dEEG(eloc,:,trl) = data;
    end
    %track progress
    prog=100*(eloc/nChan);
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
end
fprintf(1,'\n');

%finalize output
if flag==1
    EEG.data = dEEG;
    dOut = EEG;
else
    dOut = dEEG;
end

end






