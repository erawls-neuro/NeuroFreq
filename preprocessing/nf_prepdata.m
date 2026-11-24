function dOut = nf_prepdata( EEG, center, detrend, taper, analytic, verbose)
% NF_PREPDATA    Prepares data cube or EEGLAB .set for TF computation.
%
% GENERAL
% -------
%
% Prepares data matrices for TF computation. Preparation includes removing
% quadratic polynomial from time series, centering data, cosine-square
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
% 6) verbose: print verbose output? 0 or 1
%


disp('[nf_prepdata]: beginning data preparation for tf.');

if nargin<6 || isempty(verbose)
    verbose = 0;
end
if nargin<5 || isempty(analytic)
    if verbose==1
        disp('[nf_prepdata]: no argument for analytic - making signal analytic (default)');
    end
    analytic=1;
end
if nargin<4 || isempty(taper)
    if verbose==1
        disp('[nf_prepdata]: no argument for taper - tapering edges (default)');
    end
    taper=1;
end
if nargin<3 || isempty(detrend)
    if verbose==1
        disp('[nf_prepdata]: no argument for detrend - detrending signal (default)');
    end
    detrend=1;
end
if nargin<2 || isempty(center)
    if verbose==1
        disp('[nf_prepdata]: no argument for center - centering signal (default)');
    end
    center=1;
end
if nargin<1 || isempty(EEG)
    error('nf_prepdata:DataRequired','data is a required input');
end
% A standalone vector is one time series, in either orientation.
columnInput = ~isstruct(EEG) && iscolumn(EEG);
if columnInput, EEG = EEG.'; end
%get dimensions
if isstruct( EEG )
    flag = 1;
    if verbose==1
        disp('[nf_prepdata]: detected structure input - getting info from fields');
    end
    nChan = EEG.nbchan;
    nTimes = EEG.pnts;
    nTrls = EEG.trials;
    dEEG = EEG.data;
else
    flag = 0;
    if verbose==1
        disp('[nf_prepdata]: input is not a structure - determining dimensions');
    end
    dEEG = EEG;
    nChan = size(dEEG,1);
    nTimes = size(dEEG,2);
    if ndims(dEEG)==3
        nTrls = size(dEEG,3);
    else
        nTrls = 1;
    end
end

% Fractional detrending and Hilbert coefficients need floating-point storage.
if ~isfloat(dEEG), dEEG = double(dEEG); end
% Center in the original precision, exactly as the reference routine does.
% Then project in double and return to the input precision before applying the
% Hilbert transform across time, in one epoch-sized batch to bound FFT memory.
useNative = center==1 && detrend==1 && ...
    (isa(dEEG,'double') || isa(dEEG,'single')) && isreal(dEEG) && ~issparse(dEEG) && ...
    nTimes >= 3 && ndims(dEEG) <= 3 && ...
    isequal([size(dEEG,1) size(dEEG,2) size(dEEG,3)], [nChan nTimes nTrls]) && ...
    exist('hnf_prep_traces_mex','file') == 3 && ...
    (~exist('hnf_mex_available','file') || hnf_mex_available('hnf_prep_traces_mex')) && ...
    all(isfinite(dEEG(:))) && max(abs(dEEG(:))) <= realmax('double')/(4*double(nTimes));
if useNative
    weights = [];
    if taper==1, weights = tukeywin(nTimes,.1); end
    for trl=1:nTrls
        % Assigning the first analytic epoch promotes the entire output to
        % complex. The unread epochs still contain their original real input.
        epoch = real(dEEG(:,:,trl));
        epoch = bsxfun(@minus,epoch,mean(epoch,2));
        try
            epoch = hnf_prep_traces_mex(double(epoch), weights);
        catch err
            if trl ~= 1 || ~exist('hnf_mex_load_error','file') || ~hnf_mex_load_error(err)
                rethrow(err);
            end
            useNative = false;
            break
        end
        if isa(dEEG,'single'), epoch = single(epoch); end
        if analytic==1, epoch = hilbert(epoch,[],2); end
        dEEG(:,:,trl) = epoch;
    end
end
if ~useNative
    %begin preparation
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
                r = polyfit(ind,data,min(2,nTimes-1));
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
    end

end

%finalize output
if flag==1
    EEG.data = dEEG;
    dOut = EEG;
else
    dOut = dEEG;
    if columnInput, dOut = dOut.'; end
end

disp('[nf_prepdata]: data preparation for tf complete.');

end






