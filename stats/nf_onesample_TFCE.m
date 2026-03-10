function tfce = nf_onesample_TFCE(dataX,h0,chanlocs,faxis,taxis,nperm,par)
%
% GENERAL
% -------
% Uses routines from ept-TFCE to run a mass-univariate
% TFCE-corrected set of tests on input data series. Dimensions
% of input dataX MUST be subjects first. The dimensions may otherwise be
% subjectXtime, subjectsXfrequencies, subjectsXchannelsXtimes,
% subjectsXchannelsXfrequencies, or subjectsXchannelsXfrequenciesXtimes.
% The code will determine from the other arguments what the dimensions are.
%
% INPUTS
% ------
% 1) dataX - a matrix of data for statistical testing (must be subjects X
%   channels X frequencies X time, any dimensions except subject may be
%   missing)
% 2) h0 - null hypothesis value
% 3) chanlocs - channel locations (EEGLAB format; may be missing)
% 4) faxis - frequency axis, may be missing, in Hz
% 5) taxis - time points considered for testing, may be missing
% 6) nperm - number of permutations for cluster testing, recommend 10,000
% 7) par - parallelize? recommended. 0 = no, 1 = yes (default 0)
%
% OUTPUTS
% -------
% 1) tfce - a structure containing results of testing. Can be passed
% directly to nf_plottfce for plotting.
%

if nargin < 8
    varargin = {};
end

if isstruct(dataX) == 1
    if nargin < 7 || isempty(par) == 1
        par = 0;
    end

    if nargin < 6 || isempty(nperm) == 1
        nperm = 1000;
    end

    if nargin < 5 || isempty(taxis) == 1
        taxis = [];
    end

    if nargin < 4 || isempty(faxis) == 1
        faxis = [];
    end

    if nargin < 3 || isempty(chanlocs) == 1
        chanlocs = [];
    end

    if nargin < 2 || isempty(h0) == 1
        h0 = 0;
    end

    if isempty(chanlocs) == 1
        if isfield(dataX, 'chanlocs') == 1
            chanlocs = dataX.chanlocs;
        end
    end

    if isempty(faxis) == 1
        if isfield(dataX, 'freqs') == 1
            faxis = dataX.freqs;
        end
    end

    if isempty(taxis) == 1
        if isfield(dataX, 'times') == 1
            taxis = dataX.times;
        end
    end

    runner = @(X) nf_onesample_TFCE(X, h0, chanlocs, faxis, taxis, nperm, par);

    tfce = nf_tfstats_apply(dataX, runner, varargin{:});

    tfce.test = 'onesample';
    tfce.h0 = h0;
    tfce.nperm = nperm;
    tfce.par = par;

    return
end

% -----------------------------
% Numeric original behavior
% -----------------------------
if nargin < 7 || isempty(par) == 1
    par = 0;
end

if nargin < 6 || isempty(nperm) == 1
    nperm = 1000;
end

if nargin < 5 || isempty(h0) == 1
    h0 = 0;
    disp('no null hypothesis supplied; assuming 0.');
end

if nargin < 4 || isempty(taxis) == 1
    taxis = [];
    disp('no time axis supplied - assuming data have no time dimension.');
end

if nargin < 3 || isempty(faxis) == 1
    faxis = [];
    disp('no frequency axis supplied - assuming data have no frequency dimension.');
end

if nargin < 2 || isempty(chanlocs) == 1
    chanlocs = [];
    disp('no channel locations supplied - assuming single-channel data.');
end

if nargin < 1 || isempty(dataX) == 1
    error('at least data must be supplied...see help for dimensions and other inputs.');
end

if ~isempty(chanlocs)
    if par == 0
        tfc = ept_TFCE_onesample(dataX, zeros(size(dataX)) + h0, chanlocs, 'nperm', nperm);
    else
        tfc = ept_TFCE_onesample_par(dataX, zeros(size(dataX)) + h0, chanlocs, 'nperm', nperm);
    end
else
    if par == 0
        tfc = ept_TFCE_onesample(dataX, zeros(size(dataX)) + h0, [], 'nperm', nperm, 'flag_ft', 1);
    else
        tfc = ept_TFCE_onesample_par(dataX, zeros(size(dataX)) + h0, [], 'nperm', nperm, 'flag_ft', 1);
    end
end

tfce = struct();

tfce.power = dataX;
tfce.df = size(dataX, 1) - 1;

if ~isempty(taxis)
    tfce.times = taxis;
end

if ~isempty(faxis)
    tfce.freqs = faxis;
end

if ~isempty(chanlocs)
    tfce.chanlocs = chanlocs;
end

tfce.tstat = tfc.Obs;
tfce.p_vals = tfc.P_Values;
tfce.method = 'tfce';

end







