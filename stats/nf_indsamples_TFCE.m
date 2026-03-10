function tfce = nf_indsamples_TFCE(dataX,dataY,chanlocs,faxis,taxis,nperm,par)
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
% 2) dataY - a matrix of data for statistical testing (must be subjects X
%   channels X frequencies X time, any dimensions except subject may be
%   missing)
% 3) chanlocs - channel locations (EEGLAB format; may be missing)
% 4) faxis - frequency axis, may be missing, in Hz
% 5) taxis - time points considered for testing, may be missing
% 6) nperm - number of permutations for cluster testing, recommend 10,000
% 7) par - parallelize? recommended. 0 = no, 1 = yes (default 0)
%
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

    if isstruct(dataY) == 1
        runnerPair = @(XA, XB) nf_indsamples_TFCE(XA, XB, chanlocs, faxis, taxis, nperm, par);

        tfce = local_apply_pairwise_tf(dataX, dataY, runnerPair, varargin{:});

        tfce.test = 'indsamples';
        tfce.nperm = nperm;
        tfce.par = par;

        return
    end

    if isnumeric(dataY) == 1
        group = dataY(:);
        runnerSplit = @(Xobs, g) local_run_split(Xobs, g, chanlocs, faxis, taxis, nperm, par);

        tfce = local_apply_split_tf(dataX, group, runnerSplit, varargin{:});

        tfce.test = 'indsamples';
        tfce.nperm = nperm;
        tfce.par = par;

        return
    end

    error('nf_indsamples_TFCE:BadInputs', 'For struct input, dataY must be either a struct or a numeric 2-group vector.');
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

if nargin < 5 || isempty(taxis) == 1
    taxis = [];
end

if nargin < 4 || isempty(faxis) == 1
    faxis = [];
end

if nargin < 3 || isempty(chanlocs) == 1
    chanlocs = [];
end

if nargin < 2 || isempty(dataX) == 1 || isempty(dataY) == 1
    error('at least data must be supplied...see help for dimensions and other inputs.');
end

if ~isempty(chanlocs)
    if par == 0
        tfc = ept_TFCE_indsamples(dataX, dataY, chanlocs, 'nperm', nperm);
    else
        tfc = ept_TFCE_indsamples_par(dataX, dataY, chanlocs, 'nperm', nperm);
    end
else
    if par == 0
        tfc = ept_TFCE_indsamples(dataX, dataY, [], 'nperm', nperm, 'flag_ft', 1);
    else
        tfc = ept_TFCE_indsamples_par(dataX, dataY, [], 'nperm', nperm, 'flag_ft', 1);
    end
end

tfce = struct();

tfce.obs = {squeeze(mean(dataX, 1)), squeeze(mean(dataY, 1))};
tfce.sd = {squeeze(std(dataX, 0, 1)), squeeze(std(dataY, 0, 1))};
tfce.df = {size(dataX, 1) - 1, size(dataY, 1) - 1};

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

end

% -----------------------------
% Local helpers for struct mode
% -----------------------------

function out = local_run_split(Xobs, group, chanlocs, faxis, taxis, nperm, par)

u = unique(group(~isnan(group)));

if numel(u) ~= 2
    error('nf_indsamples_TFCE:BadGroup', 'Group vector must have exactly two unique values.');
end

g1 = u(1);
g2 = u(2);

X1 = Xobs(group == g1, :, :, :);
X2 = Xobs(group == g2, :, :, :);

out = nf_indsamples_TFCE(X1, X2, chanlocs, faxis, taxis, nperm, par);

end

function tfce = local_apply_split_tf(tfIn, group, runnerSplit, varargin)

runner = @(Xobs) runnerSplit(Xobs, group);

tfce = nf_tfstats_apply(tfIn, runner, varargin{:});

end

function tfce = local_apply_pairwise_tf(tfA, tfB, runnerPair, varargin)

groups = {'power', 'phase'};
slices = [];

if mod(numel(varargin), 2) ~= 0
    error('Name/value pairs required.');
end

for i = 1:2:numel(varargin)
    key = lower(char(varargin{i}));
    val = varargin{i + 1};

    if strcmp(key, 'groups') == 1
        groups = val;
    elseif strcmp(key, 'slices') == 1
        slices = val;
    end
end

tfce = struct();
tfce.results = struct();

for iG = 1:numel(groups)
    gname = groups{iG};

    if isfield(tfA, gname) == 0 || isfield(tfB, gname) == 0
        continue
    end

    if isstruct(tfA.(gname)) == 0 || isstruct(tfB.(gname)) == 0
        continue
    end

    fa = fieldnames(tfA.(gname));
    fb = fieldnames(tfB.(gname));

    if ~isequal(fa, fb)
        error('Group field mismatch between inputs.');
    end

    tfce.results.(gname) = struct();

    for iF = 1:numel(fa)
        fn = fa{iF};

        XA = tfA.(gname).(fn);
        XB = tfB.(gname).(fn);

        if ndims(XA) ~= ndims(XB)
            error('Shape mismatch between paired fields.');
        end

        sliceCount = 1;

        if ndims(XA) == 5
            sliceCount = size(XA, 5);
        end

        useSlices = 1:sliceCount;

        if isempty(slices) == 0
            useSlices = slices(:)';
        end

        res = struct();
        res.tfce = cell(numel(useSlices), 1);
        res.slice_index = useSlices(:);

        for iS = 1:numel(useSlices)
            sIdx = useSlices(iS);

            XA_s = XA;
            XB_s = XB;

            if ndims(XA) == 5
                XA_s = XA(:, :, :, :, sIdx);
                XB_s = XB(:, :, :, :, sIdx);
            end

            res.tfce{iS} = runnerPair(XA_s, XB_s);
        end

        tfce.results.(gname).(fn) = res;
    end
end

end









