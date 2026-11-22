function tfce = nf_indsamples_TFCE(dataX,dataY,chanlocs,faxis,taxis,nperm,par,varargin)
% Simple permutation tests assume independent, exchangeable observations.
% Pooled repeated trials across participants require a subject-aware design.

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



if isstruct(dataX) == 1
    hnf_stats_check_observations(dataX);
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
        hnf_stats_check_observations(dataY);
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

    error('nf_indsamples_TFCE:BadInputs', ...
        'For struct input, dataY must be either a struct or a numeric 2-group vector.');
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
    error('nf_indsamples_TFCE:DataInputRequired',...
        'at least data must be supplied...see help for dimensions and other inputs.');
end

if ~isempty(chanlocs)
    if par == 0
        tfc = hnf_ept_TFCE_indsamples(dataX, dataY, chanlocs, 'nperm', nperm);
    else
        tfc = hnf_ept_TFCE_indsamples_par(dataX, dataY, chanlocs, 'nperm', nperm);
    end
else
    if par == 0
        tfc = hnf_ept_TFCE_indsamples(dataX, dataY, [], 'nperm', nperm, 'flag_ft', 1);
    else
        tfc = hnf_ept_TFCE_indsamples_par(dataX, dataY, [], 'nperm', nperm, 'flag_ft', 1);
    end
end

tfce = struct();

tfce.obs = {squeeze(mean(dataX, 1)), squeeze(mean(dataY, 1))};
tfce.sd = {squeeze(std(dataX, 0, 1)), squeeze(std(dataY, 0, 1))};
% Preserve per-group degrees of freedom. The statistic is Welch t; the
% reported p-values are permutation maxima rather than a parametric t CDF.
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

disp('[nf_indsamples_TFCE]: tfce-corrected independent-samples t-test computation complete'); 

end

% -----------------------------
% Local helpers for struct mode
% -----------------------------

function out = local_run_split(Xobs, group, chanlocs, faxis, taxis, nperm, par)

if numel(group) ~= size(Xobs,1)
    error('nf_indsamples_TFCE:GroupMismatch','Supply one group assignment per observation.');
end
u = unique(group(~isnan(group)));

if numel(u) ~= 2
    error('nf_indsamples_TFCE:BadGroup', ...
        'Group vector must have exactly two unique values.');
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

% Use the same axis, slice and field rules as the one-sample path so
% channels or conditions cannot silently become the independent samples.
for axis = {'freqs','times'}
    name = axis{1};
    if isfield(tfA,name) && isfield(tfB,name) && ~isequal(tfA.(name)(:),tfB.(name)(:))
        error('nf_indsamples_TFCE:AxisMismatch','The two inputs have different %s coordinates.',name);
    end
end
if isfield(tfA,'chanlocs') && isfield(tfB,'chanlocs') && ...
        isfield(tfA.chanlocs,'labels') && isfield(tfB.chanlocs,'labels') && ...
        ~isequal({tfA.chanlocs.labels},{tfB.chanlocs.labels})
    error('nf_indsamples_TFCE:AxisMismatch','The two inputs have different channel labels or order.');
end
a = nf_tfstats_apply(tfA,@(x) x,varargin{:});
b = nf_tfstats_apply(tfB,@(x) x,varargin{:});
tfce = a;
groups = fieldnames(a.results);
if ~isequal(sort(groups),sort(fieldnames(b.results)))
    error('nf_indsamples_TFCE:GroupMismatch','The two inputs have different measure groups.');
end
for g = 1:numel(groups)
    group = groups{g};
    names = fieldnames(a.results.(group));
    if ~isequal(sort(names),sort(fieldnames(b.results.(group))))
        error('nf_indsamples_TFCE:GroupMismatch','The two inputs have different %s fields.',group);
    end
    for j = 1:numel(names)
        name = names{j};
        left = a.results.(group).(name);
        right = b.results.(group).(name);
        if ~isequal(left.slice_labels,right.slice_labels)
            error('nf_indsamples_TFCE:ConditionMismatch','Condition labels or ordering differ between the groups.');
        end
        for k = 1:numel(left.tfce)
            tfce.results.(group).(name).tfce{k} = runnerPair(left.tfce{k},right.tfce{k});
        end
    end
end
end
