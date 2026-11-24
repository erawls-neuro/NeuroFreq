function Results = hnf_ept_TFCE_indsamples_par(data1, data2, chanlocs, varargin)
%
% ATTENTION: THIS FUNCTION IS NOT THE ORIGINAL ept_TFCE.m
% IT HAS BEEN MODIFIED BY ERIC RAWLS 2023.
% MANY ELEMENTS OF THE CODE WERE DELETED FOR THE NEUROFREQ USE CASE.
%
% ORIGINAL DOCUMENTATION FOLLOWS:
%
% ========================================================================
%
% EEG Permutation Test (egt) using Threshold-Free Cluster Enhancement
% Copyright(C) 2012  Armand Mensen (14.12.2010)
%
% [Description]
% This tests initially computes a T-Value for each channel and sample
% These values are then enhanced or reduced depending on the size of
% the T-value and the size of the cluster the channel belongs to at
% various thresholds...
%
% TFCE is essentially the sum of the clusters to the power of 0.5 multiplied
% by the height of the current threshold squared
%
% [Input]
% Participant data should be organised into two factors
%     e.g. Dataset 1 is controls while Dataset 2 is patients
% Electrode locations file created using eeglab is required to
% calculate channel neighbours
%
% Analysis Types
% i = independent sample T-Test
%
% data1 should be a "Participants x Channels x Samples" or a
%  "Participants x Channels x Frequencies X Samples" variable
% - Channels must be in the same order as the corresponding electrodes file
% - Samples should be in chronological order
% data2 must have same dimensions as data1 but be independent
% ElecFile is not required. If not specified, must specify flag_ft = 1
%
%


%defaults
E_H        = [0.66 2]; % default parameters of E and H
nPerm      = 1000; % default number of permutations
flag_ft    = 0;
ChN = [];

% check for data
if nargin < 2
    error('at least data are required');
end

%process args
for i = 1:2:length(varargin)
    Param = varargin{i};
    Value = varargin{i+1};
    if ~ischar(Param)
        error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
        case 'e_h'
            E_H         = Value;
        case 'nperm'
            nPerm       = Value;
        case 'flag_ft'
            flag_ft     = Value;
        case 'chn'
            ChN         = Value;
        otherwise
            display (['Unknown parameter setting: ' Param])
    end
end

%set things
Data{1} = double(data1);
Data{2} = double(data2);
aData = [Data{1}; Data{2}];
e_loc = chanlocs;

%get info
nA   = size(Data{1},1);
nB   = size(Data{2},1);
nCh  = size(Data{1},2);
if min(nA,nB) < 2
    error('NeuroFreq:TooFewObservations', 'Each independent group needs at least two observations.');
end

% -- Error Checking -- %%
% Check Location File for Number of Channels
if ~flag_ft
    if ~isequal(nCh, length(e_loc))
        error ('Number of channels in data does not equal that of locations file')
    end
end

tic; % Start the timer for the entire analysis
% Calculate the channels neighbours... using the modified version ChN2
if ~flag_ft && isempty(ChN)
    disp('calculating channel neighbours')
    if nCh == 1, ChN = 1; else, ChN = ept_ChN2(e_loc); end
end

% Create all variables in loop at their maximum size to increase performance
maxTFCE = zeros(nPerm,1);




















% Calculate the actual T values of all data
disp('calculating actual differences')
% Calculate different T-values for independent and dependent groups
% Keep MATLAB's randomization order. The native kernel reads the
% selected observations without allocating two permuted group arrays.
useNative = numel(aData) >= 1024 && isreal(aData) && ~issparse(aData) && ...
    exist('hnf_mex_available','file') == 2 && hnf_mex_available('hnf_permuted_t_mex');
% Workers have independent MATLAB search paths and persistent state.
% Send the approved filename; a worker with a different or missing MEX
% computes the very same random assignment using the MATLAB formula.
workerMex = '';
if useNative, workerMex = which('hnf_permuted_t_mex'); end
T_Obs = hnf_permuted_t_map(aData, 1:(nA+nB), useNative, nA);
T_Obs = squeeze(T_Obs);
%do TFCE transformation
TFCE_Obs = hnf_tfce_transformation(T_Obs,ChN,E_H,flag_ft);

% Calculating the T value and TFCE of each different permutation
disp('calculating permutations...');
parfor i  = 1:nPerm           
    r_perm      = randperm(nA+nB); % Consider using Shuffle mex here (50-85% faster)...
    T_Perm = hnf_permuted_t_map(aData, r_perm, workerMex, nA);
    T_Perm = squeeze(T_Perm);
    % TFCE transformation...
    TFCE_Perm = hnf_tfce_transformation(T_Perm, ChN,E_H,flag_ft);
    % stores the maximum absolute value
    maxTFCE(i) = max(abs(TFCE_Perm(:)));       
end










P_Values = hnf_tfce_pvalues(TFCE_Obs, maxTFCE);

%output
Results.Obs                 = T_Obs;
Results.TFCE_Obs            = TFCE_Obs;
Results.maxTFCE             = sort(maxTFCE);
Results.P_Values            = P_Values;
%
toc
[min_P, idx] = min(Results.P_Values(:));
fprintf('Peak Welch t = %.4g, TFCE p = %.4g.\n', Results.Obs(idx), min_P);
end
