function Results = hnf_ept_TFCE_corr_par(data1, data2, chanlocs, varargin)
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
% This tests initially computes a r value for each channel and sample
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
%
% data1 should be a "Participants x Channels x Samples" or a
%  "Participants x Channels x Frequencies X Samples" variable
% - Channels must be in the same order as the corresponding electrodes file
% - Samples should be in chronological order
% data2 must be a matrix with one dimension of size "Participants".
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
e_loc = chanlocs;

%get info
nCh  = size(Data{1},2);

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
    disp('calculating channel neighbours...')
    ChN = ept_ChN2(e_loc);
end

% Create all variables in loop at their maximum size to increase performance
maxTFCE=zeros(nPerm,1);

% Calculate the actual T values of all data
disp('calculating actual differences')
%get data and regressors
D = Data{1};
regs = Data{2};    
%if TF, concatenate
siD = size(D);
siD=siD(2:end);
D = D(:,:);
%if regs are a vector, repmat them for computational efficiency
if isvector(regs)
    d=size(D); 
    d=d(2:end); 
    regs=repmat(regs(:),[1 d]);
else
    regs = regs(:,:);
end
%get t observed
D = (D - mean(D)) ./ std(D);
regs = (regs - mean(regs)) ./ std(regs);
T_Obs = sum(D .* regs) / (size(D, 1) - 1);
%reshape to size B
T_Obs = squeeze(reshape( T_Obs, [size(T_Obs,1), siD] ));
%do TFCE transformation
TFCE_Obs = tfce_transformation(T_Obs,ChN,E_H,flag_ft);
% Calculating the T value and TFCE of each different permutation
disp('calculating permutations...');
parfor i = 1:nPerm 
    %permute regs
    regs1 = regs(randperm(size(regs,1)),:);
    %get t permuted
    T_Perm = sum(D .* regs1) / (size(D,1) - 1);
    %reshape to size B
    T_Perm = reshape( T_Perm, [size(T_Perm,1), siD] );
    %do TFCE transformation
    TFCE_Perm = tfce_transformation(squeeze(T_Perm),ChN,E_H,flag_ft);
    % stores the maximum absolute value
    maxTFCE(i) = max(abs(TFCE_Perm(:))); 
end

% add observed maximum
edges = [maxTFCE;max(abs(TFCE_Obs(:)))];

%get p values
[~,bin]     = histc(abs(TFCE_Obs),sort(edges(:)));
P_Values    = 1-bin./(nPerm+2);

%output
Results.Obs                 = T_Obs;
Results.TFCE_Obs            = TFCE_Obs;
Results.maxTFCE             = sort(maxTFCE);
Results.P_Values            = P_Values;
%
toc

[min_P, idx] = min(Results.P_Values(:));
[Reg, Ch, S]      = ind2sub(size(Results.P_Values),idx);
max_Obs      = Results.Obs(idx);
display(['Peak significance found for regressor ' num2str(Reg) ', at channel ', num2str(Ch), ' at sample ', num2str(S), ': T(', num2str(size(Data{1},1)-1), ') = ', num2str(max_Obs), ', p = ', num2str(min_P)]);











end



























