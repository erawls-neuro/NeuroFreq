function TF = nf_makebehavior( TF, tab )
% NF_MAKEBEHAVIOR    Make a NeuroFreq-style 'behavior' file from a table of task stimuli and behaviors.
%
% GENERAL
% -------
% Make a NeuroFreq-style 'behavior' file from a table of task stimuli and
% behaviors. Table must have column headers and a number of rows equal to
% the number of trials in the input EEG.
%
%
% OUTPUT
% ------
% TF - TF structure including behavior.
%
% INPUT
% -----
% 1) TF - structure output by nf_tftransform.m and tf_fun functions
% 2) tab - table of task conditions, stimuli, and behaviors
%
%


behavior = readtable( tab );
vars = behavior.Properties.VariableNames;
behav=[];
if isfield(TF,'conds')
    if TF.conds ~= size(behavior,1)
        error('behavior must have the same number of entries as TF has trials');
    end
else
    if TF.ntrls ~= size(behavior,1)
        error('behavior must have the same number of entries as TF has trials');
    end
end
for h=1:size(behavior,1)
    for j=1:length(vars)
        behav(h).(vars{j}) = table2array(behavior(h,j));
    end
end
TF.behavior = behav;


end