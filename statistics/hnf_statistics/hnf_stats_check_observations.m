function hnf_stats_check_observations(tf)
% Trials pooled across participants are identifiable in nf_aggregate output.
% The simple tests cannot silently convert their null to a population null.
if ~isfield(tf,'trial_subjectindex') || ~isfield(tf,'aggtype') || ~strcmpi(tf.aggtype,'singletrial')
    return
end
subjects = unique(tf.trial_subjectindex);
if numel(subjects)>1 && numel(subjects)<numel(tf.trial_subjectindex)
    warning('NeuroFreq:PooledTrialInference', ...
        ['This analysis treats pooled trials as independent observations. For inference across participants, ' ...
         'model subject grouping or test subject-level summaries.']);
end
end
