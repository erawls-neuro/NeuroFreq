function tf = hnf_stats_layout(tf)
% Convert the stacked condition observations produced by nf_aggregate('avg')
% into subject x channel x frequency x time x condition. Use the stored row
% keys, so reordering the observations cannot change the subject pairing.
if ~isfield(tf,'aggtype') || ~strcmpi(tf.aggtype,'avg') || ...
        ~isfield(tf,'trial_subjectindex') || ~isfield(tf,'trial_condindex')
    return
end
if isfield(tf,'nf_stats_subject_first') && tf.nf_stats_subject_first, return; end
subject = tf.trial_subjectindex(:);
condition = tf.trial_condindex(:);
if isempty(subject) || numel(subject) ~= numel(condition) || ...
        any(~isfinite([subject;condition])) || any([subject;condition] < 1) || ...
        any(fix([subject;condition]) ~= [subject;condition])
    error('NeuroFreq:AverageLayout','Invalid subject or condition indices in average aggregate.');
end
n = max(subject);
k = max(condition);
[sorted,order] = sort(subject+(condition-1)*n);
if numel(sorted) ~= n*k || any(sorted ~= (1:n*k)')
    error('NeuroFreq:AverageLayout','Expected one observation per subject and condition.');
end
names = fieldnames(tf);
for j = 1:numel(names)
    name = names{j};
    if ~isempty(regexpi(name,'^(power|phase|pow_|ph_|itpc|wpli)','once')) && isnumeric(tf.(name))
        x = tf.(name);
        if size(x,4) ~= n*k
            error('NeuroFreq:AverageLayout','%s does not match trial_subjectindex.',name);
        end
        tf.(name) = permute(reshape(x(:,:,:,order), ...
            [size(x,1),size(x,2),size(x,3),n,k]),[4 1 2 3 5]);
    end
end
tf.nf_stats_subject_first = true;
end
