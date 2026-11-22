function X = hnf_stats_trial_grid(X, tf, chanlocs, freqs, times)
% Single-sensor transforms may store F x T x trials, with no explicit C.
% The fitted model uses C x F x T x trials; verify its axes against the TF
% metadata before inserting singleton dimensions.
if isempty(freqs) || isempty(times), return; end
nChan = max(1,numel(chanlocs));
n = [];
for field = {'ntrial','ntrials','ntrls'}
    if isfield(tf,field{1}) && isnumeric(tf.(field{1})) && isscalar(tf.(field{1}))
        n = tf.(field{1}); break
    end
end
if isempty(n) && isfield(tf,'trial_subjectindex'), n = numel(tf.trial_subjectindex); end
if isempty(n), n = numel(X)/(nChan*numel(freqs)*numel(times)); end
if ~isfinite(n) || n ~= fix(n) || n < 1
    error('NeuroFreq:InvalidTrialGrid','The TF axes do not match the number of observed samples.');
end
X = hnf_tf_grid(X,nChan,numel(freqs),numel(times),n,'statistical input');
end
