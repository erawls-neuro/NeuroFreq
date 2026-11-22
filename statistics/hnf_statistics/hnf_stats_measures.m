function tf = hnf_stats_measures(tf)
% Correlation/GLM wrappers have historically returned nested groups. Let
% the same code consume both nested inputs and current top-level measures.
names = fieldnames(tf);
groups = {'power','phase'};
patterns = {'^(power|pow_)','^(phase|ph_|itpc|wpli)'};
for i = 1:2
    group = groups{i};
    mapped = struct();
    if isfield(tf,group) && isstruct(tf.(group)), mapped = tf.(group); end
    for j = 1:numel(names)
        name = names{j};
        if ~isempty(regexpi(name,patterns{i},'once')) && isnumeric(tf.(name))
            if isfield(mapped,name)
                error('NeuroFreq:AmbiguousMeasure','Both %s.%s and %s are present.',group,name,name);
            end
            mapped.(name) = tf.(name);
        end
    end
    if ~isempty(fieldnames(mapped)), tf.(group) = mapped; end
end
end
