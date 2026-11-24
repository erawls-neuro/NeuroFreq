function labels = hnf_lower_labels(values)
% Case-insensitive electrode and file matching with cellstr storage.
if ~iscell(values), values = cellstr(values); end
labels = cellfun(@(label) lower(char(label)),values,'UniformOutput',false);
end
