function opts = hnf_aux_options(opts, defaults)
% Internal strict, scalar-struct options for probability/moment helpers.
if isempty(opts), opts = defaults; return; end
if ~isstruct(opts) || ~isscalar(opts)
    error('hnf:Options', 'Options must be a scalar struct.');
end
names = fieldnames(opts);
for k = 1:numel(names)
    if ~isfield(defaults, names{k})
        error('hnf:Options', 'Unknown option: %s.', names{k});
    end
    defaults.(names{k}) = opts.(names{k});
end
opts = defaults;
end
