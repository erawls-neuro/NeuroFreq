function opts = stat_options(opts, defaults)
% Private option validation; unknown fields are errors, never silent typos.
if isempty(opts), opts = struct; end
if ~isstruct(opts) || ~isscalar(opts)
    error('hnf:statistics:Options', 'opts must be a scalar struct.');
end
names = fieldnames(opts);
for k = 1:numel(names)
    if ~isfield(defaults, names{k})
        error('hnf:statistics:Option', 'Unknown option: %s.', names{k});
    end
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts, names{k}), opts.(names{k}) = defaults.(names{k}); end
end
if isfield(opts,'nonfinite')
    opts.nonfinite = stat_choice(opts.nonfinite, {'error','omit','propagate'}, 'nonfinite');
end
if isfield(opts,'overflow')
    opts.overflow = stat_choice(opts.overflow, {'error','allow','saturate'}, 'overflow');
end
end
