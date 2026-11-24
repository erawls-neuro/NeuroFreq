function o = hnf_la_opts(o, defaults)
% Internal: strict options, to catch misspelled scientific controls.
if nargin < 1 || isempty(o), o = struct; end
if ~isstruct(o) || ~isscalar(o)
    error('hnf:InvalidOptions','Options must be a scalar structure.');
end
names = fieldnames(o);
for k = 1:numel(names)
    if ~isfield(defaults,names{k})
        error('hnf:UnknownOption','Unknown option: %s.',names{k});
    end
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(o,names{k}), o.(names{k}) = defaults.(names{k}); end
end
end
