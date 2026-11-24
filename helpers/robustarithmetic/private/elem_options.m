function o = elem_options(o, supplied, id)
% Internal: strict scalar options, with a cheap empty-options path.
if isempty(supplied), return; end
if ~isstruct(supplied) || ~isscalar(supplied)
    error([id ':Options'],'Options must be a scalar structure.');
end
names = fieldnames(supplied);
for k = 1:numel(names)
    if ~isfield(o,names{k})
        error([id ':Options'],'Unknown option: %s.',names{k});
    end
    o.(names{k}) = supplied.(names{k});
end
if isfield(o,'inputPolicy')
    o.inputPolicy = elem_choice(o.inputPolicy,{'propagate','fill','error'},id);
    if ~isnumeric(o.fillValue) || ~isscalar(o.fillValue) || ~isfinite(o.fillValue)
        error([id ':Options'],'fillValue must be a finite numeric scalar.');
    end
end
if isfield(o,'overflow')
    o.overflow = elem_choice(o.overflow,{'error','saturate','allow'},id);
end
if isfield(o,'underflow')
    o.underflow = elem_choice(o.underflow,{'allow','floor','error'},id);
end
end
