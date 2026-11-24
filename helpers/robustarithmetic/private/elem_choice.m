function value = elem_choice(value, choices, id)
if isstring(value) && isscalar(value), value = char(value); end
if ~ischar(value) || ~isrow(value)
    error([id ':Options'],'Policy values must be character rows or scalar strings.');
end
value = lower(value);
if ~any(strcmp(value,choices))
    error([id ':Options'],'Unsupported policy: %s.',value);
end
end
