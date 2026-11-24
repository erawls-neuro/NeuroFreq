function value = stat_choice(value, choices, name)
if isstring(value) && isscalar(value), value = char(value); end
if ~ischar(value) || ~isrow(value) || ~any(strcmp(value,choices))
    error('hnf:statistics:OptionValue', 'Invalid %s option.', name);
end
end
