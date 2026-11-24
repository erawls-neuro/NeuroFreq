function stat_number(value, low, high, name, integer)
if nargin < 5, integer = false; end
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < low || value > high || ...
        (integer && value ~= fix(value))
    error('hnf:statistics:OptionValue', 'Invalid %s option.', name);
end
end
