function text = hnf_scalar_text(value)
% Convert scalar metadata and event codes to text without a string array.
if ischar(value)
    text = value;
elseif exist('OCTAVE_VERSION','builtin') && (isnumeric(value) || islogical(value))
    if islogical(value)
        names = {'false','true'};
        text = names{1+value};
    else
        text = num2str(value);
    end
else
    text = char(string(value));
end
end
