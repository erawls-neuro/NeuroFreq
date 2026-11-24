function id = hnf_subject_codes(subject)
% Numeric category codes in category order, with missing labels as NaN.
if ~isa(subject,'categorical'), subject = categorical(subject); end
if exist('OCTAVE_VERSION','builtin')
    % Older statistics releases do not support grp2idx on categoricals.
    names = categories(subject);
    names = names(~cellfun(@isempty,names));
    [~,id] = ismember(cellstr(subject(:)),names);
    id(id == 0 | isundefined(subject(:))) = NaN;
else
    id = double(grp2idx(subject(:)));
end
end
