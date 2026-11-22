function answer = nf_choose(question, titleText, choices, default)
% Every prompt has an explicit safe default for a dismissed dialog or batch
% session. A terminal/batch run never waits for an invisible dialog.
answer = default;
if ~usejava('desktop') || ~usejava('awt')
    fprintf('%s\n  No desktop dialog available. Choice: %s\n', question, default);
    return
end
if numel(choices) <= 3
    chosen = questdlg(question, titleText, choices{:}, default);
    if ~isempty(chosen), answer = chosen; end
else
    [selected, ok] = listdlg('Name', titleText, 'PromptString', ...
        strsplit(question, sprintf('\n')), 'SelectionMode', 'single', ...
        'ListString', choices, 'InitialValue', find(strcmp(default, choices), 1), ...
        'ListSize', [510 140]);
    if ok, answer = choices{selected}; end
end
end
