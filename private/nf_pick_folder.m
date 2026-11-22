function folder = nf_pick_folder(prompt, initial)
% uigetdir returns numeric zero when canceled. All callers receive char text.
folder = '';
if ~usejava('desktop') || ~usejava('awt'), return; end
if nargin < 2 || ~isfolder(initial), initial = pwd; end
selected = uigetdir(initial, prompt);
if ischar(selected), folder = nf_canonical(selected); end
end
