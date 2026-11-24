function yes = nf_desktop()
% Octave uses Qt dialogs whereas MATLAB's interactive installer uses AWT.
if exist('OCTAVE_VERSION','builtin')
    yes = isguirunning();
else
    yes = usejava('desktop') && usejava('awt');
end
end
