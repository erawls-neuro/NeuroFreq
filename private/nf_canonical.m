function absolute = nf_canonical(filename)
% Resolve dot components, symlinks and relative path entries before comparing.
file = javaObject('java.io.File', char(filename));
if ~file.isAbsolute(), file = javaObject('java.io.File', fullfile(pwd, char(filename))); end
absolute = char(file.getCanonicalPath());
end
