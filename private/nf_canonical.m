function absolute = nf_canonical(filename)
% Resolve dot components, symlinks and relative path entries before comparing.
if exist('OCTAVE_VERSION','builtin')
    % Resolve existing parent links even for a cache directory that has not
    % been created yet. Plain NF execution does not require Octave's JVM.
    filename = make_absolute_filename(char(filename));
    absolute = canonicalize_file_name(filename);
    if isempty(absolute)
        [parent,name,ext] = fileparts(filename);
        if strcmp(parent,filename), absolute = filename; return; end
        absolute = fullfile(nf_canonical(parent),[name ext]);
    end
    return
end
file = javaObject('java.io.File', char(filename));
if ~file.isAbsolute(), file = javaObject('java.io.File', fullfile(pwd, char(filename))); end
absolute = char(file.getCanonicalPath());
end
