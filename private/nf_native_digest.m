function hash = nf_native_digest(files, label)
% Hash the exact bytes and their names. Length prefixes distinguish file
% boundaries. Inputs are small source files, headers, or one MEX binary.
md = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
md.update(typecast(unicode2native(label,'UTF-8'),'int8'));
for i = 1:numel(files)
    file = javaObject('java.io.File',files{i});
    tag = sprintf('\n%d:%s:%d:',numel(char(file.getName())),char(file.getName()),file.length());
    md.update(typecast(unicode2native(tag,'UTF-8'),'int8'));
    bytes = javaMethod('readAllBytes','java.nio.file.Files',file.toPath());
    if ~isempty(bytes), md.update(bytes); end
end
hash = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
