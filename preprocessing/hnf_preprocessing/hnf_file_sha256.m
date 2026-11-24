function digest = hnf_file_sha256(path)
% Hash exact vendor-file bytes. An empty result denotes an unreadable file.
digest = '';
fid = fopen(path,'rb');
if fid < 0, return; end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
try
    if exist('OCTAVE_VERSION','builtin')
        digest = hash('sha256',char(fread(fid,Inf,'*uint8').'));
    else
        md = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
        while ~feof(fid)
            bytes = fread(fid,1048576,'*uint8');
            if ~isempty(bytes), md.update(typecast(bytes,'int8')); end
        end
        digest = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
    end
catch
    digest = '';
end
end
