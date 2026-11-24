function unpacked = nf_fetch_package(url, stage)
% Only extract into a new staging directory. Inspect all member names before
% extraction so even a malformed archive cannot write above that directory.
archive = fullfile(stage, 'download.zip');
unpacked = fullfile(stage, 'unpacked');
if isfile(archive) || isfolder(unpacked)
    error('NeuroFreq:StageCollision', 'The temporary download directory is already occupied.');
end
url = nf_https(url);
if exist('OCTAVE_VERSION','builtin')
    urlwrite(url,archive);
else
    websave(archive, url, weboptions('Timeout', 180));
end
zipFile = javaObject('java.util.zip.ZipFile', archive);
closeZip = onCleanup(@() zipFile.close());
entries = zipFile.entries();
while entries.hasMoreElements()
    entry = entries.nextElement();
    name = char(entry.getName());
    if ~isempty(regexp(strrep(name, '\', '/'), '(^/|^[a-zA-Z]:|(^|/)\.\.(/|$))', 'once'))
        error('NeuroFreq:InvalidArchive', 'Unsafe path in downloaded archive: %s.', name);
    end
end
clear closeZip
unzip(archive, unpacked);
end
