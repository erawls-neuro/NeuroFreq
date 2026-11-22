function url = nf_https(url)
url = regexprep(char(url), '^http://', 'https://', 'ignorecase');
if ~startsWith(lower(url), 'https://')
    error('NeuroFreq:DownloadURL', 'Expected an HTTPS download URL, got %s.', url);
end
end
