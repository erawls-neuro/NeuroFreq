function entries = nf_eeglab_releases()
% The EEGLAB extension manager gets release records (including EEGLAB
% itself) from this SCCN endpoint. "link" is the distribution ZIP. EEGLAB's
% GitHub source ZIP omits submodules such as the bundled firfilt plugin.
% Cache successful requests briefly so checking each plugin costs one fetch.
persistent records fetched
if isempty(records) || (now-fetched)*86400 > 300
    url = ['https://sccn.ucsd.edu/eeglab/plugin_uploader/' ...
           'plugin_getcountall_nowiki_json.php?type=update&upload=0'];
    fresh = webread(url,weboptions('Timeout',20,'ContentType','json'));
    if ~isstruct(fresh) || ~all(isfield(fresh,{'plugin','curversion','link'}))
        error('NeuroFreq:EEGLABCatalog','Unexpected response from the SCCN release catalog.');
    end
    records = fresh;
    fetched = now;
end
entries = records;
end
