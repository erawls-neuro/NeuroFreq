function chanlocs = hnf_topo2all(chanlocs)
% Convert EEGLAB polar geometry using its standard head radius when unknown.
% eeg_checkset can add an empty sph_radius field to every channel.
% EEGLAB 2026.1 then computes mean([]), creating all-NaN X/Y/Z. Removing
% only the wholly empty field selects EEGLAB's defined 85 mm default.
if isfield(chanlocs,'sph_radius') && all(cellfun(@isempty,{chanlocs.sph_radius}))
    chanlocs = rmfield(chanlocs,'sph_radius');
end
chanlocs = convertlocs(chanlocs,'topo2all');
end
