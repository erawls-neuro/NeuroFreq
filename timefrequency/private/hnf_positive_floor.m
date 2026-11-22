function values = hnf_positive_floor(values)
% Replace negative/zero RID values using the smallest positive value in
% the same sensor/trial. A silent epoch has exactly zero power.
for trial = 1:size(values,4)
    for sensor = 1:size(values,1)
        surface = values(sensor,:,:,trial);
        floorValue = min(surface(surface > 0));
        if isempty(floorValue), floorValue = 0; end
        surface(surface < floorValue) = floorValue;
        values(sensor,:,:,trial) = surface;
    end
end
end
