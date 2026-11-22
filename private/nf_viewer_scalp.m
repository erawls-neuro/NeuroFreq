function [Xi,Yi,Zi] = nf_viewer_scalp(x,y,values,xi,yi,angular)
% Duplicate sensor coordinates provide one interpolation site. For phase,
% interpolate unit phasors so neighbouring +/-pi values meet near pi.
[sites,~,group] = unique([x(:),y(:)],'rows');
if angular, values = exp(1i*values); end
observed = accumarray(group,values(:),[],@mean);
[Xi,Yi] = meshgrid(xi,yi);
if size(sites,1)<3 || rank([ones(size(sites,1),1) sites])<3
    Zi = nan(size(Xi));
    return
end
Zi = griddata(sites(:,1),sites(:,2),real(observed),Xi,Yi,'v4');
if angular
    Zi = Zi+1i*griddata(sites(:,1),sites(:,2),imag(observed),Xi,Yi,'v4');
    undefined = abs(Zi)<8*eps;
    Zi = angle(Zi);
    Zi(undefined) = NaN;
end
end
