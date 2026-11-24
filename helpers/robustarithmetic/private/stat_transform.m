function [raw, valid, usedCount] = stat_transform(raw, valid, method, fraction)
% Trim/winsorize in original units, then the caller chooses arithmetic scale.
count = sum(valid,1);
if strcmp(method,'trim')
    raw(~valid) = Inf;
    raw = sort(raw,1);
    discard = floor(count*fraction);
    valid = bsxfun(@gt,(1:size(raw,1))',discard) & ...
        bsxfun(@le,(1:size(raw,1))',count-discard);
else
    cut = stat_quantile(raw,valid,[fraction,1-fraction]);
    raw = bsxfun(@max,bsxfun(@min,raw,cut(2,:)),cut(1,:));
end
usedCount = sum(valid,1);
raw(~valid) = NaN;
end
